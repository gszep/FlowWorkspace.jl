function load(path::String, sample::Union{EzXML.Node,Nothing}; kwargs...)

    #######################################
    params, data = loadFCS(path; kwargs...)
    params = sample === nothing ? getMetaData(params) : findchannels(sample)

    ###################################### data with human readable channel names
    channels = dropmissing(params, "E")
    transform!(channels, ["N", "S"] => ByRow((N, S) -> ~ismissing(S) ? S != "" ? S : N : N) => "name")
    channelMap = Dict(param.N => ismissing(param.S) ? param.N : param.S for param ∈ eachrow(params))

    ###################################### biexponential transformation
    data = DataFrame(data, channels.name)
    transformationFunctions = transforms(sample; channelMap=channelMap)
    transform!(data, [channel => ByRow(func) => channel for (channel, func) ∈ transformationFunctions]...)

    ################################ load metadata from workspace
    groups = loadGroups(data, sample)
    gating = gatingGraph(sample; transforms=transformationFunctions, channelMap=channelMap)

    labels = gate(data, gating)
    return data, labels, groups, gating
end


function load(path::String; files::String=joinpath(dirname(path), "*.fcs"), cols::Symbol=:setequal, kwargs...)

    @assert(isfile(path), "no such file: $path")
    @assert(length(glob(files)) ≠ 0, "no FCS files found using pattern: $files")

    try
        workspace = root(readxml(path))
    catch error
        throw(LoadError(path, error.line, error))
    end

    workspace = root(readxml(path))
    datasets = map(dataset -> basename(dataset["uri"]), findall("//DataSet", workspace))
    @assert(length(datasets) == length(unique(datasets)), "FCS files under a workspace must have unique names")

    data, labels, groups = DataFrame(), DataFrame(), DataFrame()
    gatings = Dict()

    for path ∈ glob(files)
        @info "Loading $path"
        sample = findsample(path, workspace)

        fcs, label, group, gating = load(path, sample; kwargs...)

        append!(data, fcs, cols=cols)
        append!(labels, label, cols=:union)
        append!(groups, group, cols=:union)
        gatings[path] = gating
    end

    map(name -> replace!(labels[!, name], missing => false), names(labels))
    map(name -> replace!(groups[!, name], missing => false), names(groups))

    disallowmissing!(labels)
    disallowmissing!(groups)

    transform!(labels, AsTable(filter(x -> x ≠ "Unlabelled", names(labels))) => ByRow(x -> ~any(x)) => "Unlabelled")
    transform!(groups, AsTable(filter(x -> x ≠ "Ungrouped", names(groups))) => ByRow(x -> ~any(x)) => "Ungrouped")

    return data, labels, groups, gatings
end

function findsample(uri::String, workspace::EzXML.Node)
    uri = escapeuri(uri)

    for code ∈ [("%2F" => "/"), ("%5C" => "\\"), ("%3A" => ":"), ("%2C" => ",")]
        uri = replace(uri, code)
    end

    sample = findfirst("//DataSet[contains(@uri,'$(basename(uri))')]", workspace)
    sample === nothing && @warn("""no metadata found""")
    return sample
end

function findchannels(sample::EzXML.Node)

    params = DataFrame("keyword" => findall("..//Keywords/Keyword[starts-with(@name,'\$P') and contains(@name,'N')]", sample))
    keys = ["G", "R", "V", "S", "B", "N", "E", "AR"]

    for key ∈ keys
        transform!(params, "keyword" => ByRow(x -> getparam(chop(x["name"], head=0, tail=1) * key, sample)) => key)
    end

    transform!(params, "keyword" => ByRow(x -> parse(Int, chop(x["name"], head=2, tail=1))) => "index")
    sort!(params, "index")

    transform!(params, ["N", "S"] => ByRow((N, S) -> ~ismissing(S) ? S : stripname(N, params.N)) => "S")
    return select!(params, keys)
end

function stripname(name::String, names::AbstractArray{<:String})
    names = filter(x -> occursin(x, name) & (x ≠ name), names)
    return isempty(names) ? missing : first(names)
end

function getparam(param::String, sample::EzXML.Node)
    keyword = findfirst("""..//Keywords/Keyword[@name='$param']""", sample)
    return keyword !== nothing ? keyword["value"] != "" ? keyword["value"] : missing : missing
end


"""
    loadFCSHeader(fn::String)::Tuple{Vector{Int}, Dict{String,String}}
Efficiently extract data offsets and keyword dictionary from an FCS file.
"""
function loadFCSHeader(fn::String)::Tuple{Vector{Int},Dict{String,String}}
    open(fn) do io
        offsets = FCSFiles.parse_header(io)
        params = FCSFiles.parse_text(io, offsets[1], offsets[2])
        FCSFiles.verify_text(params)
        (offsets, params)
    end
end

"""
    getFCSSize(offsets, params)::Tuple{Int,Int}
Convert the offsets and keywords from an FCS file to cell and parameter count,
respectively.
"""
function getFCSSize(offsets, params)::Tuple{Int,Int}
    nData = parse(Int, params["\$TOT"])
    nParams = parse(Int, params["\$PAR"])

    if params["\$DATATYPE"] != "F"
        @error "Only float32 FCS files are currently supported"
        error("Unsupported FCS format")
    end

    beginData = parse(Int, params["\$BEGINDATA"])
    endData = parse(Int, params["\$ENDDATA"])

    #check that the $TOT and $PAR look okay
    if !(offsets[3] == 0 && offsets[4] == 0) && (
        (
            1 + offsets[4] - offsets[3] != nData * nParams * 4 &&
            offsets[4] - offsets[3] != nData * nParams * 4
        ) ||
        offsets[3] != beginData ||
        offsets[4] != endData
    )
        @warn "Data size mismatch, FCS is likely broken."
    end

    return (nData, nParams)
end

"""
    loadFCSSizes(fns::Vector{String})
Load cell counts in many FCS files at once. Useful as input for `slicesof`.
"""
function loadFCSSizes(fns::Vector{String})::Vector{Int}
    [(
        begin
            o, s = loadFCSHeader(fn)
            getFCSSize(o, s)[1]
        end
    ) for fn in fns]
end

"""
    loadFCS(fn::String; applyCompensation::Bool=true)::Tuple{Dict{String,String}, Matrix{Float64}}
Read a FCS file. Return a tuple that contains in order:
- dictionary of the keywords contained in the file
- raw column names
- prettified and annotated column names
- raw data matrix
If `applyCompensation` is set, the function parses and retrieves a spillover
matrix (if any valid keyword in the FCS is found that would contain it) and
applies it to compensate the data.
"""
function loadFCS(
    fn::String;
    applyCompensation::Bool = true,
)::Tuple{Dict{String,String},Matrix{Float64}}
    fcs = FileIO.load(fn)
    meta = getMetaData(fcs.params)
    data = hcat(map(x -> Vector{Float64}(fcs[x]), meta[:, :N])...)
    if applyCompensation
        spill = getSpillover(fcs.params)
        if spill != nothing
            names, mtx = spill
            cols = indexin(names, meta[:, :N])
            if any(cols .== nothing)
                @error "Unknown columns in compensation matrix" names cols
                error("Invalid compensation matrix")
            end
            compensate!(data, mtx, Vector{Int}(cols))
        end
    end
    return (fcs.params, data)
end

"""
    loadFCSSet(name::Symbol, fns::Vector{String}, pids=workers(); applyCompensation=true, postLoad=(d,i)->d)::Dinfo
This runs the FCS loading machinery in a distributed way, so that the files
`fns` (with full path) are sliced into equal parts and saved as a distributed
variable `name` on workers specified by `pids`.
`applyCompensation` is passed to loadFCS function.
See `slicesof` for description of the slicing.
`postLoad` is applied to the loaded FCS file data (and the index) -- use this
function to e.g. filter out certain columns right on loading, using `selectFCSColumns`.
The loaded dataset can be manipulated by the distributed functions, e.g.
- `dselect` for removing columns
- `dscale` for normalization
- `dtransform_asinh` (and others) for transformation
- etc.
"""
function loadFCSSet(
    name::Symbol,
    fns::Vector{String},
    pids = workers();
    applyCompensation = true,
    postLoad = (d, i) -> d,
)::Dinfo
    slices = slicesof(loadFCSSizes(fns), length(pids))
    dmap(
        slices,
        (slice) -> Base.eval(
            Main,
            :(
                begin
                    $name = vcollectSlice(
                        (i) -> last(
                            $postLoad(
                                loadFCS($fns[i]; applyCompensation = $applyCompensation),
                                i,
                            ),
                        ),
                        $slice,
                    )
                    nothing
                end
            ),
        ),
        pids,
    )
    return Dinfo(name, pids)
end

"""
    selectFCSColumns(selectColnames::Vector{String})
Return a function useful with `loadFCSSet`, which loads only the specified
(prettified) column names from the FCS files. Use `getMetaData`,
`getMarkerNames` and `cleanNames!` to retrieve the usable column names for a
FCS.
"""
function selectFCSColumns(selectColnames::Vector{String})
    ((metadata, data), idx) -> begin
        _, names = getMarkerNames(getMetaData(metadata))
        cleanNames!(names)
        colIdxs = indexin(selectColnames, names)
        if any(colIdxs .== nothing)
            @error "Some columns were not found"
            error("unknown column")
        end
        (metadata, data[:, colIdxs])
    end
end

"""
    distributeFCSFileVector(name::Symbol, fns::Vector{String}, pids=workers())::Dinfo
Distribute a vector of integers among the workers that describes which file
from `fns` the cell comes from. Useful for producing per-file statistics. The
vector is saved on workers specified by `pids` as a distributed variable
`name`.
"""
function distributeFCSFileVector(name::Symbol, fns::Vector{String}, pids = workers())::Dinfo
    sizes = loadFCSSizes(fns)
    slices = slicesof(sizes, length(pids))
    return distributeFileVector(name, sizes, slices, pids)
end

"""
    distributeFileVector(name::Symbol, sizes::Vector{Int}, slices::Vector{Tuple{Int,Int,Int,Int}}, pids=workers())::Dinfo
Generalized version of `distributeFCSFileVector` that produces the integer
vector from any `sizes` and `slices`.
"""
function distributeFileVector(
    name::Symbol,
    sizes::Vector{Int},
    slices::Vector{Tuple{Int,Int,Int,Int}},
    pids = workers(),
)::Dinfo
    dmap(
        slices,
        (slice) ->
            Base.eval(Main, :($name = collectSlice((i) -> fill(i, $sizes[i]), $slice))),
        pids,
    )
    return Dinfo(name, pids)
end

"""
    function getCSVSize(fn::String; args...)::Tuple{Int,Int}
Read the dimensions (number of rows and columns, respectively) from a CSV file
`fn`. `args` are passed to function `CSV.file`.
# Example
    getCSVSize("test.csv", header=false)
"""
function getCSVSize(fn::String; args...)::Tuple{Int,Int}
    n = 0
    k = 0
    # ideally, this will not try to load the whole CSV in the memory
    for row in CSV.File(fn, type = Float64; args...)
        n += 1
        if length(row) > k
            k = length(row)
        end
    end
    return (n, k)
end

"""
    function loadCSVSizes(fns::Vector{String}; args...)::Vector{Int}
Determine number of rows in a list of CSV files (passed as `fns`). Equivalent
to `loadFCSSizes`.
"""
function loadCSVSizes(fns::Vector{String}; args...)::Vector{Int}
    [getCSVSize(fn, type = Float64; args...)[1] for fn in fns]
end

"""
    function loadCSV(fn::String; args...)::Matrix{Float64}
CSV equivalent of `loadFCS`. The metadata (header, column names) are not
extracted. `args` are passed to `CSV.read`.
"""
function loadCSV(fn::String; args...)::Matrix{Float64}
    CSV.read(fn, DataFrame, type = Float64; args...) |> Matrix{Float64}
end

"""
    function loadCSVSet(
        name::Symbol,
        fns::Vector{String},
        pids = workers();
        postLoad = (d, i) -> d,
        csvargs...,
    )::Dinfo
CSV equivalent of `loadFCSSet`. `csvargs` are passed as keyword arguments to
CSV-loading functions.
"""
function loadCSVSet(
    name::Symbol,
    fns::Vector{String},
    pids = workers();
    postLoad = (d, i) -> d,
    csvargs...,
)::Dinfo
    slices = slicesof(loadCSVSizes(fns; csvargs...), length(pids))
    dmap(
        slices,
        (slice) -> Base.eval(
            Main,
            :(
                begin
                    $name = vcollectSlice(
                        (i) -> $postLoad(loadCSV($fns[i]; $csvargs...), i),
                        $slice,
                    )
                    nothing
                end
            ),
        ),
        pids,
    )
    return Dinfo(name, pids)
end

"""
    getMetaData(f)
Collect the meta data information in a more user friendly format.
# Arguments:
- `f`: input structure with `.params` and `.data` fields
"""
function getMetaData(meta::Dict{String,String})::DataFrame

    # declarations and initializations
    metaKeys = keys(meta)
    channel_properties = []
    defaultValue = ""

    # determine the number of channels
    pars = parse(Int, strip(join(meta["\$PAR"])))

    # determine the available channel properties
    for (key,) in meta
        if length(key) >= 2 && key[1:2] == "\$P"
            i = 3
            while i <= length(key) && isdigit(key[i])
                i += 1
            end
            if i <= length(key) && !in(key[i:end], channel_properties)
                push!(channel_properties, key[i:end])
            end
        end
    end

    # create a data frame for the results
    df = Matrix{String}(undef, pars, length(channel_properties))
    df .= defaultValue
    df = DataFrame(df, :auto)
    rename!(df, Symbol.(channel_properties))

    # collect the data from params
    for ch = 1:pars
        for p in channel_properties
            if "\$P$ch$p" in metaKeys
                df[ch, Symbol(p)] = meta["\$P$ch$p"]
            end
        end
    end

    return df
end

"""
    getMarkerNames(meta::DataFrame)::Tuple{Vector{String}, Vector{String}}
Extract suitable raw names (useful for selecting columns) and pretty readable
names (useful for humans) from FCS file metadata.
"""
function getMarkerNames(meta::DataFrame)::Tuple{Vector{String},Vector{String}}
    orig = Array{String}(meta[:, :N])
    nice = copy(orig)
    if hasproperty(meta, :S)
        for i = 1:size(meta, 1)
            if strip(meta[i, :S]) != ""
                nice[i] = meta[i, :S]
            end
        end
    end
    return (orig, nice)
end


"""
    compensate!(data::Matrix{Float64}, spillover::Matrix{Float64}, cols::Vector{Int})
Apply a compensation matrix in `spillover` (the individual columns of which
describe, in order, the spillover of `cols` in `data`) to the matrix `data`
in-place.
"""
function compensate!(data::Matrix{Float64}, spillover::Matrix{Float64}, cols::Vector{Int})
    data[:, cols] = data[:, cols] * inv(spillover)
end

"""
    parseSpillover(str::String)::Union{Tuple{Vector{String},Matrix{Float64}}, Nothing}
Parses the spillover matrix from the string from FCS parameter value.
"""
function parseSpillover(str::String)::Tuple{Vector{String},Matrix{Float64}}
    fields = split(str, ',')
    n = parse(Int, fields[1])
    if length(fields) != 1 + n + n * n
        @error "Spillover matrix of size $n expects $(1+n+n*n) fields, got $(length(fields)) instead."
        error("Invalid spillover matrix")
    end

    names = fields[2:(1+n)]
    spill = collect(transpose(reshape(parse.(Float64, fields[(2+n):length(fields)]), n, n)))

    names, spill
end

"""
    getSpillover(params::Dict{String, String})::Union{Tuple{Vector{String},Matrix{Float64}}, Nothing}
Get a spillover matrix from FCS `params`. Returns a pair with description of
columns to be applied, and with the actual spillover matrix. Returns `nothing`
in case spillover is not present.
"""
function getSpillover(
    params::Dict{String,String},
)::Union{Tuple{Vector{String},Matrix{Float64}},Nothing}
    spillNames = ["\$SPILL", "\$SPILLOVER", "SPILL", "SPILLOVER"]
    for i in spillNames
        if in(i, keys(params))
            return parseSpillover(params[i])
        end
    end
    return nothing
end

"""
    slicesof(lengths::Vector{Int}, slices::Int)::Vector{Tuple{Int,Int,Int,Int}}
Given a list of `lengths` of input arrays, compute a slicing into a specified
amount of equally-sized `slices`.
The output is a vector of 4-tuples where each specifies how to create one
slice. The i-th tuple field contains, in order:
- the index of input array at which the i-th slice begins
- first element of the i-th slice in that input array
- the index of input array with the last element of the i-th slice
- the index of the last element of the i-th slice in that array
"""
function slicesof(lengths::Vector{Int}, slices::Int)::Vector{Tuple{Int,Int,Int,Int}}
    nfiles = length(lengths)
    total = sum(lengths)
    sliceLen = repeat([div(total, slices)], slices)

    for i = 1:mod(total, slices)
        sliceLen[i] += 1
    end

    result = repeat([(0, 0, 0, 0)], slices)

    ifile = 1 # file in the window
    off = 0 # cells already taken from that window
    for i = 1:slices
        startFile = ifile
        startOff = off + 1
        avail = lengths[ifile] - off
        while avail < sliceLen[i]
            ifile += 1
            off = 0
            avail += lengths[ifile]
        end
        rest = avail - sliceLen[i]
        off = lengths[ifile] - rest
        if startOff > lengths[startFile] && startFile < ifile
            startOff = 1
            startFile += 1
        end
        result[i] = (startFile, startOff, ifile, off)
    end
    result
end

"""
    vcollectSlice(loadMtx, (startFile, startOff, finalFile, finalOff)::Tuple{Int,Int,Int,Int})::Matrix
Given a method to obtain matrix content (`loadMtx`), reconstruct a slice from
the information generated by `slicesof`.
This function is specialized for reconstructing matrices and arrays, where the
"element counts" split by `slicesof` are in fact matrix rows. The function is
therefore named _v_collect (the slicing and concatenation is _v_ertical).
The actual data content and loading method is abstracted out -- function
`loadMtx` gets the index of the input part that it is required to fetch (e.g.
index of one FCS file), and is expected to return that input part as a whole
matrix. `vcollectSlice` correctly calls this function as required and extracts
relevant portions of the matrices, so that at the end the whole slice can be
pasted together.
Example:
    # get a list of files
    filenames=["a.fcs", "b.fcs"]
    # get descriptions of 5 equally sized parts of the data
    slices = slicesof(loadFCSSizes(filenames), 5)
    # reconstruct first 3 columns of the first slice
    mySlice = vcollectSlice(
        i -> last(loadFCS(slices[i]))[:,1:3],
        slices[1])
    # (note: function loadFCS returns 4 items, the matrix is the last one)
"""
function vcollectSlice(
    loadMtx,
    (startFile, startOff, finalFile, finalOff)::Tuple{Int,Int,Int,Int},
)::Matrix
    vcat(
        [
            begin
                m = loadMtx(i)
                beginIdx = i == startFile ? startOff : 1
                endIdx = i == finalFile ? finalOff : size(m, 1)
                m[beginIdx:endIdx, :]
            end for i = startFile:finalFile
        ]...,
    )
end

"""
    collectSlice(loadVec, (startFile, startOff, finalFile, finalOff)::Tuple{Int,Int,Int,Int})::Vector
Alternative of `vcollectSlice` for 1D vectors.
"""
function collectSlice(
    loadVec,
    (startFile, startOff, finalFile, finalOff)::Tuple{Int,Int,Int,Int},
)::Vector
    vcat(
        [
            begin
                v = loadVec(i)
                beginIdx = i == startFile ? startOff : 1
                endIdx = i == finalFile ? finalOff : length(v)
                v[beginIdx:endIdx]
            end for i = startFile:finalFile
        ]...,
    )
end
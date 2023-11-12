function load(path::String, sample::Union{EzXML.Node,Nothing}; kwargs...)

    #######################################
    params, data = loadFCS(path; kwargs...)
    params = sample === nothing ? getMetaData(params) : findchannels(sample)

    ###################################### data with human readable channel names
    channels = dropmissing(params, "E")
    transform!(channels, ["N", "S"] => ByRow((N, S) -> ~ismissing(S) ? S != "" ? S : N : N) => "name")
    channelMap = Dict(param.N => ismissing(param.S) ? param.N : param.S for param ∈ eachrow(params))

    # throw error listing duplicate channel names
    duplicates = filter(x -> count(channels.name .== x) > 1, unique(channels.name))
    @assert(isempty(duplicates), "Channel names must be unique. Please resolve duplicates: $duplicates")

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

    params = DataFrame("keyword" => findall("..//Keywords/Keyword[starts-with(@name, '\$P') and substring(@name, string-length(@name)) = 'N' and translate(substring-before(substring(@name, 3), 'N'), '0123456789', '') = '']", sample))
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

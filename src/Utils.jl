function load(path::String, sample::Union{EzXML.Node,Nothing}; transform::Function=x->asinh(x/250), kwargs...)

	#######################################
    params, data = loadFCS(path; kwargs...)
	params = sample === nothing ? getMetaData(params) : findchannels(sample)
	
	###################################### data with human readable channel names
	channels = dropmissing(params,"E")
	transform!(channels, ["N","S"] => ByRow((N,S)-> ~ismissing(S) ? S != "" ? S : N : N) => "name")
	data = DataFrame(data,channels.name)

	###################################### biexponential transformation
	@. data = transform(data)

    ################################ load metadata from workspace
	groups = loadGroups(data,sample)
	gating = gatingGraph(sample; transform=transform,
		channelMap=Dict( param.N => ismissing(param.S) ? param.N : param.S for param ∈ eachrow(params) )
	)

	labels = gate(data,gating)
	return data,labels,groups,gating
end


function load(path::String; files::String=joinpath(dirname(path),"*.fcs"), transform::Function=x->asinh(x/250), cols::Symbol=:setequal, kwargs...)

	@assert( isfile(path), "no such file: $path")
	@assert( length(glob(files)) ≠ 0, "no FCS files found using pattern: $files")

	try
		workspace = root(readxml(path))
	catch error
		throw(LoadError(path,error.line,error))
	end

	workspace = root(readxml(path))
	datasets = map( dataset -> basename(dataset["uri"]), findall("//DataSet",workspace) )
	@assert( length(datasets) == length(unique(datasets)), "FCS files under a workspace must have unique names" )

	data,labels,groups = DataFrame(),DataFrame(),DataFrame()
	gatings = Dict()
	
	for path ∈ glob(files)
		@info "Loading $path"
		sample = findsample(path,workspace)

		fcs,label,group,gating = load(path,sample;transform=transform,kwargs...)
		
		append!(data,fcs,cols=cols)
		append!(labels,label,cols=:union)
		append!(groups,group,cols=:union)
		gatings[path] = gating
	end

	map( name->replace!(labels[!,name],missing=>false), names(labels) )
	map( name->replace!(groups[!,name],missing=>false), names(groups) )

	disallowmissing!(labels)
	disallowmissing!(groups)	

	transform!( labels, AsTable(filter(x->x≠"Unlabelled",names(labels))) => ByRow(x->~any(x)) => "Unlabelled" )
	transform!( groups, AsTable(filter(x->x≠"Ungrouped",names(groups))) =>  ByRow(x->~any(x)) => "Ungrouped" )

	return data,labels,groups,gatings
end

function findsample(uri::String,workspace::EzXML.Node)
	uri = escapeuri(uri)

	for code ∈ [("%2F"=>"/"), ("%5C"=>"\\"), ("%3A"=>":")]
		uri = replace(uri, code)
	end

	sample = findfirst("//DataSet[contains(@uri,'$(basename(uri))')]",workspace)
	sample === nothing && @warn("""no metadata found""")
	return sample
end

function findchannels(sample::EzXML.Node)

	params = DataFrame( "keyword" => findall("..//Keywords/Keyword[starts-with(@name,'\$P') and contains(@name,'N')]", sample))
	keys = ["G","R","V","S","B","N","E","AR"]

	for key ∈ keys
		transform!(params, "keyword" => ByRow( x->getparam(chop(x["name"],head=0,tail=1)*key, sample)) => key)
	end

	transform!(params, "keyword" => ByRow( x->parse(Int,chop(x["name"],head=2,tail=1))) => "index")
	sort!(params,"index")

	transform!(params,["N","S"]=>ByRow((N,S)-> ~ismissing(S) ? S : stripname(N,params.N))=>"S")
	return select!(params,keys)
end

function stripname(name::String,names::AbstractArray{<:String})
    names = filter(x->occursin(x,name)&(x≠name), names)
    return isempty(names) ? missing : first(names)
end

function getparam(param::String,sample::EzXML.Node)
	keyword = findfirst("""..//Keywords/Keyword[@name='$param']""", sample)
	return keyword !== nothing ? keyword["value"] != "" ? keyword["value"] : missing : missing
end
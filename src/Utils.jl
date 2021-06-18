function load(path::String; workspace::String="", transform::Function=x->asinh(x/250), channelMap::Dict=Dict(), kwargs...)

	#######################################
	path = replace(path,"\\"=>"/")
    params, data = loadFCS(path; kwargs...)
	params = getMetaData(params)
	
	###################################### data with human readable channel names
	data = DataFrame(data,params.N)
	function updateMap(laser::String,marker::String)

		if laser ∈ keys(channelMap)
			return channelMap[laser]

		elseif marker ∈ keys(channelMap)
			return channelMap[marker]

		else 
			return marker
		end
	end

	channelNames = Dict([ laser => updateMap(laser,marker)
		for (laser,marker) ∈ zip(params.N,"S" ∈ names(params) ? params.S : params.N) ])
	rename!(data,channelNames) 

	###################################### biexponential transformation
	@. data = transform(data)

    ################################ load metadata from workspace
	groups = loadGroups(path,workspace,data)
	gating = gatingGraph(path,workspace;channelMap=channelNames,transform=transform)
	labels = gate(data,gating)

	return data,labels,groups,gating
end


function load(pattern::GlobMatch; workspace::String="", transform::Function=x->asinh(x/250), channelMap::Dict=Dict(), cols::Symbol=:setequal, kwargs...)

	data,labels,groups = DataFrame(),DataFrame(),DataFrame()
	gatings = Dict()

	paths = readdir(pattern)
	@assert( length(paths) ≠ 0, "no files found using pattern $pattern" )
	
	for path ∈ paths
		@info "Loading $path"
		fcs,label,group,gating = load(path;workspace=workspace,transform=transform,channelMap=channelMap,kwargs...)
		
		append!(data,fcs,cols=cols)
		append!(labels,label,cols=:union)
		append!(groups,group,cols=:union)
		gatings[path] = gating
	end

	map( name->replace!(data[!,name],missing=>0.0), names(data) )
	map( name->replace!(labels[!,name],missing=>false), names(labels) )
	map( name->replace!(groups[!,name],missing=>false), names(groups) )

	disallowmissing!(data)
	disallowmissing!(labels)
	disallowmissing!(groups)	

	transform!( labels, AsTable(filter(x->x≠"Unlabelled",names(labels))) => ByRow(x->~any(x)) => "Unlabelled" )
	transform!( groups, AsTable(filter(x->x≠"Ungrouped",names(groups))) =>  ByRow(x->~any(x)) => "Ungrouped" )

	return data,labels,groups,gatings
end
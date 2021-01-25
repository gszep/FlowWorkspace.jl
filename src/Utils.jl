function load(path::String; workspace::String="", cofactor::Number=250, channelMap::Dict=Dict(), kwargs...)

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
	@. data = asinh(data/cofactor)
	
	if isempty(workspace)
        return data,nothing,nothing,nothing

    else ################################ load metadata from workspace

        groups = loadGroups(path,workspace,data)
		gating = gatingGraph(path,workspace;channelMap=channelNames,cofactor=cofactor)
		labels = gate(data,gating)

		return data,labels,groups,gating
    end
end


function load(pattern::GlobMatch; workspace::String="", cofactor::Number=250, channelMap::Dict=Dict(), kwargs...)

	data,labels,groups = DataFrame(),DataFrame(),DataFrame()
	gatings = Dict()
	
	for path ∈ readdir(pattern)
		fcs,label,group,gating = load(path;workspace=workspace,cofactor=cofactor,channelMap=channelMap,kwargs...)
		
		append!(data,fcs)
		gatings[path] = gating

		append!(labels,label,cols=:union)
		append!(groups,group,cols=:union)		
	end

	if isempty(workspace)
		return data,nothing,nothing,nothing
	end

	map( name->replace!(labels[!,name],missing=>false), names(labels) )
	map( name->replace!(groups[!,name],missing=>false), names(groups) )

	disallowmissing!(labels)
	disallowmissing!(groups)

	return data,labels,groups,gatings
end
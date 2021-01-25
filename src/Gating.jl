function gatingGraph(path::String, workspace::String; channelMap::Dict=Dict(), transform::Function=x->asinh(x/250))
	workspace = root(readxml(workspace))

	############################################## population names
	populations =    findall("//DataSet[contains(@uri,'$(basename(path))')]/..//Population",workspace)
	compensation = findfirst("//DataSet[contains(@uri,'$(basename(path))')]/..//transforms:spilloverMatrix",workspace)
	@assert( length(populations)>0, "gating not found in workspace for sample\n$path")

	graph = MetaDiGraph{Int64,Bool}() ############# store strategy in graph
	set_props!(graph,Dict(:sample=>path))

	for population ∈ populations
		Gate,name = findfirst("Gate",population), population["name"]
		
		id = Gate["gating:id"]
		parent_id = haskey(Gate,"gating:parent_id") ? Gate["gating:parent_id"] : nothing
		
		#################### iterate through disconnected polygons
		for gate ∈ eachelement(Gate)
			
			channels = map( dimension -> replace(dimension["data-type:name"], isnothing(compensation) ? "" : compensation["prefix"]=>"") ∈ keys(channelMap) ?
				channelMap[replace(dimension["data-type:name"], isnothing(compensation) ? "" : compensation["prefix"]=>"")] : throw("""$(dimension["data-type:name"]) not found in channels $(keys(channelMap))"""),
				findall("gating:dimension/data-type:fcs-dimension",gate) )
			
			@assert(gate.name ∈ ["PolygonGate","RectangleGate"], "$(gate.name) not supported for population label $name in sample\n$path")
			if gate.name == "PolygonGate"

				vertices = map( coordinate-> parse(Float32,coordinate["data-type:value"]),
					findall("gating:vertex/gating:coordinate",gate) )
				polygon = map( (x,y)->SVector(transform(x),transform(y)), @view(vertices[1:2:end]), @view(vertices[2:2:end]) )

			elseif gate.name == "RectangleGate"

				minima = map( dimension-> parse(Float32,dimension["gating:min"]),
					findall("gating:dimension",gate) )
				maxima = map( dimension-> parse(Float32,dimension["gating:max"]),
					findall("gating:dimension",gate) )
				
				polygon = map( (x,y)->SVector(transform(x),transform(y)),
					[first(minima),first(minima),first(maxima),first(maxima)],
					[ last(minima), last(maxima), last(maxima), last(minima)] )

			end
			push!(polygon,first(polygon))
			
			MetaGraphs.add_vertex!(graph) ################## store gate as vertex
			set_indexing_prop!( graph, MetaGraphs.nv(graph), :id,id )
			set_props!( graph, MetaGraphs.nv(graph), Dict(:channels=>channels,:polygon=>polygon,:name=>name) )
			
			################# connect parent gates to children
			~isnothing(parent_id) && MetaGraphs.add_edge!(graph,graph[parent_id,:id],MetaGraphs.nv(graph))
		end
    end
    
    return graph
end


function gate(graph::MetaDiGraph,idx::Integer;prefix::String="__gate__")
	return get_prop(graph,idx,:channels) =>
		ByRow( (x,y)->inpolygon( SVector(x,y), get_prop(graph,idx,:polygon);

			in=true, on=false, out=false)
	) => prefix*get_prop(graph,idx,:name)
end


function gate!(data::DataFrame,gating::MetaDiGraph,idx::Integer;prefix::String="__gate__")
	name = prefix*get_prop(gating,idx,:name)
	transform!( data, gate(gating,idx;prefix=prefix) )

	for parent ∈ MetaGraphs.inneighbors(gating,idx) ######### todo @gszep does not support multiple parents
		parent = prefix*get_prop(gating,parent,:name)
		transform!(data, [name,parent] => ByRow((x,y)->x&y) => name)
	end
	
	for child ∈ MetaGraphs.outneighbors(gating,idx)
		gate!(data,gating,child;prefix=prefix)
	end
end


function gate(data::DataFrame,gating::MetaDiGraph;prefix::String="__gate__")

	roots = [ idx for idx ∈ 1:MetaGraphs.nv(gating) if MetaGraphs.indegree(gating,idx)==0 ]
	names = unique([ prefix*get_prop(gating,idx,:name) for idx ∈ 1:MetaGraphs.nv(gating)]) # todo @gszep why wrap in unique()? possible bug?

	for idx ∈ roots
		gate!(data,gating,idx;prefix=prefix)
	end
	
	labels = select(data,names)
	select!(data,Not(names))

	rename!(labels,Dict([ name=>replace(name,prefix=>"") for name ∈ names ]))
	return labels
end
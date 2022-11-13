function gatingGraph(sample::EzXML.Node; channelMap::Dict=Dict(), transforms::Dict=Dict())

    ############# store strategy in graph
    name = findfirst("../SampleNode", sample)["name"]
    graph = MetaDiGraph{Int64,Bool}()
    set_props!(graph, Dict(:sample => name))

    ############################################## population names
    populations = findall("..//Population", sample)
    if length(populations) == 0

        @warn("no populations found")
        return graph
    end

    for population ∈ populations
        Gate, name = findfirst("Gate", population), population["name"]

        id = Gate["gating:id"]
        parent_id = haskey(Gate, "gating:parent_id") ? Gate["gating:parent_id"] : nothing

        #################### iterate through disconnected polygons
        for gate ∈ eachelement(Gate)

            channels = map(dimension -> dimension["data-type:name"] ∈ keys(channelMap) ? channelMap[dimension["data-type:name"]]
                                        : throw("""$(dimension["data-type:name"]) not found in channels $(keys(channelMap))"""),
                findall("gating:dimension/data-type:fcs-dimension", gate))

            transformx, transformy = transforms[first(channels)], transforms[last(channels)]

            @assert(gate.name ∈ ["PolygonGate", "RectangleGate"], "$(gate.name) not supported for population label $name in sample")
            @assert(length(channels) ≤ 2, "length($channels)>2 in $(gate.name) not supported for population label $name in sample")

            if gate.name == "PolygonGate"

                vertices = map(coordinate -> parse(Float32, coordinate["data-type:value"]),
                    findall("gating:vertex/gating:coordinate", gate))
                polygon = map((x, y) -> SVector(transformx(x), transformy(y)), @view(vertices[1:2:end]), @view(vertices[2:2:end]))

            elseif gate.name == "RectangleGate"

                minima = map(dimension -> haskey(dimension, "gating:min") ? parse(Float32, dimension["gating:min"]) : -Inf,
                    findall("gating:dimension", gate))
                maxima = map(dimension -> haskey(dimension, "gating:max") ? parse(Float32, dimension["gating:max"]) : Inf,
                    findall("gating:dimension", gate))

                polygon = map((x, y) -> SVector(transformx(x), transformy(y)),
                    [first(minima), first(minima), first(maxima), first(maxima)],
                    length(channels) == 1 ? [-Inf, Inf, Inf, -Inf] : [last(minima), last(maxima), last(maxima), last(minima)])
            end
            push!(polygon, first(polygon))

            MetaGraphs.add_vertex!(graph) ################## store gate as vertex
            set_indexing_prop!(graph, MetaGraphs.nv(graph), :id, id)
            set_props!(graph, MetaGraphs.nv(graph), Dict(:channels => channels, :polygon => polygon, :name => name))

            ################# connect parent gates to children
            (parent_id !== nothing) && MetaGraphs.add_edge!(graph, graph[parent_id, :id], MetaGraphs.nv(graph))
        end
    end

    return graph
end


function gatingGraph(sample::Nothing; kwargs...)
    return MetaDiGraph{Int64,Bool}()
end


function gate(graph::MetaDiGraph, idx::Integer; prefix::String="__gate__")
    channels, gateName = get_prop(graph, idx, :channels), prefix * get_prop(graph, idx, :name)

    return [first(channels), last(channels)] =>
        ByRow((x, y) -> inpolygon(SVector(x, y), get_prop(graph, idx, :polygon); in=true, on=false, out=false)
        ) => gateName
end


function gate!(data::DataFrame, gating::MetaDiGraph, idx::Integer; prefix::String="__gate__")
    name = prefix * get_prop(gating, idx, :name)
    transform!(data, gate(gating, idx; prefix=prefix))

    for parent ∈ MetaGraphs.inneighbors(gating, idx) ######### todo @gszep does not support multiple parents
        parent = prefix * get_prop(gating, parent, :name)
        transform!(data, [name, parent] => ByRow((x, y) -> x & y) => name)
    end

    for child ∈ MetaGraphs.outneighbors(gating, idx)
        gate!(data, gating, child; prefix=prefix)
    end
end


function gate(data::DataFrame, gating::MetaDiGraph; prefix::String="__gate__")

    roots = [idx for idx ∈ 1:MetaGraphs.nv(gating) if MetaGraphs.indegree(gating, idx) == 0]
    names = unique([prefix * get_prop(gating, idx, :name) for idx ∈ 1:MetaGraphs.nv(gating)]) # todo @gszep why wrap in unique()? possible bug?

    if length(roots) == 0
        return DataFrame("Unlabelled" => fill(true, size(data, 1)))
    end

    for idx ∈ roots
        gate!(data, gating, idx; prefix=prefix)
    end

    labels = select(data, names)
    select!(data, Not(names))

    rename!(labels, Dict([name => replace(name, prefix => "") for name ∈ names]))
    return labels
end
module FlowWorkspace

    using GigaSOM,EzXML,DataFrames
    using MetaGraphs

    using PolygonOps,StaticArrays
    using Glob: GlobMatch
    using Observables

    include("Groups.jl")
    include("Gating.jl")
    include("Utils.jl")
end

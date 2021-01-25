module FlowWorkspace

    using GigaSOM,DataFrames
    using EzXML,MetaGraphs

    using PolygonOps, StaticArrays
    using Glob: GlobMatch, @glob_str

    include("Groups.jl")
    include("Gating.jl")
    include("Utils.jl")

    export load, @glob_str
end

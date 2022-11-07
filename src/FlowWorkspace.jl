module FlowWorkspace

using GigaSOM, DataFrames
using Interpolations, EzXML, MetaGraphs

using PolygonOps, StaticArrays
using Glob: glob
using HTTP: escapeuri

include("Groups.jl")
include("Transforms.jl")
include("Gating.jl")
include("Utils.jl")

export load
end

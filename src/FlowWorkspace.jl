module FlowWorkspace

using DataFrames, FCSFiles, FileIO
using Interpolations, EzXML, MetaGraphs

using PolygonOps, StaticArrays
using Glob: glob
using URIs: escapeuri

include("Groups.jl")
include("Transforms.jl")
include("Gating.jl")
include("Utils.jl")

export load
end

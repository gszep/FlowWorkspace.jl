using FlowWorkspace
using Test

######################################################### load single file
workspace = "data/workspace.wsp"
sample = "data/fcs/101_DEN084Y5_15_E01_008_clean.fcs"

data, = load(sample)
data,labels,groups,gating = load(sample; workspace=workspace)

######################################################### multiple files
pattern = glob"data/fcs/*.fcs"
data,labels,groups,gating = load(pattern; workspace=workspace)
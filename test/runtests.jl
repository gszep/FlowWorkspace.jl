using FlowWorkspace
using Test

######################################################### load single file
workspace = "./data/den.wsp"
sample = "./data/den/101_DEN084Y5_15_E01_008_clean.fcs"

data,labels,groups,gating = load(sample)
data,labels,groups,gating = load(sample; workspace=workspace)

######################################################### multiple files
pattern = glob"./data/den/*.fcs"
data,labels,groups,gating = load(pattern)
data,labels,groups,gating = load(pattern; workspace=workspace)

######################################################### quadgates
data,labels,groups,gating = load(glob"./data/compensation/*fcs",workspace="./data/compensation.wsp")
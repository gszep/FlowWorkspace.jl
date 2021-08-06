using FlowWorkspace
using Test

######################################################### loading
load("./data/den.wsp"; files="./data/den/*.fcs")

######################################################### quadgates
load("./data/compensation.wsp"; files="./data/compensation/*fcs")
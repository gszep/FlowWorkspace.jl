# FlowWorkspace.jl
This package provides for loading and parsing of [FlowJo Workspace files](https://docs.flowjo.com/flowjo/workspaces-and-samples/ws-savinganalysis) in pure Julia. The gating strategy is parsed as a `DiMetaGraph` and group membership as `Dict`. Together with [GigaSOM.jl](https://github.com/LCSB-BioCore/GigaSOM.jl) this metadata can be used to generate event-level labels for `FCS` datasets.

[![Build Status](https://travis-ci.com/gszep/FlowWorkspace.jl.svg?branch=master)](https://travis-ci.com/gszep/FlowWorkspace.jl)
[![Coverage](https://codecov.io/gh/gszep/FlowWorkspace.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gszep/FlowWorkspace.jl)

## Installation

Type `] add FlowWorkspace` and then hit ⏎ Return at the REPL. You should see `pkg> add FlowWorkspace`.

## Basic Usage
> :warning: **FCS files under a workspace must have unique names**. This limitation will be removed in future versions

The `load` method returns a tuple of three `DataFrames` and one `Dict`:
* `data` arcsinh transformed fluorescence values for each event
* `labels` boolean telling us which events have been gated
* `groups` booleans telling us which group each event belongs to
* `gating` Dictionary of graph representations of each gating strategy

```julia
using FlowWorkspace

######################################################### load single file
workspace = "workspace.wsp"
sample = "workspace/data.fcs"
data,labels,groups,gating = load(sample; workspace=workspace)

######################################################### load multiple files with different channel names
channelMap = Dict([

    "FJComp-355 379_28-A"=>"CD3", 
    "FJComp-355 560_40-A"=>"CD8", 

    "FJComp-355 820_60-A"=>"CD4",
    "FJComp-355 670_30-A"=>"CD4",

    "FJComp-640 780_60-A"=>"CCR7",
    "FJComp-405 780_60-A"=>"CD45RA", 

    "FJComp-561 780_60-A"=>"CD127", 
    "FJComp-640 670_30-A"=>"CD25", 

    "FJComp-561 610_20-A"=>"Helios", 
    "FJComp-561 585_15-A"=>"Foxp3", 
    "Foxp3-IgM"=>"Foxp3",

    "FJComp-405 710_40-A"=>"PD-1", 
    "FJComp-640 730_35-A"=>"CXCR5", 

    "FJComp-405 670_30-A"=>"CCR6", 
    "FJComp-488 715_30-A"=>"CXCR3", 

    "FJComp-405 605_40-A"=>"CCR4", 
    "FJComp-488 525_50-A"=>"CCR10", 

    "FJComp-405 450_50-A"=>"CD103", 
    "FJComp-355 740_35-A"=>"CD69",
    "FJComp-405 515_20-A"=>"HLA-DR"
])

pattern = glob"workspace/*.fcs"
data,labels,groups,gating = load(pattern; workspace=workspace, channelMap=channelMap)
```

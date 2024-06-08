# FlowWorkspace.jl
This package provides for loading and parsing of [FlowJo Workspace files](https://docs.flowjo.com/flowjo/workspaces-and-samples/ws-savinganalysis) in pure Julia. The gating strategy is parsed as a `DiMetaGraph` and group membership as `Dict`. Together with [GigaSOM.jl](https://github.com/LCSB-BioCore/GigaSOM.jl) this metadata can be used to generate event-level labels for `FCS` datasets.

[![Build Status](https://app.travis-ci.com/gszep/FlowWorkspace.jl.svg?token=UdhYxVPspzofqoNnanY8&branch=master)](https://travis-ci.com/gszep/FlowWorkspace.jl)
[![Coverage](https://codecov.io/gh/gszep/FlowWorkspace.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gszep/FlowWorkspace.jl)

## Installation

Type `] add FlowWorkspace` and then hit ⏎ Return at the REPL. You should see `pkg> add FlowWorkspace`.

## Basic Usage
> :warning: **FCS files under a workspace must have unique names**. This limitation will be removed in future versions

The `load` method returns a tuple of three `DataFrames` and one `Dict`:
* `data` transformed fluorescence values for each event
* `labels` boolean telling us which events have been gated
* `groups` booleans telling us which group each event belongs to
* `gating` Dictionary of graph representations of each gating strategy

```julia
using FlowWorkspace
data,labels,groups,gating = load("workspace.wsp"; files="*.fcs", cols=:setequal)
```

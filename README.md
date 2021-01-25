# FlowWorkspace.jl
This package provides for loading and parsing of [FlowJo Workspace files](https://docs.flowjo.com/flowjo/workspaces-and-samples/ws-savinganalysis) in pure Julia. The gating strategy is parsed as a `DiMetaGraph` and group membership as `Dict`. Together with [GigaSOM.jl](https://github.com/LCSB-BioCore/GigaSOM.jl) this metadata can be used to generate event-level labels for `FCS` datasets.

[![Build Status](https://travis-ci.com/gszep/FlowWorkspace.jl.svg?branch=master)](https://travis-ci.com/gszep/FlowWorkspace.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/gszep/FlowWorkspace.jl?svg=true)](https://ci.appveyor.com/project/gszep/FlowWorkspace-jl)
[![Coverage](https://codecov.io/gh/gszep/FlowWorkspace.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gszep/FlowWorkspace.jl)
[![Coverage](https://coveralls.io/repos/github/gszep/FlowWorkspace.jl/badge.svg?branch=master)](https://coveralls.io/github/gszep/FlowWorkspace.jl?branch=master)

## Installation

Type `] add FlowWorkspace` and then hit ⏎ Return at the REPL. You should see `pkg> add FlowWorkspace`.

## Basic Usage

```julia
using FlowWorkspace

```

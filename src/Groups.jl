function loadGroups(path::String, workspace::String, data::DataFrame)

    workspace = root(readxml(workspace))
    sample = findfirst("//DataSet[contains(@uri,'$(basename(path))')]",workspace)
    if sample !== nothing

        return DataFrame([ (group["name"]=>fill(true,size(data,1))) for group ∈ 
            findall("""//SampleRefs/SampleRef[@sampleID='$(sample["sampleID"])']/../..""",workspace)
                if  group["name"] ≠ "All Samples" ])
    else 
        return nothing
    end
end
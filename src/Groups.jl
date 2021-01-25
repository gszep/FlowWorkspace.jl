function loadGroups(path::String, workspace::String, data::DataFrame)

    workspace = root(readxml(workspace))
    sample = findfirst("//DataSet[contains(@uri,'$(basename(path))')]",workspace)
    if ~isnothing(sample)

        return DataFrame([ (group["name"]=>fill(true,size(data,1))) for group ∈ 
            findall("""//SampleRefs/SampleRef[contains(@sampleID,'$(sample["sampleID"])')]/../..""",workspace)
                if  group["name"] ≠ "All Samples" ])
    else 
        return nothing
    end
end
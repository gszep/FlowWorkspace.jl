function loadGroups(path::String, workspace::String, data::DataFrame)
    path = basename(path).replace("%20"," ")

    workspace = root(readxml(workspace))
    sample = findfirst("//DataSet[contains(@uri,'$path')]",workspace)
    if sample !== nothing

        return DataFrame([ (group["name"]=>fill(true,size(data,1))) for group ∈ 
            findall("""//SampleRefs/SampleRef[@sampleID='$(sample["sampleID"])']/../..""",workspace)
                if  group["name"] ≠ "All Samples" ])
    else 
        return nothing
    end
end
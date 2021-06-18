function loadGroups(path::String, workspace::String, data::DataFrame)
    path = replace(basename(path),"%20"=>" ")

    try 
        workspace = root(readxml(workspace))
    catch
        @warn("""Workspace not loaded $workspace""")
        return DataFrame("Ungrouped"=>fill(true,size(data,1)))
    end

    sample = findfirst("//DataSet[contains(@uri,'$path')]",workspace)
    if sample !== nothing

        return DataFrame([ (group["name"]=>fill(true,size(data,1))) for group ∈ 
            findall("""//SampleRefs/SampleRef[@sampleID='$(sample["sampleID"])']/../..""",workspace)
                if  group["name"] ≠ "All Samples" ])
    else 
        @warn("""No metadata found for sample $path in workspace $workspace""")
        return DataFrame("Ungrouped"=>fill(true,size(data,1)))
    end
end
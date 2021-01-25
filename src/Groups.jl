function loadGroups(path::String, workspace::EzXML.Node, data::DataFrame)
    sampleID = findfirst("//DataSet[contains(@uri,'$path')]",workspace)["sampleID"]

    return DataFrame([ (group["name"]=>fill(true,size(data,1))) for group ∈ 
        findall("//SampleRefs/SampleRef[contains(@sampleID,'$sampleID')]/../..",workspace)
            if  group["name"] ≠ "All Samples" ])
end
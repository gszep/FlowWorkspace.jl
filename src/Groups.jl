function loadGroups(data::DataFrame, sample::EzXML.Node)
    return DataFrame([ (group["name"]=>fill(true,size(data,1))) for group ∈ 
        findall("""//SampleRefs/SampleRef[@sampleID='$(sample["sampleID"])']/../..""",sample)
            if  group["name"] ≠ "All Samples" ])
end

function loadGroups(data::DataFrame, sample::Nothing)
    return DataFrame("Ungrouped"=>fill(true,size(data,1)))
end
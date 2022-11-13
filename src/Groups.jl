function loadGroups(data::DataFrame, sample::EzXML.Node)
    groups = filter(group -> group["name"] ≠ "All Samples", findall("""//SampleRefs/SampleRef[@sampleID='$(sample["sampleID"])']/../..""", sample))
    return length(groups) ≠ 0 ? DataFrame([(group["name"] => fill(true, size(data, 1))) for group ∈ groups]) : loadGroups(data, nothing)
end

function loadGroups(data::DataFrame, sample::Nothing)
    return DataFrame("Ungrouped" => fill(true, size(data, 1)))
end
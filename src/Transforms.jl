function transforms(sample::EzXML.Node; channelMap::Dict=Dict())
    return Dict(map((eachelement ∘ findfirst)("../Transformations", sample)) do transform

        parameter = findfirst("data-type:parameter", transform)
        channelName = channelMap[parameter["data-type:name"]]

        transform_name = Symbol(transform.name)
        params = Dict(Symbol(a.name) => parse(Float64, nodecontent(a)) for a in attributes(transform))
        try
            return channelName => eval(:($transform_name(; $params...)))
        catch error
            throw("$transform_name transform not implemented")
        end
    end)
end

"""
FlowJo compatible linear transformation function
"""
function linear(; params...)
    return LinearInterpolation(linear_table(; params...)...; extrapolation_bc=Flat())
end

function linear_table(; minRange::Real=0, maxRange::Real=262144, gain::Real=1)
    x = range(minRange, maxRange; step=1)
    return x, gain * x
end

"""
FlowJo compatible biex transformation function
Implementation adapted from the FlowKit Python library
"""
function biex(; params...)
    return LinearInterpolation(biex_table(; params...)...; extrapolation_bc=Flat())
end

function biex_table(; width::Real=-10, neg::Real=0, pos::Real=4.418540, maxRange::Real=262144, channelRange::Real=4096, length::Real=256)
    @assert(length == 256, "length $length not supported for transforms:biex")

    width = log10(-width)
    decades = pos - (width / 2)
    extra = max(neg, 0) + (width / 2)

    zero_point = round(Int, (extra * channelRange) / (extra + decades))
    zero_point = round(Int, min(zero_point, channelRange / 2))

    if zero_point > 0
        decades = extra * channelRange / zero_point
    end

    width = width / (2 * decades)
    positive_range = Base.log(10) * decades
    negative_range = _log_root(positive_range, width)

    n_points = channelRange + 1
    values = range(0, channelRange; step=1)

    positive = exp.(values / float(n_points) * positive_range)
    negative = exp.(-values / float(n_points) * negative_range)

    s = exp((positive_range + negative_range) * (width + extra / decades))
    negative *= s
    s = positive[zero_point+1] - negative[zero_point+1]

    positive[(zero_point+1):n_points] = positive[(zero_point+1):n_points] - negative[(zero_point+1):n_points]
    positive[(zero_point+1):n_points] = maxRange / exp(positive_range) * (positive[(zero_point+1):n_points] .- s)

    neg_range = range(0, zero_point - 1; step=1)
    m = 2 * zero_point .- neg_range

    positive[neg_range.+1] = -positive[m.+1]
    return positive, values
end

function _log_root(b, w)
    x_lo = 0
    x_hi = b
    d = (x_lo + x_hi) / 2
    dx = abs(round(Int, x_lo - x_hi))
    dx_last = dx
    fb = -2 * Base.log(b) + w * b
    f = 2.0 * Base.log(d) + w * b + fb
    df = 2 / d + w

    if w == 0
        return b
    end

    for i in 0:99
        if ((((d - x_hi) * df - f) - ((d - x_lo) * df - f)) > 0) || (abs(2 * f) > abs(dx_last * df))
            dx = (x_hi - x_lo) / 2
            d = x_lo + dx
            if d == x_lo
                return d
            end
        else
            dx = f / df
            t = d
            d -= dx
            if d == t
                return d
            end
        end

        if abs(dx) < 1.0E-12
            return d
        end

        dx_last = dx
        f = 2 * Base.log(d) + w * d + fb
        df = 2 / d + w
        if f < 0
            x_lo = d
        else
            x_hi = d
        end
    end

    return d
end

"""FlowJo compatible log transformation function"""
function log(; params...)
    return LinearInterpolation(log_table(; params...)...; extrapolation_bc=Flat())
end

function log_table(; offset::Real=1, decades::Real=1, minRange::Real=0, maxRange::Real=262144)
    x = range(minRange, maxRange; step=1)
    @. y = log10((x + offset) / decades)
    return x, y
end
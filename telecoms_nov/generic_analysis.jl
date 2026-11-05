ArrayofNumbers = AbstractArray{<:Number}
ArrayofReals = AbstractArray{<:Real}

SPEED_OF_LIGHT = 299792458.0

#########
# Convenience arithmetc
#########

function dBm2mW(x)
    return 10^(x/10)
end

function mW2dBm(x)
    return 10*log10(x)
end

function dBm_sum(x...)
    P_mW = sum(Telecoms.dBm2mW.(x))
    return Telecoms.mW2dBm(P_mW)
end

function dBm_sub(x...)
    neg_mW = sum(Telecoms.dBm2mW.(x[begin+1:end]))
    pos_mW = Telecoms.dBm2mW(x[begin])
    if pos_mW < neg_mW
        # can't have a negative power so clamp it to -Inf instead of error
        return -Inf
    end
    return Telecoms.mW2dBm(pos_mW - neg_mW)
end

function ltof(l)
    return SPEED_OF_LIGHT/l
end

function ftol(f)
    return ltof(f)
end

function bw_m_Hz_approx(dl, l_0)
    return SPEED_OF_LIGHT * dl / l_0^2
end

function bw_Hz_m_approx(df, f_0)
    return bw_m_Hz_approx(df, f_0)
end

function parallel(x...)
    return sum(1 ./ x)
end

function parallel_sub(x...)
    pos_x = 1 / x[begin]
    neg_x = sum(1 ./ x[begin+1:end])

    return pos_x - neg_x
end

# these two seem to run faster than the approximate ones
function bw_m_Hz(l_high, l_low)
    return SPEED_OF_LIGHT * parallel_sub(l_low, l_high)
end

function bw_Hz_m(f_high, f_low)
    return bw_m_Hz(f_high, f_low)
end

# numerical calculus

function differentiate(x)
    dx = zeros(length(x))
    for idx in eachindex(x)
        if idx == 1
            continue
        end
        dx[idx] = x[idx] - x[idx - 1]
    end

    return dx
end

function differentiate(x, y)
    dy = differentiate(y)
    dx = x[begin+1] - x[begin]

    return dy ./ dx
end

function integral_from_module(t, x)

    problem = SciMLBase.SampledIntegralProblem(abs.(x).^2,t)
    ret = CommonSolve.solve(problem, Integrals.TrapezoidalRule())

    if ret.retcode == ReturnCode.Success
        return ret
    else
        return -1
    end
end

function integral_rect(t::T, 
    x::T) where (T<:ArrayofReals)

    # Riemann method

    if !_equal_lengths(t, x)
        return -1
    end

    dt = t[begin+1:end] .- t[begin:end-1]
    ret = sum(x[begin+1:end] .* dt)

    return ret
end
    
function integral_trap(t, x)

    if !_equal_lengths(t, x)
        return -1
    end

    dt = t[begin+1:end] .- t[begin:end-1]
    dx = x[begin+1:end] .+ x[begin:end-1]
    ret = sum(0.5 .* dx .* dt)

    return ret
end


##########
# Finding stuff
##########

function peak_search_2(x::T, 
    y::T; 
    n=1) where (T<:ArrayofReals)

    # Uses both derivatives to find points where a curve is at peak or
    # flat.

    if !_equal_lengths(x, y)
        return -1
    end

    dy = Telecoms.differentiate(x, y)
    nullos = Telecoms.find(dy, 0.0)
    
    ddy = Telecoms.differentiate(x, dy)
    n_nullos = length(nullos)
    peaks = Vector{Tuple{Float64, Float64}}()
    for i in eachindex(nullos)
        nullo_idx = nullos[i][2]
        nullo_ddy = ddy[nullo_idx]
        if nullo_ddy > 0
            n_nullos -= 1
        else
            push!(peaks, (x[nullo_idx], y[nullo_idx]))
        end
    end

    peaks_s = reverse(sort(peaks, by=x->x[2]))
    if n > length(peaks_s)
        n = length(peaks_s)
    end

    return peaks_s[begin:n]
end

function filter_flat_peaks(peaks::Vector{Tuple{Float64, Float64}}, 
    tol::Real)

    # If peak is flat, filter the rest according to a tolerance
    # This selects the leftermost side of the peak and should work
    # for any peak width as only one pass is performed.
    # ...P---.T.P---.T. etc.

    peak_idx = Vector{Int64}()
    allowed_w = Vector{Float64}()
    for p in eachindex(peaks)
        p_w = peaks[p][1]
        #if !any(x->(p_w<x+tol)&&(p_w>x-tol), allowed_w)
        if !any(x->p_w<x+tol, allowed_w)
            push!(allowed_w, p_w)
            push!(peak_idx, p)
        end
    end

    actual_peaks = typeof(peaks)()
    for idx in eachindex(peaks)
        if idx in peak_idx
            push!(actual_peaks, peaks[idx]) 
        end
    end

    return actual_peaks
end

function peak_search(x::T) where (T<:ArrayofReals)

    return n_peak_search(x, 1)[1] 
end

function peak_search(x::T,
    y::T) where (T<:ArrayofReals)

    if !_equal_lengths(x, y)
        return -1
    end

    validx = n_peak_search(x, y, 1)[1]

    return (validx[1], validx[2])
end

function n_peak_search(x::T,
    n::Integer;
    backoff::Integer = 0) where (T<:ArrayofReals)

    # Rubbish compared to peak_search_2, especially for tight spacing

    # return [(val, idx), ..., (valn, idxn)]
    # backoff is the number of indices in which one peak may be recorded 

    if n == 1
        return [findmax(x)]
    end

    # algorithm for n peak search
    points = []
    push!(points, findmax(x))
    forbidden = [points[1][2]-backoff:points[1][2]+backoff]
    for i in 1:n
        prev_max = points[end][1]
        new_max = -Inf
        idx::Int64 = 1
        max_idx::Int64 = 0

        for val in x
            # backoff region skip
            if any(x -> idx in x, forbidden)
                idx += 1
                continue 
            end

            if val > new_max && val < prev_max
                new_max = val
                max_idx = idx
                push!(forbidden, max_idx-backoff:max_idx+backoff)
            end

            idx += 1
        end

        if max_idx < 1
            # probably no more maxima so just put the maximum
            push!(points, findmax(x))
        else
            push!(points, (new_max, max_idx))
        end
    end

    return points
end

function n_peak_search(x::T,
    y::T,
    n::Integer;
    backoff::AbstractFloat = 0.0) where (T<:ArrayofReals)
    # return [(x0, y0), ..., (xn, yn)]

    if !_equal_lengths(x, y)
        return -1
    end

    # calculate backoff indices from x 
    backoff_idx = round(Int64, 
        length(x) * backoff / (maximum(x) - minimum(x)))
#    println("n_peak_search: backoff samples = " * string(backoff_idx))

    peak_validx = n_peak_search(y, n, backoff=backoff_idx)
    points = []
    for validx in peak_validx
        yn = validx[1]
        xn = x[validx[2]]

        push!(points, (xn, yn))
    end

    return points
end

function find(x::T,
    x_0::Real) where (T<:ArrayofReals)
    
    # return [(val, idx), ..., (valn, idxn)]
    # instead of using a tolerance like with findall, use inequalities
    # when error goes through 0, this is about where the point is
    points = []
    idx::Int64 = 2
    last_val = x[1]
    for val in x[2:end]
        error = val - x_0
        last_error = last_val - x_0

        if error < 0 && last_error > 0 || error > 0 && last_error < 0
            if last_error > error
                push!(points, (val, idx))
            else
                push!(points, (last_val, idx - 1))
            end
        elseif error == 0
            push!(points, (val, idx))
        end

        last_val = val
        idx += 1
    end

    return points
end

function find(x::T,
    y::T,
    y_0::Real) where (T<:ArrayofReals)

    if !_equal_lengths(x, y)
        return -1
    end

    find_validx = find(y, y_0)

    points = []
    for validx in find_validx
        yn = validx[1]
        xn = x[validx[2]]

        push!(points, (xn, yn))
    end

    return points
end

##########
# Smart slices
##########

function _max_point_limits(xdata::T,
    ydata::T, 
    xrange::AbstractFloat) where (T<:ArrayofReals)
    return [(xdata[argmax(ydata)] - xrange), 
        (xdata[argmax(ydata)] + xrange)];
end

function _min_point_limits(xdata::T,
    ydata::T, 
    xrange::AbstractFloat) where (T<:ArrayofReals)
    return [(xdata[argmin(ydata)] - xrange), 
        (xdata[argmin(ydata)] + xrange)];
end

function _data_around_limits(xdata::T, 
    ydata::T, 
    xlimits::T) where (T<:ArrayofReals)

    if xlimits == -1
        print("_data_around_limits: bad limits")
        return -1
    end

    xidx_l = find(xdata, xlimits[1]);
    xidx_h = find(xdata, xlimits[2]);

    if length(xidx_l) < 1 || length(xidx_l) < 1 
        print("_data_around_limits: no data in limits")
        return -1
    end

    xlimit = xidx_l[1][2]:xidx_h[1][2];

    slice_x = xdata[xlimit];
    slice_y = ydata[xlimit];

    return [slice_x, slice_y]
end

function data_around_peak(xdata::T, 
    ydata::T, 
    xlimits::AbstractFloat) where (T<:ArrayofReals)
    # return slice of x and y data around the y maximum within +/- xrange
    # as vector of x and y vectors

    if !_equal_lengths(xdata, ydata)
        return -1
    end

    xlimit = _max_point_limits(xdata, ydata, xlimits)

    return _data_around_limits(xdata, ydata, xlimit)
end

function data_around_trough(xdata::T, 
    ydata::T, 
    xlimits::AbstractFloat) where (T<:ArrayofReals)

    if !_equal_lengths(xdata, ydata)
        return -1
    end

    xlimit = _min_point_limits(xdata, ydata, xlimits)

    print(xlimit)

    return _data_around_limits(xdata, ydata, xlimit)
end

#########
# Checks
#########

function any_zero_nan(x)
    if any(f -> f <= 0 || f == NaN, x)
        println("Negative or undefined values in set.")
        return true
    end

    return false
end

function _equal_lengths(a::T,
    b::T) where (T<:AbstractArray)
    
    return length(a) == length(b)
end

##########
# DSP
##########

function undersample(x, 
    factor::Integer)
    if factor <= 1
        println("undersample: factor should be at least 1")
        return x
    end

    return x[1:factor:length(x)]
end

function undersample(x, 
    factor::AbstractFloat;
    tolerance=0.05)
    if factor <= 1
        println("undersample: factor should be at least 1")
        return x
    end

    println("undersample: need to make this algorithm work")

    # express factor as ratio of reasonable integers
    # factor : k 
    k = 1

    while abs(round(factor * k) - factor * k) > tolerance
        k += 1
    end

    factor = round(Int, factor * k)

    println("undersample: undersampling with ratio ", k, ":", factor)
     
    new_length = round(Int, length(x) * (k / factor))
    x_u = Vector{Float64}(undef, new_length)

    reject_c = 0
    idx_u = 1
    for v in x
        # it is still trying to place the new samples at equal spaces
        # this is inappropriate for fractional undersampling
        if reject_c >= k 
            if idx_u > length(x_u)
                break
            end
            x_u[idx_u] = v
            reject_c = 0
            idx_u += 1
        end
        reject_c += 1
    end
        
    return x_u
end

# the daq2 util.jl thing uses the same concept, but in 2D
function decimation_avg(t, x, factor::Integer)
    if factor < 1
        println("decimation_avg: factor must be over 0")
        return
    end
    if !Telecoms._equal_lengths(t, x)
        println("unequal argument lengths: ", length(t), " vs ", length(x))
        return
    end
    avg_length = div(length(x), factor)  # asuming same length
    avg_buf = Vector{typeof(x[1])}(undef, avg_length)
    idx = 1
    for i in 1:factor:length(x)-factor
        view = x[i:i+factor-1]
        avg_buf[idx] = sum(view) / factor
        idx += 1
    end
    t = t[1:avg_length] .* factor
    if length(t) != length(avg_buf)
        println("decimation_avg: dimensional arithmetic error")
        println(length(t), " vs ", length(avg_buf))
        return
    end
    return t, avg_buf
end

#####
# Analysis of spectra
#####

function OSNR_dB(w::T,
    S::T,  # logarithmic
    w_0::Real;  # signal location
    bw=0.1,  # nm - if w is in Hz, this automatically goes to Hz
    pessimistic=false,
    wavelength=true) where (T<:ArrayofReals)

    if !wavelength
        bw = Telecoms.ltof(bw)
    end

    noise_p = w_0 + bw
    noise_m = w_0 - bw
    point_p = Telecoms.find(w, noise_p)
    point_m = Telecoms.find(w, noise_m)
    if length(point_p) < 1 || length(point_m) < 1
        println("OSNR_dB: couldn't get a noise level")
        return NaN  # because this function has integer return value
    end
    idx_p = point_p[1][2]
    idx_m = point_m[1][2]
    S_p = S[idx_p]
    S_m = S[idx_m]

    point_0 = Telecoms.find(w, w_0)
    if length(point_0) < 1
        println("OSNR_dB: couldn't get signal level for some reason")
        return NaN
    end
    idx_0 = point_0[1][2]
    S_0 = S[idx_0]

    S_n = 0
    if pessimistic
        S_n = reverse(sort([S_p, S_m]))[1]
    else
        S_n = (S_p + S_m) / 2
    end
    
    OSNR = S_0 - S_n

    return OSNR
end


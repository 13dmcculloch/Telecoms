using FFTW
using DSP
using LsqFit
using Statistics
using Integrals
using CommonSolve

include("generic_analysis.jl")

ArrayofNumbers = AbstractArray{<:Number}
ArrayofReals = AbstractArray{<:Real}

###########
# FFT
###########

function fft(x::T,
    sample_rate_Hz::Real) where
    (T<:ArrayofNumbers)

    X = FFTW.fft(x)
    w = FFTW.fftfreq(length(x), sample_rate_Hz)

    return (collect(w), X)
end

function fft(i::T,
    q::T,
    sample_rate_Hz::Real) where (T<:ArrayofNumbers)
    # IQ version

    return fft(Complex.(i, q), sample_rate_Hz) 
end

function fft_norm(x::T,  
    sample_rate_Hz::Real) where
    (T<:ArrayofNumbers)
    # normalise FFT output to length of FFT
    # long FFTs conserve energy of spike hence have massive y values
    
    w, X = fft(x, sample_rate_Hz)
    inv_length = 1/length(X)

    return w, inv_length .* X
end

function fft_norm(i::T,
    q::T,
    sample_rate_Hz::Real) where (T<:ArrayofNumbers)

    return fft_norm(Complex.(i, q), sample_rate_Hz)
end

function fft_magnitude(x::T,  
    sample_rate_Hz::Real) where
    (T<:ArrayofReals)

    w, X = fft(x, sample_rate_Hz)
    
    return (w, abs.(X).^2)
end

function fft_magnitude_norm(x::T,  
    sample_rate_Hz::Real) where
    (T<:ArrayofNumbers)

    w, X = fft_norm(x, sample_rate_Hz)
    
    return (w, abs.(X).^2)
end

function fft_magnitude(i::T,
    q::T,
    sample_rate_Hz::Real) where (T<:ArrayofNumbers)

    return fft_magnitude(Complex.(i, q), sample_rate_Hz)
end


function fft_magnitude_norm_dB(x::T,  
    sample_rate_Hz::Real) where
    (T<:ArrayofNumbers)

    w, X = fft_magnitude_norm(x, sample_rate_Hz)
    
    return w, 20 .*log10.(X)
end

function fft_magnitude_dB(i::T,
    q::T,
    sample_rate_Hz::Real) where (T<:ArrayofNumbers)

    w, X = fft_magnitude(i, q, sample_rate_Hz)
    
    return w, 20 .*log10.(X)
end

function fft_magnitude_norm_dB(x::T,  
    sample_rate_Hz::Real) where
    (T<:ArrayofNumbers)

    w, X = fft_magnitude(x, sample_rate_Hz)
    
    return w, 20 .*log10.(X)
end

function fft_phase(x::T,  
    sample_rate_Hz::Real) where
    (T<:ArrayofNumbers)

    w, X = fft(x, sample_rate_Hz)
    
    return (w, angle.(X))
end

function fft_phase(i::T,
    q::T,
    sample_rate_Hz::Real) where (T<:ArrayofNumbers)

    return fft_phase(Complex.(i, q), sample_rate_Hz)
end

###########
# RFFT
###########

function rfft(x::T,  
    sample_rate_Hz::Real) where
    (T<:ArrayofReals)
    # Basic real-valued FFT (R in, C out)
    # Baseband or SSB 

    X = FFTW.rfft(x)
    w = FFTW.rfftfreq(length(x), sample_rate_Hz)

    return (collect(w), X)
end

function rfft_norm(x::T,  
    sample_rate_Hz::Real) where
    (T<:ArrayofReals)

    w, X = rfft(x, sample_rate_Hz)
    inv_length = 1/length(X)

    return w, inv_length .* X
end

function rfft_magnitude(x::T,  
    sample_rate_Hz::Real) where
    (T<:ArrayofReals)

    w, X = rfft(x, sample_rate_Hz)
    
    return (w, abs.(X).^2)
end

function rfft_magnitude_dB(x::T,  
    sample_rate_Hz::Real) where
    (T<:ArrayofReals)

    w, X = rfft_magnitude(x, sample_rate_Hz)
    
    return w, 20 .*log10.(X)
end

function rfft_magnitude_norm(x::T,  
    sample_rate_Hz::Real) where
    (T<:ArrayofReals)

    w, X = rfft_norm(x, sample_rate_Hz)
    
    return (w, abs.(X).^2)
end

function rfft_magnitude_norm_dB(x::T,  
    sample_rate_Hz::Real) where
    (T<:ArrayofReals)

    w, X = rfft_magnitude_norm(x, sample_rate_Hz)
    
    return w, 20 .*log10.(X)
end

function rfft_phase(x::T,  
    sample_rate_Hz::Real) where
    (T<:ArrayofReals)

    w, X = rfft(x, sample_rate_Hz)
    
    return (w, angle.(X))
end

##########
# PSD
##########

function psd(x::T,
    sample_rate_Hz::Real) where
    (T<:ArrayofReals)
    # real-valued version

    w, psd = psd_linear(x, sample_rate_Hz)

    return w, 10 .* log10.(psd)
end

function psd_linear(x::T,
    sample_rate_Hz::Real) where
    (T<:ArrayofReals)
    # real-valued version

    w, X = rfft_magnitude(x, sample_rate_Hz)
    
    To = sample_rate_Hz * length(x)
    PSD = 2/To .* X 

    return w, PSD
end

function psd_double(x::T,
    sample_rate_Hz::Real) where
    (T<:ArrayofNumbers)

    w, X = fft_magnitude(x, sample_rate_Hz)

    To = sample_rate_Hz * length(x)
    PSD = 1/To .* X 

    return w, 10 .* log10.(PSD)
end

function psd_double_linear(x::T,
    sample_rate_Hz::Real) where
    (T<:ArrayofNumbers)

    w, X = fft_magnitude(x, sample_rate_Hz)

    To = sample_rate_Hz * length(x)
    PSD = 1/To .* X 

    return w, PSD
end

function psd_OMA_favg(directory::AbstractString;
    peak=0  # nonzero value is the w around peak to record
    )
    # take rolling average of binary words (OMA) data PSD
    if !isdir(directory)
        println("Argument must be a directory.")
        return 0, 0
    end

    n_files = length(readdir(directory))
    if directory[end] != '/'
        directory *= '/'
    end
    
    t_temp, y_temp = load_binary_words(directory * 
        readdir(directory)[1])
    fs = 1/(t_temp[begin+1] - t_temp[begin])
    w, S_temp = psd(y_temp, fs)
    result_length = length(S_temp)
    if peak > 0
        # shift w to BB and make S_temp the width
        w = w[begin+5:end]
        S_temp = S_temp[begin+5:end]
        w_s, S_temp_s = Telecoms.data_around_peak(w, S_temp, peak)
        if w_s == -1
            println("fail")
            return -1
        end
        result_length = length(S_temp_s)
        wm = w_s[end] - w_s[begin]
        w = collect([-wm/2:fs:wm/2])
    end
    S = Vector{Float64}(undef, result_length)

    i = 1

    for filename in readdir(directory)
        print("\r                      \r" * string(i-1) * "/" * 
            string(n_files))
        x_temp, y_temp = Telecoms.load_binary_words(directory 
            * filename)
        w_tmp, S_temp = psd(y_temp, fs)
        if peak > 0
            w_tmp = w_tmp[begin+5:end]
            S_temp = S_temp[begin+5:end]
            w_s, S_temp_s = Telecoms.data_around_peak(w_tmp, S_temp, peak)
            S_temp = S_temp_s
        end
        @. S = S*(i-1)/i + S_temp/i
        i += 1

    end
    print("\r                      \r" * string(n_files) * "/" * 
        string(n_files))

    return w, S 
end

##########
# IFFT
##########

function ifft(X::T,  
    sample_rate_Hz::Real) where
    (T<:ArrayofNumbers)

    x = FFTW.ifft(X)
    
    sample_period = 1/sample_rate_Hz
    t = collect(0:sample_period:sample_period*length(x))

    return t, x
end

function ifft(i::T,
    q::T,
    sample_rate_Hz::Real) where
    (T<:ArrayofNumbers)

    return ifft(Complex.(i, q), sample_rate_Hz)
end

##########
# IRFFT
##########

function irfft(X::T,  
    sample_rate_Hz::Real) where
    (T<:ArrayofNumbers)

    # the amazing digital irfft circ- BRAAP

    transformed_length = 2*size(X,1)-1
    x = FFTW.irfft(X, transformed_length)
    
    sample_period = 1/sample_rate_Hz
    t = collect(0:sample_period:sample_period*transformed_length)

    return t, x
end

function irfft(i::T,
    q::T,
    sample_rate_Hz::Real) where
    (T<:ArrayofReals)

    return irfft(Complex.(i, q), sample_rate_Hz)
end

##########
# Analysis of distributions
##########

function _nwmm(w::T,
    X::T,
    n::Number, 
    m::Number) where (T<:ArrayofNumbers)
    # n width m maximum

    peak = maximum(X)
    hm = m * peak

    points = find(X, hm)

    if length(points) == 2
        wh = w[points[2][2]]
        wl = w[points[1][2]]
        
        return abs(wh - wl) * n
    end

    return points
end

function _nwLm(w::T,
    X::T,
    n::Number, 
    L_dB::Number) where (T<:ArrayofNumbers)
    # n width m maximum

    peak = maximum(X)
    hm = peak - L_dB

    points = find(X, hm)

    if length(points) == 2
        wh = w[points[2][2]]
        wl = w[points[1][2]]
        
        return abs(wh - wl) * n
    end

    return points
end

function fwhm(w::T,
    X::T) where (T<:ArrayofNumbers)

    return _nwmm(w, X, 1.0, 0.5)
end

function hwhm(w::T,
    X::T) where (T<:ArrayofNumbers)

    return _nwmm(w, X, 0.5, 0.5)
end

function fwhm_dB(w::T,
    X::T) where (T<:ArrayofNumbers)

    return _nwLm(w, X, 1.0, -10*log10(1/2))
end

function hwhm_dB(w::T,
    X::T) where (T<:ArrayofNumbers)

    return _nwLm(w, X, 0.5, -10*log10(1/2))
end

# none of these fit functions even work

function fit_Lorentzian(w::T, 
    X::T) where (T<:ArrayofNumbers)

    if !_equal_lengths(w, X)
        return -1
    end

    @. lorentzian_model(x, p) = p[1] / ((x - p[2])^2 + p[4]) + p[3] 
    #@. lorentzian_model(x, p) = 1 / (pi * p[1] * (1 + ((x+p[3])/p[1])^2)) + p[2];

    p0 = [1,80e6,mean(X),1]

    fit = curve_fit(lorentzian_model, w, X, p0)
    p1 = coef(fit);
    println(p1)

    lorentzian_psd = lorentzian_model(w, p1);

    return [w, lorentzian_psd]
end

function fit_Gaussian(w::T, 
    X::T) where (T<:ArrayofNumbers)
    # not tested this code

    if !_equal_lengths(w, X)
        return -1
    end

    @. gaussian_model(x, p) = p[1] * exp(-(x+p[3])^2 / (2*p[2]^2))

    # needs optimising
    p0 = [1.0, mean(X), -w[round(Int, end/2)]] 

    fit = curve_fit(gaussian_model, w, X, p0)

    p1 = coef(fit);
    gaussian_psd = gaussian_model(w, p1);

    return [w, gaussian_psd]
end

##########
# Windows
##########

# calculate the equivalent noise bandwidth of a given window
function ENBW(window::T) where (T<:ArrayofNumbers)

    #time domain after parseval
    top = length(window) * sum(window.^2)
    bottom = sum(window)^2
    enbw = top / bottom

    return enbw
end

function generate_Hanning_window(N::Integer)

    w = collect(0:1:N)

    return sin.((pi .* w) ./ N).^2;
end

function Hanning_window(x::T) where (T<:ArrayofNumbers)

    N = length(x) - 1;
    window = generate_Hanning_window(N)

    return window .* x
end

##########
# Phase and frequency analysis 
##########

function analytic(x::T) where (T<:ArrayofReals)

    return DSP.hilbert(block_DC_fft(x))
end

function unwrap_analytic(x::T) where (T<:ArrayofNumbers)

    phi_u = Vector{Float64}(undef, length(x))
    phi_u[1] = angle.(x[1]) 
    
    for idx = 2:length(phi_u)
        phi_u[idx] = phi_u[idx-1] + angle.(x[idx] * conj(x[idx-1]))
    end

    return phi_u
end

function unwrap(phi::T) where (T<:ArrayofReals)

    return DSP.unwrap(phi)
end

function frequency_spectrum(t::T,
    x::T;
    f_0=8e7) where (T<:ArrayofNumbers)

    x = Telecoms.block_DC_fft(x)
    sample_rate = 1/(t[2]-t[1])
    f_psd, S = Telecoms.psd(Telecoms.Hanning_window(x), sample_rate)
    phase = Telecoms.unwrap_analytic(Telecoms.analytic(x))
    phi_ideal = 2*pi*f_0.*t
    phi = phase .- phi_ideal
    if !Base.require_one_based_indexing(phi)
        println("Index should start at 1 for follwing operation.")
        return Nothing, Nothing, Nothing, Nothing
    end
    # differentiate phase
    df = zeros(length(phi))
    for idx in eachindex(phi)
        if idx == 1
            continue
        end
        df[idx] = phi[idx] - phi[idx - 1]
        df[idx] *= sample_rate/(2*pi)
    end
    problem = SciMLBase.SampledIntegralProblem(abs.(df).^2,t)
    ret = CommonSolve.solve(problem, Integrals.TrapezoidalRule())
    pre_hanning_power = ret.u / t[end]
    df = Telecoms.Hanning_window(df)
    f, S_f = Telecoms.psd_linear(df,
        sample_rate)
    problem = SciMLBase.SampledIntegralProblem(abs.(df).^2,t)
    ret = CommonSolve.solve(problem, Integrals.TrapezoidalRule())
    post_hanning_power = ret.u / t[end]
    S_f *= pre_hanning_power/post_hanning_power

    return f, S_f, f_psd, S
end

function frequency_spectrum_favg(filenames::AbstractArray{<:AbstractString},
    data_filename::AbstractString;
    f_0=8e7,
    quiet=false)

    t, x = Telecoms.load_binary_words(filenames[begin])
    f, S_f, f_psd, S = Telecoms.frequency_spectrum(t, x; f_0=f_0)
    S_f_avg = zeros(length(S_f))

    avg_linewidth = 0
    f_avg = 0
    lws = zeros(length(filenames))
    i = 1

    for fn in filenames
        if !quiet
            print("\r                   \r", 
                i, "/", length(filenames))
        end

        t, x = Telecoms.load_binary_words(fn)
        f, S_f, f_psd, S = Telecoms.frequency_spectrum(t, x; f_0=8e7)
        
        @. S_f_avg = S_f_avg*(i-1)/i + S_f/i
        f_avg = f

        i += 1
    end

    println()

    Telecoms.save_as_binary_floats(f_avg, S_f_avg, data_filename)

    return f_avg, S_f_avg
end

##########
# Filtering
##########

function block_DC_fft(x::ArrayofNumbers)
    israel = false
    isimaginary = false  # might be useful for stuff like DMT
    if all(y -> imag(y) == 0, x)
        israel = true
    elseif all(y -> real(y) == 0, x)
        isimaginary = true
    end

    _, X = fft(x, 10)

    X[1] = 0

    _, x = ifft(X, 10)

    if israel
        return real.(x)
    elseif isimaginary
        return imag.(x)
    end

    return x
end

function block_DC_mean(x::ArrayofNumbers)
    mu = Statistics.mean(x)

    @. x -= mu

    return x
end

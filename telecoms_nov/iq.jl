# for generation of signals etc.

# IQ functions take separate inphase and quadrature arguments

##########
# IQ wrappers
##########

function iq_amplitude(i::T,
    q::T) where (T <: AbstractArray{<:AbstractFloat})
    
    @. return abs(Complex(i, q))
end


function iq_power(i::T,
    q::T) where (T <: AbstractArray{<:AbstractFloat})
    
    @. return i^2 + q^2 
end

##########
# Generation of signals
##########

function generate_Lorentzian_signal(f_0, 
    linewidth, 
    f_s,
    N::Integer)

    sample_period = 1/f_s

    phase_increment = 2 * pi * linewidth * sample_period
    dphi_n = sqrt(phase_increment) .* randn(N)
    phi_n = cumsum(dphi_n)

    phase_term = zeros(length(phi_n))
    t = collect(0:length(phi_n)-1).*sample_period
    @. phase_term = 2*pi*f_0*t+phi_n
    x = cis.(phase_term)

    return t, x
end


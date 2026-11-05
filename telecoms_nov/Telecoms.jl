# Signal analysis module for optical fibre communications

module Telecoms

# signal analyses
include("frequency_analysis.jl")
include("generic_analysis.jl")
include("iq.jl")

# instrument control
include("oma.jl")
include("bosa.jl")
include("exfo/ltb-1.jl")

# utility
include("util.jl")

function capture_for_lw_measurement(ip::AbstractString,
    root::AbstractString,
    channel::Integer,
    n_captures::Integer,
    f_0,
    spectrum_filename::AbstractString)

    oma = OMA_connect(ip)
    OMA_setup(oma)

    data_filenames = OMA_capture_n_times(oma, n_captures, channel, root)

    close(oma)

    frequency_spectrum_favg(data_filenames, spectrum_filename;
        f_0=f_0)

    # to see stability, calculate PSD of some data
    filename = data_filenames[div(n_captures,2)]
    t, x = Telecoms.load_OMA(filename)
    fs = 1/(t[begin+1] - t[begin])
    w, psd = Telecoms.psd(x, fs)
    Telecoms.save_as_binary_floats(w[begin:end-1], psd[begin:end-1], 
        filename * "_psd")
end
end

#include("../util.jl")

function OPM_power_measurement(
    instrument::Sockets.TCPSocket,
    wavelength,
    n,
    rate_Hz)

    ready_condition = "READY"

    println(TCP_identify(instrument))
    println(Telecoms.TCP_query(instrument, "STAT?"))
    
    Telecoms.TCP_write(instrument, "UNIT:POW DBM")

    for i in 1:4
        Telecoms.TCP_write(instrument, 
            "SENS" * string(i) * ":POW:WAV " * string(wavelength) * " nm")
        Telecoms.TCP_write(instrument, 
            "SENS" * string(i) * ":POW:REF:STAT 0")
        Telecoms.TCP_write(instrument, 
            "SENS" * string(i) * ":POW:RANG:AUTO 1")
        Telecoms.TCP_write(instrument, 
            "SENS" * string(i) * ":FREQ:NCON " * string(rate_Hz))
        Telecoms.TCP_write(instrument, 
            "SENS" * string(i) * ":AVER:STAT 0")
        Telecoms.TCP_write(instrument, 
            "TRAC" * string(i) * ":POIN TRC" * string(i) * "," * string(n))
    end

    Telecoms.TCP_write(instrument, "INIT:AUTO 1, NCON")

    # busywait slopfest
    while Telecoms.TCP_query(instrument, "STAT?") != ready_condition
    end

    traces = []
    for i in 1:4
        trace_data = Telecoms.TCP_query(instrument,
            "TRAC? TRC", string(i))
        push!(traces, trace_data)
    end

end

    

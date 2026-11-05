using Sockets
using Dates

include("oma.jl")

function bosa_connect(ip::AbstractString;
    port::Integer = 10000)

    return Sockets.connect(ip,port)
end

##### From oma.jl
# disconnect
#####

function _bosa_write(bosa::Sockets.TCPSocket,
    SCPI_command::AbstractString)

    rc = write(bosa, SCPI_command)
    if rc != length(SCPI_command)
        println("_bosa_write: write error")
        return -1
    end

    return rc
end

function _bosa_query(bosa::Sockets.TCPSocket,
    SCPI_command::AbstractString;
    n_bytes::Integer = 0)
    # not specifying n_bytes results in more allocations

    if _bosa_write(bosa, SCPI_command) < 0 
        return -1
    end

    if eof(bosa)
        println("_bosa_query: EOF")
        return -1
    end

    if n_bytes < 1
        response = Vector{UInt8}()
        for c in readeach(bosa, UInt8)
            push!(response, c)
            if bytesavailable(bosa) < 1
                return response
            end
        end
    end

    response = Vector{Vector{UInt8}}()

    bytes_read = 0
    while bytes_read < n_bytes
        b = read(bosa, n_bytes)
        bytes_read += length(b)
        push!(response, b)
    end

    return [(response...)...]
end

function _bosa_query_string(bosa::Sockets.TCPSocket,
    SCPI_command::AbstractString)

    response = _bosa_query(bosa, SCPI_command)

    if length(response) < 1
        return -1
    end

    return String(response)
end

function bosa_identify(bosa::Sockets.TCPSocket)
    return _bosa_query_string(bosa, "*IDN?")
end

function _bosa_fetch(bosa::Sockets.TCPSocket;
    filename::AbstractString = "bosa_trace")
    # return trace data as raw format

    if _bosa_query_string(bosa, "format real") != "OK\r\n"
        println("_bosa_fetch: error setting format.")
    end

    n_points = _bosa_query_string(bosa, "trace:count?")
    n_points = parse(Int, n_points)
    n_bytes = n_points * 8 * 2

    start_time = time_ns()

    trace_data = _bosa_query(bosa, "trace?", n_bytes=n_bytes)

    elapsed_time = (time_ns() - start_time) / 1e9

    if length(trace_data) != n_bytes
        println("_bosa_fetch: did not fetch all content")
    end

    transfer_speed = 8 * n_bytes / elapsed_time
    println("Fetched ", n_bytes / 1e6 , " MB (", transfer_speed / 1e6, " Mb/s)")

    n_written = write(filename, trace_data)
    
    if n_written != n_bytes
        println("_bosa_fetch: did not save all content") 
    end
    if n_written != length(trace_data)
        println("_bosa_fetch: did not save all content captured") 
    end

    return trace_data
end

function _bosa_parse(bosa_raw_data::Vector{UInt8})
    # bosa raw data to w,psd

    if length(bosa_raw_data) / 16 % 1 != 0
        println("_bosa_parse: raw data n_points not multiple of 16")
        return -1
    end

    bosa_raw_data = reinterpret(Float64, bosa_raw_data)
    w = bosa_raw_data[1:2:length(bosa_raw_data)]
    psd = bosa_raw_data[2:2:length(bosa_raw_data)]

    return w, psd
end

function bosa_metadata(bosa::Sockets.TCPSocket;
    filename::AbstractString = "bosa_metadata")

    output = "BOSA trace metadata "
    output *= Dates.format(Dates.now(), ISODateTimeFormat)
    output *= ":\n"

    output *= "Measurement mode: " * _bosa_query_string(bosa, 
        "instrument:state:mode?")
    output *= "\nPolarisation: " * _bosa_query_string(bosa,
        "input:polarization?")
    output *= "\nResolution: " * _bosa_query_string(bosa,
        "sense:wavelength:resolution?")
    output *= "\nSUT power: " * _bosa_query_string(bosa,
        "input:power?")
    output *= "\nPoints: " * _bosa_query_string(bosa,
        "trace:count?")
    output *= "\nIntegral power: " * _bosa_query_string(bosa,
        "calculate:tpower?")  # only records if set on front panel
    output *= "\n"

    write(filename, output)
end

function bosa_trace(bosa::Sockets.TCPSocket)
    return _bosa_parse(_bosa_fetch(bosa))
end

function bosa_get_power(bosa::Sockets.TCPSocket)
    # return power in dBm
end

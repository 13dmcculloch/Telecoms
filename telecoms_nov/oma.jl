using Sockets

function OMA_connect(ip::AbstractString; port::Integer = 5025)
    # port for OMA is default 
    return Sockets.connect(ip,port)
end

function disconnect(OMA::Sockets.TCPSocket)
    close(OMA)
end

function twrite(OMA::Sockets.TCPSocket,
    SCPI_command::AbstractString;
    termination_characters::AbstractArray = ['\n'])

    termination_string = join(termination_characters)
    SCPI_command = join([SCPI_command, termination_string])

    rc = write(OMA, SCPI_command)
    if rc != length(SCPI_command)
        return -1
    end

    return rc
end

function query(OMA::Sockets.TCPSocket,
    SCPI_command::AbstractString;
    termination_characters::AbstractArray = ['\n'])
    # have discovered a thing where a type hinted kwarg messes up the
    # function definition if default kwarg value is not the hinted type
    
    if twrite(OMA, SCPI_command; 
        termination_characters = termination_characters) < 0
        return -1
    end

    response = readuntil(OMA, join(termination_characters))

    return response
end

function identify(instrument::Sockets.TCPSocket;
    termination_characters::AbstractArray = ['\n'])

    return query(instrument, "*IDN?"; 
        termination_characters = termination_characters)
end

function OMA_setup(OMA::Sockets.TCPSocket)
    # write a bunch of stuff to OMA to set it up for efficient transfer

    # fastest transfer is LE which is not default
    if twrite(OMA, ":WAVeform:BYTeorder LSBFirst") < 0
        return -1
    end

    # most precise data
    if twrite(OMA, ":WAVeform:FORMat WORD") < 0
        return -1
    end

    # easier to interpret packet
    if twrite(OMA, ":WAVeform:STReaming 1") < 0
        return -1
    end

    return 0
end

function _OMA_fetch(OMA::Sockets.TCPSocket;
    termination_character::AbstractChar = '\n')
    # arithmetic here is for WORD type. call OMA_setup before this
    # Old OMA: <800 Mb/s. New OMA: 947 Mb/s.

    n_points = query(OMA, ":WAVeform:POINts?")
    n_points = parse(Int, n_points)
    n_bytes = n_points * 2 + 1
    data = Vector{UInt8}(undef, n_bytes)  # allocate to speed up transfer

    if twrite(OMA, ":WAVeform:DATA?") < 0
        return -1
    end

    readuntil(OMA, '#')  # worked on old OMA
    readuntil(OMA, '0')

    byte::Vector{UInt8} = [0]
    n = 0
    start_time = time()

    read!(OMA, data)

    end_time = time()
    duration_s = end_time - start_time
    bytes_transferred = sizeof(data)
    speed = bytes_transferred / duration_s
    MB_transferred = bytes_transferred / 1e6
    println("Received ", MB_transferred, " MB (", 
        round(speed * 8 / 1e6, digits=2), " Mb/s.)")
    return data[1:end-1]  # cut out the \n
end

function _store_floats_from_words(words::Vector{UInt8},
    x::Tuple,
    y::Tuple,
    filename::AbstractString)

    if length(words) % 2 != 0
        println(length(words))
        println(words[end])
        println("Odd number of word bytes.")
        return -1
    end

    fp = open(filename, "w")
    write(fp, x[1])
    write(fp, x[2])

    word_v::Vector{UInt8} = [0, 0] 
    word::Int16 = 0
    for idx in 1:2:length(words)
        word = Int16(words[idx]) | Int16(words[idx+1]) << 8

        f = _interpret_word(word, y)

        write(fp, f)
    end

    close(fp)

    return 0
end

function load_floats_from_binary(filename::AbstractString)

    filesize = stat(filename).size
    n_bytes = convert(Int, (filesize - 4 * 8))
    if n_bytes < 1 
        println("load_floats...: not enough data")
        return -1
    elseif n_bytes % 2 != 0
        println("load_floats...: odd number of data bytes")
        return -1
    end

    n_floats = convert(Int, n_bytes / 2)
    
    floats_y = Vector{Float64}(undef, n_floats)
    floats_x = Vector{Float64}(undef, n_floats)

    fp = open(filename, "r")
    y1 = read(fp, Float64)
    y2 = read(fp, Float64)
    x1 = read(fp, Float64)
    x2 = read(fp, Float64)
    
    x = (x1, x2)
    y = (y1, y2)

    for idx in 1:length(floats_x)
        floats_x[idx] = x[2] + idx * x[1]
    end

    words::Vector{UInt8} = [0, 0] 
    word::Int16 = 0
    for idx in 1:length(floats_x)
        read!(fp, words)

        word = _pack_word(words)

        floats_y[idx] = _interpret_word(word, y)
    end

    close(fp)

    return floats_x, floats_y 
end

function load_OMA(filename::AbstractString)
    
    return load_floats_from_binary(filename)
end

function _store_ascii_floats_from_words(words::Vector{UInt8},
    x::Tuple,
    y::Tuple,
    filename::AbstractString)

    if length(words) % 2 != 0
        println(length(words))
        println(words[end])
        println("Odd number of word bytes.")
        return -1
    end

    fp = open(filename, "w")
    write(fp, string(convert(Float64, x[1])))
    write(fp, '\n')
    write(fp, string(convert(Float64, x[2])))
    write(fp, '\n')

    word_v::Vector{UInt8} = [0, 0] 
    word::Int16 = 0
    for idx in 1:2:length(words)
        word = Int16(words[idx]) | Int16(words[idx+1]) << 8

        f = _interpret_word(word, y)

        write(fp, string(convert(Float64, round(f, digits=10))))
        write(fp, '\n')
    end

    close(fp)

    return 0
end

function _pack_word(words::Vector{UInt8})

    word = Int16(words[1]) | Int16(words[2]) << 8

    return word
end

function _interpret_word(word::Int16, factor::Tuple)
    # convert word to float with given orig, increment

    return word * factor[1] + factor[2]
end

function _OMA_acquire_channel_n(OMA::Sockets.TCPSocket,
    channel::Integer)
    # p. 1607

    channel_s = join(["CHANnel", string(channel)])

    command = join([":DIGitize ", channel_s])
    if twrite(OMA, command) < 0
        return -1
    end

    command = join([":", channel_s, ":DISPlay ON"])
    if twrite(OMA, command) < 0
        return -1
    end

    return 0
end

function _OMA_acquire_all(OMA::Sockets.TCPSocket)
    # p. 294 
    for channel in 1:4
        # If you use the :DIGitize command with no parameters, 
        # the digitize operation is
        # performed on the channels that are being displayed 
        # in the Infiniium waveform viewing area.
        channel_s = join(["CHANnel", string(channel)])
        command = join([":", channel_s, ":DISPlay ON"])
        if twrite(OMA, command) < 0
            return -1
        end
    end

    command = ":DIGitize"
    if twrite(OMA, command) < 0
        return -1
    end

    return 0
end

function OMA_channel_n_data(OMA::Sockets.TCPSocket,
    channel::Integer,
    filename::AbstractString)
    # save to binary

    if _OMA_acquire_channel_n(OMA, channel) < 0
        println("Failed to acquire waveform data.")
        return -1
    end

    channel_s = join(["CHANnel", string(channel)])
    command = join([":WAVeform:SOURce ", channel_s])
    if twrite(OMA, command) < 0
        return -1
    end

    # then get the data...
    y = _OMA_voltage_conversion_factor(OMA)
    x = _OMA_time_conversion_factor(OMA)

    if length(y) < 2 || length(x) < 2
        return -1
    end

    data = _OMA_fetch(OMA)
    if length(data) < 2
        return -1
    end

    fp = open(filename, "w")

    write(fp, convert(Float64, y[1]))
    write(fp, convert(Float64, y[2]))
    write(fp, convert(Float64, x[1]))
    write(fp, convert(Float64, x[2]))
    
    write(fp, data)

    println("Written to ", filename)

    close(fp)

    return 0
end

function OMA_channel_all_data(OMA::Sockets.TCPSocket,
    filename::AbstractString)
    # waveform data for all OMA channels

    if _OMA_acquire_all(OMA) < 0
        println("Failed to acquire waveform data.")
        return -1
    end

    success = true
    for i in 1:4
        filename_s = join([filename, string(i)])
        if OMA_channel_n_data(OMA, n, filename_s) < 0
            println("Failed to transfer and save ", filename_s)
            success = false
        end
    end

    return ~success
end

function _OMA_voltage_conversion_factor(OMA::Sockets.TCPSocket)
    # voltage is y axis
    # return as (yinc, yorg)

    # y increment between integer data points
    yinc = query(OMA, ":WAVeform:YINCrement?")
    # value from which above is a deviation
    yorg = query(OMA, ":WAVeform:YORigin?")

    if length(yinc) < 2 || length(yorg) < 2 
        println("Failed to get voltage conversion factor.")
        return -1
    end

    yinc = _read_double(yinc)
    yorg = _read_double(yorg)

    return yinc, yorg
end

function _OMA_time_conversion_factor(OMA::Sockets.TCPSocket)
    # time is x axis
    # return as (xinc, xorg)

    # y increment between integer data points
    xinc = query(OMA, ":WAVeform:XINCrement?")
    # value from which above is a deviation
    xorg = query(OMA, ":WAVeform:XORigin?")

    if length(xinc) < 2 || length(xorg) < 2 
        println("Failed to get time conversion factor.")
        return -1
    end

    xinc = _read_double(xinc)
    xorg = _read_double(xorg)

    return xinc, xorg
end

function _read_double(s::AbstractString)
    # according to scanf example on p. 1606
    
    return parse(Float64, s)
end

function _floats(filename::AbstractString; n=10)
    f = open(filename, "r")

    for i in 1:n
        println(read(f, Float64))
    end

    close(f)

    return 0
end

function OMA_capture_n_times(oma::Sockets.TCPSocket,
    n_captures::Integer, 
    ch, 
    fnam::AbstractString)

    filenames = Vector{String}(undef, n_captures)

    for i in 0:n_captures-1
        output_fn = fnam * "_" * string(i) 
        filenames[i+1] = output_fn
        println("Capture ", output_fn)
        if typeof(ch) <: Integer
            if ch < 1 || ch > 4
                println("Invalid channel.")
                break
            end
            Telecoms.OMA_channel_n_data(oma, ch, output_fn)
        else
            Telecoms.OMA_channel_all_data(oma, output_fn)
        end
    end

    return filenames
end


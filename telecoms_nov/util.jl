#####
# CSV files
#####

function save_as_csv(data::Vector{<:Vector{<:Any}},
    filename::AbstractString;
    display_ratio = 7374,
    quiet = false,
    maxfilesize = 0)

    if maxfilesize > 0
        rows = length(data[1])
        cols = length(data)
        n_punctuation = rows * (cols)  # a,b,\n...
        size = n_commas
        for row in 1:rows
            for col in 1:cols
                size += sizeof(string(data[col][row]))
            end
        end

        size_MB = size / 1e6
        if maxfilesize > size_MB
            println("save_as_csv: max file size exceeded")
            return -1
        end
    end

    filename *= ".csv"  # critical for not overwriting original data...

    if ~quiet
        println("Saving to ", filename)
    end

    fp = open(filename, "w")

    rows = length(data[1])
    cols = length(data)

    for row in 1:rows 
        line = ""
        for col in 1:cols 
            line *= string(data[col][row]) * ","
        end
        line = line[1:end-1] * "\n"
        write(fp, line)
        if row % display_ratio == 0 && ~quiet
            print("\r              \r", 100 * row / rows, " %")
        end
    end

    if ~quiet
        print("\r                 \r", 100 * 1, " %\n")
    end

    close(fp)

    return 0
end

function load_csv(filename::AbstractString;
    header=false,
    delim=',',
    n_cols=0,
    display_ratio = 7374,
    quiet=false)

    fp = open(filename, "r")

    readline(fp)
    firstline = readline(fp)

    n_lines = 2
    for line in eachline(fp)
        n_lines += 1
    end
    
    if header
        n_lines -= 1
    end

    data = []
    for field in split(firstline, delim)
        if occursin(r"[+-]?(\d+([.]\d*)?([eE][+-]?\d+)?|[.]\d+([eE][+-]?\d+)?)", field)
            push!(data, Vector{Float64}(undef, n_lines))
        else
            push!(data, Vector{Float64}(undef, n_lines))
        end
    end

    seek(fp, 0)
    if header
        readline(fp)
    end

    idx = 1
    for line in eachline(fp)
        field_idx = 1
        for field in split(line, delim)
            try 
                data[field_idx][idx] = parse(Float64, field)
            catch
                data[field_idx][idx] = -90.0
            end
            field_idx += 1
        end
        idx += 1
    end

    if idx - 1 != n_lines
        println("load_csv: loaded ", idx, " rows out of ", n_lines)
    end

    close(fp)

    return data
end


#####
# My binary files
#####

function save_as_binary_floats(x::Vector{<:Number}, 
    y::Vector{<:Number}, 
    filename::AbstractString;
    maxfilesize = 0,
    quiet = false)
    # for x,y data with x axis origin and offset

    if maxfilesize > 0
        n_points = length(y) + 2  # xorg and xinc
        size = n_points * sizeof(Float64) 
        size_MB = size / 1e6
        if maxfilesize > size_MB
            println("save_as_binary_floats: max file size exceeded")
            return -1
        end
    end

    filename *= ".dat"

    if ~quiet
        println("Saving to ", filename)
    end

    xorg = x[1]
    xinc = (x[end] - x[1]) / length(x)

    fp = open(filename, "w")

    write(fp, convert(Float64, xorg))
    write(fp, convert(Float64, xinc))

    y = convert.(Float64, y)

    write(fp, y)

    close(fp)
end

function save_as_binary_floats(x::Vector{<:Number}, 
    filename::AbstractString;
    maxfilesize = 0,
    quiet = false)
    # for 1D vector 

    if maxfilesize > 0
        n_points = length(x) + 2  # xorg and xinc
        size = n_points * sizeof(Float64) 
        size_MB = size / 1e6
        if maxfilesize > size_MB
            println("save_as_binary_floats: max file size exceeded")
            return -1
        end
    end

    filename *= ".dat"

    if ~quiet
        println("Saving to ", filename)
    end

    fp = open(filename, "w")

    x = convert.(Float64, x)

    write(fp, x)

    close(fp)
end

function load_binary_floats_waveform(filename::AbstractString)
    # for x,y data with x axis origin and offset
    n_floats = stat(filename).size - 2 * 8

    fp = open(filename, "r")
    xorg = read(fp, Float64)
    xinc = read(fp, Float64)

    n_floats = convert(Int64, n_floats / 8)

    floats_y = Vector{Float64}(undef, n_floats)
    floats_x = Vector{Float64}(undef, n_floats)

    read!(fp, floats_y)

    close(fp)

    for idx in 1:n_floats
        floats_x[idx] = xorg + (idx-1) * xinc
    end

    return floats_x, floats_y
end

function load_binary_floats(filename::AbstractString)
    # 1D version
    n_floats = stat(filename).size - 2 * 8

    fp = open(filename, "r")

    n_floats = convert(Int64, n_floats / 8)

    floats_x = Vector{Float64}(undef, n_floats)

    read!(fp, floats_x)

    close(fp)

    return floats_x
end



#####
# Keysight OMA normalised integers
#####

function load_binary_words(filename::AbstractString)

    function _pack_word(words::Vector{UInt8})

        word = Int16(words[1]) | Int16(words[2]) << 8

        return word
    end

    function _interpret_word(word::Int16, factor::Tuple)
        # convert word to float with given orig, increment

        return word * factor[1] + factor[2]
    end

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

function load_OMA_words(filename::AbstractString)
    return load_binary_words(filename)
end

#####
# Agilent OSA floppy - modified input file
#####

function load_agilent_csv(filename::AbstractString)
    # Must already be processed by the shell script:
    # Loads in modified start,stop,length,trace agilent OSA csv
    data = Telecoms.load_csv(filename)
    start = data[2][1]
    stop = data[2][2]
    length = data[2][3]
    step = (stop - start) / convert(Float64, length)
    w = collect(start:step:stop-step)
    S = data[2][4:end]
    return w, S
end

#####
# Old instrument binary files
#####

function load_mathcad(filename::AbstractString)
    # output a Mathcad-type file to array of floats
    # TDS520C produces this for example

    fp = open(filename, "r")

    n_samples = parse(UInt64, readline(fp)) 
    sample_period = parse(Float64, readline(fp)) 
    readline(fp)  # not sure what this value is
    x_orig = parse(Float64, readline(fp)) 

    x = Vector{Float64}(undef, n_samples)
    y = Vector{Float64}(undef, n_samples)

    for i in 1:n_samples
        x[i] = (i-1) * sample_period + x_orig
        y[i] = parse(Float64, readline(fp)) 
    end

    close(fp)

    return (x, y)
end


#####
# interprets the data from someone's DAQ, for DAS applications
#####

function load_DAQ2(dirname::AbstractString;
    n_files=0,
    fast_range=(0, 0),  # specifically the length sample
    decimation=0,
    fast_avg=0,
    slow_avg=0,
    fastslow_avg=0,
    memory_limit=0,
    axes=false)
    # read in the sessxxx folder of ndxxxx.bin files from DAQ2

    # fast = spatial; slow = temporal
    # fast resolution is N+1 * data_resolution
    # decimation undersamples the input data by N+1
    # fast_avg takes mean of N+1 points reducing spatial resolution
    # slow_avg takes mean of N+1 points like an LPF
    # fastslow_avg takes mean of (N+1)*(N+1) blocks quickly and slowly

    # check if supplied distance range makes sense
    fast_start = fast_range[1]
    fast_end = fast_range[2]
    if fast_end < fast_start
        println("Final sample must be greater than start sample.")
        return
    end

    # get to sessxxx folder
    old_path = pwd()
    cd(dirname)
    files = readdir()
    while length(files) == 1 && ~occursin(".bin", files[1])
        cd(files[1])
        files = readdir()
    end
    if n_files > 0 && n_files <= length(files)
        files = files[1:n_files]
    end

    # set up metadata with header
    traces_count = 0
    rep_rate = 0
    sampling_freq = 0
    trace_length = 0
    for f in files
        trace_length, rep_rate, sampling_freq, n_tr = _DAQ2_header(f)
        traces_count += n_tr
    end

    # deal with ranges
    if fast_end > trace_length
        println("Specify end fast bin below trace length ", trace_length)
        cd(old_path)

        # this will help with automatically creating the correct
        # range when looping through data of diverse lengths
        return trace_length
    end
    if fast_end == 0
        memory_footprint = traces_count * trace_length * 2
    else
        memory_footprint = ((fast_end - fast_start) 
            * traces_count * 2)
    end

    if decimation > 0
        memory_footprint /= decimation + 1
    elseif fast_avg > 0
        memory_footprint /= fast_avg + 1
    elseif slow_avg > 0
        memory_footprint /= slow_avg + 1
    elseif fastslow_avg > 0
        memory_footprint /= (fastslow_avg + 1)^2
    end
    
    if memory_limit < 1
        # Julia will not GC old data unless forced to, so using 
        # free_memory is not working
        # Julia uses about 25% of memory with our libraries loaded
        memory_limit = Sys.total_memory() * 0.75
    end

    if memory_footprint > memory_limit
        println("Memory footprint of folder ", 
            memory_footprint / 1e9, " GB")
        println("Free memory ", memory_limit / 1e9, " GB")
        println("Memory limit exceeded by ", 
            (memory_footprint - memory_limit) / 1e9, " GB")
        cd(old_path)
        ret = 0
        GC.gc()
        return
    end

    if decimation > 0
        println("Decimation factor: ", decimation)
    elseif fast_avg > 0
        println("Fast averaging length factor: ", fast_avg)
    elseif slow_avg > 0
        println("Slow averaging length factor: ", slow_avg)
    elseif fastslow_avg > 0
        println("Fastslow averaging length factor: ", fastslow_avg)
    end
    
    if fast_end == 0
        println("Loading data samples from disk...")
    else
        println("Loading data samples ", fast_start, " to ", fast_end, 
            " from disk...")
    end

    if decimation > 0
        ret = _load_DAQ2bin_decimation(files,
            decimation+1,
            trace_length,
            traces_count,
            memory_footprint,
            fast_range)

    elseif fast_avg > 0
        ret = _load_DAQ2bin_fast(files,
            fast_avg+1,
            trace_length,
            traces_count,
            memory_footprint,
            fast_range)

    elseif slow_avg > 0
        ret = _load_DAQ2bin_slow(files,
            slow_avg+1,
            trace_length,
            traces_count,
            memory_footprint,
            fast_range)

    elseif fastslow_avg > 0
        ret = _load_DAQ2bin_fastslow(files,
            fastslow_avg+1,
            trace_length,
            traces_count,
            memory_footprint,
            fast_range)

    else
        ret = _load_DAQ2bin(files,
            trace_length,
            traces_count,
            memory_footprint,
            fast_range)
    end
    
    println("Loaded ", Base.summarysize(ret) / 1e9, " GB")
    cd(old_path)

    if axes
        print("Generating axes...")
        distance, time = _DAQ2_axes(ret, rep_rate, 
            sampling_freq, fast_range)

        if decimation > 0
            distance *= decimation+1
        elseif fast_avg > 0
            distance *= fast_avg+1
        elseif slow_avg > 0
            #time = time[1:slow_avg+1:end] .* (slow_avg+1)
            time *= slow_avg+1
        elseif fastslow_avg > 0
            distance, time = _DAQ2_axes(ret, rep_rate, 
                sampling_freq, fast_range; factor=fastslow_avg+1)
            time *= fastslow_avg+1
        end

        println("done")

        if size(ret, 1) != length(time) || size(ret, 2) != length(distance)
            println("Distance and time axes have different length to \
                returned sensing matrix!")
            #=
            println("Attempting adjustment...")
            # this is fucking horrific, write it properly
            if size(ret, 1) > length(time) && size(ret, 2) > length(distance)
                ret = ret[1:length(time),1:length(distance)]
            else size(ret, 1) < length(time) && size(ret, 2) < length(distance)
            =#
        end

        return ret, distance, time
    end

    return ret
end

function _DAQ2_header(first_file)
    fp = open(first_file, "r")
    head = split(readline(fp), ", ")
    head2 = split(readline(fp), ", ")
    close(fp)

    trace_length = parse(Int, 
        head2[findall(x->x=="trace_length", head)[1]])
    trig_freq = parse(Float64, 
        head2[findall(x->x=="trig_freq", head)[1]])
    # slow axis sample rate is the rep rate
    rep_rate = 1/trig_freq
    # fast axis sample rate is the sampling frequency
    sampling_freq = parse(Float64,
        head2[findall(x->x=="sampling_freq", head)[1]])
    traces_count = parse(Int,
        head2[findall(x->x=="traces_count", head)[1]])

    return trace_length, rep_rate, sampling_freq, traces_count
end

function _DAQ2_axes(ret,
    rep_rate,
    sampling_freq,
    fast_range;
    factor=0,
    c=2.99e8,
    n=1.46)

    spatial_resolution = c / n / sampling_freq / 2

    if fast_range[2] > 0
        fast_start = fast_range[1] * spatial_resolution
        fast_end = (fast_start + 
            spatial_resolution * size(ret, 2) + spatial_resolution/2)
    else
        fast_start = 0
        fast_end = size(ret,2) * spatial_resolution + spatial_resolution/2
    end

    if factor > 0
        fast_end = fast_start + size(ret,2) * spatial_resolution * factor
        distance = collect(fast_start:spatial_resolution*factor:
            fast_end-spatial_resolution*factor)
    else
        distance = collect(fast_start:spatial_resolution:
            (fast_end-spatial_resolution)*1)
    end

    time_start = 0
    time_end = size(ret, 1) * rep_rate
    time = collect(time_start:rep_rate:time_end-rep_rate/2)

    return distance, time
end

function _load_DAQ2bin(files,
    trace_length,
    traces_count,
    memory_footprint,
    fast_range)

    return _load_DAQ2bin_decimation(files, 1, trace_length,
        traces_count, memory_footprint, fast_range)
end

function _load_DAQ2bin_decimation(files, 
    decimation,
    trace_length,
    traces_count,
    memory_footprint,
    fast_range)

    # fast range use condition
    using_fast_range = false
    if fast_range[2] > 0
        using_fast_range = true
    end
    if ! using_fast_range
        decimated_trace_length = round(Int, trace_length / decimation)
        fast_start = 1
        fast_end = trace_length
    else
        decimated_trace_length = round(Int, 
            (fast_range[2] - fast_range[1]) / decimation)
        fast_start = fast_range[1] + 1
        fast_end = fast_range[2] 
    end

    ret = Matrix{Int16}(undef,
        traces_count * length(files),
        decimated_trace_length)
    processed_bytes = 0
    processed_traces = 0
    idx = 1

    for filename in files
        fp = open(filename, "r")
        header_bytes = 0
        header_bytes += length(readline(fp))
        header_bytes += length(readline(fp))
        processed_bytes += header_bytes

        for i in 1:traces_count
            buf = Vector{Int16}(undef, trace_length)
            read!(fp, buf)
            processed_traces += 1

            deci_buf = buf[fast_start:decimation:fast_end]
            processed_bytes += length(deci_buf) * 2
            
            if processed_traces % 1000 == 0
                print_loading_bar(processed_bytes, memory_footprint)
            end

            ret[idx:idx,:] = reshape(deci_buf, 1, decimated_trace_length)
            idx += 1
        end

        close(fp)
    end

    print_loading_bar(1,1)
    println()

    return ret
end

function _load_DAQ2bin_fast(files,
    factor,
    trace_length,
    traces_count,
    memory_footprint,
    fast_range)

    # fast range use condition
    using_fast_range = false
    if fast_range[2] > 0
        using_fast_range = true
    end
    if ! using_fast_range
        new_trace_length = round(Int, trace_length / factor)
        fast_start = 1
        fast_end = new_trace_length
    else
        new_trace_length = round(Int, 
            (fast_range[2] - fast_range[1]) / factor)
        fast_start = fast_range[1] + 1
        fast_end = fast_start + new_trace_length 
    end

    ret = Matrix{Int16}(undef,
        traces_count * length(files),
        new_trace_length)
    processed_bytes = 0
    processed_traces = 0
    idx = 1

    println("old trace length ", trace_length)
    println("new trace length ", new_trace_length)
    println("factor ", factor)
    println("traces to read ", traces_count)

    for filename in files
        fp = open(filename, "r")
        header_bytes = 0
        header_bytes += length(readline(fp))
        header_bytes += length(readline(fp))
        processed_bytes += header_bytes

        _, _, _, traces_count = _DAQ2_header(filename)
        for i in 1:traces_count
            buf = Vector{Int16}(undef, trace_length)
            read!(fp, buf)
            processed_traces += 1

            avg_buf = Vector{Int16}(undef, new_trace_length)
            ab_idx = 1
            for j in fast_start:factor:fast_end-factor
                view = buf[j:j+factor-1]
                if ab_idx > length(avg_buf)  # garbage bug fix
                    break
                end
                avg_buf[ab_idx] = round(Int, sum(view) / length(view))
                ab_idx += 1
            end
            processed_bytes += length(avg_buf) * 2
            
            if processed_traces % 1000 == 0
                print_loading_bar(processed_bytes, memory_footprint)
            end

            ret[idx:idx,:] = reshape(avg_buf, 1, new_trace_length)
            idx += 1
        end

        close(fp)
    end

    print_loading_bar(1,1)
    println()

    return ret
end

function _load_DAQ2bin_slow(files,
    factor,
    trace_length,
    traces_count,
    memory_footprint,
    fast_range)

    # fast range use condition
    using_fast_range = false
    if fast_range[2] > 0
        using_fast_range = true
    end
    if ! using_fast_range
        fast_start = 1
        fast_end = trace_length
        new_trace_length = trace_length
    else
        fast_start = fast_range[1] + 1
        fast_end = fast_range[2]
        new_trace_length = fast_range[2] - fast_range[1]
    end

    new_traces_count = div(traces_count, factor)

    ret = Matrix{Int16}(undef, 
        new_traces_count,
        new_trace_length)

    processed_bytes = 0
    processed_traces = 0
    idx = 1

    for filename in files
        _, _, _, traces_count = _DAQ2_header(filename)
        new_traces_count = div(traces_count, factor)

        fp = open(filename, "r")
        header_bytes = 0
        header_bytes += length(readline(fp))
        header_bytes += length(readline(fp))
        processed_bytes += header_bytes

        for i in 1:new_traces_count
            buf = Vector{Vector{Int16}}()
            for j in 1:factor
                r_buf = Vector{Int16}(undef, trace_length)
                read!(fp, r_buf)
                push!(buf, r_buf)
            end
            buf = Matrix{Int16}(transpose(stack(buf)))

            avg_buf = Vector{Int16}(undef, new_trace_length)
            ab_idx = 1
            for j in fast_start:fast_end
                view = buf[1:factor,j:j]
                if ab_idx > length(avg_buf)
                    break
                end
                avg_buf[ab_idx] = round(Int16, sum(view) / length(view))
                ab_idx += 1
            end
            
            ret[idx:idx,:] = reshape(avg_buf, 1, new_trace_length)
            processed_bytes += length(avg_buf) * 2
            idx += 1
            processed_traces += 1

            if processed_traces % 10 == 0
                print_loading_bar(processed_bytes, memory_footprint)
            end
        end

    end
    
    print_loading_bar(1,1)
    println()

    return ret
end

function _load_DAQ2bin_fastslow(files,
    factor,
    trace_length,
    traces_count,
    memory_footprint,
    fast_range)

    # fast range use condition
    using_fast_range = false
    if fast_range[2] > 0
        using_fast_range = true
    end
    if ! using_fast_range
        new_trace_length = round(Int, trace_length / factor)
        fast_start = 1
        fast_end = new_trace_length
    else
        new_trace_length = round(Int, 
            (fast_range[2] - fast_range[1]) / factor)
        fast_start = fast_range[1] + 1
        fast_end = fast_start + new_trace_length 
    end

    new_traces_count = round(Int, traces_count / factor)

    ret = Matrix{Int16}(undef, 
        new_traces_count,
        new_trace_length)

    processed_bytes = 0
    processed_traces = 0
    idx = 1

    for filename in files
        fp = open(filename, "r")
        header_bytes = 0
        header_bytes += length(readline(fp))
        header_bytes += length(readline(fp))
        processed_bytes += header_bytes

        _, _, _, traces_count = _DAQ2_header(filename)

        for i in 1:traces_count/factor
            buf = Vector{Vector{Int16}}()
            for j in 1:factor
                r_buf = Vector{Int16}(undef, trace_length)
                read!(fp, r_buf)
                push!(buf, r_buf)
            end
            buf = Matrix{Int16}(transpose(stack(buf)))

            avg_buf = Vector{Int16}(undef, new_trace_length)
            ab_idx = 1
            for j in fast_start:factor:trace_length-factor
                view = buf[1:factor,j:j+factor-1]
                if ab_idx > length(avg_buf)
                    break
                end
                avg_buf[ab_idx] = round(Int16, sum(view) / length(view))
                ab_idx += 1
            end
            
            ret[idx:idx,:] = reshape(avg_buf, 1, new_trace_length)
            processed_bytes += length(avg_buf) * 2
            idx += 1

            if processed_traces % 1000 == 0
                print_loading_bar(processed_bytes, memory_footprint)
            end

            processed_traces += 1

        end

    end
    
    print_loading_bar(1,1)
    println()

    return ret
end

function print_loading_bar(n1, n2;
    width=15)
    if n1 > n2
        pc = n2 / n1
    else
        pc = n1 / n2
    end

    buf = ""
    print("\r")
    for i in 1:width + 5
        buf *= " "
    end
    print(buf)
    print("\r")
    print(round(Int, pc * 100), "%")
    print(" [")
    equals_buf = ""
    for i in 1:round(Int, pc*width)
        equals_buf *= "="
    end
    for i in round(Int,pc*width)+1:width
        equals_buf *= " "
    end
    print(equals_buf)
    print("]")

    return pc
end


#####
# Network instrument control functions
#####

using Sockets

function TCP_connect(ip::AbstractString; port::Integer = 5025)
    return Sockets.connect(ip,port)
end

function TCP_disconnect(instrument::Sockets.TCPSocket)
    close(instrument)
end

function TCP_write(instrument::Sockets.TCPSocket,
    SCPI_command::AbstractString;
    termination_characters::AbstractArray = ['\n'])

    termination_string = join(termination_characters)
    SCPI_command = join([SCPI_command, termination_string])

    rc = write(instrument, SCPI_command)
    if rc != length(SCPI_command)
        return -1
    end

    return rc
end

function TCP_query(instrument::Sockets.TCPSocket,
    SCPI_command::AbstractString;
    termination_characters::AbstractArray = ['\n'])
    
    if twrite(instrument, SCPI_command; 
        termination_characters = termination_characters) < 0
        return -1
    end

    response = readuntil(instrument, join(termination_characters))

    return response
end

function TCP_identify(instrument::Sockets.TCPSocket;
    termination_characters::AbstractArray = ['\n'])

    return query(instrument, "*IDN?"; 
        termination_characters = termination_characters)
end


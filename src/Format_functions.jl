# format_functions.jl
#export format_intL, format_intR, format_up, format_down

# Function to format intervals for intL
function format_intL(interval)
    s, t = interval[2][2] - 1, interval[1][2] - 1
    return "<" * (s == 0 ? "̂0" : string(s)) * "' ," * (t == 0 ? "̂0" : string(t)) * "> "
end

# Function to format intervals for intR

function format_intR(interval)
    s, t = interval[1][1] - 1, interval[2][1] - 1
    return "<" * (interval[1][1] == interval[1][2] ? "̂1" : string(s) )* "," * (interval[2][1] == interval[2][2] ?  "̂1" : string(t)*"'") * "> "
end

# Function to format intervals for up
function format_up(interval)
    return "<" * string(interval[1] - 1) * "," * string(interval[2] - 1) * "> "
end

# Function to format intervals for down
function format_down(interval)
    return "<" * string(interval[1] - 1) * "', " * string(interval[2] - 1) * "'> "
end

function print_intervals(header, intervals, formatter)
    print(header)
    for interval in intervals
        print(formatter(interval))
    end
    println(" ")
end
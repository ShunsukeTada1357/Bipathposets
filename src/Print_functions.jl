# Print_functions.jl

# Function to print intervals for intL
function print_intL(interval)
    s, t = interval[2][2] - 1, interval[1][2] - 1
    return "<" * (s == 0 ? "̂0" : string(s)) * "' ," * (t == 0 ? "̂0" : string(t)) * "> "
end

# Function to print intervals for intR

function print_intR(interval)
    s, t = interval[1][1] - 1, interval[2][1] - 1
    return "<" * (interval[1][1] == interval[1][2] ? "̂1" : string(s) )* "," * (interval[2][1] == interval[2][2] ?  "̂1" : string(t)*"'") * "> "
end

# Function to print intervals for up
function print_up(interval)
    return "<" * string(interval[1] - 1) * "," * string(interval[2] - 1) * "> "
end

# Function to print intervals for down
function print_down(interval)
    return "<" * string(interval[1] - 1) * "', " * string(interval[2] - 1) * "'> "
end

function print_intervals(header, intervals, printter)
    print(header)
    for interval in intervals
        print(printter(interval))
    end
    println(" ")
end
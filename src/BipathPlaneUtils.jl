# BipathPlaneUtils.jl
# bipath plane の基本長さ
bipath_RR(n::Integer, m::Integer) = n + m + 2

# 各タイプの record / interval を bipath plane 上の点 [x,y] に写す
function bipath_left_point(up_interval::Vector{Int}, down_interval::Vector{Int}, n::Integer, m::Integer)
    RR = bipath_RR(n, m)
    a = n + m + 3 - down_interval[2]
    b = up_interval[2] + down_interval[2] - 2
    return [a, (a + b) % RR]
end

function bipath_right_point(up_interval::Vector{Int}, down_interval::Vector{Int}, n::Integer, m::Integer)
    a = up_interval[1] - 1
    b = (up_interval[2] - up_interval[1]) + (down_interval[2] - down_interval[1])
    return [a, a + b]
end

function bipath_up_point(interval::Vector{Int}, n::Integer, m::Integer)
    a = interval[1] - 1
    b = interval[2] - interval[1]
    return [a, a + b]
end

function bipath_down_point(interval::Vector{Int}, n::Integer, m::Integer)
    a = n + m + 3 - interval[2]
    b = interval[2] - interval[1]
    return [a, a + b]
end

function bipath_center_point(n::Integer, m::Integer)
    RR = bipath_RR(n, m)
    return [0, RR]
end

# 型ごとの interval データから点列をまとめて作る
function intervalstoplane(intL, intR, up, down, center, n::Integer, m::Integer)
    points = Vector{Vector{Int}}()

    for I in intL
        push!(points, bipath_left_point(I[1], I[2], n, m))
    end

    for I in intR
        push!(points, bipath_right_point(I[1], I[2], n, m))
    end

    for I in up
        push!(points, bipath_up_point(I, n, m))
    end

    for I in down
        push!(points, bipath_down_point(I, n, m))
    end

    for _ in center
        push!(points, bipath_center_point(n, m))
    end

    return points
end

# 軸ラベルも共通化したいならここに置ける
function bipath_axis_labels(n::Integer, m::Integer)
    lll = [string(i) for i in 1:n]

    push!(lll, "̂1")
    append!(lll, [string(m - i + 1) * "'" for i in 1:m])

    xaxis = vcat(lll,"̂0")   # 0,1,...,RR に対応
    yaxis = vcat("̂0", lll)

    return xaxis, yaxis
end
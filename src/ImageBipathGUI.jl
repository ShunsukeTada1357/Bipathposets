using GLMakie
GLMakie.activate!(inline=false)

function make_cell_birth_dict(FSC)
    return Dict(s[1] => s[2] for s in FSC[1])
end


function kind_color(kind)
    if kind == :left
        return :blue
    elseif kind == :right
        return :green
    elseif kind == :up
        return :orange
    elseif kind == :down
        return :purple
    elseif kind == :center
        return :red
    else
        return :black
    end
end

function representative_at_time_in_interval(basis, interval, cell_birth, i)
    if basis === nothing || interval === nothing
        return Tuple[]
    end

    if !(interval[1] <= i <= interval[2])
        return Tuple[]
    end

    return [cell for cell in basis if cell_birth[cell] <= i]
end

function disable_axis_interactions!(ax)
    for name in (:dragpan, :rectanglezoom, :scrollzoom)
        try
            GLMakie.deregister_interaction!(ax, name)
        catch
        end
    end
    return ax
end

function binary_image_for_makie(B)
    # いまの heatmap! 版と同じ向き
    return permutedims(Float32.(B), (2, 1))
end

function cubical_rep_primitives(representative)
    v_x = Float32[]
    v_y = Float32[]

    h_segments = Vector{Tuple{Vector{Float32}, Vector{Float32}}}()
    w_segments = Vector{Tuple{Vector{Float32}, Vector{Float32}}}()
    s_polys     = Vector{Tuple{Vector{Float32}, Vector{Float32}}}()

    a = 0.5f0

    for cell in representative
        tag = cell[1]

        if tag == :v
            _, x, y = cell
            push!(v_x, Float32(y) - a)
            push!(v_y, Float32(x) - a)

        elseif tag == :h
            _, x, y = cell
            push!(h_segments, ([Float32(y) - a, Float32(y) - a],
                               [Float32(x) - a, Float32(x) + a]))

        elseif tag == :w
            _, x, y = cell
            push!(w_segments, ([Float32(y) - a, Float32(y) + a],
                               [Float32(x) - a, Float32(x) - a]))

        elseif tag == :s
            _, x, y = cell
            push!(s_polys, ([Float32(y) - a, Float32(y) + 1 - a, Float32(y) + 1 - a, Float32(y) - a, Float32(y) - a],
                            [Float32(x) - a, Float32(x) - a, Float32(x) + 1 - a, Float32(x) + 1 - a, Float32(x) - a]))
        end
    end

    return (
        vertices = (v_x, v_y),
        h_segments = h_segments,
        w_segments = w_segments,
        s_polys = s_polys
    )
end

function init_binary_matrix_axis!(ax, B; title_str = "")
    m, n = size(B)
    img = binary_image_for_makie(B)

    #GLMakie.image!(ax, 1 .. n, 1 .. m, img; colormap = [:black, :white])
    GLMakie.image!(
        ax,
        0.5 .. (n + 0.5),
        0.5 .. (m + 0.5),
        img;
        colormap = [:black, :white],
        interpolate = false
    )
    GLMakie.xlims!(ax, 0.5, n + 0.5)
    GLMakie.ylims!(ax, m + 0.5, 0.5)
    ax.aspect = GLMakie.DataAspect()
    ax.title = title_str
    GLMakie.hidedecorations!(ax)
GLMakie.hidespines!(ax)
    disable_axis_interactions!(ax)

    return ax
end

function attach_rep_layers!(ax)
    vx = GLMakie.Observable(Float32[])
    vy = GLMakie.Observable(Float32[])

    GLMakie.scatter!(ax, vx, vy, color = :red, markersize = 8)

    h_lines = Tuple{GLMakie.Observable{Vector{Float32}}, GLMakie.Observable{Vector{Float32}}}[]
    w_lines = Tuple{GLMakie.Observable{Vector{Float32}}, GLMakie.Observable{Vector{Float32}}}[]
    s_lines = Tuple{GLMakie.Observable{Vector{Float32}}, GLMakie.Observable{Vector{Float32}}}[]

    return (
        vx = vx,
        vy = vy,
        h_lines = h_lines,
        w_lines = w_lines,
        s_lines = s_lines
    )
end

function ensure_line_layers!(ax, layer_store, needed_h::Int, needed_w::Int, needed_s::Int)
    while length(layer_store.h_lines) < needed_h
        ox = GLMakie.Observable(Float32[])
        oy = GLMakie.Observable(Float32[])
        GLMakie.lines!(ax, ox, oy, color = :red, linewidth = 3)
        push!(layer_store.h_lines, (ox, oy))
    end

    while length(layer_store.w_lines) < needed_w
        ox = GLMakie.Observable(Float32[])
        oy = GLMakie.Observable(Float32[])
        GLMakie.lines!(ax, ox, oy, color = :red, linewidth = 3)
        push!(layer_store.w_lines, (ox, oy))
    end

    while length(layer_store.s_lines) < needed_s
        ox = GLMakie.Observable(Float32[])
        oy = GLMakie.Observable(Float32[])
        GLMakie.lines!(ax, ox, oy, color = :blue, linewidth = 2)
        push!(layer_store.s_lines, (ox, oy))
    end
end

function update_rep_layers!(ax, layer_store, representative; title_str = "")
    prim = cubical_rep_primitives(representative)

    layer_store.vx[] = prim.vertices[1]
    layer_store.vy[] = prim.vertices[2]

    ensure_line_layers!(ax, layer_store,
        length(prim.h_segments), length(prim.w_segments), length(prim.s_polys))

    for i in eachindex(layer_store.h_lines)
        ox, oy = layer_store.h_lines[i]
        if i <= length(prim.h_segments)
            ox[] = prim.h_segments[i][1]
            oy[] = prim.h_segments[i][2]
        else
            ox[] = Float32[]
            oy[] = Float32[]
        end
    end

    for i in eachindex(layer_store.w_lines)
        ox, oy = layer_store.w_lines[i]
        if i <= length(prim.w_segments)
            ox[] = prim.w_segments[i][1]
            oy[] = prim.w_segments[i][2]
        else
            ox[] = Float32[]
            oy[] = Float32[]
        end
    end

    for i in eachindex(layer_store.s_lines)
        ox, oy = layer_store.s_lines[i]
        if i <= length(prim.s_polys)
            ox[] = prim.s_polys[i][1]
            oy[] = prim.s_polys[i][2]
        else
            ox[] = Float32[]
            oy[] = Float32[]
        end
    end

    ax.title = title_str
    return ax
end



"""
    interactive_bipath_viewer(records, B_list_a, B_list_b, FSCa, FSCb; dim=nothing, bipath_markersize=6)

Launch an interactive viewer for bipath records and their representatives.
"""
function interactive_bipath_viewer(
    records,
    B_list_a,
    B_list_b,
    FSCa,
    FSCb;
    dim = nothing,
    bipath_markersize = 6
)
    recs = isnothing(dim) ? records : [r for r in records if r.dim == dim]
    isempty(recs) && error("指定した次元の record がありません。")

    cell_birth_a = make_cell_birth_dict(FSCa)
    cell_birth_b = make_cell_birth_dict(FSCb)

    ptsx = [r.point[1] for r in recs]
    ptsy = [r.point[2] for r in recs]

    na = length(B_list_a)
    nb = length(B_list_b)
    ncols = max(na, nb)

    fig = GLMakie.Figure(size = (250 * max(ncols, 2), 700))
    n = FSCa[2] - 2
    m = FSCb[2] - 2
    RR = bipath_RR(n, m)
    xaxis, yaxis = bipath_axis_labels(n, m)

    ax_bp = GLMakie.Axis(
        fig[1, 1:ncols],
        title = isnothing(dim) ? "Bipath plane" : "Bipath plane (dim=$dim)",
        xlabel = "",
        ylabel = "",
        aspect = GLMakie.DataAspect(),
        xticks = (1:RR, xaxis),
        yticks = (0:RR-1, yaxis)
    )
    disable_axis_interactions!(ax_bp)
    GLMakie.xlims!(ax_bp, -1, RR + 1)
    GLMakie.ylims!(ax_bp, -1, RR + 1)
    GLMakie.vlines!(ax_bp, [n + 1, RR], color = :black, linestyle = :dash)
    GLMakie.hlines!(ax_bp, [0, n + 1], color = :black, linestyle = :dash)

    base_colors = [kind_color(r.kind) for r in recs]
    point_colors = GLMakie.Observable(copy(base_colors))

    bpplot = GLMakie.scatter!(
        ax_bp,
        ptsx,
        ptsy,
        color = point_colors,
        markersize = bipath_markersize
    )

    RR = (FSCa[2] - 2) + (FSCb[2] - 2) + 2
    GLMakie.xlims!(ax_bp, -1, RR + 1)
    GLMakie.ylims!(ax_bp, -1, RR + 1)

    axes_a = GLMakie.Axis[]
    for i in 1:na
        ax = GLMakie.Axis(fig[2, i], title = "A-side $i")
        disable_axis_interactions!(ax)
        push!(axes_a, ax)
    end

    axes_b = GLMakie.Axis[]
    for i in 1:nb
        ax = GLMakie.Axis(fig[3, i], title = "B-side $i")
        disable_axis_interactions!(ax)
        push!(axes_b, ax)
    end

    # 背景画像は最初に一度だけ描く
    for i in 1:na
        init_binary_matrix_axis!(axes_a[i], B_list_a[i]; title_str = "A-side $i")
    end
    for i in 1:nb
        init_binary_matrix_axis!(axes_b[i], B_list_b[i]; title_str = "B-side $i")
    end

    # representative 更新用のレイヤだけ作る
    layers_a = [attach_rep_layers!(axes_a[i]) for i in 1:na]
    layers_b = [attach_rep_layers!(axes_b[i]) for i in 1:nb]

    function draw_one_record!(rec)
        ax_bp.title = "point=$(rec.point), kind=$(rec.kind), up=$(rec.up_interval), down=$(rec.down_interval)"

        for i in 1:na
            rep_a_i = representative_at_time_in_interval(rec.up_basis, rec.up_interval, cell_birth_a, i)
            update_rep_layers!(
                axes_a[i], layers_a[i], rep_a_i;
                title_str = "A_$i | up=$(rec.up_interval)"
            )
        end

        for i in 1:nb
            rep_b_i = representative_at_time_in_interval(rec.down_basis, rec.down_interval, cell_birth_b, i)
            update_rep_layers!(
                axes_b[i], layers_b[i], rep_b_i;
                title_str = "B_$i | down=$(rec.down_interval)"
            )
        end
    end

    draw_one_record!(recs[1])

    on(GLMakie.events(fig).mousebutton, priority = 2) do event
        if event.button == GLMakie.Mouse.left && event.action == GLMakie.Mouse.press
            picked_plot, picked_idx = GLMakie.pick(ax_bp)

            if picked_plot === bpplot && picked_idx > 0 && picked_idx <= length(recs)
                rec = recs[picked_idx]

                newcols = copy(base_colors)
                newcols[picked_idx] = :black
                point_colors[] = newcols

                draw_one_record!(rec)
                return GLMakie.Consume()
            end
        end
        return
    end

    display(fig)
    return fig
end


function merge_representatives_at_time(records, side::Symbol, cell_birth, i)
    merged = Tuple[]
    seen = Set{Tuple}()

    for rec in records
        rep_i = if side == :A
            representative_at_time_in_interval(rec.up_basis, rec.up_interval, cell_birth, i)
        else
            representative_at_time_in_interval(rec.down_basis, rec.down_interval, cell_birth, i)
        end

        for cell in rep_i
            if !(cell in seen)
                push!(merged, cell)
                push!(seen, cell)
            end
        end
    end

    return merged
end


"""
    interactive_bipath_viewer_all(
        records,
        B_list_a,
        B_list_b,
        FSCa,
        FSCb;
        dim=nothing,
        bipath_markersize=10
    )

Launch an interactive viewer that displays the union of all representatives
corresponding to the clicked point in the bipath plane.
"""
function interactive_bipath_viewer_all(
    records,
    B_list_a,
    B_list_b,
    FSCa,
    FSCb;
    dim = nothing,
    bipath_markersize = 10
)
    recs = isnothing(dim) ? records : [r for r in records if r.dim == dim]
    isempty(recs) && error("指定した次元の record がありません。")

    cell_birth_a = make_cell_birth_dict(FSCa)
    cell_birth_b = make_cell_birth_dict(FSCb)

    ptsx = [r.point[1] for r in recs]
    ptsy = [r.point[2] for r in recs]

    na = length(B_list_a)
    nb = length(B_list_b)
    ncols = max(na, nb)

    fig = GLMakie.Figure(size = (250 * max(ncols, 2), 700))

    n = FSCa[2] - 2
    m = FSCb[2] - 2
    RR = bipath_RR(n, m)
    xaxis, yaxis = bipath_axis_labels(n, m)

    ax_bp = GLMakie.Axis(
        fig[1, 1:ncols],
        title = isnothing(dim) ? "Bipath plane" : "Bipath plane (dim=$dim)",
        xlabel = "",
        ylabel = "",
        xticks = (1:RR, xaxis),
        yticks = (0:RR-1, yaxis),
        aspect = GLMakie.DataAspect()
    )
    disable_axis_interactions!(ax_bp)
    GLMakie.xlims!(ax_bp, -1, RR + 1)
    GLMakie.ylims!(ax_bp, -1, RR + 1)
    GLMakie.vlines!(ax_bp, [n + 1, RR], color = :black, linestyle = :dash)
    GLMakie.hlines!(ax_bp, [0, n + 1], color = :black, linestyle = :dash)

    base_colors = [kind_color(r.kind) for r in recs]
    point_colors = GLMakie.Observable(copy(base_colors))

    bpplot = GLMakie.scatter!(
        ax_bp,
        ptsx,
        ptsy,
        color = point_colors,
        markersize = bipath_markersize
    )

    GLMakie.xlims!(ax_bp, -1, RR + 1)
    GLMakie.ylims!(ax_bp, -1, RR + 1)

    axes_a = GLMakie.Axis[]
    for i in 1:na
        ax = GLMakie.Axis(fig[2, i], title = "A-side $i")
        disable_axis_interactions!(ax)
        push!(axes_a, ax)
    end

    axes_b = GLMakie.Axis[]
    for i in 1:nb
        ax = GLMakie.Axis(fig[3, i], title = "B-side $i")
        disable_axis_interactions!(ax)
        push!(axes_b, ax)
    end

    for i in 1:na
        init_binary_matrix_axis!(axes_a[i], B_list_a[i]; title_str = "A-side $i")
    end
    for i in 1:nb
        init_binary_matrix_axis!(axes_b[i], B_list_b[i]; title_str = "B-side $i")
    end

    layers_a = [attach_rep_layers!(axes_a[i]) for i in 1:na]
    layers_b = [attach_rep_layers!(axes_b[i]) for i in 1:nb]

    function draw_record_group!(group)
        #p = group[1].point
        #ax_bp.title = "point=$([x,y]), multiplicity=$(length(group))"
        ax_bp.title = "multiplicity=$(length(group))"
        for i in 1:na
            rep_a_i = merge_representatives_at_time(group, :A, cell_birth_a, i)
            update_rep_layers!(
                axes_a[i],
                layers_a[i],
                rep_a_i;
                title_str = "A_$i "
            )
        end

        for i in 1:nb
            rep_b_i = merge_representatives_at_time(group, :B, cell_birth_b, i)
            update_rep_layers!(
                axes_b[i],
                layers_b[i],
                rep_b_i;
                title_str = "B_$i "
            )
        end
    end

    # 初期表示は最初の point に対応する全 record
    p0 = recs[1].point
    group0 = [r for r in recs if r.point == p0]
    draw_record_group!(group0)

    on(GLMakie.events(fig).mousebutton, priority = 2) do event
        if event.button == GLMakie.Mouse.left && event.action == GLMakie.Mouse.press
            picked_plot, picked_idx = GLMakie.pick(ax_bp)

            if picked_plot === bpplot && picked_idx > 0 && picked_idx <= length(recs)
                p = recs[picked_idx].point
                group = [r for r in recs if r.point == p]

                newcols = copy(base_colors)
                for j in eachindex(recs)
                    if recs[j].point == p
                        newcols[j] = :black
                    end
                end
                point_colors[] = newcols

                draw_record_group!(group)
                return GLMakie.Consume()
            end
        end
        return
    end

    display(fig)
    return fig
end
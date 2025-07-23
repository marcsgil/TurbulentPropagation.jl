using StructuredLight, CairoMakie, LinearAlgebra, Tullio, TurbulentPropagation

function get_order_indices(order)
    ((order - n, n) for n ∈ 0:order)
end

function get_indices(order)
    Iterators.flatten((get_order_indices(n) for n ∈ 0:order))
end

rs = LinRange(-12, 12, 256)

r₀ = 0.01
l₀ = 0.001
L₀ = 100
P = 3
spectrum(κ) = mvK_spectrum(κ, r₀, l₀, L₀)

ns_ft_sh = ft_sh_phase_screen(spectrum, rs, P)
normalize!(ns_ft_sh, Inf)
ns_ft_sh .*= π
##
order = 15
ψ₀ = stack(hg(rs, rs; m, n) for (m, n) ∈ get_indices(order))
ψt = stack(hg(rs, rs, 1; m, n) for (m, n) ∈ get_indices(order))

@. ψ₀ *= cis(ns_ft_sh)

ψ = free_propagation(ψ₀, rs, rs, 1, 1)

heatmap(abs2.(ψ[:, :, end]))
T = similar(ψ₀, size(ψ, 3), size(ψ, 3))

Threads.@threads for j ∈ axes(ψ₀, 3)
    for k ∈ axes(ψ₀, 3)
        T[j, k] = overlap(view(ψ, :, :, j), view(ψt, :, :, k), rs, rs)
    end
end
##
fig = Figure(; resolution=(600, 600))
ax = Axis(fig[1, 1])
ax.yreversed = true
ax.aspect = 1
heatmap!(ax, abs.(T), colormap=:jet, yreversed=true)
fig
##
F = eigen(T)
abs2.(F.vectors[:, 1]) |> maximum

lines(abs2.(F.vectors[:, 1]), color=:blue)


metric(x, L) = norm(view(x, 1:L))
vectors = sort(eachslice(F.vectors, dims=2), by=x -> metric(x, 10), rev=true)

lines(abs2.(vectors[1]), color=:blue)

abs2.(vectors[2])
##
@tullio ϕ₀[i, j] := vectors[1][k] * ψ₀[i, j, k]
##
ϕ = free_propagation(ϕ₀, rs, rs, 1)

heatmap(abs2.(ϕ))
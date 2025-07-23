using TurbulentPropagation, CairoMakie
##
r₀ = 0.01
l₀ = 0.001
L₀ = 100
P = 3
spectrum(κ) = mvK_spectrum(κ, r₀, l₀, L₀)
##
rs = LinRange(-1, 1, 256)
Rs = rs[length(rs)÷2+1:end]
ns_ft = ft_phase_screen(spectrum, rs, 512)
ns_ft_sh = ft_sh_phase_screen(spectrum, rs, P, 512)
##
heatmap(rs, rs, ns_ft_sh[:, :, end] |> Array, color=:grays)
##
Ds_ft = structure_function(ns_ft |> Array)
Ds_ft_sh = structure_function(ns_ft_sh |> Array)

plot(Rs, Ds_ft)
plot!(Rs, Ds_ft_sh)
plot!(Rs, mvK_structure_function.(Rs, r₀, l₀, L₀))
##
λ = 633e-9
k = 2π / λ

L = 1
α = 1.23 / 0.423
##
N = 128
rs = LinRange(-5e-2, 5e-2, N)
_γ, _ψ = turbulence_eigenmodes(rs, 1, 0, k)
##
ψ = reshape(_ψ[1], N, N)

heatmap(abs2.(ψ), aspect_ratio=1, size=(510, 500))
##
_γ[1] |> abs2
interp = cubic_spline_interpolation((rs, rs), ψ)
##
vizualize(ψ, ratio=4)
evolved = free_propagation(ψ, rs, rs, LinRange(0, 10, 32), k=k)

FreeParaxialPropagation.animate(evolved, ratio=4)
##
O = turbulence_eigenmodes(rs, 1, 0, k)

O[2, 3]
O[3, 2]
##
@btime turbulence_eigenmodes(rs, 1, 0, k)
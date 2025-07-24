using TurbulentPropagation, StatsBase, FFTW, CairoMakie, Bessels, SpecialFunctions

r₀ = 0.1
l₀ = 0.01
L₀ = 10

function mvK_spectrum(κ, param)
    r₀, L₀, l₀ = param
    iszero(sum(abs2, κ)) && return zero(sum(abs2, κ))
    0.49 * r₀^(-5 / 3) * exp(-sum(abs2, κ) * (l₀ / 5.92)^2) / (sum(abs2, κ) + (2π / L₀)^2)^(11 / 6)
end


function mvK_structure_function_approx(r, param)
    r₀, L₀, l₀ = param
    7.75 * r₀^(-5 / 3) * r^2 * (inv(l₀ + 2.03(r)^2)^(1 / 6) - 0.72 * (2π / L₀)^(1 / 3))
end

function mvK_structure_function(r, param)
    r₀, L₀ = param
    R = √sum(abs2, r)
    6.16 * r₀^(-5 / 3) * (3 / 5 * (L₀ / 2π)^(5 / 3) - (R * L₀ / 4π)^(5 / 6) * Bessels.besselk(5 / 6, 2π * r / L₀) / gamma(11 / 6))
end

param = (r₀, L₀, l₀)

N = 256
L = 2
δ = L / N
dest = randn(ComplexF64, N, N, 200)

plan = plan_bfft!(dest, (1,2))

dest_D = copy(dest[:, :, 1])
plan_D = plan_fft!(dest_D)
TurbulentPropagation.expected_structure_function!(dest_D, δ, δ, param, plan_D; spectrum=mvK_spectrum)

D_expected = real(dest_D[begin, begin:N÷2])
rs = (0:length(D_expected)-1) * δ
D_theoretical = [mvK_structure_function(r, param) for r in rs]

fourier_phase_screen!(dest, δ, δ, param, plan; spectrum=mvK_spectrum)
D_statistic = structure_from_statistics(real(dest))

with_theme(theme_latexfonts()) do
    fig = Figure(fontsize=24)
    ax = Axis(fig[1, 1], xlabel="r", ylabel="Structure Function")
    lines!(ax, rs, D_expected, linewidth=4)
    lines!(ax, rs, D_theoretical; linestyle=:dash, linewidth=4)
    lines!(ax, rs, 2 * D_statistic; linestyle=:dot, linewidth=4)
    fig
end
##
acorr = TurbulentPropagation.autocorrelation(real(dest))
heatmap((acorr))
##
D = structure_from_statistics(real(dest))
heatmap(D)
##
Ds = structure_from_statistics(real(dest))
Rs = (0:length(Ds)-1) * d

Ds_an = [structure_from_psd(R, hill_andrews_spectrum, param) for R in Rs]


with_theme(theme_latexfonts()) do
    fig = Figure(fontsize=20)
    ax = Axis(fig[1, 1], xlabel="r", ylabel="Structure Function")
    lines!(ax, Rs ./ r₀, Ds; label="Numerical", linewidth=3)
    lines!(ax, Rs ./ r₀, Ds_an; label="Analytical", linewidth=3, linestyle=:dash)
    axislegend(ax, fontsize=20, position=:lt)
    fig
end
##
std(real(dest))

rs = StepRangeLen(-(N ÷ 2) * d, d, N)
u = [complex(exp(-(x^2 + y^2) / 1)) for x in rs, y in rs]

plan = plan_fft!(u)
iplan = plan_ifft!(u)

cis.(real(dest))

u .*= cis.(real(dest))
heatmap(abs2.(u))

TurbulentPropagation.angular_spectrum_propagation!(u, d, d, 100, 1e-3, 2, plan, iplan)

heatmap(abs2.(u))
##
cov = mean(autocov(real(dest)), dims=2)[:, 1]
structure = 2 * (first(cov) .- cov)

rs = (0:length(structure)-1) * d
structure_hill_andrews = mvK_structure_function.(rs, r₀, l₀, L₀)

print(rs ./ r₀)

with_theme(theme_latexfonts()) do
    fig = Figure(fontsize=20)
    ax = Axis(fig[1, 1], xlabel="r", ylabel="Structure Function")
    lines!(ax, rs ./ r₀, structure; label="Structure Function", linewidth=3)
    lines!(ax, rs ./ r₀, structure_hill_andrews; label="Hill-Andrews", linewidth=3, linestyle=:dash)
    fig
end

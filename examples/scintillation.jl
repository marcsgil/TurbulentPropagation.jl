using TurbulentPropagation, StructuredLight, CairoMakie, ProgressMeter

N = 256  # Number of points in each dimension
Δx = Δy = 1.5625e-3
Lx = N * Δx  # Size of the grid
Ly = N * Δy  # Size of the grid
λ = 1e-6  # Wavelength of the wave
w = 2 * 0.0567
z = 500

rs = StepRangeLen(-Δx * (N ÷ 2), Δx, N)

nsteps = 10

Cₙ² = 1e-14
r₀ = fried_parameter(z / nsteps, Cₙ², λ)
L₀ = 1e10
l₀ = 1e-10
param = (r₀, L₀, l₀)

nsamples = 1024

u = stack(hg(rs, rs; w) for _ in 1:nsamples)

buffers = TurbulentPropagationBuffers(u)
turbulent_propagation!(buffers, hill_andrews_spectrum, param, Δx, Δx, z, λ; nsteps)

visualize(abs2.(buffers.u[:, :, 1]))
##
I = dropdims(mean(abs2, u, dims=3), dims=3)
I² = dropdims(mean(x -> abs2(x)^2, u, dims=3), dims=3)
scintillation_index = I² ./ I .^ 2 .- 1

line_scintillation = (scintillation_index[N÷2+1:end, N÷2] + scintillation_index[N÷2, N÷2+1:end]
                      + reverse(scintillation_index[1:N÷2, N÷2] + scintillation_index[N÷2, 1:N÷2])) / 4
##
@time begin
    repetitions = 10
end
##
σ² = rytov_variance(z, Cₙ², λ)

with_theme(theme_latexfonts()) do
    fig = Figure()
    ax = Axis(fig[1, 1])
    xlims!(ax, 0, 0.2)
    ylims!(ax, 0, 2σ²)
    lines!(ax, rs[N÷2+1:end], line_scintillation)
    hlines!(ax, σ²)
    fig
end
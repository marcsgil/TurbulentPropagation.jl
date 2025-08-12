using CairoMakie, TurbulentPropagation, FFTW

N = 256
L = 2
Δ = L / N
dest = Array{ComplexF64}(undef, N, N, 1000)
buffer = similar(dest, N, N)

@benchmark randn!($dest)

plan = plan_fft!(dest, (1, 2))
iplan = inv(plan)

spectrum = von_karman_spectrum

r₀ = 0.1
L₀ = 10
param = (r₀, L₀)

TurbulentPropagation.sample_fourier_phase_screen!(dest, buffer, spectrum, param, N, N, Δ, Δ, plan, iplan)


@benchmark TurbulentPropagation.sample_fourier_phase_screen!($dest, $buffer, $spectrum, $param, $N, $N, $Δ, $Δ, $plan, $iplan)


#heatmap(real(dest[:, :, 3]))
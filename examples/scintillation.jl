using TurbulentPropagation, StructuredLight, CairoMakie

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

nsamples = 128

u = stack(hg(rs, rs; w) for _ in 1:nsamples)

buffers = TurbulentPropagationBuffers(u)
turbulent_propagation!(buffers, hill_andrews_spectrum, param, Δx, Δx, z, λ; nsteps)

visualize(abs2.(buffers.u[:, :, 1]))
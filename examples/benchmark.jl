using TurbulentPropagation, StructuredLight

N = 256  # Number of points in each dimension
Δx = Δy = 1.5625e-3
Lx = N * Δx  # Size of the grid
Ly = N * Δy  # Size of the grid
λ = 1e-6  # Wavelength of the wave
w = 2 * 0.0567
z = 500

rs = StepRangeLen(-Δx * (N ÷ 2), Δx, N)

nsteps = 1

Cₙ² = 1e-14
r₀ = fried_parameter(z / nsteps, Cₙ², λ)
L₀ = 1e10
l₀ = 1e-10
param = (r₀, L₀, l₀)

nsamples = 128

u = stack(hg(rs, rs; w) for _ in 1:nsamples)

buffers = TurbulentPropagationBuffers(u)
method = AngularSpectrum()
##
@benchmark turbulent_propagation!($buffers, $hill_andrews_spectrum, $param, $Δx, $Δy, $z, $λ; nsteps=1)
##
@benchmark free_propagation!($buffers, $Δx, $Δy, $z, $λ, $method)
@benchmark sample_phase_screen!($buffers, $hill_andrews_spectrum, $param, $Δx, $Δy)
##
@benchmark TurbulentPropagation.sample_fourier_phase_screen!($buffers, $hill_andrews_spectrum, $param, $Δx, $Δy)
@benchmark TurbulentPropagation.sample_subharmonic_phase_screen!($buffers, $hill_andrews_spectrum, $param, $Δx, $Δy)
##
@benchmark TurbulentPropagation.hermitian_randn!($buffers.subharmonic_randoms)
module TurbulentPropagation
using FFTW, KernelAbstractions, Random, Integrals, StatsBase

include("free_propagation.jl")

include("spectra.jl")
export modified_von_karman_spectrum, von_karman_spectrum, kolmogorov_spectrum, hill_andrews_spectrum

include("phase_screens.jl")
export fourier_phase_screen!

include("screen_diagnostics.jl")
export structure_from_psd, structure_from_statistics

end
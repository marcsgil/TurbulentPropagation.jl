module TurbulentPropagation
using FFTW, KernelAbstractions

include("free_propagation.jl")

include("spectra.jl")
export modified_von_karman_spectrum, von_karman_spectrum, kolmogorov_spectrum, hill_andrews_spectrum
end
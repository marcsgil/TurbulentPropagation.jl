module TurbulentPropagation

using KernelAbstractions, FFTW, Random

include("utils.jl")
export TurbulentPropagationBuffers

include("free_propagation.jl")
export free_propagation!, AngularSpectrum

include("phase_masks.jl")

include("spectra.jl")
export kolmogorov_spectrum, von_karman_spectrum, modified_von_karman_spectrum, hill_andrews_spectrum

end

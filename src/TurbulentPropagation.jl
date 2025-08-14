module TurbulentPropagation

using KernelAbstractions, FFTW, Random

include("utils.jl")
export TurbulentPropagationBuffers

include("free_propagation.jl")
export free_propagation!, AngularSpectrum

include("phase_masks.jl")
export sample_phase_screen!

include("spectra.jl")
export kolmogorov_spectrum, von_karman_spectrum, modified_von_karman_spectrum, hill_andrews_spectrum

include("turbulent_propagation.jl")
export turbulent_propagation!

include("atmospheric_parameters.jl")
export fried_parameter, rytov_variance

end

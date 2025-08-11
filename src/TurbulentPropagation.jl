module TurbulentPropagation

using KernelAbstractions, FFTW

include("free_propagation.jl")
export free_propagation!

include("phase_masks.jl")

end

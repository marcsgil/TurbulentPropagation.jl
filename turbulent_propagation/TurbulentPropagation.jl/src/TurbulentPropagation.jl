module TurbulentPropagation

using Tullio,FastTransforms,Statistics,CUDA,CUDA.CUFFT,KrylovKit,LinearAlgebra
CUDA.allowscalar(false)

#include("dft_utils.jl")

include("characterization.jl")
export structure_function

include("exact_spectra.jl")
export mvK_spectrum,mvK_structure_function

include("screens.jl")
export ft_phase_screen,sh_phase_screen,ft_sh_phase_screen,sh_field

include("deterministic_bs.jl")
export deterministic_Bs_ft,deterministic_Ds_ft,deterministic_Ds_ft_sh,deterministic_Bs_sh

include("find_eigs.jl")
export turbulence_eigenmodes

end

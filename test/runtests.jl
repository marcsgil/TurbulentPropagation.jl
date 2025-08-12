using Test, TurbulentPropagation, StructuredLight, FFTW, Random, Bessels, SpecialFunctions, Statistics

Random.seed!(1234)

include("free_propagation_tests.jl")
include("structure_function_tests.jl")
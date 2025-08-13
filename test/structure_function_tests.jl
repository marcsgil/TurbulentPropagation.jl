function von_karman_structure_function(rs, param)
    r₀, L₀ = param
    r = √sum(abs2, rs)
    iszero(r) && return zero(r)
    6.16 / r₀^(5 // 3) * (3 / 5 * (L₀ / 2π)^(5 // 3) - (r * L₀ / 4π)^(5 // 6) * Bessels.besselk(5 / 6, 2π * r / L₀) / gamma(11 / 6))
end

function statistical_structure_function(Φ::AbstractVector)
    N = length(Φ)
    Φ₀ = first(Φ)
    [(Φ[m] - Φ₀)^2 for m in 1:N÷2]
end

function statistical_structure_function(Φ::Array{T,3}) where {T}
    dims = ntuple(identity, ndims(Φ) - 1) .+ 1
    slices = eachslice(Φ; dims)
    mean(statistical_structure_function, slices)
end

@testset "Structure Function Tests" begin
    N = 256
    L = 2
    Δ = L / N
    rs = (0:N÷2-1) * Δ

    r₀ = 0.1
    L₀ = 10
    param = (r₀, L₀)
    spectrum = von_karman_spectrum

    nsamples = 1000

    u = Array{ComplexF64}(undef, N, N, nsamples)
    buffers = TurbulentPropagationBuffers(u)

    TurbulentPropagation.sample_fourier_phase_screen!(buffers, spectrum, param, Δ, Δ)
    TurbulentPropagation.sample_subharmonic_phase_screen!(buffers, spectrum, param, Δ, Δ)

    D_anl = [von_karman_structure_function(r, param) for r in rs]
    D_num = statistical_structure_function(real(buffers.phase))

    @test isapprox(D_num[3:end], D_anl[3:end]; rtol=3e-2)
end
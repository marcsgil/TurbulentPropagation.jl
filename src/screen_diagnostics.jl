"""
    structure_from_psd(R, psd, param)

Compute the structure function at a separation `R` given a power spectral density `psd` with parameters `param`.

`psd` should have a signature `psd(k, param)` where `k` is the wavenumber.

The formula used is:

``
    D(R) = 8π ∫₀^∞ k² psd(k) * (1 - sinc(kR)) dk .
``
"""
function structure_from_psd(R, psd, param; reltol=1e-6, abstol=1e-6)
    f(u, param) = u^2 * psd(u / R, param) * (1 - sinc(u))
    domain = (0, Inf)
    prob = IntegralProblem(f, domain, param)
    sol = solve(prob, QuadGKJL(); reltol, abstol)
    8π * sol.u / R^3
end

function structure_from_statistics(data::AbstractVector)
    N = length(data)
    (data[begin:N÷2] .- first(data)) .^ 2
end

function structure_from_statistics(data::AbstractArray{T,N}) where {T,N}
    other_dims = ntuple(identity, N)[2:end]
    mean(structure_from_statistics, eachslice(data, dims=other_dims))
end
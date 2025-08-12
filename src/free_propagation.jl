@kernel function phase_kernel!(dest, α, qx, qy)
    J = @index(Global, NTuple)
    dest[J...] *= cis(α * (qx[J[1]]^2 + qy[J[2]]^2))
end

abstract type FreePropagationMethod end

struct AngularSpectrum <: FreePropagationMethod end

function free_propagation!(u, Δx, Δy, z, λ, plan, iplan, method::FreePropagationMethod)
    qx = fftfreq(size(u, 1), 1 / Δx)
    qy = fftfreq(size(u, 2), 1 / Δy)

    kernel! = phase_kernel!(get_backend(u))

    _free_propagation!(u, qx, qy, z, λ, plan, iplan, method, kernel!, size(u))
end

function _free_propagation!(u, qx, qy, z, λ, plan, iplan, method::AngularSpectrum, kernel!, ndrange)
    plan * u
    kernel!(u, -z * λ * π, qx, qy; ndrange)
    iplan * u
    nothing
end
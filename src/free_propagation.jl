@kernel function phase_kernel!(dest, α, qx, qy)
    J = @index(Global, NTuple)
    dest[J...] *= cis(α * (qx[J[1]]^2 + qy[J[2]]^2))
end

abstract type FreePropagationMethod end

struct AngularSpectrum <: FreePropagationMethod end

function free_propagation!(buffers::TurbulentPropagationBuffers, Δx, Δy, z, λ, method::FreePropagationMethod)
    qx = fftfreq(size(buffers.u, 1), 1 / Δx)
    qy = fftfreq(size(buffers.u, 2), 1 / Δy)

    kernel! = phase_kernel!(get_backend(buffers.u))

    _free_propagation!(buffers.u, qx, qy, z, λ, buffers.fft_plan, buffers.ifft_plan, method, kernel!, size(buffers.u))
end

function _free_propagation!(u, qx, qy, z, λ, plan, iplan, ::AngularSpectrum, kernel!, ndrange)
    plan * u
    kernel!(u, -z * λ * π, qx, qy; ndrange)
    iplan * u
    nothing
end
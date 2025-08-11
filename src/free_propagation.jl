@kernel function phase_kernel!(dest, α, qx, qy)
    J = @index(Global, NTuple)
    dest[J...] *= cis(α * (qx[J[1]]^2 + qy[J[2]]^2))
end

function free_propagation!(u, Δx, Δy, z, λ, plan, iplan)
    qx = fftfreq(size(u, 1), 1 / Δx)
    qy = fftfreq(size(u, 2), 1 / Δy)

    kernel! = phase_kernel!(get_backend(u))

    plan * u
    kernel!(u, - z * λ * π, qx, qy; ndrange=size(u))
    iplan * u
    nothing
end
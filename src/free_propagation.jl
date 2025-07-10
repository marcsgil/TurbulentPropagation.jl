@kernel function quadratic_phase_kernel!(field, qx, qy, α)
    J = @index(Global, NTuple)
    field[J...] *= cis(α * (qx[J[1]]^2 + qy[J[2]]^2))
end

function angular_spectrum_propagation!(u, Δx, Δy, z, λ, magnification, plan, iplan)
    ndrange = size(u)
    Nx, Ny = ndrange
    xs = StepRangeLen(-Δx * (Nx ÷ 2), Δx, Nx)
    ys = StepRangeLen(-Δy * (Ny ÷ 2), Δy, Ny)
    qx = fftfreq(Nx, 2π / Δx)
    qy = fftfreq(Ny, 2π / Δy)

    kernel! = quadratic_phase_kernel!(get_backend(u))
    
    u ./= magnification
    kernel!(u, xs, ys, π * (1 - magnification) / λ / z; ndrange)
    plan * u
    kernel!(u, qx, qy, -z * λ / 4π / magnification; ndrange)
    iplan * u
    kernel!(u, xs, ys, π * (magnification - 1) * magnification / λ / z; ndrange)
end
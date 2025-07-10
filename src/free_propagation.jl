@kernel function quadratic_phase_kernel!(field, qx, qy, α)
    J = @index(Global, NTuple)
    field[J...] *= cis(α * (qx[J[1]]^2 + qy[J[2]]^2))
end

"""
    angular_spectrum_propagation!(u, dx, dy, dz, λ, magnification, plan, iplan)

Propagate a field `u` using the angular spectrum method. Overwrites the result in `u`.

`dx` and `dy` are the sampling intervals in the x and y directions, respectively.
`dz` is the propagation distance.
`λ` is the wavelength of the light.
`magnification` is the magnification factor. Should be different from 1
`plan` and `iplan` are the forward and inverse FFT plans, respectively.
"""
function angular_spectrum_propagation!(u, dx, dy, dz, λ, magnification, plan, iplan)
    ndrange = size(u)
    Nx, Ny = ndrange
    xs = StepRangeLen(-dx * (Nx ÷ 2), dx, Nx)
    ys = StepRangeLen(-dy * (Ny ÷ 2), dy, Ny)
    qx = fftfreq(Nx, 2π / dx)
    qy = fftfreq(Ny, 2π / dy)

    kernel! = quadratic_phase_kernel!(get_backend(u))
    
    u ./= magnification
    kernel!(u, xs, ys, π * (1 - magnification) / λ / dz; ndrange)
    plan * u
    kernel!(u, qx, qy, -dz * λ / 4π / magnification; ndrange)
    iplan * u
    kernel!(u, xs, ys, π * (magnification - 1) * magnification / λ / dz; ndrange)
end
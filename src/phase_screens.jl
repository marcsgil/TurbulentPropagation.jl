@kernel function spectrum_kernel!(dest, qxs, qys, dqx, dqy, Cₙ², l₀, L₀, spectrum)
    J = @index(Global, NTuple)
    dest[J...] *= √(spectrum((qxs[J[1]], qys[J[2]]), Cₙ², l₀, L₀) * dqx * dqy)
end

function fourier_phase_screen!(dest, dx, dy, Cₙ², l₀, L₀, plan; spectrum=hill_andrews_spectrum)
    ndrange = size(dest)
    Nx, Ny = ndrange
    qxs = fftfreq(Nx, 2π / dx)
    qys = fftfreq(Ny, 2π / dy)

    randn!(dest)
    kernel! = spectrum_kernel!(get_backend(dest))
    kernel!(dest, qxs, qys, step(qxs), step(qys), Cₙ², l₀, L₀, spectrum; ndrange)
    plan * dest
end
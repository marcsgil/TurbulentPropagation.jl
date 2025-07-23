@kernel function spectrum_kernel!(dest, qxs, qys, dqx, dqy, param, spectrum)
    J = @index(Global, NTuple)
    dest[J...] *= √(spectrum((qxs[J[1]], qys[J[2]]), param) * dqx * dqy)
end

function fourier_phase_screen!(dest, dx, dy, param, plan; spectrum=hill_andrews_spectrum)
    ndrange = size(dest)
    Nx, Ny = ndrange
    qxs = fftfreq(Nx, 2π / dx)
    qys = fftfreq(Ny, 2π / dy)

    randn!(dest)
    kernel! = spectrum_kernel!(get_backend(dest))
    kernel!(dest, qxs, qys, step(qxs), step(qys), param, spectrum; ndrange)
    plan * dest
end

function expected_power_spectrum!(dest, dx, dy, param, plan; spectrum=hill_andrews_spectrum)
    ndrange = size(dest)
    Nx, Ny = ndrange
    qxs = fftfreq(Nx, 2π / dx)
    qys = fftfreq(Ny, 2π / dy)

    fill!(dest, one(eltype(dest)))
    kernel! = spectrum_kernel!(get_backend(dest))
    kernel!(dest, qxs, qys, step(qxs), step(qys), param, spectrum; ndrange)
    kernel!(dest, qxs, qys, step(qxs), step(qys), param, spectrum; ndrange)
    plan * dest
end

function expected_structure_function!(dest, args...; kwargs...)
    expected_power_spectrum!(dest, args...; kwargs...)
    dest .= 2 .* (first(dest) .- dest)
    #nothing
end
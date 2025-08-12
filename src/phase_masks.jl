#= @kernel function phase_screen_kernel!(dest, spectrum, param, ΔA, qx, qy)
    J = @index(Global, NTuple)
    dest[J...] *= √(spectrum((qx[J[1]], qy[J[2]]), param) * ΔA)
end =#

@kernel function evaluate_spectrum_kernel!(buffer, spectrum, param, ΔA, qx, qy)
    i, j = @index(Global, NTuple)
    buffer[i, j] = √(spectrum((qx[i], qy[j]), param) * ΔA)
end

function sample_fourier_phase_screen!(dest, buffer, spectrum, param, Nx, Ny, Δx, Δy, plan, iplan)
    qx = fftfreq(Nx, 2π / Δx)
    qy = fftfreq(Ny, 2π / Δy)
    ΔA = step(qx) * step(qy)

    kernel! = evaluate_spectrum_kernel!(get_backend(buffer))
    kernel!(buffer, spectrum, param, ΔA, qx, qy; ndrange=size(buffer))

    randn!(dest)
    #map!(x -> real(x) + imag(x), dest, dest)
    #plan * dest
    dest .*= buffer
    #iplan * dest
end
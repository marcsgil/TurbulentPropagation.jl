function idx_weight(I, J, Nx, Ny)
    i = min(I - 1, Nx - I + 1)
    j = min(J - 1, Ny - J + 1)
    (i > 1 || j > 1) && return 1 // 1
    !(iszero(i) && iszero(j)) // 2^i // 2^j
end

@kernel function evaluate_spectrum_kernel!(buffer, spectrum, param, ΔA, qx, qy)
    i, j = @index(Global, NTuple)
    buffer[i, j] = √(spectrum((qx[i], qy[j]), param) * ΔA * idx_weight(i, j, size(buffer, 1), size(buffer, 2)))
end

function hermitian_randn!(dest, plan)
    randn!(dest)
    map!(x -> real(x) + imag(x), dest, dest)
    plan * dest
    nothing
end

function hermitian_randn!(dest::AbstractArray{T,4}) where {T<:Complex}
    N = size(dest, 3)
    @assert N == size(dest, 4)
    randn!(view(dest, :, :, 1:N÷2, :))
    dest[:, :, N÷2+1:N, :] = view(dest, :, :, 1:N÷2, :)
    reverse!(view(dest, :, :, N÷2+1:N, :), dims=(3, 4))
    map!(conj, view(dest, :, :, N÷2+1:N, :), view(dest, :, :, N÷2+1:N, :))
    nothing
end

function sample_fourier_phase_screen!(dest, buffer, spectrum, param, Δx, Δy, plan, iplan)
    Nx, Ny = size(dest)
    qx = fftfreq(Nx, 2π / Δx)
    qy = fftfreq(Ny, 2π / Δy)
    ΔA = step(qx) * step(qy)

    kernel! = evaluate_spectrum_kernel!(get_backend(buffer))
    kernel!(buffer, spectrum, param, ΔA, qx, qy; ndrange=size(buffer))

    hermitian_randn!(dest, plan)
    @. dest *= buffer / √(Nx * Ny)
    iplan * dest
    nothing
end

@kernel function subharmonic_kernel!(dest, random_numbers, weight, xs, ys, qx, qy)
    i, j, k = @index(Global, NTuple)
    dest[i, j, k] += weight * random_numbers[k] * cis(qx * xs[i] + qy * ys[j])
end

function sample_subharmonic_phase_screen!(dest, buffer, spectrum, param, Δx, Δy)
    Nx, Ny = size(dest)

    xs = StepRangeLen(0, Δx, Nx)
    ys = StepRangeLen(0, Δy, Ny)

    kernel! = subharmonic_kernel!(get_backend(dest))
    ndrange = size(dest)

    N = size(buffer, 4)
    M = size(buffer, 3)

    hermitian_randn!(buffer)

    for n in axes(buffer, 4), m in axes(buffer, 3)
        (n in (N ÷ 2, N ÷ 2 + 1) && m in (M ÷ 2, M ÷ 2 + 1)) && continue
        for p in axes(buffer, 2)
            Δqx = 2π / (Nx * Δx * 3^p)
            Δqy = 2π / (Ny * Δy * 3^p)
            qx = (m - M ÷ 2 - 1 // 2) * Δqx
            qy = (n - N ÷ 2 - 1 // 2) * Δqy
            weight = sqrt(spectrum((qx, qy), param) * Δqx * Δqy)
            kernel!(dest, view(buffer, :, p, m, n), weight, xs, ys, qx, qy; ndrange)
        end
    end
end
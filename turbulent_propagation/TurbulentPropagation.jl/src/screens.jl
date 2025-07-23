function ft_phase_screen(spectrum, rs, repetitions)
    N = length(rs)
    Δk = π / last(rs)
    ks = range(-N ÷ 2, length=N) * Δk

    σs = √2 * Δk * map(k -> √spectrum(k), Iterators.product(ks, ks))

    cs = map(c -> c .* σs, eachslice(randn(ComplexF32, N, N, repetitions), dims=3)) |> stack

    plan = plan_bfft!(similar(cs, (N, N)))
    cache = similar(cs, (length(rs), length(rs)))

    for n in axes(cs, 3)
        ifftshift!(cache, view(cs, :, :, n))
        plan * cache
        fftshift!(view(cs, :, :, n), cache)
    end

    cs |> real
end

function ft_phase_screen(spectrum, rs)
    dropdims(ft_phase_screen(spectrum, rs, 1), dims=3)
end

function sh_phase_screen(spectrum, rs, P, repetitions)
    cs = √2 * randn(ComplexF64, 3, 3, P, repetitions)
    Δk = π / last(rs)

    ks = [i * Δk / 3^p for i in (-1:1), p in 1:P]

    @tullio σs[m, n, p] := Δk / 3^p * √spectrum((ks[m, p], ks[n, p]))

    @tullio ϕ[m, n, o] := cs[i, j, p, o] * σs[i, j, p] * cis(ks[i, p] * rs[m] + ks[j, p] * rs[n])

    ϕ |> real
end

function sh_phase_screen(spectrum, rs, P)
    reshape(sh_phase_screen(spectrum, rs, P, 1), length(rs), length(rs))
end

function ft_sh_phase_screen(spectrum, rs, P, repetitions)
    screens = ft_phase_screen(spectrum, rs, repetitions) + sh_phase_screen(spectrum, rs, P, repetitions)
    Threads.@threads for n in axes(screens, 3)
        screens[:, :, n] = view(screens, :, :, n) .- mean(view(screens, :, :, n))
    end
    screens
end

function ft_sh_phase_screen(spectrum, rs, P)
    dropdims(ft_sh_phase_screen(spectrum, rs, P, 1), dims=3)
end
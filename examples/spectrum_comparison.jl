using CairoMakie, TurbulentPropagation

with_theme(theme_latexfonts()) do
    fig = Figure()
    ax = Axis(fig[1, 1], xscale=log10, xlabel="Wavenumber (1/m)", ylabel="Normalized Spectrum")

    κs = logrange(1e-2, 1e6, 512)
    r₀ = 1
    L₀ = 50
    l₀ = 1e-4
    param = (r₀, L₀, l₀)

    kolmogorov = [kolmogorov_spectrum(κ, param) for κ in κs]
    von_karman = [von_karman_spectrum(κ, param) for κ in κs]
    modified_von_karman = [modified_von_karman_spectrum(κ, param) for κ in κs]
    hill_andrews = [hill_andrews_spectrum(κ, param) for κ in κs]

    lines = [von_karman ./ kolmogorov, modified_von_karman ./ kolmogorov, hill_andrews ./ kolmogorov]
    labels = ["Von Karman", "Modified Von Karman", "Hill-Andrews"]

    for (line, label) in zip(lines, labels)
        lines!(ax, κs, line, label=label)
    end

    axislegend(ax, position=:lt)

    fig
end
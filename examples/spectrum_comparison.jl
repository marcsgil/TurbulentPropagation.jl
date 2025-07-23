using TurbulentPropagation, CairoMakie

Cₙ² = 1e-14
l₀ = 1e-4
L₀ = 50

param = (Cₙ², L₀, l₀)

qs = logrange(1e-2, 1e6, 256)
kolmogorov = map(q -> kolmogorov_spectrum(q, param), qs)

with_theme(theme_latexfonts()) do
    fig = Figure(fontsize=20)
    ax = Axis(fig[1, 1], xlabel="q", ylabel="Normalized Spectrum", xscale=log10)
    linewidth = 3
    #lines!(ax, qs, ; label="Kolmogorov", linewidth)
    lines!(ax, qs, map(q-> von_karman_spectrum(q, param), qs) ./ kolmogorov; label="Von Karman", linewidth)
    lines!(ax, qs, map(q-> modified_von_karman_spectrum(q, param), qs) ./ kolmogorov; label="Modified von Karman", linewidth, linestyle=:dot)
    lines!(ax, qs, map(q-> hill_andrews_spectrum(q, param), qs) ./ kolmogorov; label="Hill-Andrews", linewidth, linestyle=:dash)
    vlines!(ax, 1/L₀, color = :black, linestyle=:dash)
    vlines!(ax, 1/l₀, color = :black, linestyle=:dash)
    axislegend(ax, position=:lt)
    text!(ax, 0.67, 0.1, text=L"l_0^{-1}", space = :relative)
    text!(ax, 0.01, 0.1, text=L"L_0^{-1}", space = :relative)
    fig
end
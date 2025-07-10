using TurbulentPropagation, CairoMakie

l₀ = 1e-4
L₀ = 50
Cₙ² = 1e-14

qs = logrange(1e-2, 1e6, 256)
kolmogorov = kolmogorov_spectrum.(qs, Cₙ²)

with_theme(theme_latexfonts()) do
    fig = Figure(fontsize=20)
    ax = Axis(fig[1, 1], xlabel="q", ylabel="Normalized Spectrum", xscale=log10)
    linewidth = 3
    #lines!(ax, qs, ; label="Kolmogorov", linewidth)
    lines!(ax, qs, von_karman_spectrum.(qs, Cₙ², L₀) ./ kolmogorov; label="Von Karman", linewidth)
    lines!(ax, qs, modified_von_karman_spectrum.(qs, Cₙ², L₀, l₀) ./ kolmogorov; label="Modified von Karman", linewidth, linestyle=:dot)
    lines!(ax, qs, hill_andrews_spectrum.(qs, Cₙ², L₀, l₀) ./ kolmogorov; label="Hill-Andrews", linewidth, linestyle=:dash)
    vlines!(ax, 1/L₀, color = :black, linestyle=:dash)
    vlines!(ax, 1/l₀, color = :black, linestyle=:dash)
    axislegend(ax, position=:lt)
    text!(ax, 0.67, 0.1, text=L"l_0^{-1}", space = :relative)
    text!(ax, 0.01, 0.1, text=L"L_0^{-1}", space = :relative)
    fig
end
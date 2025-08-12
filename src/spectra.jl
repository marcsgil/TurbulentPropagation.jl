function modified_von_karman_spectrum(κs, param)
    r₀, L₀, l₀ = param

    κ² = sum(abs2, κs)

    T = promote_type(eltype(κ²), eltype.(param)...)

    κ₀ = T(2π) / L₀
    κₘ = T(5.92) / l₀

    T(0.49) * exp(-κ² / κₘ^2) / r₀^(5 // 3) / (κ² + κ₀^2)^(11 // 6)
end

function von_karman_spectrum(κs, param)
    r₀, L₀ = param
    modified_von_karman_spectrum(κs, (r₀, L₀, false))
end

function kolmogorov_spectrum(κs, param)
    r₀ = first(param)
    von_karman_spectrum(κs, (r₀, oftype(float(r₀), Inf)))
end

function hill_andrews_spectrum(κs, param)
    r₀, L₀, l₀ = param

    κ² = sum(abs2, κs)

    T = promote_type(eltype(κ²), eltype.(param)...)

    κ₀ = T(2π) / L₀
    κₗ = T(3.3) / l₀

    K² = κ² / κₗ^2

    a₀ = T(0.49)
    a₁ = T(1.802)
    a₂ = T(0.254)

    a₀ * exp(-K²) / r₀^(5 // 3) / (κ² + κ₀^2)^(11 // 6) * (1 + a₁ * √K² - a₂ * K²^(7 // 12))
end
"""
    modified_von_karman_spectrum(q, Cₙ², L₀, l₀)

Calculate the modified von Karman spectrum for a given wavevector `q`, refractive index structure constant `Cₙ²`, inner scale `l₀`, and outer scale `L₀`.

The modified von Karman spectrum is
```math
\\Phi(q) = 0.033 Cₙ² \\frac{\\exp\\left[-|q|^2 / (5.92/l₀)^2\\right]}{\\left[q^2 + (2\\pi/L₀)^2\\right]^{11/6}}
```
"""
function modified_von_karman_spectrum(q, Cₙ², L₀, l₀)
    q² = sum(abs2, q)
    T = promote_type(eltype(q), eltype(Cₙ²), eltype(l₀), eltype(L₀))
    q²ₘ = T(5.92 / l₀)^2
    q²₀ = T(2π / L₀)^2
    T(0.033) * Cₙ² * exp(-q² / q²ₘ) / (q² + q²₀)^(11 / 6)
end

"""
    von_karman_spectrum(q, Cₙ², L₀)

Calculate the von Karman spectrum for a given wavevector `q`, refractive index structure constant `Cₙ²`, and outer scale `L₀`.

The von Karman spectrum is a special case of the modified von Karman spectrum with `l₀ = 0`:
"""
function von_karman_spectrum(q, Cₙ², L₀)
    modified_von_karman_spectrum(q, Cₙ², L₀, zero(L₀))
end

"""
    kolmogorov_spectrum(q, Cₙ²)

Calculate the Kolmogorov spectrum for a given wavevector `q` and refractive index structure constant `Cₙ²`.

The Kolmogorov spectrum is a special case of the von Karman spectrum with `L₀ = Inf`:
"""
function kolmogorov_spectrum(q, Cₙ²)
    von_karman_spectrum(q, Cₙ², oftype(Cₙ², Inf))
end

"""
    hill_andrews_spectrum(q, Cₙ², L₀, l₀)
"""
function hill_andrews_spectrum(q, Cₙ², L₀, l₀)
    q² = sum(abs2, q)
    T = promote_type(eltype(q), eltype(Cₙ²), eltype(l₀), eltype(L₀))
    qₗ² = (T(3.3) / l₀)^2
    modified_von_karman_spectrum(q, Cₙ², L₀, l₀) * (1 + T(1.802) * √(q² / qₗ²) - T(0.254) * (q² / qₗ²)^(7 // 12))
end
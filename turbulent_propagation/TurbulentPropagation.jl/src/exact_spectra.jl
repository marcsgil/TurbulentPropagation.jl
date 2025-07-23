function mvK_spectrum(κ,r₀,l₀,L₀)
    0.49 * r₀^(-5/3) * exp(-sum(abs2,κ) * (l₀ / 5.92)^2) / (sum(abs2,κ) + (2π/L₀)^2)^(11/6)
end


function mvK_structure_function(r, r₀, l₀, L₀)
    7.75 * r₀^(-5/3) * r^2 * ( inv(l₀+2.03(r)^2)^(1/6) - 0.72 * (2π/L₀)^(1/3) )
end
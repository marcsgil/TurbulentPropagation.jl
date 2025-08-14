function fried_parameter(z, Cₙ², λ)
    T = promote_type(eltype.((z, Cₙ², λ))...)
    T(0.185) * (λ^2 / Cₙ² / z)^(3 // 5)
end

function rytov_variance(z, Cₙ², λ)
    T = promote_type(eltype.((z, Cₙ², λ))...)
    T(1.23) * Cₙ² * (2π / λ)^(7 // 6) * z^(11 // 6)
end
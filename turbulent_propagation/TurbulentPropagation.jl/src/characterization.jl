function structure_function(x::AbstractArray{T,2}) where T
    m = size(x,1) ÷ 2 + 1
    Δm = size(x,1) - m
    n = size(x,2) ÷ 2 + 1
    Δn = size(x,2) - n

    x₀ = x[m, n]

    sum(( ((@view x[m:end, n]) .- x₀).^2,
          ((@view x[m:-1:m-Δm, n]) .- x₀).^2,
          ((@view x[m, n:end]) .- x₀).^2,
          ((@view x[m, n:-1:n-Δn]) .- x₀).^2, )) / 4
end

function structure_function(x::AbstractArray{T,3}) where T
    sum(structure_function(y) for y in eachslice(x,dims=3))/size(x,3)
end
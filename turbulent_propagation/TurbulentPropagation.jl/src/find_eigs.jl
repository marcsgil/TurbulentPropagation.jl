function fresnel_kernel(x,y,k,Δz)
    k / 2π / im / Δz * cis(k * (x^2+y^2) / 2 / Δz)
end

function turbulence_eigenmodes(rs,z,σ2,k)
    #=α = 1.23/.423
    r₀ = ( α * (z/k)^(5/6) / σ2 )^(5/3)
    spectrum(κ) = mvK_spectrum(κ,r₀,1e-4,100)
    screen = ft_sh_phase_screen(spectrum, rs, 3)*1e-6 |> Array=#

    w0 = .02
    zr = k*w0^2/2


    Z = z/zr
    Rs = √(1+Z^2) * rs
    @tullio _operator[μ,ν,α,β] := fresnel_kernel(Rs[μ] - rs[α],Rs[ν] - rs[β], k, z) #* cis(screen[α,β])

    N = length(rs)
    operator = reshape(_operator,N^2,N^2)

    eigsolve(operator,1)
end
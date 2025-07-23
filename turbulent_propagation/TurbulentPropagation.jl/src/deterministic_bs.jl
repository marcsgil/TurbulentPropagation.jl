function deterministic_Bs_ft(rs, r₀, l₀, L₀)

    σs = reciprocal_interval(rs) * map(k -> √mvK_spectrum(k,r₀,l₀,L₀),reciprocal_grid(rs,rs))

    #σs.^2 |> ifftshift |> bfft |> fftshift |> real

    ks = reciprocal_grid(rs)
    @tullio result[m,n] := σs[i,j]^2 * cis( ks[i]*rs[m] + ks[j]*rs[n] )
    real(result)
end

function deterministic_Bs_sh(rs, r₀, l₀, L₀,P)
    Δk = reciprocal_interval(rs)
    ks_sh = [ i * Δk/3^p for i in (-1:1), p in 1:P ]

    @tullio σs[m,n,p] := Δk / 3^p * √mvK_spectrum((ks_sh[m,p],ks_sh[n,p]),r₀,l₀,L₀)

    @tullio result[m,n] := σs[i,j,p]^2 * cis( ks_sh[i,p]*rs[m] + ks_sh[j,p]*rs[n] )
    real(result)
end

function deterministic_Ds_ft(rs,r₀,l₀,L₀)
    Bs = deterministic_Bs_ft(rs, r₀, l₀, L₀)
    2*(maximum(Bs).-Bs)
end


function deterministic_Ds_ft_sh(rs,r₀,l₀,L₀,P)
    Bs = deterministic_Bs_ft(rs, r₀, l₀, L₀) + deterministic_Bs_sh(rs, r₀, l₀, L₀,P)
    2*(maximum(Bs).-Bs)
end
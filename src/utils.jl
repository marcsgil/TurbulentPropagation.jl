struct TurbulentPropagationBuffers{T1,T2,T3,T4,T5,T6}
    u::T1
    phase::T2
    spectrum_buffer::T3
    subharmonic_randoms::T4
    fft_plan::T5
    ifft_plan::T6
    function TurbulentPropagationBuffers(u::AbstractArray, nsubharmonics=3, subharmonic_levels=6)
        phase = similar(u)
        spectrum_buffer = similar(u, size(u, 1), size(u, 2))
        subharmonic_randoms = similar(u, size(u, 3), nsubharmonics, subharmonic_levels, subharmonic_levels)
        fft_plan = plan_fft!(u, (1, 2))
        ifft_plan = inv(fft_plan)
        args = (u, phase, spectrum_buffer, subharmonic_randoms, fft_plan, ifft_plan)
        new{typeof.(args)...}(args...)
    end
end
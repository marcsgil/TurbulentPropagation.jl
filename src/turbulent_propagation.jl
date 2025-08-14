function turbulent_propagation!(buffers::TurbulentPropagationBuffers, spectrum, param, Δx, Δy, z, λ, method=AngularSpectrum(); nsteps=1)
    Δz = z / 2nsteps

    for _ in 1:nsteps
        free_propagation!(buffers, Δx, Δy, Δz, λ, method)
        sample_phase_screen!(buffers, spectrum, param, Δx, Δy)
        @. buffers.u *= cis(buffers.phase)
        free_propagation!(buffers, Δx, Δy, Δz, λ, method)
    end
end
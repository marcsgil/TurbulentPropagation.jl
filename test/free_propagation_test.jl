@testset "Free Propagation" begin
    L = 8
    N = 256
    Δ = L / N
    rs = StepRangeLen(-Δ * (N ÷ 2), Δ, N)
    zs = LinRange(0.01, 1, 32)

    #This is a gaussian mode
    u0 = lg(rs, rs)
    scalings = @. √(1 + 4 * zs^2) #Here, we introduce the scalings given by w(z)/w0

    λ = 2π

    #Now we propagate, including the scalings
    ψs = free_propagation(u0, rs, rs, zs, scalings)

    plan = plan_fft!(u0)
    iplan = inv(plan)

    buffer = similar(u0)
    for n ∈ eachindex(zs, scalings)
        copy!(buffer, u0)
        TurbulentPropagation.angular_spectrum_propagation!(buffer, Δ, Δ, zs[n], λ, scalings[n], plan, iplan)
        @test ψs[:, :, n] ≈ buffer
    end
end
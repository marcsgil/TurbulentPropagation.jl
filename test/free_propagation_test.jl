@testset "Free Propagation" begin
    L = 8
    N = 256
    d = L / N
    rs = StepRangeLen(-d * (N ÷ 2), d, N)
    zs = LinRange(0.01, 1, 32)

    #This is a gaussian mode
    u0 = [complex(exp((-x^2 - y^2))) for x in rs, y in rs]
    scalings = @. √(1 + 4 * zs^2) #Here, we introduce the scalings given by w(z)/w0

    λ = 2π

    plan = plan_fft!(u0)
    iplan = inv(plan)

    buffer = copy(u0)
    for n ∈ eachindex(zs, scalings)
        prev_scaling = n > 1 ? scalings[n-1] : 1
        current_scaling = scalings[n] / prev_scaling
        dz = n > 1 ? zs[n] - zs[n-1] : zs[n]
        TurbulentPropagation.angular_spectrum_propagation!(buffer, d, d, dz, λ, current_scaling, plan, iplan)
        d *= current_scaling
        @test isapprox(abs2.(buffer * scalings[n]), abs2.(u0), rtol=1e-6)
    end
end
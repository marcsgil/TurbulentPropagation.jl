@testset "Free Propagation" begin
    rs = LinRange(-3, 3, 256)

    u_num = Array{ComplexF64}(undef, length(rs), length(rs))
    k = 2π
    lg!(u_num, rs, rs; k)
    u_anl = similar(u_num)

    plan = plan_fft!(u_num, (1, 2))
    iplan = inv(plan)

    zs = LinRange(0, 3, 256)

    for z in zs
        free_propagation!(u_num, step(rs), step(rs), step(zs), 1, plan, iplan)
        lg!(u_anl, rs, rs, z; k)

        @test isapprox(u_num, u_anl; rtol=1e-2)
    end
end
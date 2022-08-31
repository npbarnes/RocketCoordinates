using Test
using StaticArrays
using LinearAlgebra
using RocketCoordinates

for (x,y,z) in (normalize(@SVector(randn(3))) for _ in 1:10) # Sample points uniformly on a sphere
    λ = geocentric_latitude(x,y,z)
    Ω = geocentric_longitude(x,y,z)

    N = N̂_gc(λ, Ω)
    E = Ê_gc(λ, Ω)
    D = D̂_gc(λ, Ω)

    @test N ⋅ E ≈ 0 atol=1e-12
    @test N ⋅ D ≈ 0 atol=1e-12
    @test E ⋅ D ≈ 0 atol=1e-12
    @test norm(N) ≈ 1
    @test norm(E) ≈ 1
    @test norm(D) ≈ 1
end

@test N̂_gc(0,0) ≈ SA[0,0,1]

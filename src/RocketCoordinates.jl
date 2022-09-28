module RocketCoordinates

export NED_matrix, IGRF_NED
export MZP_matrix, XYZ_matrix

using Geodesy
using SatelliteToolbox
using CSV
using DataFrames
using StaticArrays
using LinearAlgebra

using DataLoaders

const bodypositions_file = "$(pkgdir(RocketCoordinates))/data/delamereBodyPositions.txt"

geocentric_latitude(x::ECEF) = geocentric_latitude(x.x, x.y, x.z)
geocentric_latitude(x,y,z) = atan(z, √(x^2 + y^2))
geocentric_longitude(x::ECEF) = geocentric_longitude(x.x, x.y, x.z)
geocentric_longitude(x,y,z) = atan(y,x)

# Consider the (cartographic) spherical coordinate system:
# r = √(x^2 + y^2 + z^2)
# λ = atan(z, √(x^2 + y^2))
# Ω = atan(y,x)
# The functions below give unit tangent vectors to this spheical coordinate system
# analogous to r̂, θ̂, and ϕ̂ in standard spherical coordinates.
# North (N̂) is the direction of increasing latitude (not polar angle),
# East (Ê) is the direction of increasing longitude, and
# Down (D̂) is the direction of *decreasing* radius.
# in coordinates dirived from ECEF coordinate axes.
# The matrix [N̂ Ê D̂] converts a tangent vector from NED to ECEF coordinates, and
# the matrix transpose([N̂ Ê D̂]) converts a tangent vector from ECEF to NED
N̂_gc(x::ECEF) = N̂_gc(geocentric_latitude(x), geocentric_longitude(x))
Ê_gc(x::ECEF) = Ê_gc(geocentric_latitude(x), geocentric_longitude(x))
D̂_gc(x::ECEF) = D̂_gc(geocentric_latitude(x), geocentric_longitude(x))
N̂_gc(λ, Ω) = SA[-sin(λ)*cos(Ω), -sin(λ)*sin(Ω), cos(λ)]
Ê_gc(λ, Ω) = SA[-sin(Ω), cos(Ω), 0]
D̂_gc(λ, Ω) = SA[-cos(λ)*cos(Ω), -cos(λ)*sin(Ω), -sin(λ)]

"""
    NED_matrix(x::ECEF)

A matrix that rotates geocentric North, East, Down coordinates to coordinates parallel to
ECEF axes. Does not do any translation or boosting. See [`XYZ_matrix`](@ref) for a more
detailed discussion of translation, rotation, and boosting with examples.

See also [`MZP_matrix`](@ref)
"""
function NED_matrix(x::ECEF)
    λ = geocentric_latitude(x)
    Ω = geocentric_longitude(x)
    [N̂_gc(λ, Ω) Ê_gc(λ, Ω) D̂_gc(λ, Ω)]
end

"""IGRF B field in the geocentric NED coordinates described above"""
IGRF_NED(x::ECEF, year=2021.37) = igrf(year, norm(x), geocentric_latitude(x), geocentric_longitude(x), Val(:geocentric))

"""
    MZP_matrix(x::ECEF)

One possible definition of Meridional, Zonal, Parallel coordinates (not confirmed by Rob).
This matrix rotates MZP vectors to coordinates parallel to ECEF. It does not do any
translation or boosting. See [`XYZ_matrix`](@ref) for a more detailed disscussion of
translation, rotation, and boosting with examples.

The parallel direction is determined by IGRF at the point x in 2021.

See also [`NED_matrix`](@ref)
"""
function MZP_matrix(x::ECEF)
    NEDtoECEF = NED_matrix(x)
    P = normalize(NEDtoECEF * IGRF_NED(x))
    E = NEDtoECEF[:, 2] # East at x in ECEF
    Z = normalize(E - P * (E ⋅ P)) # Zonal is East projected into the perp B plane and normalized in ECEF
    M = Z × P # Meridional completes the right handed triple in ECEF
    [M Z P]
end

"""
    XYZ_matrix(x::ECEF, v::ECEF)

A matrix that rotates Hybrid Code XYZ coordinates be parallel to ECEF. `x` is the origin
of the XYZ coordinate system in ECEF (i.e. the barium release point) and `v` is the velocity
of the payload in ECEF. Note that XYZ is translated, rotated, and boosted relative to ECEF.
(We'll ignore acceleration from gravity, and roatation of B over the timescale of the
experiment.) This matrix only does the rotation; not the translation or the boost. Also, B
is IGRF at the point x which probably cooresponds with the barium canister instead of the
main payload, but the difference should be small.

Given a position x′ in XYZ coordinates we can get its ECEF coordinates by first rotating
the axes by multiplying by XYZ and then adding `x` the ECEF coordinate of the
origin of the XYZ system.
`x′_ECEF = x_ECEF + XYZ * x′_XYZ`

Given a velocity v′ in XYZ derived coordinates we can get its ECEF derived coordinates by
rotating and boosting.
`v′_ECEF = v_ECEF + XYZ * v′_XYZ`

Given a magnetic field vector we can convert from XYZ derived coordinates to ECEF derived
coordinates by `B_ECEF = XYZ * B_XYZ`. This is a special case because a non-relativistic
boost of B doesn't change it's value, and its a tangent vector so we're actually rotating
derived coordinates so we don't have to worry about translation.

Given an electric field vector we can convert from XYZ derived coordinates to ECEF derived
coordinates by `E_ECEF = v_ECEF×B_ECEF + XYZ * E_XYZ` since the non-relativistic boost
formula for E is `E′ = E + v×B`

See also [`NED_matrix`](@ref), and [`MZP_matrix`](@ref)
"""
function XYZ_matrix(x::ECEF, v::ECEF)
    MZP = MZP_matrix(x)
    v_MZP = transpose(MZP) * v
    X = normalize(MZP * SA[v_MZP[1], v_MZP[2], 0]) # X is parallel to the projection of v into the perpendicular plane
    Z = MZP[:, 3] # Z is Parallel
    Y = Z × X # Y completes the right handed coordinate system
    [X Y Z]
end

end # module

module RocketCoordinates

export NED_matrix, IGRF_NED
export MZP_matrix, XYZ_matrix, ROB_matrix
export bodypositions_file
export load_bodypositions
export df
export times, main, ba1, ba2, d1, d2, d3, d4
export R1_time, R2_time

using SatelliteToolbox
using CSV
using DataFrames
using StaticArrays
using LinearAlgebra
using Unitful

using DataLoaders

#const x = ECEF(2.2491957723032855e6, -4.980439036469659e6, 3.915796826007238e6)
const bodypositions_file = "$(pkgdir(RocketCoordinates))/data/delamereBodyPositions.txt"
const df = load_bodypositions(bodypositions_file)
const times = df.time * u"s"
const main = [SA[x,y,z]*u"m" for (x,y,z) in zip(df.main_x, df.main_y, df.main_z)]
const ba1 = [SA[x,y,z]*u"m" for (x,y,z) in zip(df.ba1_x, df.ba1_y, df.ba1_z)]
const ba2 = [SA[x,y,z]*u"m" for (x,y,z) in zip(df.ba2_x, df.ba2_y, df.ba2_z)]
const d1 = [SA[x,y,z]*u"m" for (x,y,z) in zip(df.d1_x, df.d1_y, df.d1_z)]
const d2 = [SA[x,y,z]*u"m" for (x,y,z) in zip(df.d2_x, df.d2_y, df.d2_z)]
const d3 = [SA[x,y,z]*u"m" for (x,y,z) in zip(df.d3_x, df.d3_y, df.d3_z)]
const d4 = [SA[x,y,z]*u"m" for (x,y,z) in zip(df.d4_x, df.d4_y, df.d4_z)]
#=
const mag_file = "$(pkgdir(RocketCoordinates))/data/KinetX_Release_Mag_Data.txt"
const mag = load_magdata(mag_file)
const density_file = "$(pkgdir(RocketCoordinates))/data/KDR_N.txt"
const density = load_density(density_file)
=#

# Brian Tibbetts' Release Analysis "Times in trajectory file" times of each barium release
const R1_time = 593.23u"s"
const R2_time = 627.47u"s"

geocentric_latitude(x) = geocentric_latitude(x[1], x[2], x[3])
geocentric_latitude(x,y,z) = atan(z, √(x^2 + y^2))
geocentric_longitude(x) = geocentric_longitude(x[1], x[2], x[3])
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
N̂_gc(x) = N̂_gc(geocentric_latitude(x), geocentric_longitude(x))
Ê_gc(x) = Ê_gc(geocentric_latitude(x), geocentric_longitude(x))
D̂_gc(x) = D̂_gc(geocentric_latitude(x), geocentric_longitude(x))
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
function NED_matrix(x)
    λ = geocentric_latitude(x)
    Ω = geocentric_longitude(x)
    [N̂_gc(λ, Ω) Ê_gc(λ, Ω) D̂_gc(λ, Ω)]
end

"""IGRF B field in the geocentric NED coordinates described above"""
IGRF_NED(x::AbstractArray{<:Unitful.Length}) = IGRF_NED(ustrip.(u"m", x))
IGRF_NED(x, year=2021.37) = igrf(year, norm(x), geocentric_latitude(x), geocentric_longitude(x), Val(:geocentric))

"""
    MZP_matrix(x::ECEF)

One possible definition of Meridional, Zonal, Parallel coordinates (not confirmed by Rob).
This matrix rotates MZP vectors to coordinates parallel to ECEF. It does not do any
translation or boosting. See [`XYZ_matrix`](@ref) for a more detailed disscussion of
translation, rotation, and boosting with examples.

The parallel direction is determined by IGRF at the point x in 2021.

See also [`NED_matrix`](@ref)
"""
function MZP_matrix(x)
    NEDtoECEF = NED_matrix(x)
    P = normalize(NEDtoECEF * IGRF_NED(x))
    E = NEDtoECEF[:, 2] # East at x in ECEF
    Z = normalize(E - P * (E ⋅ P)) # Zonal is East projected into the perp B plane and normalized in ECEF
    M = Z × P # Meridional completes the right handed triple in ECEF
    [M Z P]
end

"""
    ROB_matrix(x::ECEF)

Meridional, Zonal, Parallel coordinates as defined by Rob
"""
function ROB_matrix(x)
    NEDtoECEF = NED_matrix(x)
    B = normalize(NEDtoECEF * IGRF_NED(x))
    Z = normalize(B × x)
    M = normalize(Z × B)
    [M Z B]
end
"""
    XYZ_matrix(x::ECEF, v::ECEF)

A matrix that rotates Hybrid Code XYZ coordinates be parallel to ECEF. `x` is the origin
of the XYZ coordinate system in ECEF (i.e. the barium release point) and `v` is the velocity
of the payload in ECEF. Note that XYZ is translated, rotated, and boosted relative to ECEF.
(We'll ignore acceleration from gravity, and rotation of B over the timescale of the
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
function XYZ_matrix(x, v)
    MZP = MZP_matrix(x)
    v_MZP = transpose(MZP) * v
    X = normalize(MZP * SA[v_MZP[1], v_MZP[2], 0]) # X is parallel to the projection of v into the perpendicular plane
    Z = MZP[:, 3] # Z is Parallel
    Y = Z × X # Y completes the right handed coordinate system
    [X Y Z]
end

end # module

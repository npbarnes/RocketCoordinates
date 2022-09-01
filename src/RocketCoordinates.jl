module RocketCoordinates

export loaddata
export NED_matrix, IGRF_NED
export MZP_matrix, XYZ_matrix

using Geodesy
using SatelliteToolbox
using CSV
using DataFrames
using StaticArrays
using LinearAlgebra

const datafilename = "$(pkgdir(RocketCoordinates))/data/delamereBodyPositions.txt"

geocentric_latitude(x::ECEF) = geocentric_latitude(x.x, x.y, x.z)
geocentric_latitude(x,y,z) = atan(z, √(x^2 + y^2))
geocentric_longitude(x::ECEF) = geocentric_longitude(x.x, x.y, x.z)
geocentric_longitude(x,y,z) = atan(y,x)

loaddata() = CSV.read(datafilename, DataFrame;
    header = [
        :time,
        :main_x, :main_y, :main_z,
        :ba1_x, :ba1_y, :ba1_z,
        :ba2_x, :ba2_y, :ba2_z,
        :d1_x, :d1_y, :d1_z,
        :d2_x, :d2_y, :d2_z,
        :d3_x, :d3_y, :d3_z,
        :d4_x, :d4_y, :d4_z,
    ],
    delim = " ",
    ignorerepeated = true
)

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
function NED_matrix(x::ECEF)
    λ = geocentric_latitude(x)
    Ω = geocentric_longitude(x)
    [N̂_gc(λ, Ω) Ê_gc(λ, Ω) D̂_gc(λ, Ω)]
end


"""IGRF B field in the NED coordinates described above"""
IGRF_NED(x::ECEF, year=2021) = igrf(year, norm(x), geocentric_latitude(x), geocentric_longitude(x), Val(:geocentric))

"""
One possible definition of Meridional, Zonal, Parallel coordinates (not confirmed by Rob).
Multiply a tangent vector by MZP to convert from MZP to ECEF
Multiply a tangent vector by transpose(MZP) to convert from ECEF to MZP
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
Hybrid Code XYZ coordinates
Multiply a tangent vector by XYZ to convert from XYZ to ECEF
Multiply a tangent vector by transpose(XYZ) to convert from ECEF to XYZ

vec_ECEF = XYZ * vec_XYZ
vec_XYZ = transpose(XYZ) * vec_ECEF
vec_XYZ = transpose(XYZ) * MZP * vec_MZP
vec_MZP = transpose(MZP) * XYZ * vec_XYZ
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

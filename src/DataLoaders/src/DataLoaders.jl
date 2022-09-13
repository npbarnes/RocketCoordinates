module DataLoaders

export load_bodypositions, load_electricfield, load_deltamagneticfield, load_density

using CSV
using DataFrames

load_bodypositions(filename) = CSV.read(filename, DataFrame;
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

load_electricfield(filename) = CSV.read(filename, DataFrame;
    skipto = 4,
    header = [
        :time,
        :Emer, :Ezon, :E15, :E62,
        :Eperpflag, :E15flag, :E62flag,
        :Vmer, :Vzon
    ],
    delim = " ",
    ignorerepeated = true
)

load_deltamagneticfield(filename) = CSV.read(filename, DataFrame;
    skipto = 3,
    header = [
        :time,
        :dBmer, :dBzon
    ],
    delim = " ",
    ignorerepeated = true
)

load_density(filename) = CSV.read(filename, DataFrame;
    skipto = 3,
    header = [
        :time,
        :density, :flag
    ],
    delim = " ",
    ignorerepeated = true
)

end # module DataLoaders

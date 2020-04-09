module eqn

using LinearAlgebra
using DelimitedFiles

include("mhd.jl")
include("sod_shock_init.jl")
include("eigmhd.jl")
include("mhdsound.jl")

end
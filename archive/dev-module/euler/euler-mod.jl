module eqn

using LinearAlgebra
using DelimitedFiles

include("euler.jl")
include("sod_shock_init.jl")
include("eigeuler.jl")
include("eulerP.jl")

end
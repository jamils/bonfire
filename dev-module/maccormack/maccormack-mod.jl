module scheme

using LinearAlgebra
using DelimitedFiles

include("../input.jl")

# Set equation set module
if eqntype == 1
    include("../euler/euler-mod.jl")
elseif eqntype == 2
    include("../mhd/mhd-mod.jl")
end

global eqnset = eqn.eqnset;
global sound = eqn.sound;
global sod_shock_init = eqn.sod_shock_init;

include("maccormack.jl")
include("mac_solver.jl")

end
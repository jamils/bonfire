# Use Julia base modules
using LinearAlgebra
using DelimitedFiles
using PyCall

# Load input file
include("input.jl")

# Set equation set module
if eqntype == 1
    include("euler/euler-mod.jl")
elseif eqntype == 2
    include("mhd/mhd-mod.jl")
end

global eqnset = eqn.eqnset;
global sound = eqn.sound;
global sod_shock_init = eqn.sod_shock_init;

# Set solver scheme module
if schemetype == 1
    include("maccormack/maccormack-mod.jl")
elseif schemetype == 2
    include("muscl/muscl-mod.jl")
end

global solver = scheme.solver;

# Initialization
Q0, F0, neq = sod_shock_init(rho_0, p_0, u_0, dx, dt, tstop, n, γ, CFL, Bx_0, By_0, Bz_0);

# Run simulation
Q, steps = solver(dx, dt, tstop, n, γ, CFL, neq, Q0, F0);
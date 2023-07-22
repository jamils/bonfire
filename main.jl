# load dependencies
using LinearAlgebra, LoopVectorization, Documenter, OhMyREPL, Suppressor
#=
    LoopVectorization provides @turbo and 
    @tturbo macros for for loops and 
    broadcast statements. This should
    greatly increase computation time
=#

# load input file
include("input.jl")

# Select equation set
if EquationSet == 1
    printstyled("Loading ", color=:green)
    printstyled("Euler ", bold=true, color=:blue)
    print("module... \n")
    @suppress_err include("eqn/euler-module.jl")
elseif EquationSet == 2
    printstyled("Loading ", color=:green)
    printstyled("MHD ", bold=true, color=:blue)
    print("module... \n")
    @suppress_err include("eqn/mhd-module.jl")
else
    error("Invalid `EquationSet` value")
end

global eqnset = eqn.eqnset;
global sound = eqn.sound;
global eig = eqn.eig;
global SodShockInit = eqn.SodShockInit;

# Select solver
if SchemeType == 1
    @suppress_err include("scheme/maccormack-module.jl")
elseif SchemeType == 2
    @suppress_err include("scheme/muscl-module.jl")
else
    error("Invalid `SchemeType` value")
end

global solver = scheme.solver;

# Initialization
printstyled("Initilizing... \n", color=:green)
Q₀, F₀, neq = SodShockInit(ρ₀, p₀, u₀, dx, dt, tstop, n, γ, CFL, Bx₀, By₀, Bz₀);

sound(Q₀, γ);

# Run simulation
# printstyled("Running ", color=:green)
# print("main loop \n")
@time Q, steps = solver(dx, dt, tstop, n, γ, CFL, neq, Q₀, F₀);

println("\nFinished")

# include("MakePlots.jl")
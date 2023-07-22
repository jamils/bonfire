# Z-pinch equilibrium input file #
# Kolter Bradshaw #

# grab dependencies
using LinearAlgebra
using DelimitedFiles
include("z-pinch_init.jl")
include("solver.jl")
include("maccormack.jl")
include("muscl.jl")
include("mhd.jl")
include("eulereq.jl")
include("mhdsound.jl")
include("eigeuler.jl")
include("eigmhd.jl")

# choose set of equations to run
global MHD = true;
global EULER = false;

# choose numerical scheme
global RUN_MACCORMACK = false;
global RUN_MUSCL = true;

# choose whether to output data files
global OUTPUT_DATA = true;

# MUSCL parameters
if RUN_MUSCL == true
    # Limiters: 1 - MINMOD
    #           2 - SUPERBEE
    #           3 - VAN-LEER
    #           4 - MONOTONIZED
    limiter = 1;
    k = 1/3;
    eps = 0.15;
    muscl_par = [limiter, k, eps];
end

# initial properties
j_0 = [25.4648, 0];
rho_0 = [0, 0];
r_0 = 0.25;
u_0 = [0, 0];
dx = 0.001;
n = Int(1/dx);
gamma = 1.4;
dt = 1;
tstop = 0.15;
CFL = 0.75;

# initialize the problem
Q0, F0, neq = zpinch_init(j_0, rho_0, r_0, u_0, dx, dt, tstop, n, gamma, CFL);

# run fluid simulation
if RUN_MACCORMACK == true
    Q, steps = solver(dx, dt, tstop, n, gamma, CFL, neq, Q0, F0);
elseif RUN_MUSCL == true
    Q, steps = solver(dx, dt, tstop, n, gamma, CFL, neq, Q0, F0, muscl_params=muscl_par);
end
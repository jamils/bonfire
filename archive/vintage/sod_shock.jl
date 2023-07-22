# Sod shock tube input file #
# Kolter Bradshaw #

# grab dependencies
using LinearAlgebra
using DelimitedFiles
include("sod_shock_init.jl")
include("solver.jl")
include("maccormack.jl")
include("muscl.jl")
include("mhd.jl")
include("eulereq.jl")
include("mhdsound.jl")
include("eigeuler.jl")
include("eigmhd.jl")

# choose set of equations to run
global EULER = false;
global MHD = true;

# choose numerical scheme
global RUN_MACCORMACK = true;
global RUN_MUSCL = false;

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
rho_0 = [1.0, 0.125];
p_0 = [1.0, 0.1];
u_0 = [0, 0];
dx = 0.001;
n = Int(1/dx);
gamma = 1.4;
dt = 1;
tstop = 0.15;
CFL = 0.75;

# initial magnetic field
if MHD == true
    Bx_0 = [0.75, 0.75];
    By_0 = [1.0, -1.0];
    Bz_0 = [0, 0];
end

# initialize the problem
if EULER == true
    Q0, F0, neq = sod_shock_init(rho_0, p_0, u_0, dx, dt, tstop, n, gamma, CFL);
elseif MHD == true
    Q0, F0, neq = sod_shock_init(rho_0, p_0, u_0, dx, dt, tstop, n, gamma, CFL, Bx0=Bx_0, By0=By_0, Bz0=Bz_0);
end

# run fluid simulation
if RUN_MACCORMACK == true
    Q, steps = solver(dx, dt, tstop, n, gamma, CFL, neq, Q0, F0);
elseif RUN_MUSCL == true
    Q, steps = solver(dx, dt, tstop, n, gamma, CFL, neq, Q0, F0, muscl_params=muscl_par);
end
# MHD Waves input file #
# Kolter Bradshaw #

# grab dependencies
include("mhd_wave_init.jl")
include("solver.jl")
include("maccormack.jl")
include("muscl.jl")
include("mhd.jl")
include("eulereq.jl")
include("mhdsound.jl")
include("eigeuler.jl")
include("eigmhd.jl")

# initialize the problem
Q0, F0, neq = mhd_wave_init(A, Lx, rho_0, p_0, u_0, v_0, w_0, dx, dt, tstop, n, gamma, CFL, Bx_0, By_0, Bz_0);

# run fluid simulation
if RUN_MACCORMACK == true
    Q, steps = solver(dx, dt, tstop, n, gamma, CFL, neq, Q0, F0);
elseif RUN_MUSCL == true
    Q, steps = solver(dx, dt, tstop, n, gamma, CFL, neq, Q0, F0, muscl_params=muscl_par);
end
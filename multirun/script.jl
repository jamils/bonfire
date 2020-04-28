@everywhere using LinearAlgebra
@everywhere using DelimitedFiles
@everywhere using SharedArrays

@everywhere include("mhd_wave_init.jl")
@everywhere include("solver.jl")
@everywhere include("maccormack.jl")
@everywhere include("muscl.jl")
@everywhere include("mhd.jl")
@everywhere include("eulereq.jl")
@everywhere include("mhdsound.jl")
@everywhere include("eigeuler.jl")
@everywhere include("eigmhd.jl")

# choose set of equations to run
@everywhere global EULER = false;
@everywhere global MHD = true;

# choose numerical scheme
@everywhere global RUN_MACCORMACK = false;
@everywhere global RUN_MUSCL = true;

# choose whether to output data files
@everywhere global OUTPUT_DATA = true;

if MHD == true
    @everywhere global eqnstring = string("mhd");
elseif EULER == true
    @everywhere global eqnstring = string("euler");
end

if RUN_MACCORMACK == true
    @everywhere global schemestring = string("maccormack");
elseif RUN_MUSCL == true
    @everywhere global schemestring = string("muscl");
end

# MUSCL parameters
if RUN_MUSCL == true
    # Limiters: 1 - MINMOD
    #           2 - SUPERBEE
    #           3 - VAN-LEER
    #           4 - MONOTONIZED
    @everywhere limiter = 1;
    @everywhere k = 1/3;
    @everywhere eps = 0.15;
    @everywhere muscl_par = [limiter, k, eps];
end

# initial properties
@everywhere A = 0.000001;
@everywhere Lx = 1;
@everywhere rho_0 = 1.0;
@everywhere p_0 = 1.0/1.4;
@everywhere u_0 = 0;
@everywhere v_0 = 0;
@everywhere w_0 = 0;
# ---
@everywhere diff = 1;
# ---
@everywhere dx = 0.0001 * diff;
@everywhere n = Int64(10000 / diff); # Lx / dx
@everywhere gamma = 1.4;
@everywhere dt = 3.74999999e-6;
@everywhere tstop = 1.0;
@everywhere CFL = 0.75;

# initial magnetic field
if MHD == true
    @everywhere Bx_0 = 1.0;
    @everywhere By_0 = sqrt(2);
    @everywhere Bz_0 = 0.5;
end

len = Int64(10);
neq = Int64(8);

Q0all = SharedArray{Float64}(len, neq, n);
Qall  = SharedArray{Float64}(len, neq, n);

@sync @distributed for i in 1:10
    # initialize the problem
    Q0, F0, neq = mhd_wave_init(A, Lx, rho_0, p_0, u_0, v_0, w_0, dx, dt, tstop, n, gamma, CFL, Bx_0, By_0, Bz_0, i);

    # run fluid simulation
    if RUN_MACCORMACK == true
        Q, steps = solver(dx, dt, tstop, n, gamma, CFL, neq, Q0, F0, i);
    elseif RUN_MUSCL == true
        Q, steps = solver(dx, dt, tstop, n, gamma, CFL, neq, Q0, F0, i, muscl_params=muscl_par);
    end

    Q0all[i,:,:] = Q0;
    Qall[i,:,:] = Q;
end
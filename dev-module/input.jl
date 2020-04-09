# 1 == euler
# 2 == mhd
const eqntype = 1;

# 1 == maccormack
# 2 == muscl
const schemetype = 1;

# Output data from simulation into .csv files
global const OUTPUT_DATA = 1;

# initial properties
global const rho_0 = [1.0, 0.125];
global const p_0 = [1.0, 0.1];
global const u_0 = [0, 0];
global const dx = 0.001;
global const n = Int(1/dx);
global const γ = 1.4;
global const dt = 1;
global const tstop = 0.15;
global const CFL = 0.75;

# MHD initial magnetic field
global const Bx_0 = [0.75, 0.75];
global const By_0 = [1.0, -1.0];
global const Bz_0 = [0, 0];

# MUSCL Limiters: 1 - MINMOD
#           2 - SUPERBEE
#           3 - VAN-LEER
#           4 - MONOTONIZED
global const limiter = 1;
global const k = 1/3;
global const eps = 0.15;
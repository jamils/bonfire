# !! Declaring all as "const" will increase raw performance

"""
1 ⟶ Euler \n 
2 ⟶ MHD \n 
"""
EquationSet = 2;

"""
1 ⟶ MacCormack \n 
2 ⟶ MUSCL \n 
"""
SchemeType = 2;

"""
MUSCL Limiter \n 
1 ⟶ MINMOD \n 
2 ⟶ SUPERBEE \n 
3 ⟶ Van-Leer \n 
4 ⟶ Monotonized \n 
"""
Limiter = 1;

"""
Flux bias factor \n
k = -1  ⟶ Fully upwind \n
k = 0   ⟶ Upwind-bias \n
k = 1/3 ⟶ 3rd order upwind-bias \n
k = 1/2 ⟶ Leaonad's quick \n 
k = 1   ⟶ Central difference \n
"""
k = 1/3;

"""Scheme order: range 0-1"""
ϵ = 0.15;

# initial conditions

"""Density vector"""
global ρ₀ = [1.0, 0.125];

"""Pressure vector"""
global p₀ = [1.0, 0.1];

"""Velocity vector"""
global u₀ = [0, 0];

"""Grid cell size/resolution"""
global dx = 0.001;

"""Number of cells in grid"""
global n = Int(1/dx);

"""Heat capacity ratio, a constant of the material"""
global γ = 1.4;

"""Timestep"""
global dt = 1;

"""Total run-time length of the simulation"""
global tstop = 0.15;

"""Courant-Friedrichs-Lewy convergence condition"""
global CFL = 0.75;

"""Initial magnetic field for MHD - x"""
global Bx₀ = [0.75, 0.75];

"""Initial magnetic field for MHD - y"""
global By₀ = [1.0, -1.0];

"""Initial magnetic field for MHD - z"""
global Bz₀ = [0, 0];
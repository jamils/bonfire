### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 6f962a5b-9a2b-44b9-babe-4853c3eb7894
begin
	using LinearAlgebra
	using DelimitedFiles
	using Plots
	using PlutoUI
end

# ╔═╡ 95aa5d17-a75b-48bc-af41-ab371db5e0d3
md"""
# Bonfire
### Welcome to bonfire.jl! 
This is an interactive playground designed to teach elementary magnetohydrodynamics

As this is a work-in-progress, for now we will skip the necessary introduction into the theory. Let's jump into the actual code

### Code Outline
Generally, I like to draw a 'map' of any source code with which I'm working. This allows me to visualize the flow of the code, and how the individual components and functions interact with each other. I'll try to outline overall sections and the flow of the code, but please check the button below if you would like to view the direct source code inside of this notebook!
"""

# ╔═╡ 195ce57e-a21f-4872-af6e-819d17105a36
md"Show code $(@bind show_code CheckBox(default=true))"

# ╔═╡ 14369438-0b1e-4e3e-9eef-809833a37b9a
md"""
Next, we load the input paramaters. Ideally, this would be the only section of the code that a non-developer would have to interact with. Feel free to change the values of the below variables, doing so may adjust the display of the code and information herein. For this purpose, I will try to utilize some of the more *fun* features of Pluto.
"""

# ╔═╡ 83595be7-bf28-4a08-891b-99ee1dacb7fb
begin
	# 1 == euler
	# 2 == mhd
	const eqntype = 1;
	
	# 1 == maccormack
	# 2 == muscl
	const schemetype = 1;

	# MUSCL Limiters: 1 - MINMOD
	#           2 - SUPERBEE
	#           3 - VAN-LEER
	#           4 - MONOTONIZED
	global const limiter = 1;
	global const k = 1/3;
	global const ϵ = 0.15;
end;

# ╔═╡ 730195ac-480a-4d64-a3f9-a1612332e3e4
md"""
#### Equation Type
For this, there are two equation types: Euler and Magnetohydrodynamic (MHD). The euler equations are a subset of the Navier-Stokes equations and describe some basic ideal fluid flow. The MHD equations are very similar, but also include Maxwell's Equations of electromagnetism and the Lorentz Force Law. MHD essentially describes a fluid of charged particles with are evolved in changing electric and magnetic fields.
This value is denoted by the variable `eqntype` which can have a value of `1` for the Euler equations, or `2` for the MHD equations.

#### Scheme Type
The next variable `schemetype` defines with discretization scheme should be used for the equations. The first option, the *MacCormack Method*, which proceeds with a **predictor** step, and a **corrector** step. 

The second option is the *MUSCL* scheme. MUSCL stands for Monotonic Upwind Scheme for Conservation Laws - it is a type of reconstruction-evolution scheme for the finite volume method.

```math
Q_j^{n+1} = Q_j^n - \frac{\Delta t}{\Delta x}(F_{j+\frac{1}{2}}^{n*} - F_{j-\frac{1}{2}}^{n*})
```

And the flux is computed using the below functions:

```math
Q_{j+\frac{1}{2}}^{L} = Q_j + \frac{\epsilon}{4}[(1-k)\theta_{j-\frac{1}{2}}^+(Q_j - Q_{j-1}) + (1+k)\theta_{j+\frac{1}{2}}^-(Q_{j+1} - Q_j)]
```
```math
Q_{j+\frac{1}{2}}^{R} = Q_{j+1} - \frac{\epsilon}{4}[(1+k)\theta_{j+1}^-(Q_{j+1} - Q_j) + (1-k)\theta_{j+\frac{1}{2}}^+(Q_{j+2} - Q_{j+1})] 
```

The flux computations rely on two factors `k` and `ϵ`. These are also variables defined above, with the values meaning:

```math  
k = 
	 \begin{cases}
	   -1 &\quad\text{fully upwind}\\
	   0 &\quad\text{upwind-bias} \\
	   \frac{1}{3} &\quad\text{3rd order upwind-bias}\\
	   \frac{1}{2} &\quad\text{Leaonad's quick} \\
	   1 &\quad\text{central difference}\\
	 \end{cases},\,\,\,\,\,\,\,\,\,\,\,
\epsilon =     
	\begin{cases}
	   0 &\quad\text{1st order} \\
	   1 &\quad\text{2nd order or higher}\\
	 \end{cases}
```

With the input parameters loaded, the above values will now affect the code that follows.
"""

# ╔═╡ 6d03d72d-ac0a-4f68-aaac-5c902279b160
@bind n Slider(2:2:1000)

# ╔═╡ e7cd3ffb-e168-40b7-b759-c764d8b0fe1b
md"n = $n"

# ╔═╡ 9e542ab1-7ed9-490d-b4b8-7dea4db73401
begin
	# initial properties
	global const ρ₀ = [1.0, 0.125];
	global const p₀ = [1.0, 0.1];
	global const u₀ = [0, 0];
	#global const dx = 0.001;
	#global const n = Int(1/dx);
	dx = 1/n;
	global const γ = 1.4;
	global const dt = 1;
	global const tstop = 0.15;
	global const CFL = 0.75;
	
	# MHD initial magnetic field
	global const Bx₀ = [0.75, 0.75];
	global const By₀ = [1.0, -1.0];
	global const Bz₀ = [0, 0];
	
	# Make plots
	const makeplots = false;
	
	# Output data from simulation into .csv files
	global const OUTPUT_DATA = true;
end;

# ╔═╡ f4935bb2-dcbc-42a2-a564-d0d89a62b619
md"""
#### Initial Properties
This section contains the initials conditions of the simulation. Currently, the input parameters are designed for a two-dimensional system.
* ρ₀ 	- Density vector
* p₀ 	- Pressure vector
* u₀    - Velocity Vector
* dx    - Grid cell size/resolution
* n     - Number of cells in the grid
* γ     - Heat capacity ratio `` \frac{C_P}{C_V} ``, a constant of the material
* dt    - Timestep
* tstop - Total run-time length of the simulation
* CFL   - Courant-Friedrichs-Lewy coundition: convergence condition for explicity time integration of hyperbolic PDEs
The Initial Magnetic Field section contains all three directional components of the magnetic field vectors Bx₀, By₀, and Bz₀.
"""

# ╔═╡ ce5870f9-2d91-4f49-b98b-183517bbbe12
md"""
---
Once the initial conditions are set and you run the program, the first thing that bonfire does is load the discretized equation sets. This includes things like the eigenvalue and eigenvector calcuations, sound speed, and currently a default initialization method for the Sod Shock sample run. Everything is loaded into RAM as a function, and will not be analyzed until called. Once the equation sets are loaded, bonfire then needs to load the numerical scheme. The schemes are also loaded in as functions.

Once all of the necessary modules and functions are loaded, bonfire will prime and initialize itself for the Sod Shock example. This initialization runs through all of the variables, functions, and modules onces and produces the first timestep. With the Julia programming language, the initialization compiles all of the functions to low-level LLVM and assembly code. The benefit to this is that, at run time, the code can run at assembly level. 

Once the initialization is complete, bonfire will then run the solver. The solver will loop through the numerical scheme at each time step until the maximum simulation time has been reached. Depending on the values chosen, this can be anywhere between 2 to 10,000 steps. Each step forward in time depends on information calculated in the previous timestep.
"""

# ╔═╡ 220e7253-b2cc-4603-8a3b-f349d0727853
if eqntype == 1
	function eqnset(Q, γ)
		ρ = Q[1];
		ρu = Q[2];
		E = Q[3];
		P = (γ - 1)*(E - 0.5*(ρu^2)/ρ);
		F = [ρu; (ρu^2)/ρ + P; (ρu/ρ)*(E + P)];
		return F
	end
	function sound(Q, γ)
		P = @. (γ - 1)*(Q[3, :] - 0.5*(Q[2, :]^2)/Q[1, :]);
		a = @. √(γ*P/Q[1, :]);
		return a
	end
	function eig(QL, QR, γ)
		ρL = QL[1];
		ρuL = QL[2];
		EL = QL[3];
		PL = (γ - 1)*(EL - 0.5*(ρuL^2)/ρL);
		ρR = QR[1];
		ρuR = QR[2];
		ER = QR[3];
		PR = (γ - 1)*(ER - 0.5*(ρuR^2)/ρR);
		uL = ρuL/ρL;
		uR = ρuR/ρR;
		aL = √(γ*PL/ρL);
		aR = √(γ*PR/ρR);
		λL = min(uL, uR) - max(aL, aR);
		λR = min(uL, uR) + max(aL, aR);
		return λL, λR
	end
	function sod_shock_init(ρ₀, p₀, u₀, dx, dt, tstop, n, γ, CFL, Bx₀, By₀, Bz₀)
		half = Int(n/2);
		ρ = zeros(1, n);
		ρ[1:half] = ρ₀[1]*ones(1, half);
		ρ[(half + 1):n] = ρ₀[2]*ones(1, half);
	
		u = zeros(1, n);
		u[1:half] = u₀[1]*ones(1, half);
		u[(half + 1):n] = u₀[2]*ones(1, half);
	
		p = zeros(1, n);
		p[1:half] = p₀[1]*ones(1, half);
		p[(half + 1):n] = p₀[2]*ones(1, half);
	
		ρu = @. ρ*u;
	
		neq = 3;
	
		E = @. p/(γ - 1) + 0.5*(ρu^2)/ρ;
		Q = zeros(neq, n);
	
		Q[1, :] = ρ;
		Q[2, :] = ρu;
		Q[3, :] = E;
	
		F = zeros(neq, n);
	
		for i = 1:n
			F[:, i] = eqnset(Q[:, i], γ);
		end
	
		# if OUTPUT_DATA == 1
		writedlm("sod_shock_initial_conditions_euler.csv", Q, ',')
		# end
	
		return Q, F, neq
	end
	if show_code == true
	md"""
	##### Euler Equation Discretization
	"""
	end
elseif eqntype == 2
	function eqnset(Q,γ)
		ρ = Q[1];
		ρu = Q[2];
		ρv = Q[3];
		ρw = Q[4];
		Bx = Q[5];
		By = Q[6];
		Bz = Q[7];
		E = Q[8];
		B = √(Bx^2 + By^2 + Bz^2);
		P = (γ - 1)*(E - 0.5*(ρu^2 + ρv^2 + ρw^2)/ρ - B^2/2);
		F = [ρu; ρu^2/ρ - Bx^2 + P + B^2/2; ρu*ρv/ρ - Bx*By; ρu*ρw/ρ - Bx*Bz; 0; ρu*By/ρ - Bx*ρv/ρ; ρu*Bz/ρ - Bx*ρw/ρ; (E + P + B^2/2)*ρu/ρ - (Bx*ρu/ρ + By*ρv/ρ + Bz*ρw/ρ)*Bx];
		return F
	end
	function eig(QL, QR, γ)
		ρL = QL[1];
		BL = √(QL[5]^2 + QL[6]^2 + QL[7]^2);
		PL = (γ - 1)*(QL[8] - 0.5*(QL[2]^2 + QL[3]^2 + QL[4]^2)/QL[1]^2 - BL^2/2);
		ρR = QR[1];
		BR = √(QR[5]^2 + QR[6]^2 + QR[7]^2);
		PR = (γ - 1)*(QR[8] - 0.5*(QR[2]^2 + QR[3]^2 + QR[4]^2)/QR[1]^2 - BR^2/2);
	
		uL = QL[2]/ρL;
		uR = QR[2]/ρR;
		cL = √(γ*PL/ρL);
		cR = √(γ*PR/ρR);
		vL = QL[5]/√(ρL);
		vR = QR[5]/√(ρR);
		aL = √(0.5*(BL^2/ρL + cL^2 + √((BL^2/ρL + cL^2)^2 - 4*vL^2*cL^2)));
		aR = √(0.5*(BR^2/ρR + cR^2 + √((BR^2/ρR + cR^2)^2 - 4*vR^2*cR^2)));
		λL = min(uL, uR) - max(aL, aR);
		λR = min(uL, uR) + max(aL, aR);
		return λL, λR
	end
	function sound(Q, γ)
		ρ = Q[1, :];
		B = @. √(Q[5, :]^2 + Q[6, :]^2 + Q[7, :]^2);
		P = @. (γ - 1)*(Q[8, :] - 0.5*(Q[2, :]^2 + Q[3, :]^2 + Q[4, :]^2)/Q[1, :]^2 - B^2/2);
		c = @. √(γ*P/ρ);
		v = @. Q[5, :]/√(ρ);
		a = @. √(0.5*((B^2)/ρ + (c^2) + √((B^2)/ρ + (c^2))^2 - 4*(v^2)*(c^2)));
		return a
	end
	function sod_shock_init(ρ₀, p₀, u₀, dx, dt, tstop, n, γ, CFL, Bx₀, By₀, Bz₀)
		half = Int(n/2);
		ρ = zeros(1, n);
		ρ[1:half] = ρ₀[1]*ones(1, half);
		ρ[(half + 1):n] = ρ₀[2]*ones(1, half);
	
		u = zeros(1, n);
		u[1:half] = u₀[1]*ones(1, half);
		u[(half + 1):n] = u₀[2]*ones(1, half);
	
		p = zeros(1, n);
		p[1:half] = p₀[1]*ones(1, half);
		p[(half + 1):n] = p₀[2]*ones(1, half);
	
		ρu = @. ρ*u;
	
		neq = 8;
	
		Bx = zeros(1, n);
		Bx[1:half] = Bx₀[1]*ones(1, half);
		Bx[(half + 1):n] = Bx₀[2]*ones(1, half);
	
		By = zeros(1, n);
		By[1:half] = By₀[1]*ones(1, half);
		By[(half + 1):n] = By₀[2]*ones(1, half);
	
		Bz = zeros(1, n);
		Bz[1:half] = Bz₀[1]*ones(1, half);
		Bz[(half + 1):n] = Bz₀[2]*ones(1, half);
	
		B = @. √(Bx^2 + By^2 + Bz^2);
		E = @. p/(γ - 1) + 0.5*(ρu^2)/ρ + B^2/2;
		Q = zeros(neq, n);
	
		Q[1, :] = ρ;
		Q[2, :] = ρu;
		Q[3, :] = zeros(1, n);
		Q[4, :] = zeros(1, n);
		Q[5, :] = Bx;
		Q[6, :] = By;
		Q[7, :] = Bz;
		Q[8, :] = E;
	
		F = zeros(neq, n);
		for i = 1:n
			F[:, i] = eqnset(Q[:, i], γ);
		end
	
		# if OUTPUT_DATA == 1
			writedlm("../sod_shock_initial_conditions_mhd.csv", Q, ',')
		# end
	
		return Q, F, neq
	end
	if show_code == true
	md"""
	##### MHD Equation Discretization
	"""
	end
else error("Invalid value for `eqntype`")
end

# ╔═╡ 3748e928-930d-4ecd-9f3e-b681b36e9d93
if show_code == true
if eqntype == 1
	md"""
	```julia
	function eqnset(Q, γ)
		ρ = Q[1];
		ρu = Q[2];
		E = Q[3];
		P = (γ - 1)*(E - 0.5*(ρu^2)/ρ);
		F = [ρu; (ρu^2)/ρ + P; (ρu/ρ)*(E + P)];
		return F
	end
	function sound(Q, γ)
		P = @. (γ - 1)*(Q[3, :] - 0.5*(Q[2, :]^2)/Q[1, :]);
		a = @. √(γ*P/Q[1, :]);
		return a
	end
	function eig(QL, QR, γ)
		ρL = QL[1];
		ρuL = QL[2];
		EL = QL[3];
		PL = (γ - 1)*(EL - 0.5*(ρuL^2)/ρL);
		ρR = QR[1];
		ρuR = QR[2];
		ER = QR[3];
		PR = (γ - 1)*(ER - 0.5*(ρuR^2)/ρR);
		uL = ρuL/ρL;
		uR = ρuR/ρR;
		aL = √(γ*PL/ρL);
		aR = √(γ*PR/ρR);
		λL = min(uL, uR) - max(aL, aR);
		λR = min(uL, uR) + max(aL, aR);
		return λL, λR
	end
	function sod_shock_init(ρ₀, p₀, u₀, dx, dt, tstop, n, γ, CFL, Bx₀, By₀, Bz₀)
		half = Int(n/2);
		ρ = zeros(1, n);
		ρ[1:half] = ρ₀[1]*ones(1, half);
		ρ[(half + 1):n] = ρ₀[2]*ones(1, half);
	
		u = zeros(1, n);
		u[1:half] = u₀[1]*ones(1, half);
		u[(half + 1):n] = u₀[2]*ones(1, half);
	
		p = zeros(1, n);
		p[1:half] = p₀[1]*ones(1, half);
		p[(half + 1):n] = p₀[2]*ones(1, half);
	
		ρu = @. ρ*u;
	
		neq = 3;
	
		E = @. p/(γ - 1) + 0.5*(ρu^2)/ρ;
		Q = zeros(neq, n);
	
		Q[1, :] = ρ;
		Q[2, :] = ρu;
		Q[3, :] = E;
	
		F = zeros(neq, n);
	
		for i = 1:n
			F[:, i] = eqnset(Q[:, i], γ);
		end
	
		# if OUTPUT_DATA == 1
		writedlm("sod_shock_initial_conditions_euler.csv", Q, ',')
		# end
	
		return Q, F, neq
	end
	```
	"""
elseif eqntype == 2
	md"""
	```julia
	function eqnset(Q,γ)
		ρ = Q[1];
		ρu = Q[2];
		ρv = Q[3];
		ρw = Q[4];
		Bx = Q[5];
		By = Q[6];
		Bz = Q[7];
		E = Q[8];
		B = √(Bx^2 + By^2 + Bz^2);
		P = (γ - 1)*(E - 0.5*(ρu^2 + ρv^2 + ρw^2)/ρ - B^2/2);
		F = [ρu; ρu^2/ρ - Bx^2 + P + B^2/2; ρu*ρv/ρ - Bx*By; ρu*ρw/ρ - Bx*Bz; 0; ρu*By/ρ - Bx*ρv/ρ; ρu*Bz/ρ - Bx*ρw/ρ; (E + P + B^2/2)*ρu/ρ - (Bx*ρu/ρ + By*ρv/ρ + Bz*ρw/ρ)*Bx];
		return F
	end
	function eig(QL, QR, γ)
		ρL = QL[1];
		BL = √(QL[5]^2 + QL[6]^2 + QL[7]^2);
		PL = (γ - 1)*(QL[8] - 0.5*(QL[2]^2 + QL[3]^2 + QL[4]^2)/QL[1]^2 - BL^2/2);
		ρR = QR[1];
		BR = √(QR[5]^2 + QR[6]^2 + QR[7]^2);
		PR = (γ - 1)*(QR[8] - 0.5*(QR[2]^2 + QR[3]^2 + QR[4]^2)/QR[1]^2 - BR^2/2);
	
		uL = QL[2]/ρL;
		uR = QR[2]/ρR;
		cL = √(γ*PL/ρL);
		cR = √(γ*PR/ρR);
		vL = QL[5]/√(ρL);
		vR = QR[5]/√(ρR);
		aL = √(0.5*(BL^2/ρL + cL^2 + √((BL^2/ρL + cL^2)^2 - 4*vL^2*cL^2)));
		aR = √(0.5*(BR^2/ρR + cR^2 + √((BR^2/ρR + cR^2)^2 - 4*vR^2*cR^2)));
		λL = min(uL, uR) - max(aL, aR);
		λR = min(uL, uR) + max(aL, aR);
		return λL, λR
	end
	function sound(Q, γ)
		ρ = Q[1, :];
		B = @. √(Q[5, :]^2 + Q[6, :]^2 + Q[7, :]^2);
		P = @. (γ - 1)*(Q[8, :] - 0.5*(Q[2, :]^2 + Q[3, :]^2 + Q[4, :]^2)/Q[1, :]^2 - B^2/2);
		c = @. √(γ*P/ρ);
		v = @. Q[5, :]/√(ρ);
		a = @. √(0.5*((B^2)/ρ + (c^2) + √((B^2)/ρ + (c^2))^2 - 4*(v^2)*(c^2)));
		return a
	end
	function sod_shock_init(ρ₀, p₀, u₀, dx, dt, tstop, n, γ, CFL, Bx₀, By₀, Bz₀)
		half = Int(n/2);
		ρ = zeros(1, n);
		ρ[1:half] = ρ₀[1]*ones(1, half);
		ρ[(half + 1):n] = ρ₀[2]*ones(1, half);
	
		u = zeros(1, n);
		u[1:half] = u₀[1]*ones(1, half);
		u[(half + 1):n] = u₀[2]*ones(1, half);
	
		p = zeros(1, n);
		p[1:half] = p₀[1]*ones(1, half);
		p[(half + 1):n] = p₀[2]*ones(1, half);
	
		ρu = @. ρ*u;
	
		neq = 8;
	
		Bx = zeros(1, n);
		Bx[1:half] = Bx₀[1]*ones(1, half);
		Bx[(half + 1):n] = Bx₀[2]*ones(1, half);
	
		By = zeros(1, n);
		By[1:half] = By₀[1]*ones(1, half);
		By[(half + 1):n] = By₀[2]*ones(1, half);
	
		Bz = zeros(1, n);
		Bz[1:half] = Bz₀[1]*ones(1, half);
		Bz[(half + 1):n] = Bz₀[2]*ones(1, half);
	
		B = @. √(Bx^2 + By^2 + Bz^2);
		E = @. p/(γ - 1) + 0.5*(ρu^2)/ρ + B^2/2;
		Q = zeros(neq, n);
	
		Q[1, :] = ρ;
		Q[2, :] = ρu;
		Q[3, :] = zeros(1, n);
		Q[4, :] = zeros(1, n);
		Q[5, :] = Bx;
		Q[6, :] = By;
		Q[7, :] = Bz;
		Q[8, :] = E;
	
		F = zeros(neq, n);
		for i = 1:n
			F[:, i] = eqnset(Q[:, i], γ);
		end
	
		# if OUTPUT_DATA == 1
			writedlm("../sod_shock_initial_conditions_mhd.csv", Q, ',')
		# end
	
		return Q, F, neq
	end
	```
	"""
end
end

# ╔═╡ b203ba78-47bb-444b-9bd4-7eece407a3ab
if  schemetype == 1
	function maccormack(Q, F, n, neq, γ, dt, dx)
	    vis = 1; # artificial viscosity factor
	    Q̄ = zeros(neq, n);
	    Qvis = zeros(neq, n);
	    F̄ = zeros(neq, n);
	    F̄[:, 1] = F[:, 1];
	    Q̄[:, 1] = Q[:, 1];
	    Q̄[:, n] = Q[:, n];
	    Qvis[:, 1] = Q[:, 1];
	    Qvis[:, n] = Q[:, n];
	    for i = 2:(n - 1)
	        #print(source((i - 1)*dx, ti, γ, y0, i))
	        Q̄[:, i] = Q[:, i] - (dt/dx)*(F[:, i + 1] - F[:, i]);
	        F̄[:, i] = eqnset(Q̄[:, i], γ);
	        Qvis[:, i] = 0.5*(Q[:, i] + Q̄[:, i]) - dt/(2*dx)*(F̄[:, i] - F̄[:, i - 1]);
	    end
	    for i = 2:(n - 1)
	        δ₁′ = norm(Qvis[:, i + 1] - Qvis[:, i])*(Qvis[:, i + 1] - Qvis[:, i]);
	        δ₂′ = norm(Qvis[:, i] - Qvis[:, i - 1])*(Qvis[:, i] - Qvis[:, i - 1]);
	        Q[:, i] = Qvis[:, i] + (vis*dt/dx)*(δ₁′ - δ₂′);
	    end
	    return Q
	end
	function solver(dx, dt, tstop, n, γ, CFL, neq, Q, F)
	    t = 0;
	    steps = 0;
	    while t < tstop
	        a = sound(Q, γ);
	        ua = [maximum(abs.(Q[2, :]./Q[1, :] + a)), maximum(abs.(Q[2, :]./Q[1, :] - a))];
	        
	        if maximum(ua)*dt/dx > CFL
	            dt = CFL*dx/maximum(ua);
	        end
	
	        Q = maccormack(Q, F, n, neq, γ, dt, dx);
	
	        for i = 1:n
	            F[:, i] = eqnset(Q[:, i], γ);
	        end
	
	        t += dt;
	        steps += 1;
	    end
	
	    # if OUTPUT_DATA == 1
	        writedlm("sod_shock_final_maccormack.csv", Q, ',')
	    # end
	
	    return Q, steps
	end
	
	if show_code == true
	md"""
	##### MacCormack Method
	"""
	end
	
elseif schemetype == 2
	function muscl(Q, n, limiter, γ, neq, eps, k, dt, dx)
	    Q = hcat(Q[:, 1], Q[:, 1], Q[:, 1], Q[:, :], Q[:, n], Q[:, n], Q[:, n]);
	    r⁻ = zeros(neq, n + 6);
	    r⁺ = zeros(neq, n + 6);
	    θ⁻ = zeros(neq, n + 6);
	    θ⁺ = zeros(neq, n + 6);
	    QL = zeros(neq, n + 6);
	    QR = zeros(neq, n + 6);
	    FR = zeros(neq, n + 6);
	    FL = zeros(neq, n + 6);
	    F = zeros(neq, n + 6);
	    for i = 2:(n + 4)
	        for j = 1:neq
	            if Q[j, i + 1] - Q[j, i] < 10^(-6)
	                r⁻[j, i] = (Q[j, i] - Q[j, i - 1])/(10^(-6));
	                r⁺[j, i] = (Q[j, i + 2] - Q[j, i - 1])/(10^(-6)); 
	            else
	                r⁻[j, i] = (Q[j, i] - Q[j, i - 1])/(Q[j, i + 1] - Q[j, i]);
	                r⁺[j, i] = (Q[j, i + 2] - Q[j, i - 1])/(Q[j, i + 1] - Q[j, i]);
	            end
	        end
	        if limiter == 1
	            for j = 1:neq
	                θ⁻[j, i] = max(0, min(1, r⁻[j, i]));
	                θ⁺[j, i] = max(0, min(1, r⁺[j, i]));
	            end
	        elseif limiter == 2
	            for j = 1:neq
	                θ⁻[j, i] = max(0, min(1, 2*r⁻[j, i]), min(2, r⁻[j, i]));
	                θ⁺[j, i] = max(0, min(1, 2*r⁺[j, i]), min(2, r⁺[j, i]));
	            end
	        elseif limiter == 3
	            for j = 1:neq
	                θ⁻[j, i] = (r⁻[j, i] + abs(r⁻[j, i]))/(1 + abs(r⁻[j, i]));
	                θ⁺[j, i] = (r⁺[j, i] + abs(r⁺[j, i]))/(1 + abs(r⁺[j, i]));
	            end
	        elseif limiter == 4
	            for j = 1:neq
	                θ⁻[j, i] = max(0, min(2*r⁻[j, i], (1 + r⁻[j, i])/2, 2));
	                θ⁺[j, i] = max(0, min(2*r⁺[j, i], (1 + r⁺[j, i])/2, 2));
	            end
	        end
	    end
	    for i = 3:(n + 3)
	        QL[:, i] = Q[:, i] + (ϵ/4)*((1 - k)*θ⁺[:, i - 1].*(Q[:, i] - Q[:, i - 1]) + (1 + k)*θ⁻[:, i].*(Q[:, i + 1] - Q[:, i]));
	        QR[:, i] = Q[:, i + 1] - (ϵ/4)*((1 + k)*θ⁺[:, i].*(Q[:, i + 1] - Q[:, i]) + (1 - k)*θ⁻[:, i + 1].*(Q[:, i + 2] - Q[:, i + 1]));
	        FL[:, i] = eqnset(QL[:, i], γ);
	        FR[:, i] = eqnset(QR[:, i], γ);
	    end
	    for i = 3:(n + 3)
	        if eqntype == 1
	            λL, λR = eig(QL[:, i], QR[:, i], γ);
	            if λL >= 0
	                F[:, i] = FL[:, i];
				elseif λR <= 0
	                F[:, i] = FR[:, i];
	            else
	                F[:, i] = (λR*FL[:, i] - λL*FR[:, i] + λL*λR*(QR[:, i] - QL[:, i]))/(λR - λL);
	            end
	        elseif eqntype == 2
	            if i == 3 || i == (n + 3)
	                QL2 = QL[:, i];
	                QR2 = QR[:, i];
	            else
	                QL2 = QL[:, i] + 0.5*(dt/dx)*(FR[:, i - 1] - FL[:, i]);
	                QR2 = QR[:, i] + 0.5*(dt/dx)*(FR[:, i] - FL[:, i + 1]);
	            end
	            λL, λR = eig(QL2, QR2, γ);
	            if λL >= 0
	                F[:, i] = FL[:, i];
				elseif λR <= 0
	                F[:, i] = FR[:, i];
	            else
	                F[:, i] = (λR*FL[:, i] - λL*FR[:, i] + λL*λR*(QR2 - QL2))/(λR - λL);
	            end
	        end
	    end
	    
	    for i = 4:(n + 3)
	        Q[:, i] = Q[:, i] - (dt/dx)*(F[:, i] - F[:, i - 1]);
	    end
	    Q = Q[:, 4:(n + 3)];
	    return Q
	end

	function solver(dx, dt, tstop, n, γ, CFL, neq, Q, F)
	    t = 0;
	    steps = 0;
	    while t < tstop
	        a = sound(Q, γ);
	        ua = [maximum(abs.(Q[2, :]./Q[1, :] + a)), maximum(abs.(Q[2, :]./Q[1, :] - a))];
	
	        if maximum(ua)*dt/dx > CFL
	            dt = CFL*dx/maximum(ua);
	        end
	
	        Q = muscl(Q, n, limiter, γ, neq, eps, k, dt, dx);
	
	        t += dt;
	        steps += 1;
	    end
	
	    # if OUTPUT_DATA == 1
	        writedlm("sod_shock_final_muscl.csv", Q, ',')
	    # end
	
	    return Q, steps
	end
	
	if show_code == true
	md"""
	##### MUSCL Scheme
	"""
	end
	
else error("Invalid value for `schemetype`")
end

# ╔═╡ 31467804-25e6-43b8-8a9f-63738f922988
if show_code == true
if  schemetype == 1
	md"""
	```julia
	function maccormack(Q, F, n, neq, γ, dt, dx)
	    vis = 1; # artificial viscosity factor
	    Qbar = zeros(neq, n);
	    Qvis = zeros(neq, n);
	    Fbar = zeros(neq, n);
	    F̄[:, 1] = F[:, 1];
	    Q̄[:, 1] = Q[:, 1];
	    Q̄[:, n] = Q[:, n];
	    Qvis[:, 1] = Q[:, 1];
	    Qvis[:, n] = Q[:, n];
	    for i = 2:(n - 1)
	        Q̄[:, i] = Q[:, i] - (dt/dx)*(F[:, i + 1] - F[:, i]);
	        F̄[:, i] = eqnset(Qbar[:, i], γ);
	        Qvis[:, i] = 0.5*(Q[:, i] + Q̄[:, i]) - dt/(2*dx)*(F̄[:, i] - F̄[:, i - 1]);
	    end
	    for i = 2:(n - 1)
	        δ₁′ = norm(Qvis[:, i + 1] - Qvis[:, i])*(Qvis[:, i + 1] - Qvis[:, i]);
	        δ₂′ = norm(Qvis[:, i] - Qvis[:, i - 1])*(Qvis[:, i] - Qvis[:, i - 1]);
	        Q[:, i] = Qvis[:, i] + (vis*dt/dx)*(δ₁′ - δ₂′);
	    end
	    return Q
	end
	function solver(dx, dt, tstop, n, γ, CFL, neq, Q, F)
	    t = 0;
	    steps = 0;
	    while t < tstop
	        a = sound(Q, γ);
	        ua = [@. maximum(abs(Q[2, :]/Q[1, :] + a)), @. maximum(abs(Q[2, :]/Q[1, :] - a))];
	        
	        if maximum(ua)*dt/dx > CFL
	            dt = CFL*dx/maximum(ua);
	        end
	
	        Q = maccormack(Q, F, n, neq, γ, dt, dx);
	
	        for i = 1:n
	            F[:, i] = eqnset(Q[:, i], γ);
	        end
	
	        t += dt;
	        steps += 1;
	    end
	
	    # if OUTPUT_DATA == 1
	        writedlm("sod_shock_final_maccormack.csv", Q, ',')
	    # end
	
	    return Q, steps
	end
	```
	"""
	
elseif schemetype == 2
md"""
	```julia
	function muscl(Q, n, limiter, γ, neq, eps, k, dt, dx)
	    Q = hcat(Q[:, 1], Q[:, 1], Q[:, 1], Q[:, :], Q[:, n], Q[:, n], Q[:, n]);
	    r⁻ = zeros(neq, n + 6);
	    r⁺ = zeros(neq, n + 6);
	    θ⁻ = zeros(neq, n + 6);
	    θ⁺ = zeros(neq, n + 6);
	    QL = zeros(neq, n + 6);
	    QR = zeros(neq, n + 6);
	    FR = zeros(neq, n + 6);
	    FL = zeros(neq, n + 6);
	    F = zeros(neq, n + 6);
	    for i = 2:(n + 4)
	        for j = 1:neq
	            if Q[j, i + 1] - Q[j, i] < 10^(-6)
	                r⁻[j, i] = (Q[j, i] - Q[j, i - 1])/(10^(-6));
	                r⁺[j, i] = (Q[j, i + 2] - Q[j, i - 1])/(10^(-6)); 
	            else
	                r⁻[j, i] = (Q[j, i] - Q[j, i - 1])/(Q[j, i + 1] - Q[j, i]);
	                r⁺[j, i] = (Q[j, i + 2] - Q[j, i - 1])/(Q[j, i + 1] - Q[j, i]);
	            end
	        end
	        if limiter == 1
	            for j = 1:neq
	                θ⁻[j, i] = max(0, min(1, r⁻[j, i]));
	                θ⁺[j, i] = max(0, min(1, r⁺[j, i]));
	            end
	        elseif limiter == 2
	            for j = 1:neq
	                θ⁻[j, i] = max(0, min(1, 2*r⁻[j, i]), min(2, r⁻[j, i]));
	                θ⁺[j, i] = max(0, min(1, 2*r⁺[j, i]), min(2, r⁺[j, i]));
	            end
	        elseif limiter == 3
	            for j = 1:neq
	                θ⁻[j, i] = (r⁻[j, i] + abs(r⁻[j, i]))/(1 + abs(r⁻[j, i]));
	                θ⁺[j, i] = (r⁺[j, i] + abs(r⁺[j, i]))/(1 + abs(r⁺[j, i]));
	            end
	        elseif limiter == 4
	            for j = 1:neq
	                θ⁻[j, i] = max(0, min(2*r⁻[j, i], (1 + r⁻[j, i])/2, 2));
	                θ⁺[j, i] = max(0, min(2*r⁺[j, i], (1 + r⁺[j, i])/2, 2));
	            end
	        end
	    end
	    for i = 3:(n + 3)
	        QL[:, i] = Q[:, i] + (ϵ/4)*((1 - k)*θ⁺[:, i - 1].*(Q[:, i] - Q[:, i - 1]) + (1 + k)*θ⁻[:, i].*(Q[:, i + 1] - Q[:, i]));
	        QR[:, i] = Q[:, i + 1] - (ϵ/4)*((1 + k)*θ⁺[:, i].*(Q[:, i + 1] - Q[:, i]) + (1 - k)*θ⁻[:, i + 1].*(Q[:, i + 2] - Q[:, i + 1]));
	        FL[:, i] = eqnset(QL[:, i], γ);
	        FR[:, i] = eqnset(QR[:, i], γ);
	    end
	    for i = 3:(n + 3)
	        if eqntype == 1
	            λL, λR = eig(QL[:, i], QR[:, i], γ);
	            if λL >= 0
	                F[:, i] = FL[:, i];
				elseif λR <= 0
	                F[:, i] = FR[:, i];
	            else
	                F[:, i] = (λR*FL[:, i] - λL*FR[:, i] + λL*λR*(QR[:, i] - QL[:, i]))/(λR - λL);
	            end
	        elseif eqntype == 2
	            if i == 3 || i == (n + 3)
	                QL2 = QL[:, i];
	                QR2 = QR[:, i];
	            else
	                QL2 = QL[:, i] + 0.5*(dt/dx)*(FR[:, i - 1] - FL[:, i]);
	                QR2 = QR[:, i] + 0.5*(dt/dx)*(FR[:, i] - FL[:, i + 1]);
	            end
	            λL, λR = eig(QL2, QR2, γ);
	            if λL >= 0
	                F[:, i] = FL[:, i];
				elseif λR <= 0
	                F[:, i] = FR[:, i];
	            else
	                F[:, i] = (λR*FL[:, i] - λL*FR[:, i] + λL*λR*(QR2 - QL2))/(λR - λL);
	            end
	        end
	    end
	    
	    for i = 4:(n + 3)
	        Q[:, i] = Q[:, i] - (dt/dx)*(F[:, i] - F[:, i - 1]);
	    end
	    Q = Q[:, 4:(n + 3)];
	    return Q
	end

	function solver(dx, dt, tstop, n, γ, CFL, neq, Q, F)
	    t = 0;
	    steps = 0;
	    while t < tstop
	        a = sound(Q, γ);
	        ua = [maximum(abs.(Q[2, :]./Q[1, :] + a)), maximum(abs.(Q[2, :]./Q[1, :] - a))];
	
	        if maximum(ua)*dt/dx > CFL
	            dt = CFL*dx/maximum(ua);
	        end
	
	        Q = muscl(Q, n, limiter, γ, neq, eps, k, dt, dx);
	
	        t += dt;
	        steps += 1;
	    end
	
	    # if OUTPUT_DATA == 1
	        writedlm("sod_shock_final_muscl.csv", Q, ',')
	    # end
	
	    return Q, steps
	end
	```
"""
	
else error("Invalid value for `schemetype`")
end
end

# ╔═╡ 0cf039a9-0282-4ca2-91cb-927342170963
# Initialization
Q₀, F₀, neq = sod_shock_init(ρ₀, p₀, u₀, dx, dt, tstop, n, γ, CFL, Bx₀, By₀, Bz₀);

# ╔═╡ 13ea456c-9a26-4594-b6d5-deb8388e3bd4
# Run simulation
Q, steps = solver(dx, dt, tstop, n, γ, CFL, neq, Q₀, F₀);

# ╔═╡ 423bb8d9-ba8f-4f96-856a-dd6f429ce015
plot(Q[1,:])

# ╔═╡ 5baffb4c-d770-4c09-af4a-05a2124d3a5c
plot(Q[2,:])

# ╔═╡ 608a04bf-95b4-4a59-b0b2-28bfa415d40d
plot(Q[3,:])

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Plots = "~1.38.0"
PlutoUI = "~0.7.49"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.4"
manifest_format = "2.0"
project_hash = "15b52356dc35ce293a12d063329e96e4d8a9cdc1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e7ff6cadf743c098e08fca25c91103ee4303c9bb"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.6"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random", "SnoopPrecompile"]
git-tree-sha1 = "aa3edc8f8dea6cbfa176ee12f7c2fc82f0608ed3"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.20.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "00a2cccc7f098ff3b66806862d275ca3db9e6e5a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.5.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "bcc737c4c3afc86f3bbc55eb1b9fabcee4ff2d81"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.71.2"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "64ef06fa8f814ff0d09ac31454f784c488e22b29"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.71.2+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "2e13c9956c82f5ae8cbdb8335327e63badb8c4ff"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.6.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "ab9aa169d2160129beb241cb2750ca499b4e90e9"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.17"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "946607f84feb96220f480e0422d3484c49c00239"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.19"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "df6830e37943c7aaa10023471ca47fb3065cc3c4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.2"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6e9dba33f9f2c44e08a020b0caf6903be540004"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.19+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "6466e524967496866901a78fca3f2e9ea445a559"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.2"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "5b7690dd212e026bbab1860016a6601cb077ab66"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.2"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SnoopPrecompile", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "513084afca53c9af3491c94224997768b9af37e8"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.38.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eadad7b14cf046de6eb41f13c9275e5aa2711ab6"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.49"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "18c35ed630d7229c5584b945641a73ca83fb5213"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.2"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase", "SnoopPrecompile"]
git-tree-sha1 = "e974477be88cb5e3040009f3767611bc6357846f"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.11"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SnoopPrecompile]]
git-tree-sha1 = "f604441450a3c0569830946e5b33b78c928e1a85"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "e4bdc63f5c6d62e80eb1c0043fcc0360d5950ff7"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.10"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "ac00576f90d8a259f2c9d823e91d1de3fd44d348"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╟─95aa5d17-a75b-48bc-af41-ab371db5e0d3
# ╟─195ce57e-a21f-4872-af6e-819d17105a36
# ╠═6f962a5b-9a2b-44b9-babe-4853c3eb7894
# ╟─14369438-0b1e-4e3e-9eef-809833a37b9a
# ╠═83595be7-bf28-4a08-891b-99ee1dacb7fb
# ╟─730195ac-480a-4d64-a3f9-a1612332e3e4
# ╠═9e542ab1-7ed9-490d-b4b8-7dea4db73401
# ╟─f4935bb2-dcbc-42a2-a564-d0d89a62b619
# ╟─ce5870f9-2d91-4f49-b98b-183517bbbe12
# ╟─220e7253-b2cc-4603-8a3b-f349d0727853
# ╟─3748e928-930d-4ecd-9f3e-b681b36e9d93
# ╟─b203ba78-47bb-444b-9bd4-7eece407a3ab
# ╟─31467804-25e6-43b8-8a9f-63738f922988
# ╟─6d03d72d-ac0a-4f68-aaac-5c902279b160
# ╟─e7cd3ffb-e168-40b7-b759-c764d8b0fe1b
# ╠═0cf039a9-0282-4ca2-91cb-927342170963
# ╠═13ea456c-9a26-4594-b6d5-deb8388e3bd4
# ╟─423bb8d9-ba8f-4f96-856a-dd6f429ce015
# ╟─5baffb4c-d770-4c09-af4a-05a2124d3a5c
# ╟─608a04bf-95b4-4a59-b0b2-28bfa415d40d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

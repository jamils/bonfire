module eqn

using LinearAlgebra
include("../input.jl")
# ===

"""
	eqnset(Q, γ)

Euler basic equation set 
"""
function eqnset(Q, γ)
	ρ = Q[1];
	ρu = Q[2];
	E = Q[3];
	P = (γ - 1)*(E - 0.5*(ρu^2)/ρ);
	F = [ρu; (ρu^2)/ρ + P; (ρu/ρ)*(E + P)];
	return F
end

"""
	sound(Q, γ)

Sound speed ``a``
"""
function sound(Q, γ)
	P = @. (γ - 1)*(Q[3, :] - 0.5*(Q[2, :]^2)/Q[1, :]);
	a = @. √(γ*P/Q[1, :]);
	return a
end

"""
	eig(Qᴸ, Qᴿ, γ) 

Eigen vectors
"""
function eig(Qᴸ, Qᴿ, γ)
	ρᴸ = Qᴸ[1];
	ρuᴸ = Qᴸ[2];
	Eᴸ = Qᴸ[3];
	Pᴸ = (γ - 1)*(Eᴸ - 0.5*(ρuᴸ^2)/ρᴸ);
	ρᴿ = Qᴿ[1];
	ρuᴿ = Qᴿ[2];
	Eᴿ = Qᴿ[3];
	Pᴿ = (γ - 1)*(Eᴿ - 0.5*(ρuᴿ^2)/ρᴿ);
	uᴸ = ρuᴸ/ρᴸ;
	uᴿ = ρuᴿ/ρᴿ;
	aᴸ = √(γ*Pᴸ/ρᴸ);
	aᴿ = √(γ*Pᴿ/ρᴿ);
	λᴸ = min(uᴸ, uᴿ) - max(aᴸ, aᴿ);
	λᴿ = min(uᴸ, uᴿ) + max(aᴸ, aᴿ);
	return λᴸ, λᴿ
end

"""
	SodShockInit(ρ₀, p₀, u₀, dx, dt, tstop, n, γ, CFᴸ, Bx₀, By₀, Bz₀)

Initizalizes the sod shock setup problem
"""
function SodShockInit(ρ₀, p₀, u₀, dx, dt, tstop, n, γ, CFᴸ, Bx₀, By₀, Bz₀)
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

    return Q, F, neq
end

# ===
end
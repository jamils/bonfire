module eqn

using LinearAlgebra
include("../input.jl")
# ===

"""
    eqnset(Q, γ)

MHD Equations
"""
function eqnset(Q, γ)
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

"""
    eig(Qᴸ, Qᴿ, γ)

Computer eigen vectors
"""
function eig(Qᴸ, Qᴿ, γ)
    ρᴸ = Qᴸ[1];
    Bᴸ = √(Qᴸ[5]^2 + Qᴸ[6]^2 + Qᴸ[7]^2);
    Pᴸ = (γ - 1)*(Qᴸ[8] - 0.5*(Qᴸ[2]^2 + Qᴸ[3]^2 + Qᴸ[4]^2)/Qᴸ[1]^2 - Bᴸ^2/2);
    ρᴿ = Qᴿ[1];
    Bᴿ = √(Qᴿ[5]^2 + Qᴿ[6]^2 + Qᴿ[7]^2);
    Pᴿ = (γ - 1)*(Qᴿ[8] - 0.5*(Qᴿ[2]^2 + Qᴿ[3]^2 + Qᴿ[4]^2)/Qᴿ[1]^2 - Bᴿ^2/2);

    uᴸ = Qᴸ[2]/ρᴸ;
    uᴿ = Qᴿ[2]/ρᴿ;
    cᴸ = √abs(γ*Pᴸ/ρᴸ);
    cᴿ = √(γ*Pᴿ/ρᴿ);
    vᴸ = Qᴸ[5]/√(ρᴸ);
    vᴿ = Qᴿ[5]/√(ρᴿ);
    aᴸ = √(0.5*(Bᴸ^2/ρᴸ + cᴸ^2 + √((Bᴸ^2/ρᴸ + cᴸ^2)^2 - 4*vᴸ^2*cᴸ^2)));
    aᴿ = √(0.5*(Bᴿ^2/ρᴿ + cᴿ^2 + √((Bᴿ^2/ρᴿ + cᴿ^2)^2 - 4*vᴿ^2*cᴿ^2)));
    λᴸ = min(uᴸ, uᴿ) - max(aᴸ, aᴿ);
    λᴿ = min(uᴸ, uᴿ) + max(aᴸ, aᴿ);
    return λᴸ, λᴿ
end

"""
    sound(Q, γ)

Compute MHD sound speed
"""
function sound(Q, γ)
    ρ = Q[1, :];
    B = @. √(Q[5, :]^2 + Q[6, :]^2 + Q[7, :]^2);
    P = @. (γ - 1)*(Q[8, :] - 0.5*(Q[2, :]^2 + Q[3, :]^2 + Q[4, :]^2)/Q[1, :]^2 - B^2/2);
    c = @. √(γ*P/ρ);
    v = @. Q[5, :]/√(ρ);
    a = @. √abs(0.5*((B^2)/ρ + (c^2) + √((B^2)/ρ + (c^2))^2 - 4*(v^2)*(c^2)));
    return a
end

"""
    SodShockInit(ρ₀, p₀, u₀, dx, dt, tstop, n, γ, CFL, Bx₀, By₀, Bz₀)

Initizalizes the sod shock setup problem
"""
function SodShockInit(ρ₀, p₀, u₀, dx, dt, tstop, n, γ, CFL, Bx₀, By₀, Bz₀)
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

    return Q, F, neq
end

# ===
end
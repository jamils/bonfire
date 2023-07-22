module scheme

using LinearAlgebra, SharedArrays, ProgressMeter
printstyled("Loading ", color=:green)
printstyled("MUSCL ", bold=true, color=:blue)
print("module... \n")

include("../input.jl")
# Select equation set
if EquationSet == 1
    include("../eqn/euler-module.jl")
elseif EquationSet == 2
    include("../eqn/mhd-module.jl")
else
    error("Invalid `EquationSet` value")
end

global eqnset = eqn.eqnset;
global sound = eqn.sound;
global eig = eqn.eig;
global SodShockInit = eqn.SodShockInit;
# ===

"""
    muscl(Q, n, Limiter, γ, neq, ϵ, k, dt, dx)

MUSCL discretization method
"""
function muscl(Q, n, Limiter, γ, neq, ϵ, k, dt, dx) # was a "eps" value, changed to ϵ
    Q = hcat(Q[:, 1], Q[:, 1], Q[:, 1], Q[:, :], Q[:, n], Q[:, n], Q[:, n]);
    r⁻ = zeros(neq, n + 6);
    r⁺ = zeros(neq, n + 6);
    θ⁻ = zeros(neq, n + 6);
    θ⁺ = zeros(neq, n + 6);
    Qᴸ = zeros(neq, n + 6);
    Qᴿ = zeros(neq, n + 6);
    Fᴿ = zeros(neq, n + 6);
    Fᴸ = zeros(neq, n + 6);
    F = zeros(neq, n + 6);
    for i ∈ 2:(n + 4)
        for j ∈ 1:neq
            if Q[j, i + 1] - Q[j, i] < 10^(-6)
                r⁻[j, i] = (Q[j, i] - Q[j, i - 1])/(10^(-6));
                r⁺[j, i] = (Q[j, i + 2] - Q[j, i - 1])/(10^(-6)); 
            else
                r⁻[j, i] = (Q[j, i] - Q[j, i - 1])/(Q[j, i + 1] - Q[j, i]);
                r⁺[j, i] = (Q[j, i + 2] - Q[j, i - 1])/(Q[j, i + 1] - Q[j, i]);
            end
        end
        if Limiter == 1
            for j ∈ 1:neq
                θ⁻[j, i] = max(0, min(1, r⁻[j, i]));
                θ⁺[j, i] = max(0, min(1, r⁺[j, i]));
            end
        elseif Limiter == 2
            for j ∈ 1:neq
                θ⁻[j, i] = max(0, min(1, 2*r⁻[j, i]), min(2, r⁻[j, i]));
                θ⁺[j, i] = max(0, min(1, 2*r⁺[j, i]), min(2, r⁺[j, i]));
            end
        elseif Limiter == 3
            for j ∈ 1:neq
                θ⁻[j, i] = (r⁻[j, i] + abs(r⁻[j, i]))/(1 + abs(r⁻[j, i]));
                θ⁺[j, i] = (r⁺[j, i] + abs(r⁺[j, i]))/(1 + abs(r⁺[j, i]));
            end
        elseif Limiter == 4
            for j ∈ 1:neq
                θ⁻[j, i] = max(0, min(2*r⁻[j, i], (1 + r⁻[j, i])/2, 2));
                θ⁺[j, i] = max(0, min(2*r⁺[j, i], (1 + r⁺[j, i])/2, 2));
            end
        end
    end
    for i ∈ 3:(n + 3)
        Qᴸ[:, i] = Q[:, i] + (ϵ/4)*((1 - k)*θ⁺[:, i - 1].*(Q[:, i] - Q[:, i - 1]) + (1 + k)*θ⁻[:, i].*(Q[:, i + 1] - Q[:, i]));
        Qᴿ[:, i] = Q[:, i + 1] - (ϵ/4)*((1 + k)*θ⁺[:, i].*(Q[:, i + 1] - Q[:, i]) + (1 - k)*θ⁻[:, i + 1].*(Q[:, i + 2] - Q[:, i + 1]));
        Fᴸ[:, i] = eqnset(Qᴸ[:, i], γ);
        Fᴿ[:, i] = eqnset(Qᴿ[:, i], γ);
    end
    for i ∈ 3:(n + 3)
        if EquationSet == 1
            λᴸ, λᴿ = eig(Qᴸ[:, i], Qᴿ[:, i], γ);
            if λᴸ >= 0
                F[:, i] = Fᴸ[:, i];
            elseif λᴿ <= 0
                F[:, i] = Fᴿ[:, i];
            else
                F[:, i] = (λᴿ*Fᴸ[:, i] - λᴸ*Fᴿ[:, i] + λᴸ*λᴿ*(Qᴿ[:, i] - Qᴸ[:, i]))/(λᴿ - λᴸ);
            end
        elseif EquationSet == 2
            if i == 3 || i == (n + 3)
                Qᴸ2 = Qᴸ[:, i];
                Qᴿ2 = Qᴿ[:, i];
            else
                Qᴸ2 = Qᴸ[:, i] + 0.5*(dt/dx)*(Fᴿ[:, i - 1] - Fᴸ[:, i]);
                Qᴿ2 = Qᴿ[:, i] + 0.5*(dt/dx)*(Fᴿ[:, i] - Fᴸ[:, i + 1]);
            end
            λᴸ, λᴿ = eig(Qᴸ2, Qᴿ2, γ);
            if λᴸ >= 0
                F[:, i] = Fᴸ[:, i];
            elseif λᴿ <= 0
                F[:, i] = Fᴿ[:, i];
            else
                F[:, i] = (λᴿ*Fᴸ[:, i] - λᴸ*Fᴿ[:, i] + λᴸ*λᴿ*(Qᴿ2 - Qᴸ2))/(λᴿ - λᴸ);
            end
        end
    end
    
    for i ∈ 4:(n + 3)
        Q[:, i] = Q[:, i] - (dt/dx)*(F[:, i] - F[:, i - 1]);
    end

    Q = Q[:, 4:(n + 3)];
    return Q
end

"""
    solver(dx, dt, tstop, n, γ, CFL, neq, Q, F)

Actual computation function. Iterates through each
timestep and computes the values for the full grid
"""
function solver(dx, dt, tstop, n, γ, CFL, neq, Q, F)
    t = 0;
    steps = 0;

    # prog = Progress(617, 1.0)
    prog = ProgressUnknown("Running main loop ", spinner=true)
    while t < tstop
        ProgressMeter.next!(prog; showvalues=[(:t, t)], spinner="⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏")
        a = sound(Q, γ);
        ua = [maximum(abs.(Q[2, :]./Q[1, :] + a)), maximum(abs.(Q[2, :]./Q[1, :] - a))];

        if maximum(ua)*dt/dx > CFL
            dt = CFL*dx/maximum(ua);
        end
        # print("$dt, ")

        Q = muscl(Q, n, Limiter, γ, neq, ϵ, k, dt, dx);

        t += dt;
        steps += 1;
        # update!(prog, steps)
    end
    ProgressMeter.finish!(prog)

    return Q, steps
end

# ===
end
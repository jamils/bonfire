module scheme

using LinearAlgebra, LoopVectorization
printstyled("Loading ", color=:green)
printstyled("MacCormack ", bold=true, color=:blue)
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
global SodShockInit = eqn.SodShockInit;
# ===

"""
    maccormack(Q, F, n, neq, γ, dt, dx)

MacCormack stop method
"""
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
    for i ∈ 2:(n - 1)
        Q̄[:, i] = Q[:, i] - (dt/dx)*(F[:, i + 1] - F[:, i]);
        F̄[:, i] = eqnset(Q̄[:, i], γ);
        Qvis[:, i] = 0.5*(Q[:, i] + Q̄[:, i]) - dt/(2*dx)*(F̄[:, i] - F̄[:, i - 1]);
    end
    for i ∈ 2:(n - 1)
        δ₁′ = norm(Qvis[:, i + 1] - Qvis[:, i])*(Qvis[:, i + 1] - Qvis[:, i]);
        δ₂′ = norm(Qvis[:, i] - Qvis[:, i - 1])*(Qvis[:, i] - Qvis[:, i - 1]);
        Q[:, i] = Qvis[:, i] + (vis*dt/dx)*(δ₁′ - δ₂′);
    end
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

    prog = ProgressUnknown("Running main loop ", spinner=true)
    while t < tstop
        ProgressMeter.next!(prog; showvalues=[(:t, t)], spinner="⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏")
        a = sound(Q, γ);
        ua = [maximum(abs.(Q[2, :]./Q[1, :] + a)), maximum(abs.(Q[2, :]./Q[1, :] - a))];
        
        if maximum(ua)*dt/dx > CFL
            dt = CFL*dx/maximum(ua);
        end

        Q = maccormack(Q, F, n, neq, γ, dt, dx);

        for i ∈ 1:n
            F[:, i] = eqnset(Q[:, i], γ);
        end

        t += dt;
        steps += 1;
    end

    return Q, steps
end
ProgressMeter.finish!(prog)

# ===
end
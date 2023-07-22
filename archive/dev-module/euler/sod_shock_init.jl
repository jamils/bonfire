function sod_shock_init(rho0, p0, u0, dx, dt, tstop, n, γ, CFL, Bx0, By0, Bz0)
    half = Int(n/2);
    rho = zeros(1, n);
    rho[1:half] = rho0[1]*ones(1, half);
    rho[(half + 1):n] = rho0[2]*ones(1, half);

    u = zeros(1, n);
    u[1:half] = u0[1]*ones(1, half);
    u[(half + 1):n] = u0[2]*ones(1, half);

    p = zeros(1, n);
    p[1:half] = p0[1]*ones(1, half);
    p[(half + 1):n] = p0[2]*ones(1, half);

    rhou = rho.*u;

    neq = 3;

    E = p/(γ - 1) + 0.5*(rhou.^2)./rho;
    Q = zeros(neq, n);

    Q[1, :] = rho;
    Q[2, :] = rhou;
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
# Sod shock simulation code #
# Kolter Bradshaw #

function sod_shock_init(rho0, p0, u0, dx, dt, tstop, n, gamma, CFL; Bx0=nothing, By0=nothing, Bz0=nothing)

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

    if EULER == true
        neq = 3;

        E = p/(gamma - 1) + 0.5*(rhou.^2)./rho;
        Q = zeros(neq, n);

        Q[1, :] = rho;
        Q[2, :] = rhou;
        Q[3, :] = E;

        F = zeros(neq, n);
        for i = 1:n
            F[:, i] = eulereq(Q[:, i], gamma);
        end

    elseif MHD == true
        neq = 8;

        Bx = zeros(1, n);
        Bx[1:half] = Bx0[1]*ones(1, half);
        Bx[(half + 1):n] = Bx0[2]*ones(1, half);

        By = zeros(1, n);
        By[1:half] = By0[1]*ones(1, half);
        By[(half + 1):n] = By0[2]*ones(1, half);

        Bz = zeros(1, n);
        Bz[1:half] = Bz0[1]*ones(1, half);
        Bz[(half + 1):n] = Bz0[2]*ones(1, half);

        B = sqrt.(Bx.^2 + By.^2 + Bz.^2);
        E = p/(gamma - 1) + 0.5*(rhou.^2)./rho + B.^2/2;
        Q = zeros(neq, n);

        Q[1, :] = rho;
        Q[2, :] = rhou;
        Q[3, :] = zeros(1, n);
        Q[4, :] = zeros(1, n);
        Q[5, :] = Bx;
        Q[6, :] = By;
        Q[7, :] = Bz;
        Q[8, :] = E;

        F = zeros(neq, n);
        for i = 1:n
            F[:, i] = mhd(Q[:, i], gamma);
        end
    end

    if OUTPUT_DATA == true
        if EULER == true
            writedlm("sod_shock_initial_conditions_euler.csv", Q, ',')
        elseif MHD == true
            writedlm("sod_shock_initial_conditions_mhd.csv", Q, ',')
        end
    end
    return Q, F, neq
end
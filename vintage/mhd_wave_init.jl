# Sod shock simulation code #
# Kolter Bradshaw #

function mhd_wave_init(A, Lx, rho0, p0, u0, v0, w0, dx, dt, tstop, n, gamma, CFL, Bx0, By0, Bz0)

    # MHD eigenvectors
    if By0^2 + Bz0^2 == 0
        betay = 1/sqrt(2);
        betaz = 1/sqrt(2);
    else
        betay = By0/sqrt(By0^2 + Bz0^2);
        betaz = Bz0/sqrt(By0^2 + Bz0^2);
    end
    if Bx0 >= 0
        sgnBx = 1;
    else
        sgnBx = -1;
    end
    R = [0; 0; -betaz*sgnBx; betay*sgnBx; 0; betaz/sqrt(rho0); -betay/sqrt(rho0); -sgnBx*(betaz*v0 - betay*w0)];

    rho = zeros(1, n);
    for i = 1:n
        x = (i - 1)*dx;
        rho[i] = rho0 + A*R[1]*sin(2*pi*x/Lx);
    end

    p = zeros(1, n);
    for i = 1:n
        x = (i - 1)*dx;
        p[i] = p0 + A*R[8]*sin(2*pi*x/Lx);
    end

    u = zeros(1, n);
    for i = 1:n
        x = (i - 1)*dx;
        u[i] = u0 + A*R[2]*sin(2*pi*x/Lx);
    end

    v = zeros(1, n);
    for i = 1:n
        x = (i - 1)*dx;
        v[i] = v0 + A*R[3]*sin(2*pi*x/Lx);
    end

    w = zeros(1, n);
    for i = 1:n
        x = (i - 1)*dx;
        w[i] = w0 + A*R[4]*sin(2*pi*x/Lx);
    end

    rhou = rho.*u;
    rhov = rho.*v;
    rhow = rho.*w;

    neq = 8;

    Bx = zeros(1, n);
    for i = 1:n
        x = (i - 1)*dx;
        Bx[i] = Bx0 + A*R[5]*sin(2*pi*x/Lx);
    end

    By = zeros(1, n);
    for i = 1:n
        x = (i - 1)*dx;
        By[i] = By0 + A*R[6]*sin(2*pi*x/Lx);
    end

    Bz = zeros(1, n);
    for i = 1:n
        x = (i - 1)*dx;
        Bz[i] = Bz0 + A*R[7]*sin(2*pi*x/Lx);
    end

    B = sqrt.(Bx.^2 + By.^2 + Bz.^2);
    E = p/(gamma - 1) + 0.5*(rhou.^2)./rho + B.^2/2;
    Q = zeros(neq, n);

    Q[1, :] = rho;
    Q[2, :] = rhou;
    Q[3, :] = rhov;
    Q[4, :] = rhow;
    Q[5, :] = Bx;
    Q[6, :] = By;
    Q[7, :] = Bz;
    Q[8, :] = E;

    F = zeros(neq, n);
    for i = 1:n
        F[:, i] = mhd(Q[:, i], gamma);
    end

    if OUTPUT_DATA == true
        if EULER == true
            writedlm("mhd_wave_initial_conditions_euler.csv", Q, ',')
        elseif MHD == true
            writedlm("mhd_wave_initial_conditions_mhd.csv", Q, ',')
        end
    end
    return Q, F, neq
end
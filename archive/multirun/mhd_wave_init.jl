# Sod shock simulation code #
# Kolter Bradshaw #

function mhd_wave_init(A, Lx, rho0, p0, u0, v0, w0, dx, dt, tstop, n, gamma, CFL, Bx0, By0, Bz0, loopnum)

    neq = 8;

    # MHD eigenvectors
    E0 = p0/(gamma - 1) + 0.5*rho0*u0^2 + (Bx0^2 + By0^2 + Bz0^2)/2;
    a = sqrt(p0*rho0 + (gamma - 1)*rho0*p0/rho0^2);
    b = sqrt(Bx0^2 + By0^2 + Bz0^2)/sqrt(rho0);
    cA = Bx0/sqrt(rho0);
    cf = sqrt(0.5*((a^2 + b^2) + sqrt((a^2 + b^2)^2 - 4*a^2*Bx0^2/rho0)));
    cs = sqrt(0.5*((a^2 + b^2) - sqrt((a^2 + b^2)^2 - 4*a^2*Bx0^2/rho0)));
    hstar = (E0 + p0 + (Bx0^2 + By0^2 + Bz0^2)/2)/rho0;
    if Bx0 >= 0
        sgnBx = 1;
    else
        sgnBx = -1;
    end
    if By0 >= 0
        sgnBy = 1;
    else
        sgnBy = -1;
    end
    if By0^2 + Bz0^2 == 0
        betay = 1/sqrt(2);
        betaz = 1/sqrt(2);
        if cA == a
            alphaf = 1/sqrt(2);
            alphas = 1/sqrt(2);
        else
            alphaf = sqrt(a^2 - cs^2)/sqrt(cf^2 - cs^2);
            alphas = sqrt(cf^2 - a^2)/sqrt(cf^2 - cs^2);
        end
    else
        betay = By0/sqrt(By0^2 + Bz0^2);
        betaz = Bz0/sqrt(By0^2 + Bz0^2);
        alphaf = sqrt(a^2 - cs^2)/sqrt(cf^2 - cs^2);
        alphas = sqrt(cf^2 - a^2)/sqrt(cf^2 - cs^2);
    end
    if a < cA
        alpha1 = alphaf*sgnBy;
        alpha2 = alphas;
        alphabar1 = alphas*sgnBy;
        alphabar2 = alphaf;
    else
        alpha1 = alphaf;
        alpha2 = alphas*sgnBy;
        alphabar1 = alphas;
        alphabar2 = alphaf*sgnBy;
    end
    if By0^2 + Bz0^2 == 0
        betay = 1/sqrt(2);
        betaz = 1/sqrt(2);
    else
        betay = By0/sqrt(By0^2 + Bz0^2);
        betaz = Bz0/sqrt(By0^2 + Bz0^2);
    end
    c1 = -cf;
    cbar1 = -cs;
    c2 = -cs;
    cbar2 = -cf;
    if c1^2 - a^2 >= 0
        sgndiff1 = 1;
    else
        sgndiff1 = -1;
    end
    if c2^2 - a^2 >= 0
        sgndiff2 = 1;
    else
        sgndiff2 = -1;
    end
    R1 = [0; 0; -betaz*sgnBx; betay*sgnBx; 0; betaz/sqrt(rho0); -betay/sqrt(rho0); -sgnBx*(betaz*v0 - betay*w0)];
    # R1 = zeros(neq, 1);
    # R2 = [alpha1; alpha1*(u0 + c1); alpha1*v0 - alphabar1*cbar1*sgndiff1*sgnBx*betay; alpha1*w0 - alphabar1*cbar1*sgndiff1*sgnBx*betaz; 0; alphabar1*a*sgndiff1*betay/sqrt(rho0); alphabar1*a*sgndiff1*betaz/sqrt(rho0); alpha1*(hstar - a^2 - b^2 + c1^2 + u0*c1) - sgndiff1*alphabar1*cbar1*sgnBx*(v0*betay + w0*betaz)];
    R2 = zeros(neq, 1);
    # R3 = [alpha2; alpha2*(u0 + c2); alpha2*v0 - alphabar2*cbar2*sgndiff2*sgnBx*betay; alpha2*w0 - alphabar2*cbar2*sgndiff2*sgnBx*betaz; 0; alphabar2*a*sgndiff2*betay/sqrt(rho0); alphabar2*a*sgndiff2*betaz/sqrt(rho0); alpha2*(hstar - a^2 - b^2 + c2^2 + u0*c2) - sgndiff2*alphabar2*cbar2*sgnBx*(v0*betay + w0*betaz)];;
    R3 = zeros(neq, 1);

    rho = zeros(1, n);
    for i = 1:n
        x = (i - 1)*dx;
        rho[i] = rho0 + A*(R1[1] + R2[1] + R3[1])*sin(2*pi*x/Lx);
    end

    p = zeros(1, n);
    for i = 1:n
        x = (i - 1)*dx;
        p[i] = p0 + A*(R1[8] + R2[8] + R3[8])*sin(2*pi*x/Lx);
    end

    u = zeros(1, n);
    for i = 1:n
        x = (i - 1)*dx;
        u[i] = u0 + A*(R1[2] + R2[2] + R3[2])*sin(2*pi*x/Lx);
    end

    v = zeros(1, n);
    for i = 1:n
        x = (i - 1)*dx;
        v[i] = v0 + A*(R1[3] + R2[3] + R3[3])*sin(2*pi*x/Lx);
    end

    w = zeros(1, n);
    for i = 1:n
        x = (i - 1)*dx;
        w[i] = w0 + A*(R1[4] + R2[4] + R3[4])*sin(2*pi*x/Lx);
    end

    rhou = rho.*u;
    rhov = rho.*v;
    rhow = rho.*w;

    Bx = zeros(1, n);
    for i = 1:n
        x = (i - 1)*dx;
        Bx[i] = Bx0 + A*(R1[5] + R2[5] + R3[5])*sin(2*pi*x/Lx);
    end

    By = zeros(1, n);
    for i = 1:n
        x = (i - 1)*dx;
        By[i] = By0 + A*(R1[6] + R2[6] + R3[6])*sin(2*pi*x/Lx);
    end

    Bz = zeros(1, n);
    for i = 1:n
        x = (i - 1)*dx;
        Bz[i] = Bz0 + A*(R1[7] + R2[7] + R3[7])*sin(2*pi*x/Lx);
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
            writedlm("mhd_wave_initial_conditions_", eqnstring, "_", loopnum, ".csv", Q, ',')
        elseif MHD == true
            title = string("mhd_wave_initial_conditions_", eqnstring, "_", loopnum, ".csv");
            writedlm(title, Q, ',')
        end
    end
    return Q, F, neq
end
# Sod shock simulation code #
# Kolter Bradshaw #

function zpinch_init(j0, rho0, r0, u0, dx, dt, tstop, n, gamma, CFL)

    u = zeros(1, n);

    neq = 8;

    Bx = zeros(1, n);
    for i = 1:n
        x = (i - 1)*dx;
        if x <= r0
            Bx[i] =12.7324*x;
        else
            Bx[i] = 0.795775/x;
        end
    end
    By = zeros(1, n);
    Bz = zeros(1, n);

    p = zeros(1, n);
    for i = 1:n
        x = (i - 1)*dx;
        if x <= r0
            p[i] = 10.1321 - 162.114*x^2;
        else
            p[i] = 0;
        end
    end

    rho = zeros(1, n);
    for i = 1:n
        x = (i - 1)*dx;
        if x < r0
            rho[i] = p[i]/p[1];
        else
            rho[i] = 10^(-6);
        end
    end

    rhou = rho.*u;

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

    if OUTPUT_DATA == true
        if EULER == true
            writedlm("z-pinch_initial_conditions_euler.csv", Q, ',')
        elseif MHD == true
            writedlm("z-pinch_initial_conditions_mhd.csv", Q, ',')
        end
    end
    return Q, F, neq
end
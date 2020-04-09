# MacCormack numerical scheme #
# Kolter Bradshaw #

function maccormack(Q, F, n, neq, gamma, dt, dx)
    vis = 1; # artificial viscosity factor
    Qbar = zeros(neq, n);
    Qvis = zeros(neq, n);
    Fbar = zeros(neq, n);
    Fbar[:, 1] = F[:, 1];
    Qbar[:, 1] = Q[:, 1];
    Qbar[:, n] = Q[:, n];
    Qvis[:, 1] = Q[:, 1];
    Qvis[:, n] = Q[:, n];
    for i = 2:(n - 1)
        #print(source((i - 1)*dx, ti, gamma, y0, i))
        Qbar[:, i] = Q[:, i] - (dt/dx)*(F[:, i + 1] - F[:, i]);
        if EULER == true
            Fbar[:, i] = eulereq(Qbar[:, i], gamma);
        elseif MHD == true
            Fbar[:, i] = mhd(Qbar[:, i], gamma);
        end
        Qvis[:, i] = 0.5*(Q[:, i] + Qbar[:, i]) - dt/(2*dx)*(Fbar[:, i] - Fbar[:, i - 1]);
    end
    for i = 2:(n - 1)
        delprime1 = norm(Qvis[:, i + 1] - Qvis[:, i])*(Qvis[:, i + 1] - Qvis[:, i]);
        delprime2 = norm(Qvis[:, i] - Qvis[:, i - 1])*(Qvis[:, i] - Qvis[:, i - 1]);
        Q[:, i] = Qvis[:, i] + (vis*dt/dx)*(delprime1 - delprime2);
    end
    return Q
end
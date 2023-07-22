function eqnset(Q,γ)
    rho = Q[1];
    rhou = Q[2];
    rhov = Q[3];
    rhow = Q[4];
    Bx = Q[5];
    By = Q[6];
    Bz = Q[7];
    E = Q[8];
    B = sqrt(Bx^2 + By^2 + Bz^2);
    P = (γ - 1)*(E - 0.5*(rhou^2 + rhov^2 + rhow^2)/rho - B^2/2);
    F = [rhou; rhou^2/rho - Bx^2 + P + B^2/2; rhou*rhov/rho - Bx*By; rhou*rhow/rho - Bx*Bz; 0; rhou*By/rho - Bx*rhov/rho; rhou*Bz/rho - Bx*rhow/rho; (E + P + B^2/2)*rhou/rho - (Bx*rhou/rho + By*rhov/rho + Bz*rhow/rho)*Bx];
    return F
end
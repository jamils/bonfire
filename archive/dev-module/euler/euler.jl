function eqnset(Q, γ)
    rho = Q[1];
    rhou = Q[2];
    E = Q[3];
    P = (γ - 1)*(E - 0.5*(rhou^2)/rho);
    F = [rhou; (rhou^2)/rho + P; (rhou/rho)*(E + P)];
    return F
end
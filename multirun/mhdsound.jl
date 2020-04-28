# MHD sound speed #
# Kolter Bradshaw #

function mhdsound(Q, gamma)
    rho = Q[1, :];
    B = sqrt.(Q[5, :].^2 + Q[6, :].^2 + Q[7, :].^2);
    P = (gamma - 1)*(Q[8, :] - 0.5*(Q[2, :].^2 + Q[3, :].^2 + Q[4, :].^2)./Q[1, :].^2 - B.^2/2);
    c = sqrt.(gamma*P./rho);
    v = Q[5, :]./sqrt.(rho);
    a = sqrt.(0.5*((B.^2)./rho + (c.^2) + sqrt.(((B.^2)./rho + (c.^2)).^2 - 4*(v.^2).*(c.^2))));
    return a
end
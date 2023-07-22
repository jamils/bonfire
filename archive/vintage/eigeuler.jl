# Euler equation eigenvalues #
# Kolter Bradshaw #

function eigeuler(QL, QR, gamma)
    rhoL = QL[1];
    rhouL = QL[2];
    EL = QL[3];
    PL = (gamma - 1)*(EL - 0.5*(rhouL^2)/rhoL);
    rhoR = QR[1];
    rhouR = QR[2];
    ER = QR[3];
    PR = (gamma - 1)*(ER - 0.5*(rhouR^2)/rhoR);
    uL = rhouL/rhoL;
    uR = rhouR/rhoR;
    aL = sqrt(gamma*PL/rhoL);
    aR = sqrt(gamma*PR/rhoR);
    lambdaL = min(uL, uR) - max(aL, aR);
    lambdaR = min(uL, uR) + max(aL, aR);
    return lambdaL, lambdaR
end
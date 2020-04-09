function eig(QL, QR, γ)
    rhoL = QL[1];
    BL = sqrt(QL[5]^2 + QL[6]^2 + QL[7]^2);
    PL = (γ - 1)*(QL[8] - 0.5*(QL[2]^2 + QL[3]^2 + QL[4]^2)/QL[1]^2 - BL^2/2);
    rhoR = QR[1];
    BR = sqrt(QR[5]^2 + QR[6]^2 + QR[7]^2);
    PR = (γ - 1)*(QR[8] - 0.5*(QR[2]^2 + QR[3]^2 + QR[4]^2)/QR[1]^2 - BR^2/2);

    uL = QL[2]/rhoL;
    uR = QR[2]/rhoR;
    cL = sqrt(γ*PL/rhoL);
    cR = sqrt(γ*PR/rhoR);
    vL = QL[5]/sqrt(rhoL);
    vR = QR[5]/sqrt(rhoR);
    aL = sqrt(0.5*(BL^2/rhoL + cL^2 + sqrt((BL^2/rhoL + cL^2)^2 - 4*vL^2*cL^2)));
    aR = sqrt(0.5*(BR^2/rhoR + cR^2 + sqrt((BR^2/rhoR + cR^2)^2 - 4*vR^2*cR^2)));
    lambdaL = min(uL, uR) - max(aL, aR);
    lambdaR = min(uL, uR) + max(aL, aR);
    return lambdaL, lambdaR
end
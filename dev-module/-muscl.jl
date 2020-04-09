# MUSCL numerical scheme #
# Kolter Bradshaw #
function muscl(Q, n, limiter, gamma, neq, eps, k, dt, dx)
    Q = hcat(Q[:, 1], Q[:, 1], Q[:, 1], Q[:, :], Q[:, n], Q[:, n], Q[:, n]);
    rMinus = zeros(neq, n + 6);
    rPlus = zeros(neq, n + 6);
    thetaMinus = zeros(neq, n + 6);
    thetaPlus = zeros(neq, n + 6);
    QL = zeros(neq, n + 6);
    QR = zeros(neq, n + 6);
    FR = zeros(neq, n + 6);
    FL = zeros(neq, n + 6);
    F = zeros(neq, n + 6);
    for i = 2:(n + 4)
        for j = 1:neq
            if Q[j, i + 1] - Q[j, i] < 10^(-6)
                rMinus[j, i] = (Q[j, i] - Q[j, i - 1])/(10^(-6));
                rPlus[j, i] = (Q[j, i + 2] - Q[j, i - 1])/(10^(-6)); 
            else
                rMinus[j, i] = (Q[j, i] - Q[j, i - 1])/(Q[j, i + 1] - Q[j, i]);
                rPlus[j, i] = (Q[j, i + 2] - Q[j, i - 1])/(Q[j, i + 1] - Q[j, i]);
            end
        end
        if limiter == 1
            for j = 1:neq
                thetaMinus[j, i] = max(0, min(1, rMinus[j, i]));
                thetaPlus[j, i] = max(0, min(1, rPlus[j, i]));
            end
        elseif limiter == 2
            for j = 1:neq
                thetaMinus[j, i] = max(0, min(1, 2*rMinus[j,i]), min(2, rMinus[j,i]));
                thetaPlus[j, i]=max(0, min(1, 2*rPlus[j, i]), min(2, rPlus[j, i]));
            end
        elseif limiter == 3
            for j = 1:neq
                thetaMinus[j, i] = (rMinus[j, i] + abs(rMinus[j,i]))/(1 + abs(rMinus[j,i]));
                thetaPlus[j, i] = (rPlus[j, i] + abs(rPlus[j, i]))/(1 + abs(rPlus[j, i]));
            end
        elseif limiter == 4
            for j = 1:neq
                thetaMinus[j, i] = max(0, min(2*rMinus[j,i], (1 + rMinus[j,i])/2, 2));
                thetaPlus[j, i] = max(0, min(2*rPlus[j, i], (1 + rPlus[j, i])/2, 2));
            end
        end
    end
    for i = 3:(n + 3)
        QL[:, i] = Q[:, i] + (eps/4)*((1 - k)*thetaPlus[:, i - 1].*(Q[:, i] - Q[:, i - 1]) + (1 + k)*thetaMinus[:, i].*(Q[:, i + 1] - Q[:, i]));
        QR[:, i] = Q[:, i + 1] - (eps/4)*((1 + k)*thetaPlus[:, i].*(Q[:, i + 1] - Q[:, i]) + (1 - k)*thetaMinus[:, i + 1].*(Q[:, i + 2] - Q[:, i + 1]));
        if EULER == true
            FL[:, i] = eulereq(QL[:, i], gamma);
            FR[:, i] = eulereq(QR[:, i], gamma);
        elseif MHD == true
            FL[:, i] = mhd(QL[:, i], gamma);
            FR[:, i] = mhd(QR[:, i], gamma);
        end
    end
    for i = 3:(n + 3)
        if EULER == true
            lambdaL, lambdaR = eigeuler(QL[:, i], QR[:, i], gamma);
            if lambdaL >= 0
                F[:, i] = FL[:, i];
            elseif lambdaR <= 0
                F[:, i] = FR[:, i];
            else
                F[:, i] = (lambdaR*FL[:, i] - lambdaL*FR[:, i] + lambdaL*lambdaR*(QR[:, i] - QL[:, i]))/(lambdaR - lambdaL);
            end
        elseif MHD == true
            if i == 3 || i == (n + 3)
                QL2 = QL[:, i];
                QR2 = QR[:, i];
            else
                QL2 = QL[:, i] + 0.5*(dt/dx)*(FR[:, i - 1] - FL[:, i]);
                QR2 = QR[:, i] + 0.5*(dt/dx)*(FR[:, i] - FL[:, i + 1]);
            end
            lambdaL, lambdaR = eigmhd(QL2, QR2, gamma);
            if lambdaL >= 0
                F[:, i] = FL[:, i];
            elseif lambdaR <= 0
                F[:, i] = FR[:, i];
            else
                F[:, i] = (lambdaR*FL[:, i] - lambdaL*FR[:, i] + lambdaL*lambdaR*(QR2 - QL2))/(lambdaR - lambdaL);
            end
        end
    end
    
    for i = 4:(n + 3)
        Q[:, i] = Q[:, i] - (dt/dx)*(F[:, i] - F[:, i - 1]);
    end
    Q = Q[:, 4:(n + 3)];
    return Q
end
function sound(Q, γ)
    P = (γ - 1)*(Q[3, :] - 0.5*(Q[2, :].^2)./Q[1, :]);
    a = sqrt.(γ*P./Q[1, :]);

    return a
end
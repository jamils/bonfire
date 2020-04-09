function solver(dx, dt, tstop, n, γ, CFL, neq, Q, F)
    t = 0;
    steps = 0;
    while t < tstop
        a = sound(Q, γ);
        ua = [maximum(abs.(Q[2, :]./Q[1, :] + a)), maximum(abs.(Q[2, :]./Q[1, :] - a))];

        if maximum(ua)*dt/dx > CFL
            dt = CFL*dx/maximum(ua);
        end

        Q = muscl(Q, n, limiter, γ, neq, eps, k, dt, dx);

        t += dt;
        steps += 1;
    end

    # if OUTPUT_DATA == 1
        writedlm("sod_shock_final_muscl.csv", Q, ',')
    # end

    return Q, steps
end
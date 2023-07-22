function solver(dx, dt, tstop, n, gamma, CFL, neq, Q, F; muscl_params=nothing)
    t = 0;
    steps = 0;
    while t < tstop
        if EULER == true
            P = @. (gamma - 1)*(Q[3, :] - 0.5*(Q[2, :]^2)/Q[1, :]);
            a = @. sqrt(gamma*P/Q[1, :]);
        elseif MHD == true
            a = mhdsound(Q, gamma);
        end
        ua1 = @. abs(Q[2, :]/Q[1, :] + a);
        ua2 = @. abs(Q[2, :]/Q[1, :] - a);
        ua = [maximum(ua1), maximum(ua2)];
        ua_max::Float64 = maximum(ua);
        if ua_max*dt/dx > CFL
            dt = CFL*dx/ua_max;
        end
        #---
        if steps%1000 == 0
            # println("dt = ", dt)
            # print_time = (tstop - t);
            # println("Time = ", print_time)
        end
        #---
        if RUN_MACCORMACK == true
            Q = maccormack(Q, F, n, neq, gamma, dt, dx);
            for i = 1:n
                if EULER == true
                    F[:, i] = eulereq(Q[:, i], gamma);
                elseif MHD == true
                    F[:, i] = mhd(Q[:, i], gamma);
                end
            end
        elseif RUN_MUSCL == true
            limiter = muscl_params[1];
            k = muscl_params[2];
            eps = muscl_params[3];
            Q = muscl(Q, n, limiter, gamma, neq, eps, k, dt, dx);
        end
        t += dt;
        steps += 1;
    end
    if OUTPUT_DATA == true
        if EULER == true
            if RUN_MACCORMACK == true
                writedlm("sod_shock_final_euler_maccormack.csv", Q, ',')
            elseif RUN_MUSCL == true
                writedlm("sod_shock_final_euler_muscl.csv", Q, ',')
            end
        elseif MHD == true
            if RUN_MACCORMACK == true
                writedlm("sod_shock_final_mhd_maccormack.csv", Q, ',')
            elseif RUN_MUSCL == true
                writedlm("sod_shock_final_mhd_muscl_1e4.csv", Q, ',')
            end
        end
    end
    return Q, steps
end
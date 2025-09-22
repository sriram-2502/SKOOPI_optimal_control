function u = get_swing_up_control(lqr_params, x_op, x_desired)

    % parse info
    K_lqr = lqr_params.K_lqr;

    if (abs(x_desired(2)-x_op(2)))*(180/pi) <= 30 %
            u_volt = -K_lqr*(x_op'-x_desired');
            u_volt = saturate_fun(u_volt,12,-12);
            x_dot  = x_op(3);
            u_lqr  = volts_to_force(x_dot,u_volt);
            u      = u_lqr;
        else
            u_energy = get_energy_based(x_op);
            u = u_energy;  
    end
end
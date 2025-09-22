function force_out = volts_to_force(x_dot, volts)

    % motor parameters for cart pole
    
    Rm  = 12.5;        % Motor Resistance
    kb  = 0.031;       % Motor back emf constant
    kt  = 0.031;       % Motor Torque constant
    r   = 0.006;       % radius of motor shaft
    
    force_out = (kt*volts*r - kt*kb*x_dot)/(Rm*(r^2));
end
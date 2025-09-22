function joint_ID = ID_planarRR(q_des,q_dot_des,q_ddot_des,robot_params,sys_params,navigation_params,euler_params)

    q = navigation_params.x_ini'; q_dot = [0;0];
    joint_ID.angles = q';
    joint_ID.vel = q_dot';
    joint_ID.control = [0;0]';
    dt = euler_params.step_size; N = euler_params.n_steps;
    
    % wrap around SE2 torous for q1 and q2
    if(q(1) > 2*pi)
        disp('wrap:0');
        q(1) = q(1,1) - 2*pi; 
    end
    if(q(2,1) > 2*pi)
        disp('wrap:0')
        q(2) = q(2,1) - 2*pi; 
    end
    if(q(1) < 0)
        disp('wrap:0')
        q(1) = 2*pi + q(1,1);
    end
    if(q(2) < 0)
        disp('wrap:0')
        q(2) = 2*pi + q(2,1);
    end

    t = 0;
    for i = 1:N-2
        x = [q;q_dot];  
        % get M C G matrices
        sys_params.add_noise =  true;
        sys_info = planarRR_info(t, x, [0;0], robot_params, sys_params);
        M = sys_info.M; C = sys_info.C; G = sys_info.G; D = sys_info.D;

        % vanilla inverse dynamics
        Kp  = 5; Kv = 10;
        e = q - q_des(i,:)'; e_dot = q_dot - q_dot_des(i,:)';
%       u_id = M*(q_ddot_des(i,:)') + C + G + D*q_dot + M*(-Kp*e -Kv.*e_dot);
        u_id = M*(q_ddot_des(i,:)') + C + G + M*(-Kp*e -Kv.*e_dot);

        % apply control to simulate using euler
        dxdt = dynamics_planarRR(x, u_id, sys_info);
        x_update = x + dt.*dxdt;
        
        % update states
        q = x_update(1:2); q_dot = x_update(3:4);
        
         % wrap around SE2 torous for q1 and q2
        if(q(1) > 2*pi)
            disp('wrap:'); disp(i);
            q(1) = q(1,1) - 2*pi; 
        end
        if(q(2) > 2*pi)
            disp('wrap:'); disp(i);
            q(2) = q(2,1) - 2*pi;
        end
        if(q(1) < 0)
            disp('wrap:'); disp(i);
            q(1) = 2*pi + q(1,1);
        end
         if(q(2) < 0)
            disp('wrap:'); disp(i);
            q(2) = 2*pi + q(2,1);
        end 

        joint_ID.angles = [joint_ID.angles; q'];
        joint_ID.vel = [joint_ID.vel; q_dot'];
        joint_ID.control = [joint_ID.control; u_id'];
        t = t + dt;
    end
end
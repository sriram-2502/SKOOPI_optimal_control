function joint_SKOOPI = ID_planarRR_SKOOPI(q_des,q_dot_des,q_ddot_des,robot_params,sys_params,navigation_params,euler_params,lqr_params_transformed,sys_info,L_list,s_x_list1)

    q = navigation_params.x_ini'; q_dot = [0;0];
    joint_SKOOPI.angles = q';
    joint_SKOOPI.vel = q_dot';
    joint_SKOOPI.control = [0;0]';
    dt = euler_params.step_size; N = euler_params.n_steps;
    dynamics = @dynamics_planarRR;
    R = lqr_params_transformed.R;
    Lambda = sys_info.eig_vals;
    W = sys_info.eig_vectors;
    A = sys_info.A;
    
    %% wrap around SE2 torous for q1 and q2
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

    %% simulation loop
    t = 0;
    for i = 1:N-2
        x = [q;q_dot];
        sys_info = planarRR_info(t, x, [0;0], robot_params, sys_params);
        sys_info.eig_vals = Lambda;
        sys_info.eig_vectors = W;
        sys_info.A = A;
        g = sys_info.dynamics_g;
        M = sys_info.M; C = sys_info.C; G = sys_info.G; D = sys_info.D;

        % get SKOOPI control
        % ------ compute eigfn based control ------
        phi = compute_path_integrals(x, dynamics, sys_info);
        phi_x_op            = phi.phi_x_op';
        phi_x_op = [real(phi_x_op(1)); imag(phi_x_op(1));real(phi_x_op(3)); imag(phi_x_op(3));];
        grad_phi_x_op       = compute_gradients(phi);
        grad_phi_x_op = [
            real(grad_phi_x_op(1,:));  % real part of row 1
            imag(grad_phi_x_op(1,:));  % imag part of row 1
            real(grad_phi_x_op(3,:));  % real part of row 3
            imag(grad_phi_x_op(3,:));  % imag part of row 3
        ];

        % get SCOOPI controller
        L      = L_list{i};
        s_x    = s_x_list1{i};
        K1     = 2*phi_x_op'*L*grad_phi_x_op;
        K2     = 2*s_x'*grad_phi_x_op;
%         scale = [1;1];
%         scale = [4;10];
          scale = [5;12];
%         scale = 1e-20*[1;1];
        % u_klqr = M*(q_ddot_des(i,:)') + C + G + D*q_dot + scale*M*(-inv(R)*g(x)'*(K1+K2)');
%         u_klqr = scale.*(-inv(R)*g(x)'*(K1+K2)');
        u_klqr = M*(q_ddot_des(i,:)') + C + G + D*q_dot + sys_info.tau_d + scale.*M*(-g(x)'*(K1+K2)');
        % u_klqr = scale.*(g(x)'*(K1+K2)');

        if(norm(u_klqr) > 100)
            disp('!!! ---- SKOOPI: high control values -----')
            u_klqr
        end

        % apply control to simulate using euler
        dxdt = dynamics_planarRR(x, u_klqr, sys_info);
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

        joint_SKOOPI.angles = [joint_SKOOPI.angles; q'];
        joint_SKOOPI.vel = [joint_SKOOPI.vel; q_dot'];
        joint_SKOOPI.control = [joint_SKOOPI.control; u_klqr'];
        t = t + dt;
    end
end
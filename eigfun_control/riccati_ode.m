function P_dot = riccati_ode(t,P,A,B,Q,R)
    % function to output Differetial Ricciati ODE 
    P = reshape(P,size(A));
    P_dot = -(A'*P + P*A - P*B*inv(R)*B'*P + Q);
    P_dot = P_dot(:); % vectorize for ode solving
end
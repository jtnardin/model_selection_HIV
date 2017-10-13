function vprime = ODE_model_1_w_sens(t,v,q)

    
    % q = [c,delta, N,n,T0,V0]
    % v = [T,V,S11,S12,S21,S22]
    
    k = q(1)/(q(3)*q(5));
    
    A = zeros(6);

    %T
    A(1,1:2) = [-q(2) (1-q(4))*k*q(5)];
    %V
    A(2,1:2) = [q(3)*q(2) -q(1)];
    %S_11 = dTdc
    A(3,:) = [0 (1-q(4))/q(3) -q(2) 0 (1-q(4))*k*q(5) 0];
    %S_12 = dTddelta
    A(4,:) = [-1 0 0 -q(2) 0 (1-q(4))*k*q(5)];
    %S_21 = dVdc
    A(5,:) = [0 -1 q(3)*q(2) 0 -q(1) 0];
    %S_22 = dvddelta
    A(6,:) = [q(3) 0 0 q(3)*q(2) 0 -q(1)];
    
    
    vprime = A*v;

end
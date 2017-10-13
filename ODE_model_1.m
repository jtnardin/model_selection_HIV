function vprime = ODE_model_1(t,v,q)

    
    % q = [c,delta, N,n,T0,V0]

    k = q(1)/(q(3)*q(5));

    A = [-q(2) (1-q(4))*k*q(5) ; q(3)*q(2) -q(1)];

    vprime = A*v;

end
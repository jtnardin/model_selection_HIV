function vprime = ODE_model_2(t,v,q)

    % v = [T,Ts,V]
    
    % q = [c,delta,n,p,N,dT,f,k2,Tmax,T0,V0,s]
    %      1 2     3 4 5 6  7 8  9    10 11 12 
   
    % k = c/(T0*N)
    k = q(1)/(q(10)*q(5));
    % s = dT*T0 + cV0/N
    s = q(6)*q(10)+q(1)*q(11)/q(5);
    
    

    A = [q(4)-q(6) 0 0 ; 0 -q(2) 0 ; 0 q(5)*q(2) -q(1)];
    
    nonlin = [s + q(4)/q(9)*v(1)*(v(1)+v(2)) - (1-q(3))*k*v(1)*v(3) ;...
        (1-q(3))*k*v(1)*v(3) ; 0];

    vprime = A*v + nonlin;

end
function vprime = ODE_model_2_w_sens(t,v,q)

    % v = [T,Ts,V,Tc,Tsc,Vc,Td,Tsd,Vd]
    %      1 2  3 4  5   6  7  8   9
    
    % q = [c,delta,n,p,N,dT,f,k2,Tmax,T0,V0,s]
    %      1 2     3 4 5 6  7 8  9    10 11 12 
   
    % k = c/(T0*N)
    k = q(1)/(q(10)*q(5));
    % s = dT*T0 + cV0/N
    s = q(6)*q(10)+q(1)*q(11)/q(5);
    
    

    A = [q(4)-q(6) 0 0 ; 0 -q(2) 0 ; 0 q(5)*q(2) -q(1)];
    Ac = [0 0 0 ; 0 0 0 ; 0 0 -1];
    Ad = [0 0 0 ; 0 -1 0 ; 0 q(5) 0];
    
    zero3 = zeros(3);
    
    A_all = [A zero3 zero3 ; Ac A zero3 ; Ad zero3 A];
    
    nonlin = [s + q(4)/q(9)*v(1)*(v(1)+v(2)) - (1-q(3))*k*v(1)*v(3) ;...
        (1-q(3))*k*v(1)*v(3) ; 0];
    
    nonlinc = [q(11)/q(5) + q(4)/q(9)*(2*v(1)*v(4) + v(1)*v(5) + v(4)*v(2)) - (1-q(3))/(q(5)*q(10))*(v(1)*v(3) +...
        q(1)*v(4)*v(3) + q(1)*v(1)*v(6)) ; (1-q(3))/(q(5)*q(10))*(v(1)*v(3) + q(1)*v(4)*v(3) + q(1)*v(1)*v(6)) ; 0];
    
    nonlind = [q(4)/q(9)*(2*v(1)*v(7) + v(1)*v(8) + v(7)*v(2)) - (1-q(3))*k*(v(7)*v(3) + v(1)*v(9)) ; (1-q(3))*k*(v(7)*v(3) + v(1)*v(9)) ; 0 ];

    nonlin_all = [nonlin ; nonlinc ; nonlind];
    
    
    vprime = A_all*v + nonlin_all;

end
function J = cost_function_model2(q,qknown,tdata,ydata,ode_solve,h)

    

    % q = [c,delta]
    % qknown = [n,p,N,dT,f,k2,Tmax,T0,V0]
    
    %calculate initial condition for Tstar = cV0/(N*delta)
    T0s = q(1)*qknown(9)/(qknown(3)*q(2));
    
    
    switch ode_solve
    
        case 'forward_euler'
    
            t = 0:h:tdata(end);
            
            model = euler_sim_model2([q;qknown],t,[qknown(8);T0s;qknown(9)],h);
            
                      
        case 'rk4'
            
            options = odeset('AbsTol',h);
            
            [t,model] = ode45(@(t,model) ODE_model_2(t,model,[q;qknown]),[tdata(1) tdata(end)]...
                ,[qknown(8);T0s;qknown(9)],options);
    
                        
        otherwise
        
            error('invalid ODE solver specified')
            
    end

    model_interp = interp1(t,model(:,3),tdata);
    

    res = (log(model_interp) - log(ydata))/log(10);
    J = 1/length(ydata)*(res'*res);
    

end
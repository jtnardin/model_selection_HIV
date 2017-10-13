function J = cost_function(q,qknown,tdata,ydata,ode_solve,h)

    t = 0:h:tdata(end);

    % q = [c,delta]
    % qknown = [N,n,T0,V0]
    
    %calculate initial condition for Tstar
    T0s = q(1)*qknown(4)/(qknown(1)*q(2));
    
    
    switch ode_solve
    
        case 'forward_euler'
        
            model = euler_sim_model([q;qknown],t,h);
            
            model_interp = interp1(t,model(2,:),tdata);
            
        case 'rk4'
            
            options = odeset('AbsTol',h);
            
            [t,model] = ode45(@(t,model) ODE_model_1(t,model,[q;qknown]),t...
                ,[T0s;qknown(4)],options);
    
            model_interp = interp1(t,model(:,2),tdata);
            
        otherwise
        
            error('invalid ODE solver specified')
            
    end

    

    res = (log(model_interp) - log(ydata))/log(10);
    J = 1/length(ydata)*(res'*res);
    

end
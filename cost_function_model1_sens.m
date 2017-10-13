function [J,grad] = cost_function_model1_sens(q,qknown,tdata,ydata,ode_solve,h)

    

    % q = [c,delta]
    % qknown = [N,n,T0,V0]
    
    %calculate initial condition for Tstar
    T0s = q(1)*qknown(4)/(qknown(1)*q(2));
    S110 = qknown(4)/(qknown(1)*q(2));
    S120 = -q(1)*qknown(4)/(qknown(1)*q(2)^2);
    
    model0 = [T0s;qknown(4);S110;S120;0;0];
    
    switch ode_solve
    
        case 'forward_euler'
            
            t = 0:h:tdata(end);
        
            model = euler_sim_model1_sens([q;qknown],t,[T0s;qknown(4);...
                zeros(4,1)],h);
            
        case 'rk4'
            
            options = odeset('AbsTol',h);
            
            [t,model] = ode45(@(t,model) ODE_model_1_w_sens(t,model,[q;qknown]),...
                [tdata(1) tdata(end)],model0,options);
    
                       
        otherwise
        
            error('invalid ODE solver specified')
            
    end
    
    model_interp = interp1(t,model(:,[2 5 6]),tdata);

    

    %%%%%% abs of V for now -- be careful of this
    res = (log(abs(model_interp(:,1))) - log(ydata))/log(10);
    J = 1/length(ydata)*(res'*res);
    
%     res_sens = res.*(1./model_interp(:,1) - 1./ydata)/log(10).*(model_interp(:,2:3));
    res_sens = res.*(1./model_interp(:,1))/log(10).*(model_interp(:,2:3));
    grad = 2/length(ydata)*sum(res_sens,1)';
    

end
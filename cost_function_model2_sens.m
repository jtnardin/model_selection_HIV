function [J,grad] = cost_function_model2_sens(q,qknown,tdata,ydata,ode_solve,h)

    
    % q = [c,delta]
%          1 2  
    % qknown = [n,p,N,dT,f,k2,Tmax,T0,V0]
%               1 2 3 4  5 6  7    8  9 
    %calculate initial condition for Tstar = cV0/(N*delta)
    T0s = q(1)*qknown(9)/(qknown(3)*q(2));
    
    %T0sc = V0/(Ndelta)
    T0sc = qknown(9)/(qknown(3)*q(2));
    T0sd = -q(1)*qknown(9)/(qknown(3)*q(2)^2);
    
    model0 = [qknown(8);T0s;qknown(9);0;T0sc;0;0;T0sd;0];
    
    switch ode_solve
    
        case 'forward_euler'
            
            t = 0:h:tdata(end);
        
            model = euler_sim_model2_sens([q;qknown],t,model0,h);
            
        case 'rk4'
            
            options = odeset('AbsTol',h);
            
            [t,model] = ode45(@(t,model) ODE_model_2_w_sens(t,model,[q;qknown]),...
                [tdata(1) tdata(end)],model0,options);
    
                       
        otherwise
        
            error('invalid ODE solver specified')
            
    end
    
    model_interp = interp1(t,model(:,[3 6 9]),tdata);

    

    %%%%%% abs of V for now -- be careful of this
    res = (log(abs(model_interp(:,1))) - log(ydata))/log(10);
    J = 1/length(ydata)*(res'*res);
    
%     res_sens = res.*(1./model_interp(:,1) - 1./ydata)/log(10).*(model_interp(:,2:3));
    res_sens = res.*(1./model_interp(:,1))/log(10).*(model_interp(:,2:3));
    grad = 2/length(ydata)*sum(res_sens,1)';
    

end
% model 1

qall = [0.638  0.638 ; 0.648 0.648 ; .649 .649 ; .649 .649 ]';

patient = 1;
stepsz = 1;
q = qall(:,stepsz);

load('Viral_data.mat')
%how do we want to solve ODE;
ode_solve = 'forward_euler';

%number of patients
pnum = length(T);

%data time vector
tdata = t(1:end-1);
N = length(tdata);

%known parameters
N=480;
n=0.8;

%last point always a little off
tdata = t(1:end-1);
y = y(1:end-1,:);

ydata = y(:,patient);


V0 = ydata(1);
T0 = T(patient);

qknown = [N,n,T0,V0]';

switch ode_solve
    
    case 'forward_euler'
        
        %step sizes considered
        h = [1e-1 1e-2 1e-3 1e-4];
        
    case 'rk4'
        
        %tolerances considered
        h = [1e-3 1e-6 1e-9 ];
        
    otherwise
        
        error('invalid solver specified')
        
end


  t = 0:h:tdata(end);

    % q = [c,delta]
    % qknown = [N,n,T0,V0]
    
    %calculate initial condition for Tstar
    T0s = q(1)*qknown(4)/(qknown(1)*q(2));
    
    
    switch ode_solve
    
        case 'forward_euler'
        
            model = euler_sim_model1([q;qknown],t,[T0s;qknown(4)],h);
            
        case 'rk4'
            
            options = odeset('AbsTol',h);
            
            [t,model] = ode45(@(t,model) ODE_model_1(t,model,[q;qknown]),t...
                ,[T0s;qknown(4)],options);
    
                       
        otherwise
        
            error('invalid ODE solver specified')
            
    end
 
 model_interp = interp1(t,model(:,2),tdata);
 
 
 semilogy(tdata,ydata,'k.')
 hold on
 semilogy(t,model(:,2))


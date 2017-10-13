%%% likelihood_compute_all_patients_make_table.m written 10-5-17 by JTN to
%%% fit simple viral model to all patient data using various step sizes. We
%%% will then check how the choice of step size influences the likelihood
%%% computation


%display when simulation began
c = clock;
disp(['time is ' date ' ' num2str(c(4)) ':' num2str(c(5)) ])


load('Viral_data.mat')
%how do we want to solve ODE;
ode_solve = 'rk4';

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

% 
% %table to be made in latex
% A = struct;
% A.data = zeros(pnum,3*length(h));


%for parallelization
A = cell(pnum,1);

parfor i = 1
    
    A{i} = zeros(1,3*length(h));

    ydata = y(:,i);
    V0 = ydata(1);
    T0 = T(i);
    
    for j = 1:length(h)
        
       
        qknown = [N,n,T0,V0]';
        
        options = optimset;
        tic
        [q,J] = fminsearch(@(q) cost_function_model1(q,qknown,tdata,ydata,ode_solve,h(j)),[1 1]',options);
        toc
        
        
        A{i}(3*(j-1)+1:3*(j-1)+2) = q';
        A{i}(3*j) = -tN/2*(1+log(2*pi)+log(J));
        
    end
end
   
save('/scratch/summit/jona8898/model_selection/model_1_tol.mat','A')
%%% likelihood_compute_all_patients_make_table.m written 10-5-17 by JTN to
%%% fit simple viral model to all patient data using various step sizes. We
%%% will then check how the choice of step size influences the likelihood
%%% computation

load('Viral_data.mat')
%how do we want to solve ODE;
ode_solve = 'rk4';

%number of patients
pnum = length(T);

%data time vector
tdata = t(1:end-1);
N = length(tdata);


%known parameters
n=0.8;
p = 0.03;
N=480;
dT = 0.02;
f = 0.03;
k2 = 0.01;
Tmax = 1500;


%last point always a little off
tdata = t(1:end-1);


switch ode_solve
    
    case 'forward_euler'
        
        %step sizes considered
        h = [1e-1 1e-2 1e-3 1e-4];
        
    case 'rk4'
        
        %tolerances considered
        h = [1e-3 1e-6 1e-9];
        
    otherwise
        
        error('invalid solver specified')
        
end


%table to be made in latex
A = struct;
A.data = zeros(pnum,3*length(h));


tic

for i = 1:pnum
   
         
        ydata = y(1:end-1,i);
        V0 = ydata(1);
        T0 = T(i);
   
    
    for j = 1:length(h)
        
         
        qknown = [n,p,N,dT,f,k2,Tmax,T0,V0]';
        
        options = optimoptions('fmincon','CheckGradients',false,'SpecifyObjectiveGradient',true,'algorithm','sqp','display','iter');%optimset('algorithm','sqp');%('display','iter');
        tic
        [q,J] = fmincon(@(q) cost_function_model2_sens(q,qknown,tdata,ydata,ode_solve,h(j))...
            ,[2 1]',[],[],[],[],[],[],[],options);
        toc

        
        A.data(i,3*(j-1)+1:3*(j-1)+2) = q';
        A.data(i,3*j) = -N/2*(1+log(2*pi)+log(J));
        j
    end
    i
end
   
toc

A.dataFormat = {'%.3f'};
A.makeCompleteLatexDocument = 1;
latexTable(A);
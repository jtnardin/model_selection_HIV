%%% likelihood_compute_all_patients_make_table.m written 10-5-17 by JTN to
%%% fit simple viral model to all patient data using various step sizes. We
%%% will then check how the choice of step size influences the likelihood
%%% computation

patient = 8;
tol_ind = 3;

save_fig = 0;

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



     

ydata = y(1:end-1,patient);
V0 = ydata(1);
T0 = T(patient);

qknown = [n,p,N,dT,f,k2,Tmax,T0,V0]';


c = linspace(0,40,20);
delta = linspace(0,1,20);

J = zeros(length(c),length(delta));

tic

for i = 1:length(c)
    for j = 1:length(delta)

        q = [c(i) ; delta(j)];

        J(i,j) = cost_function_model2(q,qknown,tdata,ydata,ode_solve,h(tol_ind));
    end
end

toc

figure

[C,D] = meshgrid(c,delta);

[Jcont,g] = contour(C,D,J');
clabel(Jcont,g)


view(2)

xlabel('c')
ylabel('$\delta$','interpreter','latex')
if patient == 10
    title(['Cost function, patient 2' num2str(patient)])
else
    title(['Cost function, patient 20' num2str(patient)])
end


if save_fig == 1
    exportfig(gcf,['cost_contour_2_pat_' num2str(patient) '_' ode_solve '_'...
        num2str(tol_ind) '.eps'],'color','rgb','fontsize',1.5)
    saveas(gcf,['cost_contour_2_pat_' num2str(patient) '_' ode_solve '_'...
        num2str(tol_ind) '.fig'])
end
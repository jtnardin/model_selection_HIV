load('Viral_data.mat')

patient = 1;

tdata = t(1:end-2);
ydata = y(1:end-2,patient);
T0 = T(patient);


N=480;
n=0.8;
V0 = ydata(1);

qknown = [N,n,T0,V0]';

options = optimset('display','iter');

h = 1e-1;

[q,J] = fminsearch(@(q) cost_function(q,qknown,tdata,ydata,h),[1 1]',options);

t = 0:h:tdata(end);
model = euler_sim_model([q;qknown],t,h);
semilogy(t,model(2,:))
hold on
semilogy(tdata,ydata,'k.')
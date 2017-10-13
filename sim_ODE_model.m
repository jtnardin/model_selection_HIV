t = linspace(0,14,100);

c = 13.7;
delta = 0.35;
N=480;
n=0.8;
V0 = 1e4;

q = [c,delta,N,n,V0]';

y0 = [V0,q(1)*q(5)/(q(3)*q(2))];

[t,y] = ode15s(@(t,y) ODE_model(t,y,q),t,y0);

semilogy(t,y(:,2),'b')
hold on

tdata = [0:14];

ydataperf = interp1(t,y(:,2),tdata,'cubic');


var1 = .2;
ydatanoise1 = exp(sqrt(var1)*randn(size(ydataperf))).*ydataperf;


var2 = 1;
ydatanoise2 = exp(sqrt(var2)*randn(size(ydataperf))).*ydataperf;


var3 = 5;
ydatanoise3 = exp(sqrt(var3)*randn(size(ydataperf))).*ydataperf;


semilogy(tdata,ydatanoise1,'k.','markersize',10)
semilogy(tdata,ydatanoise2,'r.','markersize',10)
semilogy(tdata,ydatanoise3,'m.','markersize',10)

legend('model',num2str(var1),num2str(var2),num2str(var3))

xlabel('Time')
ylabel('ln(v)')
title('Simulated data sets')

% exportfig(gcf,'viral_data_sim.eps','color','rgb','fontsize',1.5)
% 
% save('sim_data_viral','tdata','ydataperf','var1','var2','var3',...
%     'ydataperf','ydatanoise1','ydatanoise2','ydatanoise3','C','delta','N','n','V0')
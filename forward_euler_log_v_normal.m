%written 9-27-17 to use forward Euler to simulate virion model

%in this version, we simulate log(y) instead of y
    %parameters

clear all; clc


c = 13.7;
delta = 0.35;
N=480;
n=0.8;
V0 = 1e4;

q = [c,delta,N,n,V0]';

T0s = q(1)*q(5)/(q(3)*q(2));

k = q(1)/(q(3)*T0s);

A = [-q(2) (1-q(4))*k*T0s ; q(3)*q(2) -q(1)];

y0 = [V0;q(1)*q(5)/(q(3)*q(2))];
ylog0 = log([V0;q(1)*q(5)/(q(3)*q(2))]);

h = [0.1 0.05 0.01 0.001 .0001];

sigmas = zeros(length(h)+1,6);

for j = 1:length(h)

    tfinal = 14;

    t = 0:h(j):tfinal;

    %number points
    tn = length(t);

    y = zeros(2,tn);
    y(:,1) = y0;

    ylog = zeros(2,tn);
    ylog(:,1) = ylog0;


    for i = 2:tn

        ylog(:,i) = ylog(:,i-1)+h(j)*[(1-q(4))*k*T0s*exp(ylog(2,i-1)-ylog(1,i-1))...
            - q(2) ; q(3)*q(2)*exp(ylog(1,i-1)-ylog(2,i-1)) - q(1)];

        y(:,i) = y(:,i-1)+h(j)*A*y(:,i-1);

    end

    load('sim_data_viral.mat','tdata','ydatanoise1','ydatanoise2','ydatanoise3','ydataperf')

    yinterp = interp1(t,y(2,:),tdata,'pchip');
    yloginterp = interp1(t,ylog(2,:),tdata,'pchip');




    res1 = log(yinterp) - log(ydatanoise1);
    sigmas(j+1,1) = 1/length(tdata)*res1*res1';

    logres1 = yloginterp - log(ydatanoise1);
    sigmas(j+1,2) = 1/length(tdata)*logres1*logres1';


    res2 = log(yinterp) - log(ydatanoise2);
    sigmas(j+1,3) = 1/length(tdata)*res2*res2';

    logres2 = yloginterp - log(ydatanoise2);
    sigmas(j+1,4) = 1/length(tdata)*logres2*logres2';


    res3 = log(yinterp) - log(ydatanoise3);
    sigmas(j+1,5) = 1/length(tdata)*res3*res3';

    logres3 = yloginterp - log(ydatanoise3);
    sigmas(j+1,6) = 1/length(tdata)*logres3*logres3';

end

restrue1 = log(ydataperf) - log(ydatanoise1);
restrue2 = log(ydataperf) - log(ydatanoise2);
restrue3 = log(ydataperf) - log(ydatanoise3);

Jtrue1 = 1/length(tdata)*restrue1*restrue1';
Jtrue2 = 1/length(tdata)*restrue2*restrue2';
Jtrue3 = 1/length(tdata)*restrue3*restrue3';

sigmas(1,:) = [Jtrue1 Jtrue1 Jtrue2 Jtrue2 Jtrue3 Jtrue3];

sigmas

%written 9-27-17 to use forward Euler to simulate virion model

%in this version, we simulate log(y) instead of y
    %parameters

    
    
    c = 13.7;
    delta = 0.35;
    N=480;
    n=0.8;
    V0 = 1e4;

    q = [c,delta,N,n,V0]';
    
    T0s = q(1)*q(5)/(q(3)*q(2));
    
    k = q(1)/(q(3)*T0s);

    y0 = log([V0;q(1)*q(5)/(q(3)*q(2))]);
    
    h = 0.0001;
    
    tfinal = 14;
    
    t = 0:h:tfinal;
    
    %number points
    tn = length(t);
    
    y = zeros(2,tn);
    y(:,1) = y0;
    
    
    for i = 2:tn
       
        y(:,i) = y(:,i-1)+h*[(1-q(4))*k*T0s*exp(y(2,i-1)-y(1,i-1)) - q(2) ;...
            q(3)*q(2)*exp(y(1,i-1)-y(2,i-1)) - q(1)];
        
    end
      
    load('sim_data_viral.mat','tdata','ydatanoise1','ydatanoise2','ydatanoise3','ydataperf')
    
    yinterp = interp1(t,y(2,:),tdata,'pchip');
    
    restrue1 = log(ydataperf) - log(ydatanoise1);
    Jtrue1 = 1/length(tdata)*restrue1*restrue1'
    
    res1 = yinterp - log(ydatanoise1);
    J1 = 1/length(tdata)*res1*res1'
    
    
    restrue2 = log(ydataperf) - log(ydatanoise2);
    Jtrue2 = 1/length(tdata)*restrue2*restrue2'
    
    
    res2 = yinterp - log(ydatanoise2);
    J2 = 1/length(tdata)*res2*res2'
    
    
    restrue3 = log(ydataperf) - log(ydatanoise3);
    Jtrue3 = 1/length(tdata)*restrue3*restrue3'
    
    
    res3 = yinterp - log(ydatanoise3);
    J3 = 1/length(tdata)*res3*res3'
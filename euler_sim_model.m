%parameters
function y = euler_sim_model(q,t,h) 
    
    % q = [c, delta, N,n,T0,V0];

    T0s = q(1)*q(6)/(q(3)*q(2));

    k = q(1)/(q(3)*q(5));

    A = [-q(2) (1-q(4))*k*q(5) ; q(3)*q(2) -q(1)];


    y0 = [T0s;q(6)];


    %number points
    tn = length(t);

    y = zeros(2,tn);
    y(:,1) = y0;


    for i = 2:tn

        y(:,i) = y(:,i-1)+h*A*y(:,i-1);

    end

end
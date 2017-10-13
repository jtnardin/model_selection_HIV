%parameters
function v = euler_sim_model2(q,t,v0,h) 
    
    %number points
    tn = length(t);

    v = zeros(3,tn);
    v(:,1) = v0;


    for i = 2:tn

        v(:,i) = v(:,i-1)+h*ODE_model_2(t(i),v(:,i-1),q);

    end
    
    v = v';

end
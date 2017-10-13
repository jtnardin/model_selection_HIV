%parameters
function v = euler_sim_model2_sens(q,t,v0,h) 
    
    %number points
    tn = length(t);

    v = zeros(length(v0),tn);
    v(:,1) = v0;


    for i = 2:tn

        v(:,i) = v(:,i-1)+h*ODE_model_2_w_sens(t(i),v(:,i-1),q);

    end
    
    v = v';

end
%%% likelihood_compute_all_patients_make_table.m written 10-5-17 by JTN to
%%% fit simple viral model to all patient data using various step sizes. We
%%% will then check how the choice of step size influences the likelihood
%%% computation

clear all; clc

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


    switch ode_solve

        case 'forward_euler'

            %step sizes considered
            h = [1e-1 1e-2 1e-3 1e-4];

        case 'rk4'

            %tolerances considered
            h = [1e-3 1e-6 1e-9 1e-10];

        otherwise

            error('invalid solver specified')

    end


    %table to be made in latex
    A = struct;
    A.data = zeros(pnum,3*length(h));

    tic

    LB = [0 0]';

    for i = 1:pnum
        for j = 1:length(h)

            ydata = y(1:end-1,i);
            V0 = ydata(1);
            T0 = T(i);

            qknown = [N,n,T0,V0]';

            options = optimoptions('fmincon','CheckGradients',false,'SpecifyObjectiveGradient',true,'algorithm','sqp','display','iter');%optimset('algorithm','sqp');%('display','iter');

            [q,J] = fmincon(@(q) cost_function_model1_sens(q,qknown,tdata,ydata,ode_solve,h(j))...
                ,[13,2]',[],[],[],[],LB,[],[],options);

            A.data(i,3*(j-1)+1:3*(j-1)+2) = q';
            A.data(i,3*j) = -N/2*(1+log(2*pi)+log(J));

        end
        i
    end

    toc

    A.dataFormat = {'%.3f'};
    A.makeCompleteLatexDocument = 1;
    latexTable(A);
    

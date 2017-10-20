%%% likelihood_compute_all_patients_make_table.m written 10-5-17 by JTN to
%%% fit simple viral model to all patient data using various step sizes. We
%%% will then check how the choice of step size influences the likelihood
%%% computation

clear all; clc

load('Viral_data.mat')
%how do we want to solve ODE;
% ode_solve = 'rk4';

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


for k = 1
    for l = 1

        %which routine are we using
        switch k

            case 1

                ode_solve = 'forward_euler';
                %step sizes considered
                h = [1e-2 1e-3 1e-4];

            case 2

                ode_solve = 'rk4';
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
        q0 = [28 0.25]';

        for i = 1:pnum


                ydata = y(1:end-1,i);
                V0 = ydata(1);
                T0 = T(i);


            for j = 1:length(h)


                qknown = [n,p,N,dT,f,k2,Tmax,T0,V0]';

                
                switch l
                    
                    case 1
                        options = optimoptions('fmincon','SpecifyObjectiveGradient',false,'algorithm','sqp','display','iter','CheckGradients',false);%optimset('algorithm','sqp');%('display','iter');
                        tic
                        [q,J] = fmincon(@(q) cost_function_model2(q,qknown,tdata,ydata,ode_solve,h(j))...
                            ,q0,[],[],[],[],LB,[],[],options);
                        toc
                        
                        sens_id = '';

                    case 2
                        options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'algorithm','sqp','display','iter','CheckGradients',false);%optimset('algorithm','sqp');%('display','iter');
                        tic
                        [q,J] = fmincon(@(q) cost_function_model2_sens(q,qknown,tdata,ydata,ode_solve,h(j))...
                            ,q0,[],[],[],[],LB,[],[],options);
                        toc
                        
                        sens_id = 'sens_';

                    otherwise
                        
                        error('Invalid l value')
                        
                end
                        
                A.data(i,3*(j-1)+1:3*(j-1)+2) = q';
                A.data(i,3*j) = -N/2*(1+log(2*pi)+log(J));

            end
            i
        end

        toc

        A.dataFormat = {'%.3f'};
        A.makeCompleteLatexDocument = 1;
        text = latexTable(A);
% 
        q_text = ' ';

        for i = 1:length(text)
            q_text = [q_text text{i}];
        end

        fileID = fopen(['model2_' sens_id  'sqp_' ode_solve '2.tex'],'w');
        fprintf(fileID,'%s',q_text);
        fclose(fileID)

        
    end
end
clear all; clc

load('fit_model2_2.mat')


%currently, A = fe_no_grad
%           B = rk4_no_grad
%           C = fe_grad
%           D = rk4_grad    

h_fe = 10.^[-2 -3 -4];
h_rk4 = 10.^[-3 -6 -9 -10];



ellA = A(:,3:3:end);
ellB = B(:,3:3:end);

ellC = C(:,3:3:end);
ellD = D(:,3:3:end);

figure('units','normalized','outerposition',[0 0 1 1])

for i = 1:10
    
    subplot(2,5,i)
    
    semilogx(h_fe,ellA(i,:),'k-','linewidth',1)
    hold on
    semilogx(h_rk4,ellB(i,:),'r--','linewidth',1)
%     hold on
    semilogx(h_fe,ellC(i,:),'b.-','linewidth',1)
    semilogx(h_rk4,ellD(i,:),'m-','linewidth',.5)
    
    semilogx(h_fe,ellA(i,:),'k.','markersize',15)
    semilogx(h_rk4,ellB(i,:),'r.','markersize',15)
    semilogx(h_fe,ellC(i,:),'b.','markersize',15)
    semilogx(h_rk4,ellD(i,:),'m.','markersize',15)
        
    
    if i == 1 || i ==6
        ylabel('$\ell(q^{MLE})$','interpreter','latex')
    end
    
    if i > 5
        xlabel('h')
    end
    
    if i == 5
       legend('rk4','feuler w/ grad','rk4 w/ grad','location','west') 
    end
    
    title(['Patient ' num2str(i)])
    
end

% exportfig(gcf,'likelihood_compute.eps','fontsize',1.5,'color','rgb')
% saveas(gcf,'likelihood_compute.fig')

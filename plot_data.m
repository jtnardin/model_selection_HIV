figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:10
    subplot(2,5,i)
    semilogy(t,y(:,i),'k.','markersize',20)
    axis([0 t(end) 1e0 1e5])
    if i == 10
        title(['Patient 2' num2str(i)])
    else
        title(['Patient 20' num2str(i)])
    end
    xlabel('Time (days)')
    ylabel('Viral Load')
end
exportfig(gcf,'v_data.eps','color','rgb','fontsize',1.5)
saveas(gcf,'v_data.fig')
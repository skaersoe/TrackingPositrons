function [cpu,gpu1,gpu2] = benchmark_energy_bremsstrahlung_heap()
    data = csvread('/Users/johannes/Dropbox/Bachelor/TrackingPositrons/Simulation/Data/benchmark_energy_bremsstrahlung',1);
    cpu = [];
    gpu1 = [];
    gpu2 = [];
    comparison = [];
    for i=1:length(data)
        if data(i,1) == 0
            cpu = [cpu;data(i,2:3)];
        elseif data(i,1) == 1
            gpu1 = [gpu1;data(i,2:3)];
        else
            gpu2 = [gpu2;data(i,2:3)];
        end
    end
    for i=1:length(gpu2)
       comparison(i,:) = [cpu(i,1),cpu(i,2)/gpu2(i,2)];
    end
    fontsize=24;
    f = figure;
    %cpu(:,2) = cpu(:,2) / max(data(:,3));
    %gpu1(:,2) = gpu1(:,2) / max(data(:,3));
    %gpu2(:,2) = gpu2(:,2) / max(data(:,3));
    hold off
    plot(cpu(1:length(gpu1),1),cpu(1:length(gpu1),2),'-ob','LineWidth',2);
    hold on
    plot(gpu1(:,1),gpu1(:,2),'-o','Color',[1 0 0],'LineWidth',2);
    plot(gpu2(:,1),gpu2(:,2),'-o','Color',[0 0.75 0],'LineWidth',2);
    % plot(comparison(:,1),comparison(:,2),'--m','LineWidth',2);
    title('Variable versus static heap for bremsstrahlung','FontSize',fontsize);
    xlabel('Energy of incident electrons [MeV]','FontSize',fontsize);
    ylabel('Runtime [s]','FontSize',fontsize);
    l1=legend('CPU','GPU variable heap','GPU max heap','Location','SouthEast');
    set(l1,'FontSize',fontsize);
    set(gca,'FontSize',fontsize,'YScale','log','XScale','log','XTick',[1e1,1e2,1e3,1e4,1e5,5e5,1e6],'YLim',[1e-2 1e2],'YTick',[1e-2 1e-1 1e0 1e1 1e2]);
    % axis([min(gpu1(:,1)) max(gpu1(:,1)) -1 10*(ceil(0.1*max(data(:,3))))]);
    set(f,'Position',[800,600,800,600],'PaperPositionMode','auto');
    print(f,'-depsc','../Plots/benchmark_energy_bremsstrahlung_heap.eps');
end
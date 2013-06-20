function [cpu,gpu1,gpu2] = benchmark_energy_bremsstrahlung()
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
    cpu(:,2) = cpu(:,2) / max(data(:,3));
    gpu1(:,2) = gpu1(:,2) / max(data(:,3));
    gpu2(:,2) = gpu2(:,2) / max(data(:,3));
    hold off
    plot(cpu(1:length(gpu1),1),cpu(1:length(gpu1),2),'-ob','LineWidth',2);
    hold on
    plot(gpu1(:,1),gpu1(:,2),'-o','Color',[0 0.75 0],'LineWidth',2);
    %plot(gpu2(:,1),gpu2(:,2),'-o','Color',[0 0.75 0],'LineWidth',2);
    plot(comparison(:,1),comparison(:,2),'--m','LineWidth',2);
    title('Bremsstrahlung benchmark for 32 incident electrons','FontSize',20);
    xlabel('Energy of incident electrons','FontSize',20);
    ylabel('Normalized runtime [1/max]','FontSize',20);
    l1=legend('CPU','GPU variable heap','CPU time / GPU Time','Location','NorthWest');
    set(l1,'FontSize',20);
    set(gca,'FontSize',20,'XScale','log','YScale','log');
    axis([min(gpu1(:,1)) max(gpu1(:,1)) 0 max(gpu1(:,2))]);
end
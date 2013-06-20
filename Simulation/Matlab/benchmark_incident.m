function [cpu,gpu] = benchmark_incident()
    data = csvread('/Users/johannes/Dropbox/Bachelor/TrackingPositrons/Simulation/Data/benchmark_incident',1);
    data(:,3) = data(:,3) / max(data(:,3));
    cpu = [];
    gpu = [];
    comparison = [];
    for i=1:length(data)
        if data(i,1) == 0
            cpu = [cpu;data(i,2:3)];
        else
            gpu = [gpu;data(i,2:3)];
        end
    end
    for i=1:length(cpu)
       comparison(i,:) = [cpu(i,1),gpu(i,2)/cpu(i,2)];
    end
    hold off
    plot(cpu(:,1),cpu(:,2),'-ob','LineWidth',2);
    hold on
    plot(gpu(:,1),gpu(:,2),'-o','Color',[0 0.75 0],'LineWidth',2);
    plot(comparison(:,1),comparison(:,2),'--r','LineWidth',2);
    title('Benchmark for number of incident particles','FontSize',20);
    xlabel('Number of incident particles','FontSize',20);
    ylabel('Normalized runtime [1/max]','FontSize',20);
    l1=legend('CPU','GPU','GPU time / CPU time','Location','NorthEast');
    set(l1,'FontSize',17);
    set(gca,'FontSize',20,'XScale','log','XTick',[2^7,2^8,2^9,2^10,2^11,2^12,2^13,2^14],'YTick',[0,0.25,0.5,0.75,1,1.25,1.5,1.75,2]);
    a = [2^7 max(data(:,2)) 0 2];
    axis(a);
    plot([a(1),a(2)],[1 1],'--k');
end
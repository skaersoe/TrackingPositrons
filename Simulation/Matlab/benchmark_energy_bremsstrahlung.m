function [cpu,gpu] = benchmark_energy_bremsstrahlung()
    data = csvread('/Users/johannes/Dropbox/Bachelor/TrackingPositrons/Simulation/Data/benchmark_energy_bremsstrahlung',1);
    % data(:,3) = data(:,3) / max(data(:,3));
    cpu = [];
    gpu = [];
    comparison = [];
    for i=1:length(data)
        if data(i,1) == 0
            cpu = [cpu;data(i,2:3)];
        elseif data(i,1) == 1
            gpu = [gpu;data(i,2:3)];
        end
    end
    for i=1:length(cpu)
       comparison(i,:) = [cpu(i,1),gpu(i,2)/cpu(i,2)];
    end
    fontsize=24;
    f = figure;
    hold off
    x = comparison(:,1);
    [haxes,hline1,hline2]=plotyy(x,cpu(:,2),x,comparison(:,2));
    set(haxes(1),'Box','off')
    set(haxes(2),'Box','off')
    set(haxes(2), 'XTickLabel','','XAxisLocation','Top') 
    axes(haxes(1));
    y_tick = [1e-2 1e-1 1e0 1e1 1e2];
    set(gca,'FontSize',fontsize,'XScale','log','YScale','log','ycolor',[0 0 0],'YLim',[min(y_tick) max(y_tick)],'YTick',y_tick);
    xlabel(haxes(1),'Energy of incident electrons [MeV]','FontSize',fontsize);
    ylabel(haxes(1),'Runtime [s]','FontSize',fontsize);
    axes(haxes(2));
    set(gca,'FontSize',fontsize,'XScale','log','YLim',[0 max(comparison(:,2))],'YTick',0:1:max(comparison(:,2)),'ycolor',[1 0 0]);
    ylabel(haxes(2),'GPU / CPU runtime','FontSize',fontsize);
    set(hline1,'LineStyle','-','Marker','o','Color',[0 0 1],'LineWidth',2);
    set(hline2,'LineStyle','--','Color',[1 0 0],'LineWidth',2);
    hold(haxes(1));
    % hold(haxes(2));
    plot(haxes(1),gpu(:,1),gpu(:,2),'-o','Color',[0 0.75 0],'LineWidth',2);
    %plot(comparison(:,1),comparison(:,2),'--r','LineWidth',2);
    title('Bremsstrahlung benchmark for 64 incident electrons','FontSize',fontsize);
    l1=legend('GPU time / CPU time','CPU','GPU','Location','North');
    set(l1,'FontSize',fontsize);
    %set(gca,'FontSize',fontsize,'XScale','log');
    %axis([min(data(:,2)) max(data(:,2)) 0 1]);
    set(f,'Position',[800,600,880,660],'PaperPositionMode','auto');
    print(f,'-depsc','../Plots/benchmark_energy_bremsstrahlung.eps');
end
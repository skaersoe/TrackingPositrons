function [cpu,gpu] = benchmark_incident()
    data = csvread('/Users/johannes/Dropbox/Bachelor/TrackingPositrons/Simulation/Data/benchmark_incident');
    % data(:,3) = data(:,3) / max(data(:,3));
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
    fontsize=24;
    f = figure;
    hold off
    x = comparison(:,1);
    [haxes,hline1,hline2]=plotyy(x,cpu(:,2),x,comparison(:,2));
    set(haxes(1),'Box','off')
    set(haxes(2),'Box','off')
    set(haxes(2), 'XTickLabel','','XAxisLocation','Top') 
    axes(haxes(1));
    y_tick = [1e-1 1e0 1e1 1e2 1e3];
    set(gca,'FontSize',fontsize,'XScale','log','YScale','log','ycolor',[0 0 0],'YLim',[5*min(y_tick) max(y_tick)],'YTick',y_tick,'XTick',[2^6,2^7,2^8,2^9,2^10,2^11,2^12,2^13,2^14]);
    xlabel(haxes(1),'Number of incident electrons at 10 MeV','FontSize',fontsize);
    ylabel(haxes(1),'Runtime [s]','FontSize',fontsize);
    axes(haxes(2));
    set(gca,'FontSize',fontsize,'XScale','log','YLim',[0 2],'YTick',0:0.25:2,'ycolor',[1 0 0]);
    ylabel(haxes(2),'GPU / CPU runtime','FontSize',fontsize);
    set(hline1,'LineStyle','-','Marker','o','Color',[0 0 1],'LineWidth',2);
    set(hline2,'LineStyle','--','Color',[1 0 0],'LineWidth',2);
    hold(haxes(1));
    plot(haxes(1),gpu(:,1),gpu(:,2),'-o','Color',[0 0.75 0],'LineWidth',2);
    hold(haxes(2));
    a = [2^6 max(data(:,2)) 0 2];
    plot(haxes(2),[a(1),a(2)],[1 1],'--k');
    %plot(comparison(:,1),comparison(:,2),'--r','LineWidth',2);
    title('Incident particle count benchmark for 10 MeV electrons','FontSize',fontsize);
    l1=legend('CPU','GPU','GPU time / CPU time','Location','NorthEast');
    set(l1,'FontSize',fontsize);
    %set(gca,'FontSize',fontsize,'XScale','log');
    %axis([min(data(:,2)) max(data(:,2)) 0 1]);
    set(f,'Position',[800,600,800,600],'PaperPositionMode','auto');
    print(f,'-depsc','../Plots/benchmark_incident.eps');
end

function [cpu,gpu] = old()
    data = csvread('/Users/johannes/Dropbox/Bachelor/TrackingPositrons/Simulation/Data/benchmark_incident',1);
    % data(:,3) = data(:,3) / max(data(:,3));
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
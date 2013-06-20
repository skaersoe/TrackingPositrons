function [] = stepsperlaunch_bremsstrahlung()
  subplot(1,2,2);
%   data = csvread('/Users/johannes/Dropbox/Bachelor/TrackingPositrons/Simulation/Data/bremsstrahlung_benchmark',1);
%   data(:,2) = data(:,2)/max(data(:,2));
%   hold off
%   plot(data(:,1),data(:,2),'-og','LineWidth',2);
%   hold on
  data = csvread('/Users/johannes/Dropbox/Bachelor/TrackingPositrons/Simulation/Data/benchmark_steps_bremsstrahlung',1);
  data(:,2) = data(:,2)/max(data(:,2));
  hold off
  plot(data(1:length(data)-1,1),data(1:length(data)-1,2),'-o','Color',[0 0.75 0],'LineWidth',2);
  hold on
  cpu = [min(data(:,1)),max(data(:,1));
         data(length(data),2),data(length(data),2)];
  lowestgpu = min(data(1:length(data)-1,2));
  lowestgpu_line = [cpu(1,:);lowestgpu,lowestgpu];
  plot(cpu(1,:),cpu(2,:),'--b','LineWidth',3);
  plot(lowestgpu_line(1,:),lowestgpu_line(2,:),'--','Color',[0 0.75 0],'LineWidth',1);
  axis([min(data(:,1)),max(data(:,1)),min(data(:,2))-0.04,1.04]);
  title('B) Energy loss and bremsstrahlung (high complexity, small dimensions)','FontSize',14);
  xlabel('Steps per kernel launch','FontSize',16);
  ylabel('Runtime [1/max]','FontSize',16);
  l1=legend('GPU','CPU','Fastest GPU configuration','Location','North');
  set(l1,'FontSize',14);
  set(gca,'FontSize',14);
  % Energy loss benchmark
  subplot(1,2,1);
  data = csvread('/Users/johannes/Dropbox/Bachelor/TrackingPositrons/Simulation/Data/benchmark_steps',1);
  data(:,2) = data(:,2)/max(data(:,2));
  hold off
  plot(data(1:length(data)-1,1),data(1:length(data)-1,2),'-o','Color',[0 0.75 0],'LineWidth',2);
  hold on
  cpu = [min(data(:,1)),max(data(:,1));
         data(length(data),2),data(length(data),2)];
  lowestgpu = min(data(1:length(data)-1,2));
  lowestgpu_line = [cpu(1,:);lowestgpu,lowestgpu];
  plot(cpu(1,:),cpu(2,:),'--b','LineWidth',3);
  plot(lowestgpu_line(1,:),lowestgpu_line(2,:),'--','Color',[0 0.75 0],'LineWidth',1);
  title('A) Energy loss (low complexity, large dimensions)','FontSize',14);
  axis([min(data(:,1)),max(data(:,1)),min(data(:,2))-0.04,1.04]);
  xlabel('Steps per kernel launch','FontSize',16);
  ylabel('Runtime [1/max]','FontSize',16);
  l2=legend('GPU','CPU','Fastest GPU configuration','Location','East');
  set(l2,'FontSize',14);
  set(gca,'FontSize',14);
end
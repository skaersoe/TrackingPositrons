function [] = stepsperlaunchdiverse_benchmark()
  data = csvread('/Users/johannes/Dropbox/Bachelor/TrackingPositrons/Simulation/Data/stepsperlaunchdiverse_benchmark',1);
  plot(data(:,1),data(:,2),'-xg');
  title('Steps per kernel launch benchmark (diverse geometry)');
  xlabel('Steps per kernel launch');
  ylabel('Runtime [s]');
end
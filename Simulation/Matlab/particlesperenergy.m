function [x,y,eb] = particlesperenergy()
  data = csvread('/Users/johannes/Dropbox/Bachelor/TrackingPositrons/Simulation/Data/particles_per_energy');
  x = [];
  y = [];
  eb = [];
  current = 0;
  start = 1;
  minimum_y = min(data(:,2));
  for i=1:length(data)
    if data(i,1) ~= current
       if current ~= 0
         x = [x;current];
         this_y = mean(data(start:i-1,2));
         y = [y;this_y];
         this_eb = std(data(start:i-1,2));
         if this_y - this_eb < 0
             this_eb = this_y - minimum_y;
         end
         eb = [eb;this_eb];
       end
       current = data(i,1);
       start = i;
    end
  end
  errorbar(x,y,eb,'-xr');
  set(gca,'XScale','log','YScale','log','FontSize',16);
  xlabel('Energy [MeV]');
  ylabel('Average generated particles');
  title('Particles generated per electron energy');
end
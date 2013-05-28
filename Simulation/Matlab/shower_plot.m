function [] = shower_plot()
    data = csvread('/Users/johannes/Dropbox/Bachelor/TrackingPositrons/Simulation/Data/bremsstrahlung_shower');
    % electron = [];
    % photon = [];
    n_electrons = 0;
    n_photons = 0;
    hold all;
    axis([0 4e2 -2e2 2e2]);
    xlabel('x [cm]');
    ylabel('y [cm]');
    for i=1:size(data,1)
        if data(i,1) == 11
            % electron = [electron;data(i,:)];
            h_e = plot([data(i,2);data(i,6)],[data(i,3);data(i,7)],'-r');
            n_electrons = n_electrons + 1;
        else
            % photon = [photon;data(i,:)];
            h_p = plot([data(i,2);data(i,6)],[data(i,3);data(i,7)],'-b');
            if n_photons == 0
                legend([h_e,h_p],{'Electron/positron tracks','Photon tracks'},'Location','NorthWest');
            end
            n_photons = n_photons + 1;
        end
    end
end
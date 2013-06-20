function [] = shower_plot()
    data = csvread('/Users/johannes/Dropbox/Bachelor/TrackingPositrons/Simulation/Data/bremsstrahlung_shower');
    disp('Plotting for data of size: ');
    disp(size(data));
    dim_plot(1,data,1,2);
    % dim_plot(2,data,1,3);
    dim_plot(2,data,2,3);
end

function [] = dim_plot(sub,data,dim1,dim2)
    subplot(1,2,sub);
    n_electrons = 0;
    n_photons = 0;
    hold all;
    axis1 = [min([data(:,1+dim1);data(:,5+dim1)]) max([data(:,1+dim1);data(:,5+dim1)])];
    axis2 = [min([data(:,1+dim2);data(:,5+dim2)]) max([data(:,1+dim2);data(:,5+dim2)])];
    max1 = max(abs(axis1));
    max2 = max(abs(axis2));
    axis1 = [-max1 max1];
    axis2 = [-max2 max2];
    title('Start and end vertices for bremsstrahlung system','FontSize',16);
    if dim1 == 1
        xlabel('x [cm]','FontSize',16);
        if axis1(1,1) < 0
            axis1(1,1) = 0;
        end
%         if axis1(1,2) > 1000
%             axis1(1,2) = 1000;
%         end
        axis1(1,1) = 80;
    elseif dim1 == 2
        xlabel('y [cm]','FontSize',16);
    else
        xlabel('z [cm]','FontSize',16);
    end
    if dim2 == 1
        ylabel('x [cm]','FontSize',16);
        if axis2(1,1) < 0
            axis2(1,1) = 0;
        end
%         if axis2(1,2) > 1000
%             axis2(1,2) = 1000;
%         end
    elseif dim2 == 2
        ylabel('y [cm]','FontSize',16);
    else
        ylabel('z [cm]','FontSize',16);
    end
    set(gca,'FontSize',14);
    for i=1:size(data,1)
        if data(i,1) == 11
            % electron = [electron;data(i,:)];
            h_e = plot([data(i,1+dim1);data(i,5+dim1)],[data(i,1+dim2);data(i,5+dim2)],'-r');
            %n_electrons = n_electrons + 1;
        else
            % photon = [photon;data(i,:)];
            h_p = plot([data(i,1+dim1);data(i,5+dim1)],[data(i,1+dim2);data(i,5+dim2)],'-b');
            %n_photons = n_photons + 1;
        end
    end
    axis([axis1,axis2]);
    if dim1 == 1
        l=legend([h_e,h_p],{'Electron/positron tracks','Photon tracks'},'Location','SouthWest');
    else
        l=legend([h_e,h_p],{'Electron/positron tracks','Photon tracks'},'Location','NorthWest');
    end
    set(l,'FontSize',16);
end
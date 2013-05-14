h = scatter3(B(:,3),B(:,2),B(:,1),100,B(:,4),'.');
axis([min(B(:,3)) max(B(:,3)) -1 max(B(:,2)) min(B(:,1)) max(B(:,1))]);
set(gca,'YDir','Reverse');
xlabel('z [cm]');
ylabel('y [cm]');
zlabel('x [cm]');
colorbar;
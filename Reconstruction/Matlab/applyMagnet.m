function [di,dj] = applyMagnet(i,j)

B = [0 0 1];

z = 0; % dummy for the cross product.

di = (j * B(3) - z * B(2)) * dt;
dj = (z * B(1) - i * B(3)) * dt;

return;
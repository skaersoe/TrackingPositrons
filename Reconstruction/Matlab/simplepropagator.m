% A simple particle propagator

% Global variables

global dt = 0.001;
    
global MAX_X = 1000;

global MIN_X = 0;

global MAX_Y = 500;

global MIN_Y = -500;

% Initialize the state, assume unity mass and charge

x = 0;

y = 0;

spread = 0.01;

p = normrnd(100,1);

j = spread*rand()*p;

i = (1-j)*p;

r = [x y i j p]';

% Kalman variables

%Propagator code

while (outOfBounds(r(1),r(2)) == FALSE)
    
    

axis tight;
plot(realpos(:,1),realpos(:,2),'-',measuredpos(:,1),measuredpos(:,2),'.');
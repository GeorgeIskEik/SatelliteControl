function [ydot] = Omega(t,y)
%OMEGA Summary of this function goes here
%   Detailed explanation goes here

% Define the angular velocity as a function of time 
omega = [sin(0.1*t) 0.01 cos(0.1*t)]'*20*pi/180;

% write the B transformatio matrix
B = [-y(2) -y(3) -y(4);
    y(1) -y(4) y(3);
    y(4) y(1) -y(2);
    -y(3) y(2) y(1)];
% quat rates
ydot = (0.5.*B)*omega;
end


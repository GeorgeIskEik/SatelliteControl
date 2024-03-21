
clear all
% Initial Euler angles (yaw, pitch, roll) in degrees
psi = deg2rad(40);   % yaw
theta = deg2rad(30); % pitch
phi = deg2rad(80);   % roll

% Time parameters
t0 = 0;      % Initial time
tf = 60;     % Final time (1 minute)
dt = 0.01;   % Time step size

% Initialize time array
t = t0:dt:tf;

% Initialize arrays to store Euler angles
psi_array = zeros(size(t));
theta_array = zeros(size(t));
phi_array = zeros(size(t));

% Body angular velocity components in the body frame (B frame)
B_omega = @(t) [sin(0.1*t); 0.01; cos(0.1*t)] * deg2rad(20); % Function handle for B_omega

% Numerical integration using Euler's method
for i = 1:length(t)
    % Store current Euler angles
    psi_array(i) = psi;
    theta_array(i) = theta;
    phi_array(i) = phi;
    
    % Compute time derivative of Euler angles using angular velocity
    psi_dot = B_omega(t(i))' * [1; 0; -sin(theta)]; % Derivative of yaw
    theta_dot = B_omega(t(i))' * [0; cos(phi); sin(phi)*cos(theta)]; % Derivative of pitch
    phi_dot = B_omega(t(i))' * [0; -sin(phi); cos(phi)*cos(theta)]; % Derivative of roll
    
    % Update Euler angles using Euler's method
    psi = psi + psi_dot * dt;
    theta = theta + theta_dot * dt;
    phi = phi + phi_dot * dt;
end

% Calculate Euler angle norm at time step 42s
psi_42 = psi_array(4201); % Index corresponding to 42s
theta_42 = theta_array(4201);
phi_42 = phi_array(4201);
euler_angle_norm = sqrt(psi_42^2 + theta_42^2 + phi_42^2);

% Display Euler angle norm at time step 42s
disp(['Euler angle norm at time step 42s: ', num2str(euler_angle_norm)]);
%%
% Initial EP vector
beta_0 = [0.408248; 0; 0.408248; 0.816497];

% Time parameters
t0 = 0;      % Initial time
tf = 60;     % Final time (1 minute)
dt = 0.01;   % Time step size

% Initialize time array
t = t0:dt:tf;

% Initialize array to store the EP vector components
beta_array = zeros(4, length(t));
beta_array(:,1) = beta_0;

% Body angular velocity components in the body frame (B frame)
B_omega = @(t) [sin(0.1*t); 0.01; cos(0.1*t)] * deg2rad(20); % Function handle for B_omega

% Numerical integration using Euler's method
for i = 2:length(t)
    % Compute time derivative of EP vector using angular velocity
    beta_dot = 0.5 * [0, -B_omega(t(i)); B_omega(t(i)), -skew_symmetric(beta_array(1:3,i-1))] * beta_array(:,i-1);
    
    % Update EP vector using Euler's method
    beta_array(:,i) = beta_array(:,i-1) + beta_dot * dt;
    
    % Normalize EP vector to ensure unit quaternion
    beta_array(:,i) = beta_array(:,i) / norm(beta_array(:,i));
end

% Find the index corresponding to 42 seconds
index_42s = find(t == 42);

% Calculate the norm of the EP vector component at 42 seconds
beta_norm_squared = sqrt(beta_array(1,index_42s)^2 + beta_array(2,index_42s)^2 + beta_array(3,index_42s)^2);

% Display the norm of the EP vector component at 42 seconds
disp(['Norm squared of the EP vector component at 42 seconds: ', num2str(beta_norm_squared)]);
%%
% Define the body angular velocity function
B_omega = @(t) [sin(0.1*t); 0.01; cos(0.1*t)] * deg2rad(20); % B frame components

% Define the EP differential kinematic equations
matrix = [   -beta(2),        -beta(3),        -beta(4);
               0,               beta(4),         -beta(3);
             -beta(4),         0,               beta(2);
               beta(3),        -beta(2),         0];
ep_dynamics = @(t, beta) 0.5 * matrix * B_omega(t);

% Initial EP vector
beta_init = [0.408248; 0; 0.408248; 0.816497]; % Initial EP vector

% Time span
tspan = [0, 60]; % Time span from 0 to 60 seconds

% Solve the differential equations using ode23
[t, beta] = ode23(ep_dynamics, tspan, beta_init);

% Find the index corresponding to 42 seconds
index_42s = find(t >= 42, 1);

% Calculate the norm of the EP vector component at 42 seconds
beta_norm_squared = beta(index_42s, 1)^2 + beta(index_42s, 2)^2 + beta(index_42s, 3)^2;

% Display the norm of the EP vector component at 42 seconds
disp(['Norm squared of the EP vector component at 42 seconds: ', num2str(beta_norm_squared)]);
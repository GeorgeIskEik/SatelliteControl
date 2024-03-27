% Numerical time integrator 

xinit = [deg2rad(40);deg2rad(30);deg2rad(80)];           %initial angle values 
tfinal = 60;                                             %(seconds)
xn = xinit;                                              %current step

% vectors for storring the values of each angle rate 
storage1 = zeros(1,tfinal);
storage2 = zeros(1,tfinal);
storage3 = zeros(1,tfinal);

% integration
for t = 1: 1 : tfinal
    f = Function(0.5236, 1.3963,t);

    x_nxt = xn+f;

    xn = x_nxt
    if t==42
        norm = sqrt(xn(1)^2+xn(2)^2+xn(3)^2)
    end
    storage1(t) = xn(1);
    storage2(t) = xn(2);
    storage3(t) = xn(3);
end

% plotting
t=1:tfinal;
figure,
subplot(311)
plot(t,storage1)
subplot(312)
plot(t,storage2)
subplot(313)
plot(t,storage3)



%%
clear all
% Initial Euler Parameters
q_init = [0.408248,0.,0.408248,0.816497];

% Time parameters
t0 = 0;      % Initial time
tf = 70;     % Final time (1 minute)
dt = 0.01;   % Time step size

% Initialize time array
t = t0:dt:tf;

% Initialize arrays to store Euler angles
beta_array = zeros(4,length(t));
beta_array(:,1) = q_init';
% Body angular velocity components in the body frame (B frame)
B_omega = @(t) [sin(0.1*t); 0.01; cos(0.1*t)] * deg2rad(20); % Function handle for B_omega

% Numerical integration using Euler's method
for i = 1:length(t)
    % Store current Euler angles
    if i == 1
        b0_array(i) = beta_array(1,i);
        b1_array(i) = beta_array(2,i);
        b2_array(i) = beta_array(3,i);
        b3_array(i) = beta_array(4,i);
    else
        b0_array(i) = beta_array(1,i-1);
        b1_array(i) = beta_array(2,i-1);
        b2_array(i) = beta_array(3,i-1);
        b3_array(i) = beta_array(4,i-1);
    end
    updateMatrix = [ -b1_array(i), -b2_array(i), -b3_array(i);
                    b0_array(i), -b3_array(i), b2_array(i);
                    b3_array(i), b0_array(i), -b1_array(i);
                    -b2_array(i), b1_array(i), b0_array(i)];
    % Compute time derivative of Euler angles using angular velocity
    beta_dot = 0.5*updateMatrix*B_omega(t(i));
    % Update parameters using Euler's method
    beta_array(:,i) = beta_array(:,i) + beta_dot * dt;
    % Normalize EP vector to ensure unit quaternion
    beta_array(:,i) = beta_array(:,i) / norm(beta_array(:,i));
end

figure 
plot(t,beta_array(2,:))
% Find the index corresponding to 42 seconds
index_42s = find(t == 42);

% Calculate the norm of the EP vector component at 42 seconds
beta_norm = sqrt(beta_array(2,index_42s)^2 + beta_array(3,index_42s)^2 + beta_array(4,index_42s)^2)

%% Ode45 integration EP
clear all
y0 = [0.408248,0.,0.408248,0.816497];
% set the absolute and relative tolerance accordingly
options = odeset('RelTol', 1e-20, 'AbsTol', 1e-22);
[t,y] = ode45(@Omega,[0 60],y0,options);
beta_0 = y(:,1);
beta_1 = y(:,2);
beta_2 = y(:,3);
beta_3 = y(:,4);
theta = (2*acos(beta_0))*180/pi;
figure
subplot(4,1,1)
plot(t,theta,'linewidth',2)
legend("b_0")
subplot(4,1,2)
plot(t,beta_1,'linewidth',2)
legend("b_1")
subplot(4,1,3)
plot(t,beta_2,'linewidth',2)
legend("b_2")
subplot(4,1,4)
plot(t,beta_3,'linewidth',2)
legend("b_3")
norm = sqrt(0.396537^2+0.587116^2+0.413106^2)
% correct for 0.8201

%% Ode45 integration CRP
clear all
% initial CRP's 
y0 = [0.4,0.2,-0.1];

% set the absolute and relative tolerance accordingly
options = odeset('RelTol', 1e-20, 'AbsTol', 1e-22);
[t,y] = ode45(@OmegaCRP,[0 60],y0,options);

q_1 = y(:,1);
q_2 = y(:,2);
q_3 = y(:,3);

% plot results.
figure
subplot(3,1,1)
plot(t,q_1,'linewidth',2)
legend("q_1")
subplot(3,1,2)
plot(t,q_2,'linewidth',2)
legend("q_2")
subplot(3,1,3)
plot(t,q_3,'linewidth',2)
legend("b_3")

% find the closest index to the value of 42.
[~, index] = min(abs(t - 42));
% calculate norm @ specified value.
norm = sqrt(y(index,1)^2+y(index,2)^2+y(index,3)^2)
% correct for 1.2.

%%
%FUNCTION Mapps the measured body angular velocities to the angle rates 
function rates = Function(theta,Phi,t)

Omega = [sin(0.1*t); 0.01; cos(0.1*t)]*0.3491;
matrix = [0, sin(Phi), cos(Phi);
          0, cos(Phi)*cos(theta), -cos(theta)*sin(Phi);
          cos(theta), sin(Phi)*sin(theta), cos(Phi)*sin(theta)];
rates = 1/cos(theta)*matrix*Omega;

end
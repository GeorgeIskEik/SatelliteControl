load 'Final model.mat'

mu = 398600.5;
re = 6378137e-3;
f = 1/298.257223563;
rp = re - f*re;

ra = 400;
rc = (re + rp)/2 + ra;
w_0 = sqrt(mu/rc^3)

%%
i_x = 0.0021; i_y = 0.0023; i_z = 0.0024;
Ib = diag([0.0021; 0.0023;0.0024])
I_s = 1.464*10^-5;
%%
global startup;
startup = false;
global startup1;
startup1 = false;

wsinit = zeros(3,1);
ws1 = wsinit(1);
ws2 = wsinit(2);
ws3 = wsinit(3);
J = Ib - A3*I_s*A3.'

J_inv = inv(J)

Initq_bob = [1,0,0,0];
n = Initq_bob(:,1);
e1 = Initq_bob(:,2);
e2 = Initq_bob(:,3);
e3 = Initq_bob(:,4);

R_bo_init = subs(R_b_o,[n,e_1,e_2,e_3],Initq_bob);

c_2 = R_bo_init(:,2);
c_3 = R_bo_init(:,3);

taua = [0; 0; 0];

theta_initd = [60; 60; 60];
theta_initr = deg2rad(theta_initd)
q = eul2quat(theta_initr')

Wbob = [1.4; 1; 1.5];
Ws = wsinit;
k_e = 0.3;
D = eye(3).*k_e;
e = q(:,2:4)
k2 = 2
%%
V_dot = w_0.*Wbob.*cross(c_2,(I_s.*A3*(A3*Wbob+Ws)))-Wbob.*k_e.*e'+D*Wbob+w_0.*cross(c_2,(I_s.*A3*(A3*Wbob+Ws)))+k2.*Wbob.*e'
%%

% q = eul2quat(theta_initr.')
% e = q(:,2:4)
% e1 = e(1);
% e2 = e(2);
% e3 = e(3);
% dot_n = q(1)
% R_bo = [1-2*e2^2-2*e3^2 2*e1*e2-2*e3*dot_n 2*e1*e3+2*e2*dot_n;
%          2*e1*e2+2*e3*dot_n 1-2*e1^2-2*e3^2 2*e2*e3-2*e1*dot_n;
%          2*e1*e3-2*e2*dot_n 2*e2*e3+2*e1*dot_n 1-2*e1^2-2*e2^2];
% 
% c_2 = R_bo(:,2)
% c_3 = R_bo(:,3)
% Wbob = [14; 10; 15];
% Ws = wsinit;
% k_e = 10;
% D = eye(3).*k_e;
% taua = k_e.*e.' + D*Wbob + w_0.*cross(c_2,A3*I_s*(A3.'*Wbob+Ws)) - w_0^2*cross(c_2,Ib*c_2)+3*w_0^2*cross(c_3,Ib*c_3)

function [dotWbob,dotWs,quart_bo,R_bo] = SystemNonlinearModel(taua)
% Function that calculates the system states  
%   This system function receives initializing input parameters and control
%   signals and calculates the components of the ordinal differential
%   equations of each state.

A = eye(3);
global Ib J 
M1 = eye(3);


wbib = wbob - w_0*c_2;

% The angular velocity of the satellite ode

fhi = J\(-cross(wbib,(Ib*wbib+A(:,1)*I_s.*ws1+A(:,2)*I_s.*ws2+A(:,3)*I_s.*ws3)));%inertial component equation
fht = -J\A*taua3;                     %input torque equation
fhg = J\(3*w0^2*cross(c_3,Ib*c_3));     %gravitational component equation
dotc = -cross(wbob,c_2);                
fha = dotc*(w0);                      %additional (orbital velocity component)

dotWbob = fhi+fht+fhg+dotc+fha;
% The reaction wheel angular velocity ode

fsi = A.'*(J\(cross(wbib,(Ib*wbib+A(:,1)*I_s.*ws1+A(:,2)*I_s.*ws2+A(:,3)*I_s.*ws3))));%mod to accomodate ws
fst = (A.'*(J\A)+1/I_s)*taua;
fsg = -A.'*(J\(3*w0^2*cross(c_3,Ib*c_3)));

dotWs = fsi+fst+fsg;
ws1 = dotWs(1);
ws2 = dotWs(2);
ws3 = dotWs(3);

% update quaternion positioning

dot_n = -1/2*e.'*dotWbob;
dot_e = 1/2*(n*M1*dotWbob + cross(e,dotWbob));

quart_bo = [dot_n; dot_e];
e1 = dot_e(1,:);
e2 = dot_e(2,:);
e3 = dot_e(3,:);

% update rotation matrix

R_bo = [1-2*e2^2-2*e3^2 2*e1*e2-2*e3*dot_n 2*e1*e3+2*e2*dot_n;
         2*e1*e2+2*e3*dot_n 1-2*e1^2-2*e3^2 2*e2*e3-2*e1*dot_n;
         2*e1*e3-2*e2*dot_n 2*e2*e3+2*e1*dot_n 1-2*e1^2-2*e2^2];
c_2 = R_bo(:,2);
c_3 = R_bo(:,3);
end


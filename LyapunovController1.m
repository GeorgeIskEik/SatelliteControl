function [taua] = LyapunovController1(Reference,R_bo,w_0,Feedback)
%Lyapunov controller function 
%   Control strategy based on a Lyapunov function candidate. The core idea
%   is that the reaction wheel angular velocity is interpreted as an
%   external control input and is therefore not included in the statespace
%   of the system.

dotWbob = Feedback(1);
dotWs = Feedback(2);
quart_bo = Feedback(3);

% Angular error rate
deltaW = dotWbob-refW;
deltaWdot = -Ib^-1*deltaW*cross(refW,Ib)-corss(deltaW,(Ib.*refW))-I^-1*cross(deltaW,(Ib.*deltaW))+Ib^-1*u-refWdot-Ib^-1*cross(refW,(Ib.*dotWbob));
%
% Quaternion error
qref = Initq_bob;
qref = eul2quat(Reference);
refinv = inv(qref);
deltaq = quatmultiply(quart_bo,refinv);
deltaq_4 = deltaq(1);                   %deltaq_4 is in fact the first element but is named for math model referencing
deltaq_v = deltaq(2:4);

deltaq_vdot = (1/2).*(-cross(deltaW,deltaq_v)-2.*cross(deltaW,deltaq_v)+deltaq_4*deltaW);
deltaq_4dot = -(1/2).*deltaW.'.*deltaq_v;
deltaq_dot = [deltaq_4dot,deltaq_vdot];

%



dot_n = deltaq_dot(1);
e1 = deltaq_dot(2);
e2 = deltaq_dot(3);
e3 = deltaq_dot(4);

% calculate error (REFERENCE-FEEDBACK)

% Look into quaternion error

% Update DCM with the calculated error parameters
R_bo = [1-2*e2^2-2*e3^2 2*e1*e2-2*e3*dot_n 2*e1*e3+2*e2*dot_n;
         2*e1*e2+2*e3*dot_n 1-2*e1^2-2*e3^2 2*e2*e3-2*e1*dot_n;
         2*e1*e3-2*e2*dot_n 2*e2*e3+2*e1*dot_n 1-2*e1^2-2*e2^2];

c_2 = R_bo(:,2);
c_3 = R_bo(:,3);

taua = k_e.*deltaq_vdot.' + D*deltaWdot + w_0.*cross(c_2,A3*I_s*(A3.'*deltaWdot+dotWs)) - w_0^2*cross(c_2,Ib*c_2)+3*w_0^2*cross(c_3,Ib*c_3);
end


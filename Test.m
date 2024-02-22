eul = [60;60;60]
Initq_bob = [1 0 0 0]
trans = deg2rad(eul)
Reference = trans
Ws = zeros(3,1);
wbob_init=[1.4; 1; 1.5];
% quaternions = eul2quat(Reference')
% refinv = quatinv(quaternions)
% deltaq = quatmultiply(quart_bo,refinv)
refW = [12; 13; 5]
k_e = 0.3;
D = eye(3).*k_e; 
startup = true;
quaternions = zeros(1,4);
AngVel = zeros(3,1);
WAngVel = zeros(3,1);

%% Controller (Reference, R_bo, w_0, Feedback, Initq_bob, Ib, refW, k_e, D, A3, I_s, startup1)
u = zeros(3,1);
M = ones(3,1);

% Angular error rate
refWdot = -Ib^-1*cross(refW,(Ib*refW))+Ib^-1*u;
deltaW = dotWbob-refW;
deltaWdot = -Ib^-1*deltaW.*cross(refW,Ib*M)-cross(deltaW,(Ib*refW))-Ib^-1*cross(deltaW,(Ib*deltaW))+Ib^-1*u-refWdot-Ib^-1*cross(refW,(Ib*dotWbob));
%
% Quaternion error

qref = eul2quat(Reference');
refinv = quatinv(qref);
deltaq = quatmultiply(quart_bo',refinv)
deltaq_4 = deltaq(1);                   %deltaq_4 is in fact the first element but is named for math model analogy
deltaq_v = deltaq(2:4);

deltaq_vdot = (1/2).*(-cross(deltaW,deltaq_v')-2.*cross(deltaW,deltaq_v')+deltaq_4'*deltaW)
deltaq_4dot = -(1/2).*deltaW'*deltaq_v'
deltaq_dot = [deltaq_4dot deltaq_vdot']
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

taua = k_e.*deltaq_vdot + D*deltaWdot + w_0.*cross(c_2,A3*I_s*(A3'*deltaWdot+dotWs)) - w_0^2*cross(c_2,Ib*c_2)+3*w_0^2*cross(c_3,Ib*c_3)
u = taua;

%% System (taua,Init,startup,w_0,J,A3,Ib,I_s)


% Decide whether to use initial values or the one calculated in the
% previous itteration
if(startup == true)
   
    n = Initq_bob(1);
    e1 = Initq_bob(2);
    e2 = Initq_bob(3);
    e3 = Initq_bob(4);

    R_bo = [1-2*e2^2-2*e3^2 2*e1*e2-2*e3*n 2*e1*e3+2*e2*n;
        2*e1*e2+2*e3*n 1-2*e1^2-2*e3^2 2*e2*e3-2*e1*n;
        2*e1*e3-2*e2*n 2*e2*e3+2*e1*n 1-2*e1^2-2*e2^2];
    c_2 = R_bo(:,2);
    c_3 = R_bo(:,3);

    ws1 = Ws(1);
    ws2 = Ws(2);
    ws3 = Ws(3);
    wbob = wbob_init;


else
    q_bob = deltaq_dot;
    n = q_bob(1);
    e1 = q_bob(2);
    e2 = q_bob(3);
    e3 = q_bob(4);

    R_bo = [1-2*e2^2-2*e3^2 2*e1*e2-2*e3*n 2*e1*e3+2*e2*n;
        2*e1*e2+2*e3*n 1-2*e1^2-2*e3^2 2*e2*e3-2*e1*n;
        2*e1*e3-2*e2*n 2*e2*e3+2*e1*n 1-2*e1^2-2*e2^2]
    c_2 = R_bo(:,2);
    c_3 = R_bo(:,3);

    ws1 = WAngVel(1);
    ws2 = WAngVel(2);
    ws3 = WAngVel(3);
    wbob = AngVel;
end
%EndInit

M1 = eye(3);
wbib = wbob - w_0*c_2;

% The angular velocity of the satellite ode

fhi = J\(-cross(wbib,(Ib*wbib+A3(:,1)*I_s.*ws1+A3(:,2)*I_s.*ws2+A3(:,3)*I_s.*ws3)));%inertial component equation
fht = -J\A3*taua;                     %input torque equation
fhg = J\(3*w_0^2*cross(c_3,Ib*c_3));     %gravitational component equation
dotc = -cross(wbob,c_2)
fha = dotc*(w_0);                      %additional (orbital velocity component)

dotWbob = fhi+fht+fhg+dotc+fha
% The reaction wheel angular velocity ode

fsi = A3'*(J\(cross(wbib,(Ib*wbib+(A3(:,1)*I_s.*ws1+A3(:,2)*I_s.*ws2+A3(:,3)*I_s.*ws3)))))%mod to accomodate ws
fst = (A3'*(J\A3)+1/I_s)*taua;
fsg = -A3'*(J\(3*w_0^2*cross(c_3,Ib*c_3)));

dotWs = fsi+fst+fsg;
ws1 = dotWs(1);
ws2 = dotWs(2);
ws3 = dotWs(3);

% update quaternion positioning
e = [e1; e2; e3]
dot_n = -1/2*e'*dotWbob;
dot_e = 1/2*(n*M1*dotWbob + cross(e,dotWbob));

quart_bo = [dot_n; dot_e];
e1 = dot_e(1);
e2 = dot_e(2);
e3 = dot_e(3);

% update rotation matrix

R_bo = [1-2*e2^2-2*e3^2 2*e1*e2-2*e3*dot_n 2*e1*e3+2*e2*dot_n;
    2*e1*e2+2*e3*dot_n 1-2*e1^2-2*e3^2 2*e2*e3-2*e1*dot_n;
    2*e1*e3-2*e2*dot_n 2*e2*e3+2*e1*dot_n 1-2*e1^2-2*e2^2];

% update calculated data

 quaternions = quart_bo;
 AngVel = dotWbob;
 WAngVel = dotWs;

%update flag
startup = false;

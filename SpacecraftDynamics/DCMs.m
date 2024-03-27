%%
b1 = 1/3*[1;2;-2]
b2 = 1/sqrt(2)*[0;1;1]
b3 = 1/(3*sqrt(2))*[4;-1;1]

f1 = 1/4*[3;-2;sqrt(3)]
f2 = 1/2*[-1;0;sqrt(3)]
f3 = -1/4*[sqrt(3);2*sqrt(3);1]

rad2deg(acos(b1))
rad2deg(acos(b2))
rad2deg(acos(b3))
rad2deg(acos(f1))
rad2deg(acos(f2))
rad2deg(acos(f3))

BN = [b1'; b2'; b3']
FN = [f1'; f2'; f3']
BF = BN*FN'
%% Init 
% Angle set 1
ang1 = deg2rad(20);
ang2 = deg2rad(10);
ang3 = deg2rad(-10);
matrix1 = M1(ang3);
matrix2 = M2(ang2);
matrix3 = M3(ang1);

matrix1s = M1s();
matrix2s = M2s();
matrix3s = M3s();
%% Angle set 2
alpha1 = deg2rad(-5);
alpha2 = deg2rad(5);
alpha3 = deg2rad(5);

inter1 = M1(alpha3);
inter2 = M2(alpha2);
inter3 = M3(alpha1);
%% Compute the 3-2-1 sequence DCM
ThTwO = matrix1*matrix2*matrix3
ThTwOs = matrix1s*matrix2s*matrix3s
OTwTh = matrix3*matrix2*matrix1;
OTwThs = matrix3s*matrix2s*matrix1s;
ThOThs = matrix3s*matrix1s*matrix3s;
ThOTh = matrix3*matrix1*matrix3;
TwThTws = matrix2s*matrix3s*matrix2s
%%
RN = inter1*inter2*inter3

% operation BR = BN*RN'
BR = ThTwO*RN'
%% Extract angles for 3-1-3 sequence
%trans = ThTwO';
PSI = atan2(BR(1,2),BR(1,1))
THETA = -asin(BR(1,3))
PHI = atan2(BR(2,3),BR(3,3))
norm1 = sqrt(PSI^2+THETA^2+PHI^2)
rad2deg(PSI)
rad2deg(THETA)
rad2deg(PHI)


%

Omega = atan2(ThTwO(3,1),-ThTwO(3,2))
i = acos(ThTwO(3,3))
omega = atan2(ThTwO(1,3),ThTwO(2,3))
%norm2 = sqrt(Omega^2+i^2+omega^2)
%angles = rotm2eul(BR,'')
rad2deg(Omega)
rad2deg(i)
rad2deg(omega)
%%
syms L r fi
theta = -(L*fi)/r
BN = M1s()*M3s()
Nb = BN.'
NB = [cos(fi), -sin(fi)*cos(theta), sin(fi)*sin(theta);
        sin(fi), cos(fi)*cos(theta), -cos(fi)*sin(theta);
        0, sin(theta), cos(theta)]
bv = [2, 1, 1]
nv = NB*bv'
%%

rot = [0.925417, 0.336824, 0.173648;
       0.0296956, -0.521281, 0.852869;
       0.377786, -0.784102, -0.492404]
theta = acos(1/2*(ThTwO(1,1)+ThTwO(2,2)+ThTwO(3,3)-1))
trans = [ThTwO(2,3)-ThTwO(3,2); ThTwO(3,1)-ThTwO(1,3); ThTwO(1,2)-ThTwO(2,1)];
e_inter = 1/2*(sin(theta))*trans;
e = e_inter/norm(e_inter)  

%% form DCM from quat sequence
qsec = [0.235702,0.471405,-0.471405,0.707107];
b0 = qsec(1);
b1 = qsec(2);
b2 = qsec(3);
b3 = qsec(4);
DCM = [b0^2+b1^2-b2^2-b3^2, 2*(b1*b2+b0*b3), 2*(b1*b3-b0*b2);
    2*(b1*b2-b0*b3), b0^2-b1^2+b2^2-b3^2, 2*(b2*b3+b0*b1);
    2*(b1*b3+b0*b2), 2*(b2*b3-b0*b1), b0^2-b1^2-b2^2+b3^2]

dcm = quat2dcm(qsec)
%% Given rotation 3-2-1 EA(20,10,-10) waht is the eq EP set?

ang1 = deg2rad(20);
ang2 = deg2rad(10);
ang3 = deg2rad(-10);
matrix1 = M1(ang3);
matrix2 = M2(ang2);
matrix3 = M3(ang1);

matrix1s = M1s();
matrix2s = M2s();
matrix3s = M3s();

ThTwO = matrix1*matrix2*matrix3;
ThTwOs = matrix1s*matrix2s*matrix3s
quat = DCM2EP(ThTwO)

%% Adding euler parameters 
clear all
% quaternion1 B relative to N
BN = [0.774597; 0.258199; 0.516398; 0.258199];

% quaternion1 F relative to B
FB = [0.359211; 0.898027; 0.179605; 0.179605];

FN = quatAdd(FB)*BN
% q = quatmultiply(FB',BN')
% rad2deg(2*acos(q(1)))
% rad2deg(2*acos(-q(1)))
if rad2deg(2*acos(FN(1)))<0
    FN(1) = -FN(1)
end
%possible answers :0.709738, 0.528365, -0.146008, -0.373204
%%
clear all

FN = [0.359211,0.898027,0.179605,0.179605];
BN = [-0.377964,0.755929,0.377964,0.377964];
NB = [-0.377964;-0.755929;-0.377964;-0.377964]
FB = quatAdd(FN)*NB

rad2deg(2*acos(FB(1)))
rad2deg(2*acos(-FB(1)))
  
%%
FN = [0.359211,0.898027,0.179605,0.179605];
rho = FN(2)/FN(1)
%%
clear all
q_bn = [0.1; 0.2; 0.3];
dcm = crp2dcm(q_bn)

%% DCM2CRP : DCM->EP->CRP
clear all

DCM = [0.333333 -0.666667 0.666667;
       0.871795 0.487179 0.0512821;
       -0.358974 0.564103 0.74359];
% transform DCM to quaternions
quat = DCM2EP(DCM)
% extract CRP using the quat-CRP relationship
q1 = quat(2)/quat(1);
q2 = quat(3)/quat(1);
q3 = quat(4)/quat(1);
q = [q1,q2,q3]

%% Subtraction FB
clear all
q_fn = [0.1; 0.2; 0.3];
q_bn = [-0.3; 0.3; 0.1];
q_bf = (q_bn-q_fn+cross(q_bn,q_fn))/(1+q_bn'*q_fn)

%%
q_nf = -q_fn
dcm_f = crp2dcm(q_nf);
dcm_i = crp2dcm(q_bn);
dcm_bf = dcm_i*dcm_f;
quat = DCM2EP(dcm_bf)

q1 = quat(2)/quat(1);
q2 = quat(3)/quat(1);
q3 = quat(4)/quat(1);
q = [q1,q2,q3]
%%
sig = [0.1, 0.2, 0.3];

% flip to shadow set
sig_s = -sig/norm(sig'*sig)

%% Converting from MRP to DCM
clear all
sig =  [0.1; 0.2; 0.3];
cdm = mrp2dcm(sig)
% verified
%% Map the DCM to the equivalent MRP set (short route)
clear all
DCM = [0.763314, 0.0946746, -0.639053
        -0.568047, -0.372781, -0.733728
        -0.307692, 0.923077, -0.230769];

% Shepphard's method to go to quaternions
quat = DCM2EP(DCM);
sig = quat(2:4)/(1+quat(1));
sig_s = -sig/norm(sig'*sig);
% verified

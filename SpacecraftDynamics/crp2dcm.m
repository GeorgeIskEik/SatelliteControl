function [DCM] = crp2dcm(q)
%CRP2DCM Summary of this function goes here
%   Detailed explanation goes here
q1 = q(1);
q2 = q(2);
q3 = q(3);
matrix = [1+q1^2-q2^2-q3^2 2*(q1*q2+q3) 2*(q1*q3-q2);
            2*(q2*q1-q3) 1-q1^2+q2^2-q3^2 2*(q2*q3+q1);
            2*(q3*q1+q2) 2*(q3*q2-q1) 1-q1^2-q2^2+q3^2];
DCM = 1/(1+q'*q)*matrix;
end


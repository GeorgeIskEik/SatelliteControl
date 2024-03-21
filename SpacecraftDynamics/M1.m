function matrix1 = M1(ang)
%M1 Funtion that calculates the first rotation matrix
%   Takes the angle of rotation about the first axis and computes the
%   matrix for the first rotation in the DCM

matrix1 = [1 0 0;
     0  cos(ang) sin(ang);
     0 -sin(ang) cos(ang)];
end


function matrix1s = M1s()
%M1s Funtion that prints the symbolic first rotation matrix
%   Takes the angle of rotation about the first axis and prints the
%   symbolic matrix for the first rotation in the DCM
syms phi
matrix1s = [1 0 0;
     0  cos(phi) sin(phi);
     0 -sin(phi) cos(phi)];
end


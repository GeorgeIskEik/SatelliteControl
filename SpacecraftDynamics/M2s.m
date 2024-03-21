function matrix2s = M2s()
%M2s Funtion that prints the symbolic second rotation matrix
%   Takes the angle of rotation about the second axis and prints the
%   symbolic matrix for the second rotation in the DCM
syms theta
matrix2s = [cos(theta), 0, -sin(theta);
            0, 1, 0;
        sin(theta), 0, cos(theta)];
end


function matrix2 = M2(ang)
%M1 Funtion that calculates the second rotation matrix
%   Takes the angle of rotation about the second axis and computes the
%   matrix for the second rotation in the DCM
matrix2 = [cos(ang), 0, -sin(ang);
            0, 1, 0;
        sin(ang), 0, cos(ang)];
end


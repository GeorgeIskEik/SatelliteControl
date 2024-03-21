function matrix3 = M3(ang)
%M1 Funtion that calculates the third rotation matrix
%   Takes the angle of rotation about the third axis and computes the
%   matrix for the third rotation in the DCM

matrix3 = [cos(ang), sin(ang), 0;
      -sin(ang), cos(ang), 0;
      0,  0, 1];
end


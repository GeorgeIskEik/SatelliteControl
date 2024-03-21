function matrix3s = M3s()
%M3S Funtion that prints the symbolic third rotation matrix
%   Takes the angle of rotation about the third axis and prints the
%   symbolic matrix for the third rotation in the DCM
syms Psi
matrix3s = [cos(Psi), sin(Psi), 0;
      -sin(Psi), cos(Psi), 0;
      0,  0, 1];
end


function [matrix] = skew(vect,order)
%SKEW Summary of this function goes here
%   Detailed explanation goes here
switch order
    case 3
        e1 = vect(1);
        e2 = vect(2);
        e3 = vect(3);
        matrix = [0 -e3 e2;
                  e3 0 -e1;
                  -e2 e1 0];
    case 4
        e1 = vect(1);
        e2 = vect(2);
        e3 = vect(3);
        e4 = vect(4);
        matrix = [0 -e1 -e2 -e3;
                  e1 0 -e4 e2;
                  e2 e4 0 -e1;
                  e3 -e2 e1 0];
    otherwise
        disp('Not supported!')
end


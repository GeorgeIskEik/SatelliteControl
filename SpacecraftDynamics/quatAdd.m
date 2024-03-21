function [matrix] = quatAdd(q)
b0 = q(1);
b1 = q(2);
b2 = q(3);
b3 = q(4);
matrix = [b0, -b1, -b2, -b3;
            b1, b0, b3, -b2;
            b2, -b3, b0, b1;
            b3, b2, -b1, b0];
end


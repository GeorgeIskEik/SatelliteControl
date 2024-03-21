function [q] = DCM2EP(DCM)
max = 0;
b0c = sqrt(1/4*(1+trace(DCM)));
b1c = sqrt(1/4*(1+2*DCM(1,1)-trace(DCM)));
b2c = sqrt(1/4*(1+2*DCM(2,2)-trace(DCM)));
b3c = sqrt(1/4*(1+2*DCM(3,3)-trace(DCM)));
storage = [b0c,b1c,b2c,b3c];
for i = 1:length(storage)
    if storage(i) > max
        max = storage(i);
    end
end

switch max
    case b0c
        b0 = b0c;
        b1 = (DCM(2,3)-DCM(3,2))/(b0*4);
        b2 = (DCM(3,1)-DCM(1,3))/(b0*4);
        b3 = (DCM(1,2)-DCM(2,1))/(b0*4);
    case b1c
        b1 = b1c;
        b0 = (DCM(2,3)-DCM(3,2))/(b1*4);
        b2 = (DCM(1,2)+DCM(2,1))/(b1*4);
        b3 = (DCM(3,1)+DCM(1,3))/(b1*4);
    case b2c
        b2 = b2c;
        b0 = (DCM(3,1)-DCM(1,3))/(b2*4);
        b1 = (DCM(1,2)+DCM(2,1))/(b2*4);
        b3 = (DCM(2,3)+DCM(3,2))/(b2*4);
    case b3c
        b3 = b3c;
        b0 = (DCM(1,2)-DCM(2,1))/(b3*4);
        b1 = (DCM(3,1)+DCM(1,3))/(b3*4);
        b2 = (DCM(2,3)+DCM(3,2))/(b3*4);
end
q = [b0,b1,b2,b3];

if (pi-(2*acos(b0)))>0
    disp("shortest route");
else
    disp("Longest route");
    if (pi-(2*acos(b0))) == 0
        disp("ambiguous");
    end
end
end


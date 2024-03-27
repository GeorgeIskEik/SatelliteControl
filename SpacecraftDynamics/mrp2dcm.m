function [matrix] = mrp2dcm(sig)

%MRP2DCM converts a set of Modified Rodrigues Parameters into a direction
% cosine matrix form.

% matrix = eye(3)+(8*norm(skew(sig,3))^2-4*(1-sig'*sig)*norm(skew(sig,3)))/(1+sig'*sig).^2;
s1 = sig(1);
s2 = sig(2);
s3 = sig(3);

formator = [4*(s1^2-s2^2-s3^2) + (1-sig'*sig).^2, 8*s1*s2 + 4*s3*(1-sig'*sig), 8*s1*s3 - 4*s2*(1-sig'*sig);
            8*s2*s1 - 4*s3*(1-sig'*sig), 4*(-s1^2+s2^2-s3^2) + (1-sig'*sig).^2, 8*s2*s3 + 4*s1*(1-sig'*sig);
            8*s3*s1 + 4*s2*(1-sig'*sig), 8*s3*s2 - 4*s1*(1-sig'*sig), 4*(-s1^2-s2^2+s3^2) + (1-sig'*sig).^2];

matrix = ((1+sig'*sig).^2)^-1*formator;
end


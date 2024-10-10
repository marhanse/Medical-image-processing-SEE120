function f = dist(a1,a2,b1,b2)
%FTUNDIST Summary of this function goes here
%   Detailed explanation goes here
 load('Project5-3.mat')
 dx = a1*tmat + a2*tmat^2;
 dy = b1*tmat + b2*tmat^2;
 u = uu;
 v = vv; 
 tang1 = 2.*pi*(v.*dy + u.*dx);
 E = exp(i.*tang1);
 ft3 = ft.*E;
 imshow(ifft2(fftshift(ft3)), [0 400])
end


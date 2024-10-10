function [restoredIm] = imRestor(a1,a2,b1,b2,ftDist)
load Project5-3.mat
dx =a1*tmat + a2*tmat^2;
dy =tmat*b1 + b2*tmat^2;
tang2=-2*pi*(uu.*dx+vv.*dy);
E=exp(i*tang2);
restoredIm = ifft2(fftshift(ftDist./E));
end
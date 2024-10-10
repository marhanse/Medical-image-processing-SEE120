function f = undist(a1,a2,b1,b2,ftofimage)
%undist: Function to remove the distortion 

% a1 = constant for linear motion in x-direction 
% a2 = constant for acc in x - direction 
% b1 = constant for linear motion in y-direction 
% b2 = constant for acc in y - direction 
% ftofimage = ft of the distorted image

 load('Project5-3.mat')
 
 dx = a1*tmat + a2*tmat^2;
 dy = b1*tmat + b2*tmat^2;
 u = uu;
 v = vv; 
 tang1 = 2.*pi*(v.*dy + u.*dx);
 E = exp(-i.*tang1);
 f = (ftofimage./E); 
 %imshow(ifft2(fftshift(ft3)),[0 400]);
 title('Undistorted')
end
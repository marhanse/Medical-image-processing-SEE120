%% info 


%img – An undistorted image (size 512,512). Use imshow(img.[0, 400]) to see fine structure in this image.
%ft – A 512x512 matrix containing the measured Fourier components of the above undistorted image. This Fourier transform is centred with the zero spatial frequency measurement in pixel 257,257. This FT can be duplicated from the true image img using MATLAB commands fft2 and fftshift.
%uu – A matrix showing the u values of each of the 512x512 spatial frequencies measured in the above ft file. 
% vv – A matrix showing the v values of each of the 512x512 spatial frequencies measured in the above ft file.
%tmat – A matrix containing the time each Fourier component was observed during the scan from time -0.5 to +0.5 and with time equal to zero at central point of the matrix 257,257.
%tang – Assuming the subject moves 2cm at constant velocity in the vertical direction during the scan this matrix gives the phase error (in radians) introduced at each Fourier sample due to subject motion.
%img2 – This is the resulting distorted image due to the above image motion. This image is obtained following equation 3 by applying the phase errors in file tang to the true values of the FT of the image in file ft, then inverse Fourier transforming and taking the real part. Again use imshow(img2,[0, 400]) to see fine structure in the image; this time observing the ripples introduced by image motion.


%% Load image 
clear
clc

load Project5-3.mat
load Project5-challenges.mat

one = Proj5_Chal1_uv;
two = Proj5_Chal2_uv;
three = Proj5_Chal3_uv; 

disp('Enter (1) for linear y - motion (2) for linear x and y - motion (3) for linnear and acc (4) for test image img2');
Challenge=input('This image will be undistorted ');

if Challenge == 1 
    ftDist = one;
elseif Challenge == 2 
    ftDist = two;
elseif Challenge == 3 
    ftDist = three;
else
    ftDist = ft2; 
end
imDist = ifft2(fftshift(ftDist));

%[imDist,ftDist]=imDistort(7,3,10,2,ft); %create a distorted image  imDistort(a1,a2,b1,b2,myImg,ft); 


[numRows,numCols] = size(ftDist);%size of the image
len=min([numRows,numCols])/2;


a2=0;% set accelaration to zero for approx linear motion
b2=0;


% b1 
u=0;
phiV=[];
for v = 0:len-1
    ev=angle(ftDist(257-v,257-u))+angle(ftDist(257+v,257+u));% Angle of every point along v-axis 
    phiV = [phiV, ev];

end

DDV= conv2(phiV,0.5*[-1 2 -1],"same");%takes the seconde derivative 
Kb1=median(DDV);% Estimation of k. Ty phi proportional to kv^2 (11)
dv =vv(1,1)-vv(512,1); %change in dv
dt= tmat(1,1)-tmat(512,1); %change in dt
b=(Kb1*dv/(dt*2*pi)); % Formula gained from eqv (6), (10), (11)
b1=b*(numRows*numCols)/(100*2);% Converting from m/pixels^2 to cm/s 



optMethod=input('Optimal method for estimating a1?(Yes=1) (No=0) '); %Ask if the user wants to use the optimum method 

if optMethod==1
    
    %a1+b1 = a1b1 optimum method for singnal to noise ratio. 
    v=0;
    phiU=[];
    for u = 0:len-1
        eu=angle(ftDist(257-u,257+u))+angle(ftDist(257+u,257-u)); % Angle of every point along diagonal axis 
        phiU = [phiU, eu];
    end  
    DDU= conv2(phiU,0.5*[-1 2 -1],"same");
    Kab=median(DDU);% Estimation of k. Ty phi proportional to kv^2 (11)
    du =uu(1,512)-uu(512,1);
    dv =vv(1,512)-vv(512,1);
    dt=tmat(1,512)-tmat(512,1); %Change in dt
    ab=(Kab*(du+dv)*0.5/(dt*2*pi));% Formula gained from eqv (6), (10), (11)
    a1b1=ab*(numRows*numCols)/(100*2);% Converting from m/pixels^2 to cm/s 
    
    a1=a1b1-b1;
else
    % b1 
    v=0;
    phiU=[];
    for u = 0:len-1
        eu=angle(ftDist(257-v,257-u))+angle(ftDist(257+v,257+u));% Angle of every point along v-axis 
        phiU = [phiU, eu];
    
    end
    
    DDU= conv2(phiU,0.5*[-1 2 -1],"same");%takes the seconde derivative 
    Ka1=median(DDU);% Estimation of k. Ty phi proportional to ku^2 (11)
    du =uu(1,512)-uu(1,1); %change in du
    dt= tmat(1,512)-tmat(1,1); %change in dt
    a=(Ka1*du/(dt*2*pi)); % Formula gained from eqv (6), (10), (11)
    a1=a*(numRows*numCols)/(100*2);% Converting from m/pixels^2 to cm/s 
        
end
%% 


% test to restor image linear motion
dx =a1*tmat + a2*tmat^2;
dy =tmat*b1 + b2*tmat^2;
tang2=-2*pi*(uu.*dx+vv.*dy);
E=exp(i*tang2);


im_diff = ifft2(fftshift(ft)) - ifft2(fftshift(ftDist./E));
%std2(im_diff);
%% restor acceleration 

if std2(im_diff)<10^-3
    disp("This restoration is okey")
    disp(" ")
else % try to restore the image when guesing the acceleration of the movment
    disp("Linear restoration NOT, try with accelarton restoration")


    a2=0;% set a2=0, Because it given in the instructions
    b2=10; %Guess that it is 0 to try to calculate a2 first 
    okStd2=0.01; % when standardivation is lower then this the loop stops
    maxRun=40;%maximum of loops 
    step=2;
    stop=0;
    %b2
    disp("Calculate b2")
    valb=[0;0 ];
    while std2(ifft2(fftshift(ft)) - imRestor(a1,a2,b1,b2,ftDist))>okStd2 && length(valb)<maxRun && stop==0 %try to change so it stops when standard varitaion gets bigger again
        b2=b2-step;
        standDev1= std2(ifft2(fftshift(ft)) - imRestor(a1,a2,b1,b2,ftDist));
        standDev2= std2(ifft2(fftshift(ft)) - imRestor(a1,a2,b1,b2+step,ftDist));

         if standDev2<standDev1 %if the standard deviation do increase insted of decreasing
             b2=b2+step;
             step=step/2; %makes the step smaller to be more accurate
             b2=b2+step;
            if std2(ifft2(fftshift(ft)) - imRestor(a1,a2,b1,b2,ftDist))<okStd2 || length(valb)>maxRun
                stop=1;
                disp("break")
            end
             
         end

         valb=[valb(1,:),b2;valb(2,:), standDev1];      
         fprintf("StandDev of diff %0.6f ", std2(img - imRestor(a1,a2,b1,b2,ftDist)))
         disp(" ")


    end
end




%% ploting the result
dx =a1*tmat + a2*tmat^2;
dy =tmat*b1 + b2*tmat^2;

tang2=-2*pi*(uu.*dx+vv.*dy);
E=exp(i*tang2);

restoredIm = ifft2(fftshift(ftDist./E));

figure(2)
subplot(2,2,1)
imshow(real(img),[0 400]);
title('Original image') 

subplot(2,2,2)
imshow(real(imDist),[0 400]);
title('Distorted image') 

subplot(2,2,3)
imshow(real(restoredIm),[0 400]);
title('Image restoration')

subplot(2,2,4)
imshow(real(img-restoredIm),[]) % Displays the error image (true - restored) 
title('Error image')

Max_pixel_error = max(max(real(img-restoredIm))); % Max pixel error 
fprintf('The maximum pixel error is = %.5f',Max_pixel_error)
disp(" ")


standDev= std2(ifft2(fftshift(ft)) - restoredIm);
fprintf('Std2 of the diffrenc bettwen real and restored = %.5f',standDev)
disp("  ")
fprintf('a1= %.2f a2= %.2f b1= %.2f b2= %.2f', a1,a2,b1,b2);
disp(" ")

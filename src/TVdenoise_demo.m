%TVdenoise_demo  Demo script for TV denoising algorithm applied to EAR
%
% This script illustrates the use of the TV denoising algorithm
% implemented in the function TVdenoise; it produces Fig. 2 in the
% accompanying paper. The user can easily modify the script for
% other images and noise levels.
%
% The script loads a clean image and adds Gaussian noise, and then
% computes two TV reconstructions for two different values of tau.  Note
% that the computation of each TV reconstruction takes several seconds.

% J. Dahl^1, P.C. Hansen^2, S.H. Jensen^1 & T.L. Jensen^1
% CSI project: (1) Aalborg University, (2) Technical University of Denmark
% April 28, 2009.

clear, clc
close all
disp('Starting TVdenoise_demo')

% Set parameters.
im_name   = 'Pirate.tif';   % Name of clean image.
noise_std = 25;             % Standard deviation of image noise.
tau       = 0.85;           % Factor in residual bound.

% Load the clearn image and add Gaussian noise; make sure the pixels of the
% noisy image are in the range 0,...,255.
Xc = double(imread(im_name));
B = Xc + noise_std*randn(size(Xc));
B(B<0) = 0; B(B>255) = 255;

figure(1), clf, colormap(gray)
subplot(2,2,1)
  imagesc(Xc); axis image off
  title('Original clean image')
subplot(2,2,2)
  imagesc(B), axis image off
  title('Noisy image');
    
% Set the residual bound delta, and compute reconstruction.
disp(['Computing TV reconstruction for tau = ',num2str(tau),...
      ', this takes several seconds ...'])
delta = tau*sqrt(numel(B))*noise_std;
[X,info1] = TVdenoise(B,delta);

subplot(2,2,3)
  imagesc(X); axis image off
  title(['TV denoised image,  \tau = ',num2str(tau)])
      
% Also try a larger value of tau.
tau = 0.9;
disp(['Computing TV reconstruction for tau = ',num2str(tau),...
      ', this takes even longer ...'])
delta =  tau*sqrt(numel(B))*noise_std;
[XX,info2] = TVdenoise(B,delta);

subplot(2,2,4)
  imagesc(XX); axis image off
  title(['TV denoised image,  \tau = ',num2str(tau)])
% Median filtered and TV denoised client script based on MxTV toolbox.
%
% Reference: 
% J. Dahl^1, P.C. Hansen^2, S.H. Jensen^1 & T.L. Jensen^1
% CSI project: (1) Aalborg University, (2) Technical University of Denmark
% April 28, 2009.


% clear, clc
close all
disp('Starting TVdenoise_demo')

% Set parameters.
im_name   = 'text2.tif';   % Name of clean image.
noise_std = 25;             % Standard deviation of image noise.
tau       = 0.1;           % Factor in residual bound.

% Load the clearn image and add Gaussian noise; make sure the pixels of the
% noisy image are in the range 0,...,255.
% X0=imread(im_name);
% Xc = double(imread(im_name));

%%B = Xc + noise_std*randn(size(Xc));
% B = Xc;
% B(B<0) = 0; B(B>255) = 255;
Xc = proj2_n;
  
xc_max = max(Xc(:));
xc_min = min(Xc(:));
B = 255 * (Xc -xc_min)/(xc_max - xc_min );
figure(1), clf, colormap(gray)
subplot(2,2,1)
  imagesc(Xc); axis image off
  title('Original image')
% subplot(2,2,2)
%   imagesc(B), axis image off
%   title('Noisy image');
 
% Set the residual bound delta, and compute reconstruction.
disp(['Computing TV reconstruction for tau = ',num2str(tau),...
      ', this takes several seconds ...'])
delta = tau*sqrt(numel(B))*noise_std;
[X,info1] = TVdenoise(B,delta);

subplot(2,2,2)
  imagesc(X); axis image off
  title(['TV denoised image,  \tau = ',num2str(tau)])
      
% Also try a larger value of tau.
% tau = 1.2;
% disp(['Computing TV reconstruction for tau = ',num2str(tau),...
%       ', this takes even longer ...'])
% delta =  tau*sqrt(numel(B))*noise_std;
% [XX,info2] = TVdenoise(B,delta);

Xm=medfilt2(B,[3 3]);
subplot(2,2,3)
  imagesc(Xm); axis image off
  title('Median filter denoised image');


error=[B(:,250) Xm(:,250) X(:,250)];
figure
% error=[Xm(:,250) X(:,250)];
% subplot(2,2,4);
plot(error);
legend('original noisy image','Median denoise','Tv-denoise');
% legend('Median denoise','Tv-denoise');



proj_TV = X/255*(xc_max - xc_min )+xc_min;
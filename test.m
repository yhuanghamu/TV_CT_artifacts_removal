%% Physical vaule to grayscale transformation and invert-transformation
clear;
data_path = 'C:\Users\HYC\Documents\GitHub\TV_CT_artifacts_removal\data\120214';
file_proj = fullfile(data_path,'proj_800.mat');
load(file_proj);

[m,n] = size(proj_800);
proj = zeros(m,n);
% 0.45 <->3.65605786E-07
for i = 1:n
    proj(:,i) = proj_800(:,i) - mean(proj_800(:,i));
end

min_v = min(proj(:));
max_v = max(proj(:));

proc_gs = (proj - min_v)/(max_v - min_v)*255; % transform to gray scale value.

% Set parameters.
noise_std = 3;             % Standard deviation of image noise.
tau       = 0.85;           % Factor in residual bound.

delta = tau*sqrt(numel(proc_gs))*noise_std;
[proj_800_tv,info1] = TVdenoise(proc_gs,delta);

save(fullfile(data_path,'proj_800_tv.mat'),'proj_800_tv');

recon_800 = recon(proj_800_tv);
save(fullfile(data_path,'recon_800.mat'),'recon_800');



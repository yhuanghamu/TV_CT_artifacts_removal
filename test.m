%% Physical vaule to grayscale transformation and invert-transformation
clear;
proj_name = 'proj_800.mat';
data_path = '/Users/hyc/Documents/Github/TV_CT_artifacts_removal/data/150516';
% file_proj = fullfile(data_path,name);
% load(file_proj);
load(proj_name);
proj = proj_800;
[m,n] = size(proj);
% proj_n = zeros(m,n);
% 0.45 <->3.65605786E-07
% for i = 1:n
%     proj_n(:,i) = proj(:,i) - mean(proj(:,i));
% end

min_v = min(proj(:));
max_v = max(proj(:));
proc_gs = (proj - min_v)/(max_v - min_v)*255; % transform to gray scale value.

% Set parameters.
noise_std = 4;             % Standard deviation of image noise.
tau       = 0.85;           % Factor in residual bound.
delta = tau*sqrt(numel(proc_gs))*noise_std;
[proj_tv,info1] = TVdenoise(proc_gs,delta);

proj_tv = proj_tv /255 * (max_v - min_v) + min_v;
varname = strcat('tv_',proj_name);
save(fullfile(data_path,varname),'proj_tv');


recon_tv = recon(proj_tv);
varname = strcat('recon_tv_',proj_name);
save(fullfile(data_path,varname),'recon_tv');

recon_normal = recon(proj);
varname = strcat('recon_normal_',proj_name);
save(fullfile(data_path,varname),'recon_normal');



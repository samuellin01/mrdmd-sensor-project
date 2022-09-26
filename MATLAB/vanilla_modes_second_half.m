%% Script to create maps of POD modes for 2005-16 PM2.5 data

clear; close all; clc;
datpath = '../DATA/';
figpath = '../FIGURES/';

% Xcoord = ncread([datpath,'2007_June_PM25.nc'],'x');
% Ycoord = ncread([datpath,'2007_June_PM25.nc'],'y');
Band = ncread([datpath,'Downsampled_Annual_2000to2016_PM25_nocompression.nc'], 'Band1');
mask = ncread([datpath,'PM25_mask_array_coarse2_final.nc'],'Band1');

% Latter 4096 days of data, excluding corrupted dates
Band = Band(:,:,setdiff((2082:6180),[3291,5689,5690]));

[m,n,p] = size(Band); %p is the time length

N = m*n;
X = zeros(N,p);

M = length(mask(mask==1));
Y = zeros(M, p);

for i=1:ceil(p/2)
    snapshot = reshape(Band(:,:,i),N,1);
    Y(:,i) = snapshot(mask==1);
end

% train on entirety of data
Iord = 1:p;
Itrain = Iord(1:end);
%Itest = Iord(~ismember(Iord,Itrain));

Train = Y(:,Itrain);
timeavg = mean(Train,2);
Train = bsxfun(@minus,Train,timeavg);

[Psi,S,V] = svd(Train,'econ');
[m,n] = size(Train);
sing = diag(S);
thresh = optimal_SVHT_coef(n/m,0)*median(sing);
r_opt = length(sing(sing>=thresh))

R = [1 2 r_opt 10];

%% print incremental eigenmode to file

for r=r_opt
    
    display_fig(Psi(:,r),mask,[],[]);
    fig_name = strcat(figpath, 'FIG_SST_R=',string(r_opt),'.fig');
    savefig(fig_name);
end

%% plot singular values
close all;
plot(sing,'.','color',.7*[1 1 1]);
hold on
plot(sing(1:r_opt),'b.');

plot(R, sing(R),'ro');
plot(R, sing(R),'r.');
grid on
savefig([figpath,'FIG_SINGULAR_VALUES.fig']);

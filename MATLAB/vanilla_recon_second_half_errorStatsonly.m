%% Script to reconstruct PM2.5 data using POD/QR for 2000-2011

clear; close all; clc;
datpath = '../DATA/';
figpath = '../FIGURES/';

PRINT_FIG = false;

Band = ncread([datpath,'Downsampled_Annual_2000to2016_PM25_nocompression.nc'], 'Band1');
mask = ncread([datpath,'PM25_mask_array_coarse2_final.nc'],'Band1');

% Latter 4096 days of data, excluding corrupted dates
Band = Band(:,:,setdiff((2082:6180),[3291,5689,5690]));


[m,n,p] = size(Band); %p is the time length

%% Code

N = m*n;
X = zeros(N,p);

M = length(mask(mask==1));
Y = zeros(M,p);
for i=1:p
    snapshot = reshape(Band(:,:,i),N,1);
    Y(:,i) = snapshot(mask==1);
end

% train on entirety of the data
Iord = 1:p;

Itrain = Iord(1:end);
%Itest = Iord(~ismember(Iord,Itrain));

% Test dates for reconstruction
%Itest = [1065 1521 2373 3196]; % 12/1/02, 03/1/04, 07/1/06, 10/1/08
%indt = 1; % Index of Itest to run reconstruction (1 to 4 here)
rng(10)
Itest = randperm(2082, 400) %extrapolate
%Itest = randperm(4000, 400) %interpolate


Train = Y(:,Itrain);
timeavg = mean(Train,2);
Train = bsxfun(@minus,Train,timeavg);

[Psi,S,V] = svd(Train,'econ');
[m,n] = size(Train);
sing = diag(S);
thresh = optimal_SVHT_coef(n/m,0)*median(sing);
r_opt = length(sing(sing>=thresh))

% select validation snapshot
%x = Y(:,Itest(indt))-timeavg; 
% bounds = [min(x+timeavg) max(x+timeavg)];
%display_fig(x+timeavg,mask,[],[-1.8 35.6]); 
%if PRINT_FIG
%    file_name = strcat(figpath, 'FIG_SNAPSHOT.fig');
%    savefig(file_name);
%end
%printFormattedEPS(gca,gcf,[figpath,'FIG_enso_true.eps'],[]);

% seed random number generator for reproducibility
rng(729);

% target ranks
R = [10, 25, 50, 100, 250, 500, 750, 1000, r_opt, 1500, 2000];
%R = [10, 25, 50]
rmse_output = zeros(400);



%extrapolate error
Band = ncread([datpath,'Downsampled_Annual_2000to2016_PM25_nocompression.nc'], 'Band1');
mask = ncread([datpath,'PM25_mask_array_coarse2_final.nc'],'Band1');

% First 4096 days of data, excluding corrupted dates
Band = ncread([datpath,'Downsampled_Annual_2000to2016_PM25_nocompression.nc'], 'Band1');
Band = Band(:,:,setdiff((1:4097),[3291]));

Y = zeros(M,p);
for i=1:p
    snapshot = reshape(Band(:,:,i),N,1);
    Y(:,i) = snapshot(mask==1);
end


%% Reconstruction with random vs. qr sensors
for r = R(1:end-1)

    %% POD approximation with r eigenmodes
    %disp(r)
    %xproj = Psi(:,1:r)*(Psi(:,1:r)'*x);
    %[rmse, mpe] = rmse_mpe(xproj+timeavg,Band(:,:,Itest(indt)),mask);
    %close all; display_fig(xproj+timeavg,mask,[],[-1.8 35.6]);
    %if PRINT_FIG
    %    file_name = strcat(figpath, 'FIG_POD_PROJ_R_OPT=',num2str(r),'_RMSE=',num2str(rmse),'_MPE=',num2str(mpe),'.fig');
    %    savefig(file_name);
    %    png_name = strcat(figpath, 'FIG_POD_PROJ_R_OPT=',num2str(r),'_RMSE=',num2str(rmse),'_MPE=',num2str(mpe),'.png');
    %    saveas(gcf,png_name);
    %end
    %printFormattedEPS(gca,gcf,[figpath,'FIG_enso_proj',num2str(r),'.eps'],[]);
    
    %% Random reconstruction with r sensors
    % Algorithm takes entire 871x413 grid and traverses each row one pixel
    % at a time, placing a sensor every sparsity_value number of pixels.
    % Mask is then applied.
    
    %%sensors = randperm(m,r);
    %sparsity_value = 143;
    
    %sensors = [];
    %find_sensors = zeros(871*413,1);
    %for i=1:size(find_sensors,1)
    %    if mod(i,sparsity_value) == 0
    %        find_sensors(i) = 1;
    %    end
    %end
    %disp(strcat('num: ',num2str(ceil(0.97*m/r))));
    %find_sensors = find_sensors(mask==1);
    %for i=1:size(find_sensors,1)
    %    if find_sensors(i)==1
    %        sensors = [sensors i];
    %    end
    %end

    %close all; display_sensors(x+timeavg,mask,sensors);
    %if PRINT_FIG
    %    file_name = strcat(figpath, 'FIG_RAND_SENSORS_sensors=', num2str(size(sensors,2)),'.fig');
    %    savefig(file_name);
    %end
    %printFormattedEPS(gca,gcf,[figpath,'FIG_enso_randmask',num2str(r),'.eps'],[]);
    
    %xls = Psi(:,1:r)*(Psi(sensors,1:r)\x(sensors));
    %[rmse, mpe] = rmse_mpe(xls+timeavg, Band(:,:,Itest(indt)),mask);
    %close all; display_fig(xls+timeavg,mask,[],[-1.8 35.6]);
    %if PRINT_FIG
    %    file_name = strcat(figpath, 'FIG_RAND_PROJECTION_R_OPT=',num2str(r),'_RMSE=',num2str(rmse),'_MPE=',num2str(mpe),'.fig');
    %    savefig(file_name);
    %    png_name = strcat(figpath, 'FIG_RAND_PROJECTION_R_OPT=',num2str(r),'_RMSE=',num2str(rmse),'_MPE=',num2str(mpe),'.png');
    %    saveas(gcf,png_name);
    %end
    %%printFormattedEPS(gca,gcf,[figpath,'FIG_enso_rand_',num2str(r),'.eps'],[]);
    
    %% QDEIM with r QR sensors
    
    [~,~,pivot] = qr(Psi(:,1:r)','vector');
    sensors = pivot(1:r);
    
    %close all; display_sensors(x+timeavg,mask,sensors);
    %if PRINT_FIG
    %    file_name = strcat(figpath, 'FIG_QR_SENSORS_R_OPT=',num2str(r),'.fig');
    %    savefig(file_name);
    %end
    %printFormattedEPS(gca,gcf,[figpath,'FIG_enso_opt_mask',num2str(r),'.eps'],[]);


    fprintf("starting loop")
    disp(r)
    for i=1:400%size(Xrecon,2)
        x = Y(:,Itest(i))-timeavg;  
        xls = Psi(:,1:r)*(Psi(sensors,1:r)\x(sensors));
        [rmse, mpe] = rmse_mpe(xls+timeavg,Band(:,:,Itest(i)),mask);
        rmse_output(i) = rmse;
    end 
   
    fprintf("means")
    mean(real(rmse_output))
    median(real(rmse_output))
    rmse_output = zeros(400);

    
    %close all; display_fig(xls+timeavg,mask,[],[-1.8 35.6]);
    %if PRINT_FIG
    %    file_name = strcat(figpath, 'FIG_QR_RECON_R_OPT=',num2str(r),'_RMSE=',num2str(rmse),'_MPE=',num2str(mpe),'.fig');
    %    savefig(file_name);
    %    png_name = strcat(figpath, 'FIG_QR_RECON_R_OPT=',num2str(r),'_RMSE=',num2str(rmse),'_MPE=',num2str(mpe),'.png');
    %    saveas(gcf,png_name);
    %end
    %printFormattedEPS(gca,gcf,[figpath,'FIG_enso_opt_',num2str(r),'.eps'],[]);
end



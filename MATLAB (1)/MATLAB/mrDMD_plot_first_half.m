%% Script to display mrDMD mode amplitude maps of PM2.5 data, 2000-2011

clear; close all; clc

datpath = '../DATA/';
figpath = '../FIGURES/';

% Set PRINT_FIG=true to export figures
PRINT_FIG = true;
set(0,'defaultaxescolororder',parula(5));

% Import data
dat = ncread([datpath,'Downsampled_Annual_2000to2016_PM25_nocompression.nc'], 'Band1'); 
mask = ncread([datpath,'PM25_mask_array_coarse2_final.nc'],'Band1');

% First 4096 days of the data, not counting corrupted date
%dat = dat(:,:,setdiff((1:6180),[3291,5689,5690]));
dat = dat(:,:,setdiff((1:4097),[3291]));

numyears = 11;
nweeks = numyears*52;

Y = zeros(length(mask(mask==1)),size(dat,3));
for i=1:size(dat,3)
    Band = dat(:,:,i);
    Y(:,i) = Band(mask==1);
end

[N,M] = size(Y);
dt = 1; %time interval of 1

% First 4096 days of data
time = setdiff((1:4097),[3291])';

% Edit r, max_cyc, and levels here
% max_cyc = 1, L=13
tree = mrDMD_fb(Y(:,1:length(time)),dt,10,1,13,true);

%% plot amplitudes of the mrDMD
figure;
[ptree, map, low_f_cutoff] = mrDMD_map(tree);
[L,J] = size(tree);

imagesc(-sqrt(map));
shortInt = tree{L,J}.T;
T = datetime(2000,1,1,0,0,0) + days(time(1:shortInt:end));
set(gca, 'YTick', 0.5:(L+0.5), 'YTickLabel', low_f_cutoff*365);
set(gca,'XTick', (0:J-1) + 0.5);
axis xy; axis tight; box on; grid on
set(gca,'XTickLabel', strcat(num2str(month(T)),'/',num2str(year(T))));
set(gca, 'XTickLabelRotation',45);

colormap(gca,'pink');shading flat

if PRINT_FIG
    file_name = strcat(figpath, 'FIG_MRDMD_MAP.fig');
    savefig(file_name);
end

%% collect and display unique mrDMD modes

lmap = []; jmap = [];
Lib = []; Lib2 = [];
Omega = []; Omega2 = [];
Amp = []; Amp2 = [];
Periods = [];
Lambda = []; Lambda2 = [];

tol = 1e-2; % 1e-2

for l=1:L
    for j=1:2^(l-1)
        Lib = [Lib tree{l,j}.Phi(1:N,:)];        
        Omega = [Omega; tree{l,j}.omega*365]; % yearly (!)
        Amp = [Amp; tree{l,j}.P];
        Periods = [Periods; tree{l,j}.T];
        Lambda = [Lambda; tree{l,j}.lambda];
        
        ind = find(imag(tree{l,j}.omega*365)>tol);
%       Code from mrDMD paper to produce maps of the modes; however, the 
%       tol value may not be calibrated correctly, so a lot of modes are
%       printed and the program takes too long to run. Use
%       mrDMD_specific_modes_first_half.m
% 
%         for kk=1:length(ind)
%             figure;
%             display_fig(abs(tree{l,j}.Phi(1:N,ind(kk))),mask,[],[]);
%             str = ['Phi97' num2str(l) ',' num2str(j) '_' ...
%                 num2str(round(imag(tree{l,j}.omega(ind(kk))*365),2))];
%             if PRINT_FIG
%                 file_name = strcat(figpath, 'FIG_MRDMD_MODE_L=',string(l),'_J=',string(j),'.fig');
%                 savefig(file_name);
%             end
%         end
            
%         Lib2 = [Lib2 tree2{l,j}.Phi(1:N,:)];
%         Omega2 = [Omega2; tree2{l,j}.omega*52]; % yearly
%         Amp2 = [Amp2; tree2{l,j}.P];
%         Lambda2 = [Lambda2; tree2{l,j}.lambda];
%         
%         ind = find(imag(tree2{l,j}.omega*52)>tol);
%         for kk=1:length(ind)
%             figure;
%             display_fig(abs(tree2{l,j}.Phi(1:N,ind(kk))),mask,[],[]);
%             str = ['Phi16' num2str(l) ',' num2str(j) '_' ...
%                 num2str(round(imag(tree2{l,j}.omega(ind(kk))*52),2))];
%             if PRINT_FIG
%                 export_fig([figpath str],'-pdf');
%             end
%         end
    end
end

%% Filter library to keep only oscillatory and background mode
ind = find(abs(imag(Omega))>tol);
%ind = [1 ;ind];

figure;
hold on
idx = find(abs(imag(Omega))>.1);
scatter(real(Omega),imag(Omega),150,0.5*[1 1 1])
scatter(real(Omega(idx)),imag(Omega(idx)),150,'r','filled');

grid on; box on

if PRINT_FIG
    file_name = strcat(figpath, 'FIG_MRDMD_FREQ.fig');
    savefig(file_name);
end

%% display QDEIM sensor locations and corresponding time series

r = 4; % Change number of sensors here (for timeseries)
cmap = lines(r);
[Phi,S,V] = svd(Y(:,1:length(time)),'econ');
Phi = Phi(:,1:r);
[~,~,piv] = qr(Phi','vector');
piv = piv(1:r);

figure;
display_mrdmd_sensors(mask,piv,cmap);box on
if PRINT_FIG
    file_name = strcat(figpath, 'FIG_POD_SENSORS_R=',string(r),'.fig');
    savefig(file_name);
end

figure;
plot(Y(piv,1:length(time))','LineWidth',2); axis tight
shortInt = tree{L,J}.T;
T = datetime(2000,1,1,0,0,0) + days(time(1:shortInt:end));
set(gca,'XTick', (0:J-1)*shortInt + 0.5);
axis tight; box on;
set(gca,'XTickLabel', strcat(num2str(month(T)),'/',num2str(year(T))));
set(gca, 'XTickLabelRotation',45);
set(gca,'XTick',[]);set(gcf,'Position',[380   320   540   200]);

if PRINT_FIG
    file_name = strcat(figpath, 'FIG_POD_TIMESERIES_R=',string(r),'.fig');
    savefig(file_name);
end

%% display multiscale sensor locations and corresponding time series
cmap = lines(r);

Phi = Lib(1:N,idx);
[~,~,piv] = qr(Phi','vector');
piv = piv(1:r);

figure;
display_mrdmd_sensors(mask,piv,cmap); box on
if PRINT_FIG
    file_name = strcat(figpath, 'FIG_MRDMD_SENSORS_R=',string(r),'.fig');
    savefig(file_name);
end

figure;
plot(Y(piv,1:length(time))','LineWidth',2); axis tight; box on

set(gca,'XTick', (0:J-1)*shortInt + 0.5);
axis tight; box on;
set(gca,'XTickLabel', strcat(num2str(month(T)),'/',num2str(year(T))));
set(gca, 'XTickLabelRotation',45);

set(gca,'XTick',[]);set(gcf,'Position',[380   320   540   200]);

if PRINT_FIG
    file_name = strcat(figpath, 'FIG_MRDMD_TIMESERIES_R=',string(r),'.fig');
    savefig(file_name);
end

%% display POD energy spectrum 

sig = diag(S)/sum(diag(S));

figure;
hold on
semilogy(sig,'.','Color',.5*[1 1 1]); 
scatter(1:r,sig(1:r),150,'r','filled');

set(gca,'yscale','log'); axis tight
grid on, box on

if PRINT_FIG
    file_name = strcat(figpath, 'FIG_POD_SPECTRUM.fig');
    savefig(file_name);
end



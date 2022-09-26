%% Script to display particular modes of an mrDMD map produced from mrDMD_plot script

clear; close all; clc

datpath = '../DATA/';
figpath = '../FIGURES/';

% Set PRINT_FIG=true to export figures
PRINT_FIG = true;
set(0,'defaultaxescolororder',parula(5));

% Import data
dat = ncread([datpath,'Downsampled_Annual_2000to2016_PM25_nocompression.nc'], 'Band1'); 
mask = ncread([datpath,'PM25_mask_array_coarse2_final.nc'],'Band1');

%dat = dat(:,:,setdiff((1:4097),[3291]));
dat = dat(:,:,setdiff((2082:6180),[3291,5689,5690]));

% These two variables are now obsolete
numyears = 11;
nweeks = numyears*52;

Y = zeros(length(mask(mask==1)),size(dat,3));
for i=1:size(dat,3)
    Band = dat(:,:,i);
    Y(:,i) = Band(mask==1);
end

[N,M] = size(Y);
dt = 1; %time interval of 1

%time = setdiff((1:4097),[3291])';
time = setdiff((2082:6180),[3291,5689,5690])'; 

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

%% Plot background mode (l=1,j=1) and save image
for l=1
    for j=1
        Lib = [Lib tree{l,j}.Phi(1:N,:)];        
        Omega = [Omega; tree{l,j}.omega*365]; % yearly (!)
        Amp = [Amp; tree{l,j}.P];
        Periods = [Periods; tree{l,j}.T];
        Lambda = [Lambda; tree{l,j}.lambda];

        figure;
        display_fig(abs(tree{l,j}.Phi(1:N,1)),mask,[],[]);
        colorbar;
        if PRINT_FIG
            file_name = strcat(figpath, 'FIG_MRDMD_MODE_L=',string(l),'_J=',string(j),'.fig');
            savefig(file_name);
            png_name = strcat(figpath, 'FIG_MRDMD_MODE_L=',string(l),'_J=',string(j),'.png');
            saveas(gcf,png_name);
        end
    end
end

%% Plot other mrDMD modes

% Custom input l and j values, such that the ith mode is on level number
% l_vals(i) and is the j_vals(i)-th time window in that level
l_vals = [5,6,7,9,9,9,10];
j_vals = [1,10,40,37,42,171,451];

for i=1:size(l_vals,2)
    l = l_vals(i);
    j = j_vals(i);
    Lib = [Lib tree{l,j}.Phi(1:N,:)];        
    Omega = [Omega; tree{l,j}.omega*365]; % yearly (!)
    Amp = [Amp; tree{l,j}.P];
    Periods = [Periods; tree{l,j}.T];
    Lambda = [Lambda; tree{l,j}.lambda];
    
    ind = find(imag(tree{l,j}.omega*365)>tol);

    for kk=1:length(ind)
        figure;
        display_fig(abs(tree{l,j}.Phi(1:N,ind(kk))),mask,[],[]);
        colorbar;
        str = ['Phi97' num2str(l) ',' num2str(j) '_' ...
            num2str(round(imag(tree{l,j}.omega(ind(kk))*365),2))];
        if PRINT_FIG
            file_name = strcat(figpath, 'FIG_MRDMD_MODE_L=',string(l),'_J=',string(j),'.fig');
            savefig(file_name);
            png_name = strcat(figpath, 'FIG_MRDMD_MODE_L=',string(l),'_J=',string(j),'.png');
            saveas(gcf,png_name);
        end
    end
end

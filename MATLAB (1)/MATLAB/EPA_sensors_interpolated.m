%% Script to interpolate EPA sensor locations (1706 on coarser 871x413) to full scale map

clear; close all; clc
datpath = '../DATA/';
figpath = '../FIGURES/FIRES_RECONS/EPA/2007/EPAsensors_2007SUMMER_';

% Set PRINT_FIG=true to export figures
PRINT_FIG = true;
set(0, 'defaultfigurecolor', 'w');

% Read in the data
dat = ncread([datpath,'Downsampled_Annual_2000to2016_PM25_nocompression.nc'], 'Band1'); 
mask = ncread([datpath,'PM25_mask_array_coarse2_final.nc'],'Band1');
epa = ncread([datpath, 'EPA_locations_coarse.nc'],'EPA_BAND');

% First 4096 days of data, excluding corrupted dates
dat = dat(:,:,setdiff((1:4097),[3291]));
%dat = dat(:,:,setdiff((2082:6180),[3291,5689,5690]));

%indt = linspace(912,943+31,63) %July-August 2002 fire
%indt = linspace(2739,2739+62,63) %July-August 2007 fire
%indt = linspace(2483,2483+122,123) %July-Oct 2012 fire
%indt = linspace(720,720+45,46) %January 2002 winter
%indt = linspace(2557-45,2557,46) %January 2007 winter
%indt = linspace(2301-40,2301+40,81) %Jan 2012 winter
%indt = linspace(2483-45,2483,46) %July-Oct 2012 summer
indt = linspace(2739-62,2739,63) %July-August 2007 summer

time = setdiff((1:4097),[3291])';
%time = setdiff((2082:6180),[3291,5689,5690]);
times = datetime(2000,1,1,0,0,0) + days(time(indt));


for i=1:size(indt,2)
    t = times(i);

    %get data only from EPA locations
    epa_locations = dat(:,:,indt(i)).*epa;
    %Set all other grids to NANs
    epa_locations(epa_locations==0) = NaN; 
    %Interpolate using only EPA data
    epa_interpolate = inpaint_nans(epa_locations); %add this package to the main file!

    %Plot EPA interpolated data with mask
    epa_interpolate_mask = zeros(length(mask(mask==1)),1);
    epa_interpolate_mask = epa_interpolate(mask==1);
    
    [rmse, mpe] = rmse_mpe(epa_interpolate_mask,dat(:,:,indt(i)),mask);
    
    figure;
    display_fig(epa_interpolate_mask,mask,[],[-1.8 35.6]);
    if PRINT_FIG
        file_name = strcat(figpath, 'FIG_EPA_RECON_',num2str([month(t),day(t),year(t)],'%d_%d_%d'),'_RMSE=',num2str(rmse),'_MPE=',num2str(mpe),'.fig');
        savefig(file_name);
        png_name = strcat(figpath, 'FIG_EPA_RECON_=',num2str([month(t),day(t),year(t)],'%d_%d_%d'),'_RMSE=',num2str(rmse),'_MPE=',num2str(mpe),'.png');
        saveas(gcf,png_name);
    end
    
    [AE, FE] = ae_fe_map(epa_interpolate_mask,dat(:,:,indt(i)),mask);
    
    figure;
    display_fig_hot(AE,mask,[],[0 10]);
    if PRINT_FIG
        file_name = strcat(figpath, 'FIG_EPA_AE_',num2str([month(t),day(t),year(t)],'%d_%d_%d'),'_RMSE=',num2str(rmse),'_MPE=',num2str(mpe),'.fig');
        savefig(file_name);
        png_name = strcat(figpath, 'FIG_EPA_AE_=',num2str([month(t),day(t),year(t)],'%d_%d_%d'),'_RMSE=',num2str(rmse),'_MPE=',num2str(mpe),'.png');
        saveas(gcf,png_name);
    end
    
    figure;
    display_fig_bwr(FE,mask,[],[-100 100]);
    if PRINT_FIG
        file_name = strcat(figpath, 'FIG_EPA_FE_',num2str([month(t),day(t),year(t)],'%d_%d_%d'),'_RMSE=',num2str(rmse),'_MPE=',num2str(mpe),'.fig');
        savefig(file_name);
        png_name = strcat(figpath, 'FIG_EPA_FE_=',num2str([month(t),day(t),year(t)],'%d_%d_%d'),'_RMSE=',num2str(rmse),'_MPE=',num2str(mpe),'.png');
        saveas(gcf,png_name);
    end
    
end




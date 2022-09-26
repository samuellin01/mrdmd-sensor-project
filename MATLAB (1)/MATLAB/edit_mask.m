% Code to filter out locations and overwrite an existing mask
% based on dataset input

clear; close all; clc;
datpath = '../DATA/';
figpath = '../FIGURES/';

% Dataset
Band = ncread([datpath,'Downsampled_Annual_2000to2016_PM25_nocompression.nc'], 'Band1');
% Existing mask
mask = ncread([datpath,'PM25_mask_array_coarse2_copy.nc'],'Band1');

% % Code for retrieving data at a particular time snapshot
% disp(size(Band));
% 
% figure;
% imagesc(Band(:,:,3291));
% colormap(jet);
% set(gca,'view',[-90,90]);
% 
% file_name = strcat(figpath, 'FIG_BAND_3291.fig');
% savefig(file_name);

% 2000-16 data has corrupted data that needs to be removed
% Band = Band(:,:,setdiff((1:6180),[3291,5689,5690]));

[m,n,p] = size(Band); %p is the time length

%% Fix conflicting values in mask
% Store all points of discrepancy between mask and Band in (x,y,t) triplets
discrepancies = [];
for i=1:p
    if nnz(isnan(Band(:,:,i))) ~= nnz(mask==0)
        [band_row,band_col] = find(isnan(Band(:,:,i)));
        [mask_row,mask_col] = find(mask==0);
        [discrepancies,loc] = find_discrepancy(band_row(1:size(mask_row,1)),band_col(1:size(mask_row,1)),mask_row,mask_col,discrepancies,i);
        while loc~=0
            band_row(loc) = [];
            band_col(loc) = [];
            [discrepancies,loc] = find_discrepancy(band_row(1:size(mask_row,1)),band_col(1:size(mask_row,1)),mask_row,mask_col,discrepancies,i);
        end
    end
end

% Sort and remove duplicates from discrepancies
discrepancies = sortrows(unique(discrepancies(:,1:3),'rows'));
% Fix the discrepancies in the mask
for i=1:size(discrepancies,1)
     mask(discrepancies(i,1),discrepancies(i,2)) = 0;
end

% Overwrite mask with new mask array
ncwrite('../DATA/PM25_mask_array_coarse2_copy.nc','Band1',mask)

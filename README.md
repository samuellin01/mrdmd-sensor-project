# Modal Decomposition Techniques on PM2.5 data

This repository contains Matlab code that applies matrix QR pivoting and multiresolution dynamic mode decomposition (mrDMD) to optimize sensor placement for PM
2.5 data. 

Due to GitHub file size limitations, the following analyzed dataset is linked [online](https://zenodo.org/record/5809461#.YzEJFnbMLtU).

## Repository organization

The repository is organized into three directories: /DATA, /FIGURES, and /MATLAB. All runnable scripts are located in the /MATLAB directory. The /DATA directory stores the PM2.5 data files and the mask. The /FIGURES directory is used to store output from the MATLAB scripts. At the top of each executable MATLAB script is a line that sets figpath to the name of the output folder. All figpaths are currently set to FIGURES but should be changed if additional folders are created to store output. 

The folder sensor_project_results contains results from past work and does not affect the code in any way.

Runnable MATLAB scripts come in pairs, one for the first 4096 days (2000-2011) and one for the latter 4096 days (2005-2016). At present, the mask will not work unless data at three dates are excluded (already programmed). The 2000-16 data file contains data for 6180 days, and the excluded dates are at indices 3291, 5689, and 5690.

## POD/QR 

**vanilla_modes_first_half.m** and **vanilla_modes_second_half.m** perform SVD and produce two figures: the eigenmode at the optimal rank, and a graph of the singular values.

**vanilla_recon_first_half.m** and **vanilla_recon_second_half.m** produce reconstructions at selected dates. In particular, they produce 1) a POD approximation (non-sensor) with the optimal number of eigenmodes, 2) a grid-based random reconstruction with that number of sensors, and 3) a QR reconstruction with the same number of sensors. To configure dates for reconstruction, change the values in the array Itest (line 36). See additional notes on configuring dates.

Unlike the mrDMD recon scripts, these vanilla recon scripts can only produce reconstructions for one date at a time. Many indices can be listed on line 36. The single index to run reconstruction is configured on line 37 and should be changed accordingly.

RMSE and MPE values are stored in the names of outputted files.

## mrDMD

**mrDMD_plot_first_half.m** and **mrDMD_plot_second_half.m** perform mrDMD on their respective 4096-day periods. They produce an mrDMD tree, a plot of the mrDMD frequencies, QR and mrDMD timeseries and map for most important r sensors, and a figure of the POD singular values. 
1. Parameters for the mrDMD algorithm, particularly max_cyc and levels, can be changed on line 37. 
2. The original version based on the Manohar paper prints out mrDMD modes using this script. This script is unable to do that. See documentation on the following scripts.
3. The timeseries for the most important r sensors is created in two sections that start on line 131. The value of r can be changed on line 133.

**mrDMD_specific_modes_first_half.m** and **mrDMD_specific_modes_second_half.m** produce maps of specific modes from an mrDMD tree (see above). Modes are defined by the variables l and j, where l is the level number (counting from bottom to top) and j denotes the location of the mode at that level. To elaborate, if the lth level is divided into n subsets, then the jth subset will be printed. At present, the l and j values must be determined manually, unfortunately.
1. Parameters for the mrDMD algorithm, particularly max_cyc and levels, can be changed on line 37.
2. The background mode (l=1,j=1) will always be printed.
3. Custom l and j values can be configued on lines 96-97 in the arrays l_vals and j_vals. As an example, to print the 10th mode on the 6th level, then l_vals should contain the value 6, and the corresponding j_vals value is 10.

**mrDMD_recon_first_half.m** and **mrDMD_recon_second_half.m** produce mrDMD reconstructions at selected dates. Specifically, they will print out an mrDMD tree, reconstructions at selected dates, maps of the actual data at those same dates, and a map of the locations of mrDMD sensors.

1. Parameters for the mrDMD algorithm, particularly max_cyc and levels, can be changed on line 36. 
2. Dates selected for reconstruction can be configured with the indt array on line 44. The scripts will create reconstructions for all of these dates in one execution (unlike vanilla recon). See additional notes on configuring dates.
3. RMSE and MPE values are stored in the names of outputted files. Imaginary values are to be expected here and can be ignored.

## Configuring dates for reconstruction

Both vanilla and mrDMD recon produce reconstructions at dates specified by the user. Dates are specified as indices of the Band array. To determine the index for a particular date ENDDATE, calculate the number of days between 1/1/2000 and ENDDATE, not including end date. (timeanddate.com/date/duration.html) Some adjustments may be necessary because of the corrupted data at some dates.

If running on first 4096 days of data:
That number is the correct index.

If running on latter 4096 days of data:
Subtract 2082 from that number.


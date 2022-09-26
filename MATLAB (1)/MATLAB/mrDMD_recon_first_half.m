%% Script to reconstruct snapshots of PM2.5 data (2000-11) using mrDMD algorithm

clear; close all; clc
datpath = '../DATA/';
figpath = '../FIGURES/250sensors';

% Set PRINT_FIG=true to export figures
PRINT_FIG = true;
set(0, 'defaultfigurecolor', 'w');

% Read in the data
dat = ncread([datpath,'Downsampled_Annual_2000to2016_PM25_nocompression.nc'], 'Band1'); 
mask = ncread([datpath,'PM25_mask_array_coarse2_final.nc'],'Band1');

% First 4096 days of data, excluding corrupted dates
dat = dat(:,:,setdiff((1:4097),[3291]));

% Obsolete variables
numyears = 11;
nweeks = numyears*52;

Y = zeros(length(mask(mask==1)),size(dat,3));
for i=1:size(dat,3)
    Band = dat(:,:,i);
    Y(:,i) = Band(mask==1);
end

[N,M] = size(Y);
dt = 1;

% mrDMD, 4096-day period beginning in 2000
time = setdiff((1:4097),[3291])';

% Edit r, max_cyc, levels here
% max_cyc = 1 and L=10
tree = mrDMD_fb(Y(:,1:length(time)),dt,10,1,10,true);

[U,S,V] = svd(Y(:,1:length(time)),'econ');
[~,~,piv] = qr(U(:,1:30)',0);
qdeim = piv(1:30);

%% plot amplitudes of mrDMD

indt = [1065 1521 2373 3196]; % 12/1/02, 03/1/04, 07/1/06, 10/1/08

[ptree, map, low_f_cutoff] = mrDMD_map(tree);
[L, J] = size(tree);

imagesc(-sqrt(map));

shortInt = tree{L,J}.T;

T = datetime(2000,1,1,0,0,0) + days(time(1:shortInt:end));
set(gca, 'YTick', 0.5:(L+0.5), 'YTickLabel', low_f_cutoff*365);
set(gca,'XTick', (0:J-1) + 0.5,'XTickLabel',num2str([month(T),year(T)],'%d/%d'));
axis xy; axis tight; box on; grid on

ylim = get(gca,'ylim');

times = year(datetime(2000,1,1,0,0,0) + days(time(indt)));
hold on
for i=1:length(indt)
    stem(1+indt(i)/shortInt, ylim(2),'k.','LineWidth',1.5)
    
end

set(gca, 'XTickLabelRotation',45,'ylim',ylim,'FontSize',12);
colormap(gca,'pink');shading flat

if PRINT_FIG
    file_name = strcat(figpath, 'FIG_MRDMD_RECON_TIMES.fig');
    savefig(file_name);
end

%% collect unique mrDMD modes

lmap = []; jmap = [];
Lib = [];
Omega = [];
Amp = [];
Periods = [];

for l=1:L
    for j=1:2^(l-1)
        Lib = [Lib tree{l,j}.Phi(1:N,:)];
        Omega = [Omega; tree{l,j}.omega*365]; % yearly
        Amp = [Amp; tree{l,j}.P];
        Periods = [Periods; tree{l,j}.T];
    end
end

[~,~,piv] = qr(Lib.','vector');
%sens = piv(1:size(Lib,2));
sens = piv(1:250);

Phi = cell(1,L);
omega = cell(L,1);
b = cell(L,1);
t0 = cell(L,1);

PhiJ = cell(1,J); win =1;
lab = [];
xL = zeros(N,L); %prediction at each level
update = false;
Lev_coeff = [];

Xrecon = [];

%% compute reconstruction of PM25 data from multiscale sensor measurements
for i=1:length(time)
    x = Y(:,i);
    
    %xavg = L(:,i); % filtered data for sensors
    
    if (i>2 && i < length(time)-1)
        xavg = mean(Y(:,i-2:i+2),2);
    else 
        xavg = x;
    end


    % at each level update t-t0
    for l=1:L
        Int = tree{l,1}.T;
        
        if (mod(i,Int)==0 && i< length(time))
            t0{l} = time(i+1);
            ind = floor((i+1)/Int)+1;
            update = true;
        elseif (i==1)
            ind = 1;
            t0{l} = time(1);
            update = true;
        end
        
        if (update) 
            %disp([i ind]) % update modes when time windows change

            Phi{l} = tree{l,ind}.Phi(1:N,:);
            omega{l} = tree{l,ind}.omega;
            b{l} = tree{l,ind}.P;
            %sub{l} = ceil(1/8/pi/tree{l,ind}.rho/dt);
            

            Modes = cell2mat(Phi);
            
            ti = (time(i)-t0{l}); % time in weeks
            xL(:,l) = Phi{l}*(exp(2*pi*omega{l}*ti).*b{l});
            if (i==1)
                x0 = x;
            else
                x0 = Y(:,i+1);
            end
            lev_coeff = xL\x0;           
            
            if (l==L) % store J sets of modes by last window
                PhiJ{win} = cell2mat(Phi);
                lab = [lab; win*ones(size(PhiJ{win},2),1)];
                Lev_coeff = [Lev_coeff lev_coeff];
                win = win+1;
            end
            
            update = false;
            r = size(Modes,2);              
            disp([size(Modes,2) length(sens)]);
        end                
            
        ti = (time(i)-t0{l}); % time in weeks
        xL(:,l) = Phi{l}*(exp(2*pi*omega{l}*ti).*b{l});
    end
    
    
    % weighted mrDMD reconstruction
    % least-squares fit of by-level approximation to true snapshot
    xpred = xL*lev_coeff;    
    
    % collect Phi at all levels and approximate
    xhat = Modes*(Modes(sens,:)\x(sens));
    xhat2 = Lib*(Lib(sens,:)\x(sens));  
    
    xpod = U(:,1:r)*(U(qdeim,1:r)\x(qdeim));

    if (ismember(i, indt))
        Xrecon = [Xrecon xhat];
    end
    errpred(i) = norm(x-xpred)/norm(x);
    errsens(i) = norm(x-xhat)/norm(x);
    errsens2(i) = norm(x-xhat2)/norm(x);
    errpod(i) = norm(x-xpod)/norm(x);
    
end

%% display reconstructions and actual maps of PM25 data from multiscale sensors

time = setdiff((1:4097),[3291]);
mask = ncread([datpath,'PM25_mask_array_coarse2_final.nc'],'Band1');
times = datetime(2000,1,1,0,0,0) + days(time(indt));
disp(size(Xrecon));

for i=1:size(Xrecon,2)
    figure; 
    t = times(i);
    display_fig(Xrecon(:,i),mask,[],[-1.8 35.6]);
    % RMSE and MPE calculations
    [rmse, mpe] = rmse_mpe(Xrecon(:,i),dat(:,:,indt(i)),mask);
    title(num2str([month(t),day(t),year(t)],'%d/%d/%d'));
    if PRINT_FIG
        % Save reconstruction image
        file_name = strcat(figpath, 'FIG_RECON_',num2str([month(t),day(t),year(t)],'%d_%d_%d'),'_RMSE=',num2str(rmse),'_MPE=',num2str(mpe),'.fig');
        savefig(file_name);
        png_name = strcat(figpath, 'FIG_RECON_',num2str([month(t),day(t),year(t)],'%d_%d_%d'),'_RMSE=',num2str(rmse),'_MPE=',num2str(mpe),'.png');
        saveas(gcf,png_name);
        % Save actual image
        figure;
        C = dat(:,:,indt(i));
        b = imagesc(C, [-1.8 35.6]);
        set(b,'AlphaData',(~isnan(C)));
        set(gca,'view',[-90,90]);
        colormap(jet);
        axis off;
        title(num2str([month(t),day(t),year(t)],'%d/%d/%d'));
        file_name = strcat(figpath, 'FIG_RECON_ACTUAL_indt=',num2str(indt(i)),'_',num2str([month(t),day(t),year(t)],'%d_%d_%d'));
        savefig(file_name);
    end
end

%% display sensors
figure;

% 'r' parameter makes sensor points red, currently overwritten in
% display_mrdmd_sensors file
display_mrdmd_sensors(mask,sens,'r');

if PRINT_FIG
    file_name = strcat(figpath, 'FIG_RECON_SENSORS_',num2str(length(sens)),'.fig');
    savefig(file_name);
end

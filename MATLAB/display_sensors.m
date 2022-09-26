function [ ] = display_sensors( x,mask, pivot )
%   display_sensors Displays PM2.5 QR sensors on map
%   
%   INPUTS
%   x: data matrix of actual PM2.5 (obsolete for current implementation)
%   mask: binary mask for displaying continental US
%   pivot: sensor location indices
%

    
    sensors = zeros(871,413);
    P = zeros(size(x)); P(pivot)=1:length(pivot);    
    sensors(mask==1) = P;
       
    %mask(mask==1)=x;
    b = imagesc(mask');
    set(gca,'ydir','reverse');
    axis xy;
    shading flat, colormap(gray), drawnow
    alpha 0.3    
    hold on
    
    %scale x for colormap
    x = x+abs(min(x));
    x = x/max(x);
    
    % Color gradient that ranks the sensors by importance 
    first_half = size(pivot,2)/2;
    second_half = size(pivot,2) - first_half;
    % Amount of red decreases differently for 1st and 2nd half
    r1 = linspace(1,0.983,first_half)';
    r2 = linspace(0.982,0.431,second_half)';
    r = vertcat(r1,r2);
    % Amount of green decreases exponentially
    g = linspace(1,0.544,size(pivot,2))';
    g = exp(g) - 1.7183;
    % Amount of blue
    b = linspace(1,0.09803,size(pivot,2))';
    cmap = horzcat(r,g,b);
    % Flip so that white is least important, dark red most important
    cmap = flip(cmap); 

    S = reshape(real(sensors)',871*413,1);
    [C,IC,~] = unique(S);
    
    % align Ilin with pivot somehow
    
    [I,J] = ind2sub(size(sensors'),IC(2:end));
    
    scatter(J,I,25,cmap,'filled','MarkerEdgeColor','k');
        
    
    axis off
end


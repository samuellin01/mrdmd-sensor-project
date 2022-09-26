function display_mrdmd_sensors( mask, indS, cmap )
%   display_mrdmd_sensors Displays PM2.5 mrDMD sensors on map
%   
%   INPUTS
%   mask: binary mask for displaying continental US
%   indS: sensor location indices
%   cmap: colormap for sensors
%

    sensors = zeros(871,413);
    
    % number sensors by sensor order (if desired, color by order can be
    % implemented)
    % This is necessary for unique() to work properly
    P = zeros(size(mask(mask==1))); 
    P(indS)=1:length(indS);    
    
    sensors(mask==1) = P;
    b = imagesc(mask');
    set(gca,'ydir','reverse');
    axis xy;
    shading flat, colormap(gray), drawnow
    alpha 0.3
    hold on

    % Color gradient that ranks the sensors by importance 
    first_half = ceil(size(indS,2)/2);
    second_half = size(indS,2) - first_half;
    % Amount of red decreases differently for 1st and 2nd half
    r1 = linspace(1,0.983,first_half)';
    r2 = linspace(0.982,0.431,second_half)';
    r = vertcat(r1,r2);
    disp(size(r));
    % Amount of green decreases exponentially
    g = linspace(1,0.544,size(indS,2))';
    g = exp(g) - 1.7183;
    disp(size(g));
    % Amount of blue
    b = linspace(1,0.09803,size(indS,2))';
    disp(size(b));
    cmap = horzcat(r,g,b);
    % Flip so that white is least important, dark red most important
    cmap = flip(cmap);
    
    S = reshape(sensors',871*413,1);
    [~,IC,~] = unique(S);
    
    % IC(2:end) contains linear indices of sensor locations
    [I,J] = ind2sub(size(sensors'),IC(2:end));
    
    % set sensor markersize to 15
    scatter(J,I,30,cmap,'filled','MarkerEdgeColor','k');
    %alpha 0.1
    
    axis off; hold off
end


function [ output_args ] = display_fig( x, mask, pivot, cbounds )
%   display_fig Displays PM2.5 data on white map of continental US 
%   with gray areas excluded in mask

    snapshot = NaN*zeros(871*413,1);
    snapshot(mask==1) = x;
    
    sensors = zeros(871,413);
    P = zeros(size(x)); P(pivot)=ones(size(P(pivot)));
    sensors(mask==1) = P;
    
    % cbounds determines the range of the colorbar
    C = reshape(real(snapshot),871,413)';
    if (~isempty(cbounds))
        b = imagesc(C,cbounds);
    else 
        b = imagesc(C);%,[-1.799 30.77999]);
    end
    set(gca,'ydir','reverse');
    axis xy;
    shading interp, colormap(jet), drawnow;
    set(b,'AlphaData',(~isnan(C)));
    if (~isempty(sensors))
        hold on
        S = reshape(sensors,871,413)';
        [I,J] = find(S>0);
        
        scatter(J,I,'b.');
        
    end
    axis off
end


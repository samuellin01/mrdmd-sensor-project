function [AE,FE] = ae_fe_map(recon,actual,mask)
% Calculates Root Mean Squared Error and 
% Mean Percentage Error of reconstruction with actual data
% recon is 213804 x 1, actual is 871 x 413 with NaN values

    % Apply mask on actual and reshape into vector
    X = zeros(length(mask(mask==1)),1);
    X = actual(mask==1);
    
    % No modifications needed for recon
    Y = recon; 
    %X = actual;

    % Compute RMSE 
    AE = abs(Y(:)-X(:));
    
    % Compute MPE 
    FE = (100 * 2 * (Y-X) ./ (X + Y)); 
end


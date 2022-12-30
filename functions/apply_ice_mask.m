function [X_temp] = apply_ice_mask(X_temp,SIC)
    [~,~,dep] = size(X_temp);
    ice_mask = X_temp(:,:,1) > SIC;
    for i = 1:dep
        % Step 1: Apply ice mask
        X_temp2 = X_temp(:,:,i);
        X_temp2(~ice_mask) = NaN;
        X_temp(:,:,i) = X_temp2;
    end

end

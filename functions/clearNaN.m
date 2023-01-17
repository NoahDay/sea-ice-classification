function [X_out,row_idx,X_new] = clearNaN(X_temp)
    % Remove all the NaNs from the data 
    % Noah Day, UofA, Jan 2023

    [~,~,dep] = size(X_temp);
    X_new = [];
    X_temp_out = [];
    for i = 1:dep
        clear idx1 idx2 X_new
        % Step 1: removing NaN and taking SIC > 0.15
        idx1 = logical(isnan(X_temp(:,:,i))); % NaN = 1
        idx2= logical(X_temp(:,1,i) > 0.15);
        idx1(:,1) = ~idx2; % SIC < 0.15 or NaN = 1
        idx = logical(idx1) ; % Remove cells with either NaN or SIC < 0.15
        [len, wid] = size(X_temp(:,:,i));
        
        wid = wid - 2; % Take off lon and lat
        count = 1;
        for j = 1:len
            if sum(idx(j,:)) == 0 % If conditions are met store the data, else delete it
                X_new(count,:) = X_temp(j,:,i);
                count  = count + 1;
            end
        end
        
        % Calculate the row index of each layer so we can extract the
        % individual files from the 2D array later
        if isempty(X_new)
            row_idx(i) = 0;
        else
            [row_idx(i),~] = size(X_new(:,:));
        end
        
        % Append on the data for the current layer
        if ~isempty(X_new)
            X_temp_out = [X_temp_out; X_new(:,:)];
        end
    end

    X_out = X_temp_out;

end
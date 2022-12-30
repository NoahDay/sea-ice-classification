function [X_out,row_idx] = clearNaN(X_temp)
    [~,~,dep] = size(X_temp);
    X_new = [];
    for i = 1:dep

        % Step 1: removing NaN
        idx1 = isnan(X_temp(:,:,i));
        idx = logical(idx1);
        [len, wid] = size(X_temp(:,:,i));
        
        wid = wid - 2;
        count = 1;
        for j = 1:len
            if sum(idx(j,:)) == 0
                X_new(count,:,i) = X_temp(j,:,i);
                count  = count + 1;
            end
        end
        %X_new_unstandard = X_new;
        %clear X_new_unstandard
        if isempty(X_new)
            row_idx(i) = 0;
        else
            [row_idx(i),~] = size(X_new(:,:,i));
        end
    end

    X_temp = [];
    for i = 1:dep
        if ~isempty(X_new)
            X_temp = [X_temp; X_new(:,:,i)];
        end
    end
    X_out = X_temp;

end
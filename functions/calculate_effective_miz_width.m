function [MIZ_width, class,MIZ] = calculate_effective_miz_width(filename,sector,k_means,miz_def)
    % To measure the width of the identified marginal ice zone we calculate
    % the 'effective penetration' \citep{wadhams1975}, $x_e$, by 
    % integrating over the sea ice concentration $P_z$ along the meridional 
    % line $x$ such that, x_e = \int_0 ^x P_z dz.
    % The quantity $x_e$ represents the width of the MIZ assuming the ice 
    % cover is completely compact (i.e., has a SIC of 100%).


    aice = data_format_sector(filename,'aice',sector);
    [tlat,tlon,~,ulat,ulon] = grid_read('om2');
    ice_mask =  aice > 0.15;
    % Assign all data which isn't in the ice mask to NaN
    k_means(~ice_mask) = NaN; 
    idx_miz = isnan(k_means(1,:));
    [len,~] = size(aice);

    % Loop around the south pole, each i is a line of constant longitude
    for i = 1:len
        temp = k_means(i,~idx_miz);
        [len1,wid1] = size(temp);
        if len1 == 0
            edge_class(i) = 0;
        elseif wid1 == 0
            edge_class(i) = 0;
        else
            % Take the last cell that > 0.15 SIC
            edge_class(i) = temp(end); 
        end
    end

    idx = isnan(edge_class);
    class = mode(edge_class(~idx));
    MIZ = k_means == miz_def;%class;  
    clear MIZ_width
    
    for i = 1:len
        ulat_vec_temp = ulat(i,:);
        ulon_vec_temp = ulon(i,:);
        MIZ_vec_temp = MIZ(i,:);
        if sum(MIZ_vec_temp) == 0
            MIZ_width(i,1:2) = 0;
        else
            %south_point = find(MIZ_vec_temp, 1, 'first');
            %north_point = find(MIZ_vec_temp, 1, 'last');
            f = find(diff([0,MIZ_vec_temp,0]==1));
            p = f(1:2:end-1);  % Start indices
            y = f(2:2:end)-p;  % Consecutive ones counts
            max_y = find(y == max(y));
            p = p(max_y);
            y = y(max_y);
            south_point = p;
            north_point = p+y;
            % Use u-grid as it is positioned in the north east corner of
            % the cell
            MIZ_distance = pathdist([ulat_vec_temp(south_point-1),ulat_vec_temp(north_point)],[ulon_vec_temp(south_point),ulon_vec_temp(north_point)],'km');
            MIZ_width(i,1:2) = MIZ_distance(1:2);

            for j = south_point:north_point
                % Calculate the histogram P_z dz.
                temp(j,1:2) = aice(i,j).*pathdist([ulat_vec_temp(j-1),ulat_vec_temp(j)],[ulon_vec_temp(j-1),ulon_vec_temp(j)],'km');
            end
            % Integrate over, x_e = \int_0 ^x P_z dz.
            MIZ_width_corrected(i,1) = sum(temp(:,2));

        end
    end
    %MIZ_width = MIZ_width(:,2);
    MIZ_width = MIZ_width_corrected;
end


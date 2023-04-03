clc
clear all
cd /Users/noahday/GitHub/sea-ice-classification/
%%
cd /Users/noahday/GitHub/sea-ice-classification/
T = readtable('kmeans.csv');
% '/Users/noahday/Gadi/ia40/waves-10/history/kmeans.csv'
%
T.k;
lon = T.longitude;
lat = T.latitude;
k = T.k;
date = T.date;
ai = T.aice;
r = T.fsdrad;

clear T
%
date_idx = date == datetime('2002-04-01');
lon = lon(date_idx);
lat = lat(date_idx);
k = k(date_idx);
ai = ai(date_idx);
r = r(date_idx);
[LT,LN,~,ulat,ulon] = grid_read('om2');
k_means = ones(size(LN))*NaN;
aice = ones(size(LN))*NaN;
fsdrad = ones(size(LN))*NaN;
%

for i = 1:sum(date_idx)
    [lon_pos,lat_pos,~] = near2(LN,LT,lon(i),lat(i));
    k_means(lon_pos,lat_pos) = k(i)+1; 
    aice(lon_pos,lat_pos) = ai(i); 
    fsdrad(lon_pos,lat_pos) = r(i); 
    %k_means(~sector_mask) = NaN;
end
cd /Users/noahday/GitHub/sea-ice-classification/
%%
close all
kmeans_plot = ones(size(LN))*NaN;
kmeans_plot(aice>0.15) = k_means(aice>0.15);
figure
w = worldmap('world');
        axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
        setm(w, 'Origin', [-90 0 0]);
        setm(w, 'maplatlimit', [-90,-65]);
        setm(w, 'grid', 'off');
        setm(w, 'frame', 'off');
        setm(w, "FontColor",[0.5, 0.5, 0.5])
        setm(w, 'labelrotation', 'on')
        %setm(w, 'meridianlabel', 'on','FontSize',font_size)
        %setm(w, 'parallellabel', 'off','FontSize',font_size)
        setm(w, 'mlabellocation', 60);
        setm(w, 'plabellocation', 5);
        pcolorm(LT,LN,kmeans_plot,'FaceAlpha',0.99)
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(w, land, 'FaceColor', [0.5 0.5 0.5],'FaceAlpha',.5)
        colorbar; %cmocean('deep');
         scalebar('length',500,'location','se')
%%
figure
w = worldmap('world');
        axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
        setm(w, 'Origin', [-90 0 0]);
        setm(w, 'maplatlimit', [-90,-50]);
        setm(w, 'grid', 'off');
        setm(w, 'frame', 'off');
        setm(w, "FontColor",[0.5, 0.5, 0.5])
        setm(w, 'labelrotation', 'on')
        %setm(w, 'meridianlabel', 'on','FontSize',font_size)
        %setm(w, 'parallellabel', 'off','FontSize',font_size)
        setm(w, 'mlabellocation', 60);
        setm(w, 'plabellocation', 5);
        pcolorm(LT,LN,aice,'FaceAlpha',0.99)
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(w, land, 'FaceColor', [0.5 0.5 0.5],'FaceAlpha',.5)
        colorbar; %cmocean('deep');

figure
w = worldmap('world');
        axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
        setm(w, 'Origin', [-90 0 0]);
        setm(w, 'maplatlimit', [-90,-50]);
        setm(w, 'grid', 'off');
        setm(w, 'frame', 'off');
        setm(w, "FontColor",[0.5, 0.5, 0.5])
        setm(w, 'labelrotation', 'on')
        %setm(w, 'meridianlabel', 'on','FontSize',font_size)
        %setm(w, 'parallellabel', 'off','FontSize',font_size)
        setm(w, 'mlabellocation', 60);
        setm(w, 'plabellocation', 5);
        pcolorm(LT,LN,fsdrad,'FaceAlpha',0.99)
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(w, land, 'FaceColor', [0.5 0.5 0.5],'FaceAlpha',.5)
        colorbar; %cmocean('deep');
%%
pathname = '/Users/noahday/Gadi/ia40/kmeans_2009.nc';
ncdisp(pathname)
kmeans_data = ncread(pathname,'k');
lon_k = ncread(pathname,'LON');
lat_k = ncread(pathname,'LAT');
size(kmeans_data)
size(kmeans_data)

close all
figure
w = worldmap('world');
        axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
        setm(w, 'Origin', [-90 0 0]);
        setm(w, 'maplatlimit', [-90,-50]);
        setm(w, 'grid', 'off');
        setm(w, 'frame', 'off');
        setm(w, "FontColor",[0.5, 0.5, 0.5])
        setm(w, 'labelrotation', 'on')
        %setm(w, 'meridianlabel', 'on','FontSize',font_size)
        %setm(w, 'parallellabel', 'off','FontSize',font_size)
        setm(w, 'mlabellocation', 60);
        setm(w, 'plabellocation', 5);
        pcolorm(lat_k,lon_k,kmeans_data(:,:,180),'FaceAlpha',0.99)
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(w, land, 'FaceColor', [0.5 0.5 0.5],'FaceAlpha',.5)
        colorbar; %cmocean('deep');
%%
addpath functions
addpath /Users/noahday/Github/cice-plotting-tools/functions
addpath /Users/noahday/GitHub/CICE-analyser/processing
kmeans_data = ncread("/Users/noahday/Gadi/ia40/kmean_2009.nc","k");
for i = 1:365
   miz_def = 0;
  [MIZ_width(i,:)] = miz_width(kmeans_data(:,:,i),miz_def);
end
%%
size(MIZ_width)

close all
figure

plot(1:365, mean(MIZ_width'))
%%
filename1 = '/Users/noahday/Gadi/ia40/waves-10/iceh.2005-09-01.nc';
filename2 = '/Volumes/NoahDay5TB/WIM_on/history/fullyears/iceh.2015-09-01.nc';

lat1 = ncread(filename1,"TLAT");
lat2 = ncread(filename2,"TLAT");

lon1 = ncread(filename1,"TLON");
lon2 = ncread(filename2,"TLON");

fsdrad1 = ncread(filename1,"fsdrad");
fsdrad2 = ncread(filename2,"fsdrad");

close all
figure
w = worldmap('world');
        axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
        setm(w, 'Origin', [-90 0 0]);
        setm(w, 'maplatlimit', [-90,-50]);
        setm(w, 'grid', 'off');
        setm(w, 'frame', 'off');
        setm(w, "FontColor",[0.5, 0.5, 0.5])
        setm(w, 'labelrotation', 'on')
        %setm(w, 'meridianlabel', 'on','FontSize',font_size)
        %setm(w, 'parallellabel', 'off','FontSize',font_size)
        setm(w, 'mlabellocation', 60);
        setm(w, 'plabellocation', 5);
        pcolorm(lat1,lon1,fsdrad1-fsdrad2,'FaceAlpha',0.99)
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(w, land, 'FaceColor', [0.5 0.5 0.5],'FaceAlpha',.5)
        colorbar; %cmocean('deep');
%%


 %%
x = [12 24 36];
cpu_time = [1.66 1.5 2.2];
wall_time = [10.66 4 3.78];

close all
conFigure(6)
f = figure;
plot(x, cpu_time,'-o')
ylim([0,3])
ylabel('CPU time [minutes]')
xlabel("Number of CPUs")
exportgraphics(f,'cputime.pdf','ContentType','vector')


conFigure(6)
f = figure;
plot(x, wall_time,'-o')
ylim([0,12])
ylabel('Wall time [minutes]')
xlabel("Number of CPUs")
exportgraphics(f,'walltime.pdf','ContentType','vector')



function [MIZ_width] = miz_width(k_means,miz_def)

% Loop around the south pole, each i is a line of constant longitude
[len, wid, ~] = size(k_means);
    %for i = 1:len
    %    temp = k_means(i,~idx_miz);
    %    [len1,wid1] = size(temp);
    %    if len1 == 0
    %        edge_class(i) = 0;
    %    elseif wid1 == 0
    %        edge_class(i) = 0;
    %    else
            % Take the last cell that > 0.15 SIC
    %        edge_class(i) = temp(end); 
    %    end
    %end

    %idx = isnan(edge_class);
    %class = mode(edge_class(~idx));
    MIZ = k_means == miz_def;%class;  
    clear MIZ_width
    [LT,LN,~,ulat,ulon] = grid_read('om2');
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

            %for j = south_point:north_point
            %    % Calculate the histogram P_z dz.
            %    temp(j,1:2) = aice(i,j).*pathdist([ulat_vec_temp(j-1),ulat_vec_temp(j)],[ulon_vec_temp(j-1),ulon_vec_temp(j)],'km');
            %end
            % Integrate over, x_e = \int_0 ^x P_z dz.
            %MIZ_width_corrected(i,1) = sum(temp(:,2));

        end
    end
    MIZ_width = MIZ_width(:,2);
    %MIZ_width = MIZ_width_corrected;

end

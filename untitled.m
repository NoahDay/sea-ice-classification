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
pathname = '/Users/noahday/Gadi/ia40/waves-025/kmeans_2023-03-30.nc';
ncdisp(pathname)
kmeans_data = ncread(pathname,'k');
lon_k = ncread(pathname,'TLON');
size(kmeans_data)
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



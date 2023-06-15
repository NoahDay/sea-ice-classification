%% Plotting tools for python formatted files
close all
clear all
cd /Users/noahday/GitHub/sea-ice-classification/
Cmap = [0.9805, 0.5000, 0.4453; 0.4416, 0.7490, 0.4322; 0.3639, 0.5755, 0.748];
    
    
%% World map
filename = "/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_nofsd_2019.nc";
ncdisp(filename)
k = ncread(filename,"k");

close all
C1 = linspecer(3);
%Cmap = C1([2,3,1],:);
%Cmap(1,:) = [251,128,114]/256;
Cmap = [0.9805, 0.5000, 0.4453; 0.4416, 0.7490, 0.4322; 0.3639, 0.5755, 0.748];
Cmap = Cmap([2,3,1],:);

lat = ncread(filename,"LAT");
lon = ncread(filename,"LON");


plot_lon = lon;
plot_lon(end+1,:) = (lon(1,:));% + lon(end,:))/2;
plot_lat = lat;
plot_lat(end+1,:) = (lat(1,:) + lat(end,:))/2;
plot_data = k(:,:,90);
%plot_data = k(:,:,273);
plot_data(end+1,:) = plot_data(1,:);
landmask = ncread(filename,'tmask');
landmask(end+1,:) = landmask(1,:);

close all
f = figure;
w = worldmap('world');
    axesm eqaazim;
    setm(w, 'Origin', [-90 0 0]);
    %setm(w, 'maplatlimit', [-90,-53]);
    setm(w, 'maplatlimit', [-90,-63]);
    setm(w, 'maplonlimit', [0,360]);
    setm(w, 'grid', 'off');
    setm(w, 'frame', 'off')
    pcolorm(plot_lat,plot_lon,plot_data);
    colormap(Cmap)
    bedmap2('patchgl');
    %outlineashelf('all','color','k');
    %scalebar('color','k','units','km','location','se')
    scalebar('length',1000,...
        'units','km',...
        'color','k','location','sw',...
        'fontangle','italic','FontSize',14)
exportgraphics(f,'worldmap_2019-03-01_nofsd.pdf', 'ContentType', 'vector','BackgroundColor','none')
%% Mean variables in each class

% Read in CSV
var_list = ["aice", "hi", "hs", "fsdrad", "iage", "alvl", "longitude", "latitude", "date", "k"];
filename = "/Volumes/NoahDay1TB/sea-ice-classification/data/kmeans_2019.csv";
standard_table = readtable(filename);


filename = "/Volumes/NoahDay1TB/sea-ice-classification/data/raw_2019.csv";
raw_table = readtable(filename);

raw_table.k = standard_table.k;
%%
%mean(standard_table.k==0)

groupsummary(raw_table,"k","mean","aice")
%% Distribution of SIC, ITD, FSD for each class
nIDs = 3;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); % {'(a)','(b)','(c)','(d)'}
%text(0.025,0.95,charlbl{1},'Units','normalized')

font_size = 6;
%var_vec = [1,2,4];
%label_plot = label_vec(var_vec);
label_plot = {'SIC [-]', 'Thickness [m]', 'Floe size [m]'};
addpath /Users/noahday/GitHub/CICE-plotting-tools/functions/
addpath /Users/noahday/Documents/MATLAB/matlab2tikz/src/

clear A

bins = linspace(min(raw_table.aice), max(raw_table.aice), 101);

A.MIZ = histcounts(raw_table.aice(raw_table.k == 0),bins);
A.Thin = histcounts(raw_table.aice(raw_table.k == 1),bins);
A.Thick = histcounts(raw_table.aice(raw_table.k == 2),bins);

A.MIZ = A.MIZ./n_points;
A.Thin = A.Thin./n_points;
A.Thick = A.Thick./n_points;

A.MIZ(A.MIZ < 1e-4) = 1e-4;
A.Thin(A.Thin < 1e-4) = 1e-4;
A.Thick(A.Thick < 1e-4) = 1e-4;

line_width = 2;
n_points = length(raw_table.aice);
close all
conFigure(11,3)

f = figure;
set(gcf, 'Position',  [0, 0, 30, 10])

tiledlayout(1,3)
%nexttile_vec = [1 2 5 6];
nexttile
%nhist(A,'noerror')
histogram('BinCounts',A.MIZ, 'BinEdges', bins,'DisplayStyle','stairs','EdgeColor',Cmap(1,:),'LineWidth',2)
hold on
histogram('BinCounts',A.Thin, 'BinEdges', bins,'DisplayStyle','stairs','EdgeColor',Cmap(2,:),'LineWidth',2)
histogram('BinCounts',A.Thick, 'BinEdges', bins,'DisplayStyle','stairs','EdgeColor',Cmap(3,:),'LineWidth',2)
xlabel('Sea ice concentration [-]')
ylabel('Density [-]')
set(gca,'YScale','log')
set(gca,'XScale','linear')
%ylim([1e-4,1e0])
%legend("MIZ", "Thin", "Thick")
text(0.025,0.9,charlbl{1},'Units','normalized')
%exportgraphics(f,'kmean_hist_sic_total.pdf','ContentType','vector','BackgroundColor','none')
%matlab2tikz('kmean_hist_sic_total.tex', 'standalone', true);


clear A

bins = linspace(min(raw_table.hi), max(raw_table.hi), 101);

A.MIZ = histcounts(raw_table.hi(raw_table.k == 0),bins);
A.Thin = histcounts(raw_table.hi(raw_table.k == 1),bins);
A.Thick = histcounts(raw_table.hi(raw_table.k == 2),bins);

% A.MIZ = A.MIZ./n_points;
% A.Thin = A.Thin./n_points;
% A.Thick = A.Thick./n_points;
% 
% A.MIZ(A.MIZ < 1e-4) = 1e-4;
% A.Thin(A.Thin < 1e-4) = 1e-4;
% A.Thick(A.Thick < 1e-4) = 1e-4;

line_width = 2;
n_points = length(raw_table.hi);

nexttile
%nhist(A,'noerror')
histogram('BinCounts',A.MIZ./n_points, 'BinEdges', bins,'DisplayStyle','stairs','EdgeColor',Cmap(1,:),'LineWidth',2)
hold on
histogram('BinCounts',A.Thin./n_points, 'BinEdges', bins,'DisplayStyle','stairs','EdgeColor',Cmap(2,:),'LineWidth',2)
histogram('BinCounts',A.Thick./n_points, 'BinEdges', bins,'DisplayStyle','stairs','EdgeColor',Cmap(3,:),'LineWidth',2)
xlabel('Ice thickness [m]')
ylabel('Density [-]')
set(gca,'YScale','linear')
set(gca,'XScale','linear')
%ylim([1e-4,1e0])
legend("1", "2", "3")
text(0.025,0.9,charlbl{2},'Units','normalized')
%exportgraphics(f,'kmean_hist_sic_total.pdf','ContentType','vector','BackgroundColor','none')
%matlab2tikz('kmean_hist_sic_total.tex', 'standalone', true);


clear A

bins = linspace(min(raw_table.fsdrad), max(raw_table.fsdrad), 101);

A.MIZ = histcounts(raw_table.fsdrad(raw_table.k == 0),bins);
A.Thin = histcounts(raw_table.fsdrad(raw_table.k == 1),bins);
A.Thick = histcounts(raw_table.fsdrad(raw_table.k == 2),bins);

A.MIZ = A.MIZ./n_points;
A.Thin = A.Thin./n_points;
A.Thick = A.Thick./n_points;
 
A.MIZ(A.MIZ < 1e-4) = 1e-4;
A.Thin(A.Thin < 1e-4) = 1e-4;
A.Thick(A.Thick < 1e-4) = 1e-4;


line_width = 2;
n_points = length(raw_table.fsdrad);

nexttile
%nhist(A,'noerror')
histogram('BinCounts',A.MIZ, 'BinEdges', bins,'DisplayStyle','stairs','EdgeColor',Cmap(1,:),'LineWidth',2)
hold on
histogram('BinCounts',A.Thin, 'BinEdges', bins,'DisplayStyle','stairs','EdgeColor',Cmap(2,:),'LineWidth',2)
histogram('BinCounts',A.Thick, 'BinEdges', bins,'DisplayStyle','stairs','EdgeColor',Cmap(3,:),'LineWidth',2)
xlabel('Floe radius [m]')
ylabel('Density [-]')
set(gca,'YScale','log')
set(gca,'XScale','linear')
ylim([1e-4,1e0])
%legend("1", "2", "3")
text(0.025,0.9,charlbl{3},'Units','normalized')

fontsize(gcf,scale=2)
exportgraphics(f,'kmean_hist_all_total.pdf','ContentType','vector','BackgroundColor','none')
matlab2tikz('kmean_hist_all_total.tex', 'standalone', true);



%% Time series
clc
[a n] = unique(datenum(raw_table.date));
D1 = sortrows([a n],2);
unique_dates = datestr(D1(:,1));
unique_dates = datetime(unique_dates);
for i = 1:length(unique_dates)
    date_lp = unique_dates(i);
    mean_data(1,i) = mean(raw_table.fsdrad(raw_table.date==date_lp & raw_table.k==0));
    mean_data(2,i) = mean(raw_table.fsdrad(raw_table.date==date_lp & raw_table.k==1));
    mean_data(3,i) = mean(raw_table.fsdrad(raw_table.date==date_lp & raw_table.k==2));
    mean_data(4,i) = mean(raw_table.fsdrad(raw_table.date==date_lp));

    std_data(1,i) = std(raw_table.fsdrad(raw_table.date==date_lp & raw_table.k==0));
    std_data(2,i) = std(raw_table.fsdrad(raw_table.date==date_lp & raw_table.k==1));
    std_data(3,i) = std(raw_table.fsdrad(raw_table.date==date_lp & raw_table.k==2));
end
disp('Done')



%% Plot the ts
%% Floe size
close all
conFigure(11,2)

%ts_dates_CI  = [unique_dates', fliplr(unique_dates')];
%ts_data_CI = [mean_data-std_data, fliplr(mean_data+std_data)];

f = figure;
for i = 1:3
    plot(unique_dates, mean_data(i,:),"Color",Cmap(i,:))
    hold on
%     fill(ts_dates_CI, ts_data_CI(i,:) , 1,....
%         'facecolor',Cmap(i,:), ...
%         'edgecolor','none', ...
%         'facealpha', 0.2);
end
plot(unique_dates, mean_data(4,:))
ylabel('Floe size [m]')
set(gca,'YScale','linear')
exportgraphics(f,'ts_floesize.pdf','ContentType','vector','BackgroundColor','none')

%% Age
close all
conFigure(11,2)

ts_dates_CI  = [unique_dates', fliplr(unique_dates')];
ts_data_CI = [mean_data-std_data, fliplr(mean_data+std_data)];

f = figure;
for i = 1:3
    plot(unique_dates, mean_data(i,:),"Color",Cmap(i,:))
    hold on
    fill(ts_dates_CI, ts_data_CI(i,:) , 1,....
        'facecolor',Cmap(i,:), ...
        'edgecolor','none', ...
        'facealpha', 0.2);
end
ylabel('Ice age [years]')
set(gca,'YScale','linear')
exportgraphics(f,'ts_age.pdf','ContentType','vector','BackgroundColor','none')

%% Thickness
close all
conFigure(11,2)

ts_dates_CI  = [unique_dates', fliplr(unique_dates')];
ts_data_CI = [mean_data-std_data, fliplr(mean_data+std_data)];

f = figure;
for i = 1:3
    plot(unique_dates, mean_data(i,:),"Color",Cmap(i,:))
    hold on
    fill(ts_dates_CI, ts_data_CI(i,:) , 1,....
        'facecolor',Cmap(i,:), ...
        'edgecolor','none', ...
        'facealpha', 0.2);
end
ylabel('Ice thickness [m]')
set(gca,'YScale','linear')
exportgraphics(f,'ts_itd.pdf','ContentType','vector','BackgroundColor','none')

%% Alvl
close all
conFigure(11,2)

ts_dates_CI  = [unique_dates', fliplr(unique_dates')];
ts_data_CI = [mean_data-std_data, fliplr(mean_data+std_data)];

f = figure;
for i = 1:3
    plot(unique_dates, mean_data(i,:),"Color",Cmap(i,:))
    hold on
    fill(ts_dates_CI, ts_data_CI(i,:) , 1,....
        'facecolor',Cmap(i,:), ...
        'edgecolor','none', ...
        'facealpha', 0.2);
end
ylabel('Level ice [-]')
set(gca,'YScale','linear')
exportgraphics(f,'ts_sice.pdf','ContentType','vector','BackgroundColor','none')


%% Mean variables in each class

% Read in CSV
var_list = ['daidtt','daidtd','Tsfc','shear','divu','strength','frazil',...
    'congel','Tair','trsig','uvel','vvel','strairx','strairy','strocnx',...
    'strocny','strintx','strinty','strcorx','strcory','wave_sig_ht',...
    'peak_period','sst','frzmlt', "longitude", "latitude", "date", "k"];
filename = "/Volumes/NoahDay1TB/sea-ice-classification/data/analysis_raw_2019.csv";
analysis_table = readtable(filename);
%
filename = "/Volumes/NoahDay1TB/sea-ice-classification/data/kmeans_2019.csv";
standard_table = readtable(filename);

analysis_table.k = standard_table.k;

%% Time series
clc
[a n] = unique(datenum(analysis_table.date));
D1 = sortrows([a n],2);
unique_dates = datestr(D1(:,1));
unique_dates = datetime(unique_dates);
for i = 1:length(unique_dates)
    date_lp = unique_dates(i);
    for class = 0:2
        %mean_data(class+1,i) = mean(analysis_table.fsdrad(analysis_table.date==date_lp & analysis_table.k==class));
        %std_data(class+1,i) = std(analysis_table.fsdrad(analysis_table.date==date_lp & analysis_table.k==class));
        %mean_data(2,i) = mean(raw_table.daidtt(raw_table.date==date_lp & raw_table.k==1));
        %mean_data(3,i) = mean(raw_table.daidtt(raw_table.date==date_lp & raw_table.k==2));
        %std_data(2,i) = std(raw_table.daidtt(raw_table.date==date_lp & raw_table.k==1));
        %std_data(3,i) = std(raw_table.daidtt(raw_table.date==date_lp & raw_table.k==2));
    end

end
disp('Done')
%% V vel and daidtd
[a n] = unique(datenum(analysis_table.date));
D1 = sortrows([a n],2);
unique_dates = datestr(D1(:,1));
unique_dates = datetime(unique_dates);
for i = 1:length(unique_dates)
    date_lp = unique_dates(i);
    for class = 0
        mean_data_vvel(class+1,i) = mean(analysis_table.vvel(analysis_table.date==date_lp & analysis_table.k==class));
        std_data_vvel(class+1,i) = std(analysis_table.vvel(analysis_table.date==date_lp & analysis_table.k==class));

        mean_data_daidtd(class+1,i) = mean(analysis_table.daidtd(analysis_table.date==date_lp & analysis_table.k==class));
        std_data_daidtd(class+1,i) = std(analysis_table.daidtd(analysis_table.date==date_lp & analysis_table.k==class));

    end

end
disp('Done')
%%
close all
conFigure(11,1.5)

ts_dates_CI  = [unique_dates', fliplr(unique_dates')];
ts_data_CI = [mean_data-std_data, fliplr(mean_data+std_data)];

f = figure;
%for i = 1:3
yyaxis left
plot(unique_dates, movmean(mean_data_vvel(:),15))%,"Color",Cmap(1,:))
ylabel("Zonal velocity [m/s]")
hold on
yyaxis right
plot(unique_dates, -movmean(mean_data_daidtd(:),15))%,"Color",Cmap(2,:))
ylabel("Loss of SIC from dynamics [$\%$/day]")
%      fill(ts_dates_CI, ts_data_CI(i,:) , 1,....
%          'facecolor',Cmap(i,:), ...
%          'edgecolor','none', ...
%          'facealpha', 0.2);
%end
%ylim([0,1])
%ylabel('Sig. wave height [$m$]')
set(gca,'YScale','linear')
exportgraphics(f,'ts_vvel_daidtd.pdf','ContentType','vector','BackgroundColor','none')

%%

close all
conFigure(11,1.5)

ts_dates_CI  = [unique_dates', fliplr(unique_dates')];
ts_data_CI = [mean_data-std_data, fliplr(mean_data+std_data)];

f = figure;
for i = 1:3
    plot(unique_dates, movmean(mean_data(i,:),15),"Color",Cmap(i,:))
    hold on
%      fill(ts_dates_CI, ts_data_CI(i,:) , 1,....
%          'facecolor',Cmap(i,:), ...
%          'edgecolor','none', ...
%          'facealpha', 0.2);
end
%ylim([0,1])
ylabel('Sig. wave height [$m$]')
set(gca,'YScale','linear')
exportgraphics(f,'ts_wave_sig_ht.pdf','ContentType','vector','BackgroundColor','none')
%%
close all
conFigure(11,1.5)

ts_dates_CI  = [unique_dates', fliplr(unique_dates')];
ts_data_CI = [mean_data-std_data, fliplr(mean_data+std_data)];

f = figure;
for i = 1:3
    plot(unique_dates, movmean(mean_data(i,:),30),"Color",Cmap(i,:))
    hold on
%      fill(ts_dates_CI, ts_data_CI(i,:) , 1,....
%          'facecolor',Cmap(i,:), ...
%          'edgecolor','none', ...
%          'facealpha', 0.2);
end
%ylim([0,1])
ylabel('Meridional velocity [m/s]')
set(gca,'YScale','linear')
exportgraphics(f,'ts_frzmlt.pdf','ContentType','vector','BackgroundColor','none')


%% Time series velocities
clc
[a n] = unique(datenum(analysis_table.date));
D1 = sortrows([a n],2);
unique_dates = datestr(D1(:,1));
unique_dates = datetime(unique_dates);
for i = 1:length(unique_dates)
    date_lp = unique_dates(i);
    for class = 0:2
        temp_datax = analysis_table.uvel(analysis_table.date==date_lp & analysis_table.k==class);
        temp_datay = analysis_table.vvel(analysis_table.date==date_lp & analysis_table.k==class);
        mean_data(class+1,i) = mean(sqrt(temp_datax.^2 + temp_datay.^2));
        std_data(class+1,i) = std(sqrt(temp_datax.^2 + temp_datay.^2));
    end

end
disp('Done')

%%
close all
conFigure(11,1.5)

ts_dates_CI  = [unique_dates', fliplr(unique_dates')];
ts_data_CI = [mean_data-std_data, fliplr(mean_data+std_data)];

f = figure;
for i = 1:3
    plot(unique_dates, movmean(mean_data(i,:),10),"Color",Cmap(i,:))
    hold on
%     fill(ts_dates_CI, ts_data_CI(i,:) , 1,....
%         'facecolor',Cmap(i,:), ...
%         'edgecolor','none', ...
%         'facealpha', 0.2);
end
ylim([0,0.3])
ylabel('Ice speed [m/s]')
set(gca,'YScale','linear')
exportgraphics(f,'ts_speed.pdf','ContentType','vector','BackgroundColor','none')




%% Changes in FSD

% Read in CSV
%var_list = ['daidtt','daidtd','Tsfc','shear','divu','strength','frazil',...
%    'congel','Tair','trsig','uvel','vvel','strairx','strairy','strocnx',...
%    'strocny','strintx','strinty','strcorx','strcory','wave_sig_ht',...
%    'peak_period','sst','frzmlt', "longitude", "latitude", "date", "k"];
filename = "/Volumes/NoahDay1TB/sea-ice-classification/data/kmeans_2018.csv";
standard_table = readtable(filename);

filename = "/Volumes/NoahDay1TB/sea-ice-classification/data/analysis_fsd_raw_2018.csv";
analysis_fsd_table = readtable(filename);
analysis_fsd_table.k = standard_table.k;

%% Heatmap
close all
miz_dt = mean(diff(effective,1));
%reshape 
% find(lon == 0.5)   281
plot_data = (diff(effective(:,:),1));
plot_data = [plot_data(281:end,:); plot_data(1:280,:)];

plot_lon = [lon(281:end); lon(1:280)];

conFigure(30)
figure
h = pcolor(ts_dates(1:end),plot_lon(1:end-1), movmean(plot_data,20,1));
 %h = pcolor();
 set(h, 'EdgeColor', 'none');
cmocean('balance','pivot');
%clim([-200,200])
cb = colorbar;
cb.Label.String = 'Change in MIZ width [m/day]';
xlabel 'Days of the year'
ylabel 'Longitude'


%% Time series
clc
count = 1;
n_year = length(2011:2019);
for year_lp = 2019:2019
    disp(num2str(year_lp))
    filename = strcat('/Volumes/NoahDay1TB/sea-ice-classification/data/kmeans_',num2str(year_lp),'.csv');
    standard_table = readtable(filename);
    filename = strcat('/Volumes/NoahDay1TB/sea-ice-classification/data/analysis_fsd_raw_',num2str(year_lp),'.csv');
    analysis_fsd_table = readtable(filename);
    analysis_fsd_table.k = standard_table.k;

    [a n] = unique(datenum(analysis_fsd_table.date));
    D1 = sortrows([a n],2);
    unique_dates = datestr(D1(:,1));
    unique_dates = datetime(unique_dates);
    for i = 1:length(unique_dates)
        date_lp = unique_dates(i);
        for class = 0:2
            mean_data_wave(class+1,i,count) = sum(analysis_fsd_table.dafsd_wave(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
            std_data_wave(class+1,i,count) = std(analysis_fsd_table.dafsd_wave(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
    
            mean_data_newi(class+1,i,count) = sum(analysis_fsd_table.dafsd_newi(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
            std_data_newi(class+1,i,count) = std(analysis_fsd_table.dafsd_newi(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
    
            mean_data_weld(class+1,i,count) = sum(analysis_fsd_table.dafsd_weld(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
            std_data_weld(class+1,i,count) = std(analysis_fsd_table.dafsd_weld(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
        end
    
    end
    count = count + 1;
end
disp('Done')

%%
close all
%addpath functions
conFigure(11,0.8)

%ts_dates_CI  = [unique_dates', fliplr(unique_dates')];
%ts_data_CI = [mean_data-std_data, fliplr(mean_data+std_data)];
plot_wave = mean(mean_data_wave,3);
plot_weld = mean(mean_data_weld,3);
plot_newi = mean(mean_data_newi,3);

f = figure;
t = tiledlayout(3,1,'TileSpacing',"tight");
nexttile
for i = 1:3
    plot(unique_dates, -movmean(plot_wave(i,:),30),"Color",Cmap(i,:),'LineWidth',1.5)
     hold on
%      fill(ts_dates_CI, ts_data_CI(i,:) , 1,....
%          'facecolor',Cmap(i,:), ...
%          'edgecolor','none', ...
%          'facealpha', 0.2);
end
ylabel('Decrease in floe size [m/day]')
xticklabels('')
title('Wave induced breakup')
set(gca,'YScale','linear')
%ylim([0.01,300])

nexttile
for i = 1:3
    plot(unique_dates, movmean(plot_weld(i,:),30),"Color",Cmap(i,:),'LineWidth',1.5)
     hold on
%      fill(ts_dates_CI, ts_data_CI(i,:) , 1,....
%          'facecolor',Cmap(i,:), ...
%          'edgecolor','none', ...
%          'facealpha', 0.2);
end
ylabel('Change in floe size [m/day]')
xticklabels('')
title('New ice formation')
set(gca,'YScale','linear')
%ylim([0.01,300])

nexttile
for i = 1:3
    plot(unique_dates, movmean(plot_newi(i,:),30),"Color",Cmap(i,:),'LineWidth',1.5)
     hold on
%      fill(ts_dates_CI, ts_data_CI(i,:) , 1,....
%          'facecolor',Cmap(i,:), ...
%          'edgecolor','none', ...
%          'facealpha', 0.2);
end
ylabel('Change in floe size [m/day]')
title('Welding')

%ylim([0.01,300])
fontsize(gcf,scale=1)
set(gca,'YScale','linear')
exportgraphics(f,'ts_dafsd_all.pdf','ContentType','vector','BackgroundColor','none')
%exportgraphics(f,'ts_dafsd_all.png','ContentType','image')

%%
close all
cmap_dafsd = [252,141,89; 120,198,121; 8,104,172]/256;
conFigure(11,0.5)
f = figure;
t = tiledlayout(3,1,'TileSpacing',"tight");

title_string = {"MIZ","Thin interior","Thick interior"};
for i = 1:3
    nexttile
    p = plot(unique_dates, movmean(plot_wave(i,:),30),"Color",cmap_dafsd(1,:),'LineWidth',1.5);
    %area(unique_dates, movmean(plot_wave(i,:),30) , 1,....
    %      'facecolor',cmap_dafsd(1,:), ...
    %      'edgecolor','none', ...
    %      'facealpha', 0.2);
    ylabel("Change in floe size [m/day]")
    xticklabels({[]})
    title(title_string{i})
    hold on
    plot(unique_dates, movmean(plot_weld(i,:),30),"Color",cmap_dafsd(2,:),'LineWidth',1.5)
    ylabel("Change in floe size [m/day]")
     plot(unique_dates, movmean(plot_newi(i,:),30),"Color",cmap_dafsd(3,:),'LineWidth',1.5)
    ylim([-3,3])
    ylabel("Change in floe size [m/day]")
    if i == 1
        legend("Waves","Welding","New ice",'Orientation','horizontal')
    end
    if i == 3
        xticklabels('auto')
    end
end
%%
close all
area_plot_data = cumsum([movmean(plot_newi(i,:),30);  movmean(plot_weld(i,:),30)])';

figure
t = tiledlayout(3,1,'TileSpacing',"tight");
for i = 1:3
    nexttile
    %p = plot(unique_dates, movmean(plot_wave(i,:),30),"Color",cmap_dafsd(1,:),'LineWidth',1.5);
    area(unique_dates, movmean(plot_wave(i,:),30))
    hold on
    area(unique_dates, area_plot_data)% , 1,....
        % 'facecolor',cmap_dafsd(1,:), ...
        % 'edgecolor','none', ...
        % 'facealpha', 0.2);
    %hold on
    %area(unique_dates, area_plot_data , 1,....
         %'facecolor',cmap_dafsd(:,:), ...
         %'edgecolor','none', ...
         %'facealpha', 0.2);
    %ylabel("Change in floe size [m/day]")
    xticklabels({[]})
    title(title_string{i})
    hold on
    %plot(unique_dates, movmean(plot_weld(i,:),30),"Color",cmap_dafsd(2,:),'LineWidth',1.5)
    %ylabel("Change in floe size [m/day]")
    % plot(unique_dates, movmean(plot_newi(i,:),30),"Color",cmap_dafsd(3,:),'LineWidth',1.5)
    ylim([-3,3])
    ylabel("Change in floe size [m/day]")
    if i == 1
        legend("Waves","New ice","Welding",'Orientation','horizontal')
    end
    if i == 3
        xticklabels('auto')
    end
end
%%
figure
Y = [1, 5, 3;
     3, 2, 7;
     1, 5, 3;
     2, 6, 1];
area(Y)
grid on
colormap summer
set(gca,'Layer','top')
title 'Stacked Area Plot'


%% Time series comparing MIZ with and without floe size
clc
count = 1;
n_year = length(2011:2018);
for year_lp = 2019:2019
    disp(num2str(year_lp))
    filename = strcat('/Volumes/NoahDay1TB/sea-ice-classification/data/kmeans_',num2str(year_lp),'.csv');
    standard_table = readtable(filename);
    filename = strcat('/Volumes/NoahDay1TB/sea-ice-classification/data/analysis_fsd_raw_',num2str(year_lp),'.csv');
    analysis_fsd_table = readtable(filename);
    analysis_fsd_table.k = standard_table.k;

    [a n] = unique(datenum(analysis_fsd_table.date));
    D1 = sortrows([a n],2);
    unique_dates = datestr(D1(:,1));
    unique_dates = datetime(unique_dates);
    for i = 1:length(unique_dates)
        date_lp = unique_dates(i);
        for class = 0:2
            mean_data_wave(class+1,i,count) = mean(analysis_fsd_table.dafsd_wave(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
            std_data_wave(class+1,i,count) = std(analysis_fsd_table.dafsd_wave(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
    
            mean_data_newi(class+1,i,count) = mean(analysis_fsd_table.dafsd_newi(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
            std_data_newi(class+1,i,count) = std(analysis_fsd_table.dafsd_newi(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
    
            mean_data_weld(class+1,i) = mean(analysis_fsd_table.dafsd_weld(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
            std_data_weld(class+1,i,count) = std(analysis_fsd_table.dafsd_weld(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
        end
    
    end
    % No_Fsd
    filename = strcat('/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_nofsd_',num2str(year_lp),'.csv');
    standard_table = readtable(filename);
    filename = strcat('/Volumes/NoahDay1TB/sea-ice-classification/data/analysis_fsd_raw_',num2str(year_lp),'.csv');
    analysis_fsd_table = readtable(filename);
    analysis_fsd_table.k = standard_table.k;
    for i = 1:length(unique_dates)
        date_lp = unique_dates(i);
        for class = 0:2
            mean_data_wave_no_fsd(class+1,i,count) = mean(analysis_fsd_table.dafsd_wave(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
            std_data_wave_no_fsd(class+1,i,count) = std(analysis_fsd_table.dafsd_wave(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
    
            mean_data_newi_no_fsd(class+1,i,count) = mean(analysis_fsd_table.dafsd_newi(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
            std_data_newi_no_fsd(class+1,i,count) = std(analysis_fsd_table.dafsd_newi(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
    
            mean_data_weld_no_fsd(class+1,i) = mean(analysis_fsd_table.dafsd_weld(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
            std_data_weld_no_fsd(class+1,i,count) = std(analysis_fsd_table.dafsd_weld(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
        end
    
    end
    count = count + 1;
end
disp('Done')

%%

close all
figure
pcolor(1:365, 1:3, mean_data_newi(:,:))
colorbar

%%
cmap_dafsd = [252,141,89; 120,198,121; 8,104,172]/256;

close all
%addpath functions
conFigure(11,0.8)

%ts_dates_CI  = [unique_dates', fliplr(unique_dates')];
%ts_data_CI = [mean_data-std_data, fliplr(mean_data+std_data)];

f = figure;
t = tiledlayout(3,1,'TileSpacing',"tight");
nexttile

plot(unique_dates, -movmean(mean_data_wave(1,:),30),"Color",cmap_dafsd(1,:),'LineWidth',1.5)
hold on
plot(unique_dates, -movmean(mean_data_wave_no_fsd(3,:),30),"Color",cmap_dafsd(1,:)*0.5,'LineWidth',1.5)

%      fill(ts_dates_CI, ts_data_CI(i,:) , 1,....
%          'facecolor',Cmap(i,:), ...
%          'edgecolor','none', ...
%          'facealpha', 0.2);

ylabel('Decrease in $r_a$ [m/day]')
xticklabels('')
title('Wave induced breakup')
legend('With floe size','Without floe size')
%
nexttile
plot(unique_dates, movmean(mean_data_newi(1,:),30),"Color",cmap_dafsd(2,:),'LineWidth',1.5)
hold on
plot(unique_dates, movmean(mean_data_newi_no_fsd(3,:),30),"Color",cmap_dafsd(2,:)*0.5,'LineWidth',1.5)

ylabel('Change in $r_a$ [m/day]')
xticklabels('')
title('New ice formation')

nexttile
plot(unique_dates, movmean(mean_data_weld(1,:),30),"Color",cmap_dafsd(3,:),'LineWidth',1.5)
hold on
plot(unique_dates, movmean(mean_data_weld_no_fsd(3,:),30),"Color",cmap_dafsd(3,:)*0.5,'LineWidth',1.5)

ylabel('Change in $r_a$ [m/day]')
title('Welding')

%ylim([0,500])
fontsize(gcf,scale=1)
set(gca,'YScale','linear')
exportgraphics(f,'ts_dafsd_miz_comp.pdf','ContentType','vector','BackgroundColor','none')
%exportgraphics(f,'ts_dafsd_all.png','ContentType','image')
%
%%
close all
%addpath functions
conFigure(11,1.5)

%ts_dates_CI  = [unique_dates', fliplr(unique_dates')];
%ts_data_CI = [mean_data-std_data, fliplr(mean_data+std_data)];

f = figure;
plot(unique_dates, movmean(mean_data_wave(1,:),30),"Color",cmap_dafsd(1,:),'LineWidth',1.5)
hold on
area(unique_dates, movmean(mean_data_wave(1,:),30) , 1,....
          'facecolor',cmap_dafsd(1,:), ...
          'edgecolor','none', ...
          'facealpha', 0.2);
plot(unique_dates, movmean(mean_data_newi(1,:),30),"Color",cmap_dafsd(2,:),'LineWidth',1.5)
area(unique_dates, movmean(mean_data_newi(1,:),30) , 1,....
          'facecolor',cmap_dafsd(2,:), ...
          'edgecolor','none', ...
          'facealpha', 0.2);
plot(unique_dates, movmean(mean_data_weld(1,:),30),"Color",cmap_dafsd(3,:),'LineWidth',1.5)
area(unique_dates, movmean(mean_data_weld(1,:),30) , 1,....
          'facecolor',cmap_dafsd(3,:), ...
          'edgecolor','none', ...
          'facealpha', 0.2);
legend('Wave breakup','New floe formation','Welding of floes')
ylabel('Change in mean floe radius [m/day]')
exportgraphics(f,'ts_dafsd_miz_allproc.pdf','ContentType','vector','BackgroundColor','none')
%%

clc
count = 1;
n_year = length(2019:2019);

for i = 1:length(unique_dates)
    date_lp = unique_dates(i);
    for class = 0:2
        mean_data_wave(class+1,i,count) = mean(analysis_fsd_table.dafsd_wave(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
        std_data_wave(class+1,i,count) = std(analysis_fsd_table.dafsd_wave(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));

        mean_data_newi(class+1,i,count) = mean(analysis_fsd_table.dafsd_newi(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
        std_data_newi(class+1,i,count) = std(analysis_fsd_table.dafsd_newi(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));

        mean_data_weld(class+1,i) = mean(analysis_fsd_table.dafsd_weld(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
        std_data_weld(class+1,i,count) = std(analysis_fsd_table.dafsd_weld(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));

        mean_data_latg(class+1,i) = mean(analysis_fsd_table.dafsd_latg(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
        std_data_latg(class+1,i,count) = std(analysis_fsd_table.dafsd_latg(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));

        mean_data_latm(class+1,i) = mean(analysis_fsd_table.dafsd_latm(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
        std_data_latm(class+1,i,count) = std(analysis_fsd_table.dafsd_latm(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
    end

end

%%

area_plot = [movmean(mean_data_wave(1,:),30); movmean(mean_data_latm(1,:),30); movmean(mean_data_newi(1,:),30); movmean(mean_data_weld(1,:),30); movmean(mean_data_latg(1,:),30);];
f = figure;
area(unique_dates, area_plot(1:2,:)')
hold on
area(unique_dates, area_plot(3:5,:)')% , 1,....

          %'edgecolor','none', ...
          %'facealpha', 0.2);
legend('Wave breakup', 'Lateral melt', 'New floe formation','Welding of floes', 'Lateral growth')
ylabel('Change in mean floe radius [m/day]')
colormap(cmap_dafsd)
exportgraphics(f,'ts_dafsd_miz_allproc.pdf','ContentType','vector','BackgroundColor','none')



%% Spider plot

% Read in CSV
var_list = ["SIC [$\%$]", "Thickness [m]", "Snow thickness [m]", "Floe size [m]", "Age [years]", "Level ice [$\%$]"];

temp_P = zeros(3,6);
count = 1;
for year_lp = 2011:2018
    disp(num2str(year_lp))

    filename = strcat("/Volumes/NoahDay1TB/sea-ice-classification/data/kmeans_",num2str(year_lp),".csv");
    standard_table = readtable(filename);
    
    disp("Standard read")
    filename = strcat("/Volumes/NoahDay1TB/sea-ice-classification/data/raw_",num2str(year_lp),".csv");
    raw_table = readtable(filename);
    disp("Raw read")
    raw_table.k = standard_table.k;
    D = groupsummary(raw_table,"k","mean");
    temp_P = temp_P + table2array(D(:,4:end-3))./length(2011:2018);
    
end
disp("Done!")
%% plot it
addpath /Users/noahday/GitHub/CICE-plotting-tools/functions/
addpath /Users/noahday/Documents/MATLAB/matlab2tikz/src/

temp_P = [0.7490    0.3541    0.0385   25.1202 0.2710    0.6655;
          0.9850    0.6547    0.1262  250.2048 0.2393    0.8598; 
          0.9635    1.3164    0.2137  728.3942 0.4047    0.7247];


P = temp_P;
P(:,1) = P(:,1)*100;
P(:,6) = P(:,6)*100;
P = P(end:-1:1,:);
var_switch = [1, 2, 6, 3, 5, 4];
P = P(:,var_switch);


% Initialize data points
%D1 = [5 3 9 1 2];
%D2 = [5 8 7 2 9];
%D3 = [8 2 1 4 6];
%P = [D1; D2; D3];
% Delete variable in workspace if exists

% Spider plot
close all 
conFigure(11,1.7)

f = figure;
set(gcf, 'Position',  [0, 0, 15, 10])
if exist('s', 'var')
    delete(s);
end
s = spider_plot_class(P);
% Legend properties
 s.LegendLabels = {'Cluster 1', 'Cluster 2', 'Cluster 3'};
 s.LegendVisible = 'off';
 s.LegendHandle.Location = 'northeastoutside';
 axes_labels = cellstr(var_list(1:6));
 s.AxesLabels = axes_labels(var_switch);
 s.FillOption = {'on', 'off', 'off'};
 s.Color = Cmap(end:-1:1,:);
 s.LineWidth = 2.5;
 s.Marker = 'none';
 s.AxesFontSize = 12;
 s.LabelFontSize = 12;
 s.FillOption = {'on', 'on', 'on'};
 s.FillTransparency = [0.1, 0.1, 0.1];
 s.AxesColor = [0.8, 0.8, 0.8];
 axes_limits = [0, 0, 0, 0, 0, 0; 100, 2, 0.4, 800, 0.7, 100];
 s.AxesLimits = axes_limits(:,var_switch); %[1, 1, 0, 1, 1, 1; 100, 100, 100, 100, 100, 100];
 s.AxesPrecision = [0, 1, 1, 1, 1, 1];
 s.AxesLabelsEdge = 'none';
 s.AxesRadial = 'off';
 s.AxesPrecision = 1;
 s.AxesInterpreter = 'latex';
 s.AxesScaling = 'linear';
 s.AxesInterval = 2;
 %s.AxesWebType = 'circular';

%s.AxesLimits = [1, 2, 1, 1, 1, 1; 10, 8, 9, 5, 10, 10];
%s.AxesFont = 'Times New Roman';
%s.LabelFont = 'Times New Roman';

matlab2tikz('spider_plot.tex', 'standalone', true);
exportgraphics(f,'spider_plot.pdf','ContentType','vector','BackgroundColor','none')
%%
clc
filename = "/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_nofsd_2019.csv";
standard_table = readtable(filename);

disp("Standard read")
filename = strcat("/Volumes/NoahDay1TB/sea-ice-classification/data/raw_2019.csv");
raw_table = readtable(filename);
disp("Raw read")
k_vec = standard_table.k;
thick_mask = k_vec == 0;
thin_mask = k_vec == 1;
miz_mask = k_vec == 2;
k_vec(thin_mask) = 2;
k_vec(thick_mask) = 1;
k_vec(miz_mask) = 0;
raw_table.k = k_vec;
D = groupsummary(raw_table,"k","mean");

temp_P = table2array(D(:,4:end-3));


%%
P = temp_P;
P(:,1) = P(:,1)*100;
P(:,6) = P(:,6)*100;
P = P(end:-1:1,:);
var_switch = [1, 2, 6, 3, 5, 4];
P = P(:,var_switch);


% Initialize data points
%D1 = [5 3 9 1 2];
%D2 = [5 8 7 2 9];
%D3 = [8 2 1 4 6];
%P = [D1; D2; D3];
% Delete variable in workspace if exists

% Spider plot
close all 
conFigure(11,1.7)

f = figure;
set(gcf, 'Position',  [0, 0, 15, 10])
if exist('s', 'var')
    delete(s);
end
s = spider_plot_class(P);
% Legend properties
 s.LegendLabels = {'Cluster 1', 'Cluster 2', 'Cluster 3'};
 s.LegendVisible = 'off';
 s.LegendHandle.Location = 'northeastoutside';
 axes_labels = cellstr(var_list(1:6));
 s.AxesLabels = axes_labels(var_switch);
 s.FillOption = {'on', 'off', 'off'};
 s.Color = Cmap(end:-1:1,:);
 s.LineWidth = 2.5;
 s.Marker = 'none';
 s.AxesFontSize = 12;
 s.LabelFontSize = 12;
 s.FillOption = {'on', 'on', 'on'};
 s.FillTransparency = [0.1, 0.1, 0.1];
 s.AxesColor = [0.8, 0.8, 0.8];
 axes_limits = [0, 0, 0, 0, 0, 0; 100, 2, 0.4, 800, 0.7, 100];
 s.AxesLimits = axes_limits(:,var_switch); %[1, 1, 0, 1, 1, 1; 100, 100, 100, 100, 100, 100];
 s.AxesPrecision = [0, 1, 1, 1, 1, 1];
 s.AxesLabelsEdge = 'none';
 s.AxesRadial = 'off';
 s.AxesPrecision = 1;
 s.AxesInterpreter = 'latex';
 s.AxesScaling = 'linear';
 s.AxesInterval = 2;
 %s.AxesWebType = 'circular';

%s.AxesLimits = [1, 2, 1, 1, 1, 1; 10, 8, 9, 5, 10, 10];
%s.AxesFont = 'Times New Roman';
%s.LabelFont = 'Times New Roman';

%matlab2tikz('spider_plot.tex', 'standalone', true);
exportgraphics(f,'spider_plot_nofsd.pdf','ContentType','vector','BackgroundColor','none')

%%
%filename = 
tarea = ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_2011.nc','tarea');
count = 0;
for year_lp = 11:19
    k = ncread(strcat('/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_20',num2str(year_lp),'.nc'),'k');
    aice = ncread(strcat('/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_20',num2str(year_lp),'.nc'),'aice');
    count = count + 1;
    for i = 1:365
        for class = 0:2
            k_day = k(:,:,i);
            idx = k_day == class;
            areas_year(class+1,i,count) = sum(nansum(tarea(idx)));
        end
    end
end
areas = mean(areas_year,3);

%% Sea ice extent

t1 = datetime(2019,1,1,0,0,0);
t2 = datetime(2019,12,31,0,0,0);
t = t1:t2;
conFigure(11,2)
f = figure;
for i = 1:3
    plot(t,areas(i,:)./(10^12),'color',Cmap(i,:),'linewidth',1.5)
    hold on
end
%Returns handles to the patch and line objects
%chi=get(gca, 'Children');
%Reverse the stacking order so that the patch overlays the line
%set(gca, 'Children',flipud(chi))

%colormap(Cmap)
ylabel('Sea ice extent [$10^6$ km$^2$]')
legend('Cluster 1','Cluster 2','Cluster 3','location','northwest');
legend boxoff  
set(gca,'box','off')

exportgraphics(f,'cluster-SIE.pdf', 'ContentType', 'vector','BackgroundColor','none')
addpath /Users/noahday/Documents/MATLAB/matlab2tikz/src/
matlab2tikz('cluster-SIE.tex', 'standalone', true);

%% MIZ % of ice cover
close all
oct_start = find(t == datetime(2019,10,1,0,0,0));

changed_dates = t;
changed_dates([oct_start:365]) = changed_dates([oct_start:365]) -365;
changed_dates = changed_dates([oct_start:365, 1:oct_start-1]);
conFigure(11,2)
f = figure;
plot(changed_dates,movmean(areas(1,[oct_start:365, 1:oct_start-1])./(sum(areas(:,[oct_start:365, 1:oct_start-1]))),30),'color',Cmap(1,:),'linewidth',1.5)

%Returns handles to the patch and line objects
%chi=get(gca, 'Children');
%Reverse the stacking order so that the patch overlays the line
%set(gca, 'Children',flipud(chi))

%colormap(Cmap)
ylabel('MIZ fraction')
legend('Cluster 1','Cluster 2','Cluster 3','location','northeast');
legend boxoff  
set(gca,'box','off')
ylim([0,1])
exportgraphics(f,'cluster-SIE-fraction.pdf', 'ContentType', 'vector','BackgroundColor','none')
addpath /Users/noahday/Documents/MATLAB/matlab2tikz/src/
matlab2tikz('cluster-SIE-fraction.tex', 'standalone', true);

%% Mapping the amount of pancake ice

filename = "/Volumes/NoahDay1TB/sea-ice-classification/data/analysis_fsd_raw_2019.csv";
analysis_table = readtable(filename);

%[len, wid] = size(analysis_table);
%%
[a n] = unique(datenum(analysis_fsd_table.date));

D1 = sortrows([a n],2);
unique_dates = datestr(D1(:,1));
unique_dates = datetime(unique_dates);
[LAT,LON] = grid_read('om2');

[len,wid] = size(LAT);
pancake_map = zeros(len, wid, length(a));

for day_lp = 1:length(unique_dates)
    date_lp = unique_dates(day_lp);
    temp_table = analysis_fsd_table(analysis_fsd_table.date==date_lp,:);
    [len, wid] = size(temp_table);

    for row_lp = 1:len
        [loc1, loc2] = near2(LON,LAT,analysis_table.longitude(row_lp),analysis_table.latitude(row_lp));
        pancake_map(loc1,loc2,day_lp) = analysis_table.pancake(row_lp);
    end
end

%%
clear areas pancake_area
tarea = ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_2019.nc','tarea');

for year_lp = 19:19
    k = ncread(strcat('/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_20',num2str(year_lp),'.nc'),'k');
    aice = ncread(strcat('/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_20',num2str(year_lp),'.nc'),'aice');
    fsdrad = ncread(strcat('/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_20',num2str(year_lp),'.nc'),'fsdrad');
    hi = ncread(strcat('/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_20',num2str(year_lp),'.nc'),'hi');
    iage = ncread(strcat('/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_20',num2str(year_lp),'.nc'),'iage');
    count = count + 1;
    for i = 1:365
        iage_day = iage(:,:,i);
        iage_mask = iage_day < 1/6;
        %k_day = k(:,:,i);
        %k_mask = k_day == 0;
        idx = logical(iage_mask);
        pancake_temp = pancake_map(1:360,:,i); 
        pancake_area(i) = sum(sum(tarea(idx).*pancake_temp(idx)));
%         for class = 0:2
%             k_day = k(:,:,i);
%             k_mask = k_day == class;
%             %fsdrad_day = fsdrad(:,:,i);
%             %fsdrad_mask = fsdrad_day < 20;
%             %hi_day = hi(:,:,i);
%             %hi_mask = hi_day < 0.2;
%             %iage_day = iage(:,:,i);
%             %iage_mask = iage_day < 1/4;
%             %idx = k_mask.*iage_mask;%k_mask.*fsdrad_mask.*hi_mask.*iage_mask;
%             %pancake_area(class+1,i) = sum(sum(tarea(logical(idx))));
%             
% 
%             k_day = k(:,:,i);
%             idx = k_day == class;
%             SIE_area(class+1,i) = sum(nansum(tarea(idx)));
%         end
    end
end
count = 0;
for year_lp = 15:19
    k = ncread(strcat('/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_20',num2str(year_lp),'.nc'),'k');
    %aice = ncread(strcat('/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_20',num2str(year_lp),'.nc'),'aice');
    fsdrad = ncread(strcat('/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_20',num2str(year_lp),'.nc'),'fsdrad');
    %hi = ncread(strcat('/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_20',num2str(year_lp),'.nc'),'hi');
    %iage = ncread(strcat('/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_20',num2str(year_lp),'.nc'),'iage');
    disp(num2str(year_lp))
    count = count + 1;
    for i = 1:365
        fsdrad_temp = squeeze(fsdrad(:,:,i));
        for class = 0:2
            k_day = k(:,:,i);
            k_mask = k_day == class;
            k_day = k(:,:,i);
            idx = k_day == class;
            SIE_area(class+1,i,count) = sum(nansum(tarea(idx)));
            fsdrad_vec(class+1,i,count) = nanmean(nanmean(fsdrad_temp(idx)));
        end
    end
end
%%
count = 1;
n_year = length(2011:2019);
for year_lp = 2018:2019
    disp(num2str(year_lp))
    filename = strcat('/Volumes/NoahDay1TB/sea-ice-classification/data/kmeans_',num2str(year_lp),'.csv');
    standard_table = readtable(filename);
    filename = strcat('/Volumes/NoahDay1TB/sea-ice-classification/data/analysis_fsd_raw_',num2str(year_lp),'.csv');
    analysis_fsd_table = readtable(filename);
    analysis_fsd_table.k = standard_table.k;

    [a n] = unique(datenum(analysis_fsd_table.date));
    D1 = sortrows([a n],2);
    unique_dates = datestr(D1(:,1));
    unique_dates = datetime(unique_dates);
    for i = 1:length(unique_dates)
        date_lp = unique_dates(i);
        for class = 0:2
            mean_data_wave(class+1,i,count) = sum(analysis_fsd_table.dafsd_wave(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
            std_data_wave(class+1,i,count) = std(analysis_fsd_table.dafsd_wave(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
    
            mean_data_newi(class+1,i,count) = sum(analysis_fsd_table.dafsd_newi(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
            std_data_newi(class+1,i,count) = std(analysis_fsd_table.dafsd_newi(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
    
            mean_data_weld(class+1,i,count) = sum(analysis_fsd_table.dafsd_weld(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
            std_data_weld(class+1,i,count) = std(analysis_fsd_table.dafsd_weld(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
        end
    
    end
    count = count + 1;
end
disp('Done')
%areas = mean(areas_year,3);
%%
close all
t1 = datetime(2019,1,1,0,0,0);
t2 = datetime(2019,12,31,0,0,0);
t = t1:t2;

% SIE 
SIE_area_temp = SIE_area;
SIE_area = mean(SIE_area_temp,3);
conFigure(11,3)

%f = figure;
f = figure('Position',[0, 0, 40, 10]);
conFigure(11,4)


plot_data = [movmean(SIE_area(1,:)./10^(12),1); movmean(SIE_area(2,:)./10^(12),1); movmean(SIE_area(3,:)./10^(12),1); movmean(pancake_area./10^(12),1)];

ax = plot(t, [plot_data(1:4,:)],'linewidth',3);

colormap(Cmap)

legend('Unsupervised MIZ','Cluster 2','Cluster 3','Pancake ice','Rate of wave breakup','location','northwest','orientation','vertical');
legend boxoff  
set(gca,'box','off')
ax(1).Color = Cmap(1,:);
ax(2).Color = Cmap(2,:);
ax(3).Color = Cmap(3,:);
ax(4).Color = [136,65,157]/256;

ylabel 'Sea ice extent [$10^6$ km$^2$]';
% 
% hAx(2).YColor = [247,104,161]/256;
% hAx(2).YLabel.String = 'Decrease in floe radius from wave breakup [m/day]';
% 
% hLine2(1).Color = [247,104,161]/256;
%ylim([0,1])

exportgraphics(f,'SIE2.pdf', 'ContentType', 'vector','BackgroundColor','none')
addpath /Users/noahday/Documents/MATLAB/matlab2tikz/src/
matlab2tikz('SIE2.tex', 'standalone', true);


%%
close all
conFigure(11,3)
f = figure('Position',[0, 0, 40, 10]);
conFigure(11,4)
set(gca,'FontSize',24) 

%f = figure;
plot_data = [movmean(mean(SIE_area(1,:,:),3)./10^(12),1); movmean(pancake_area./10^(12),1)];
%plot_data = [movmean(nanmean(fsdrad_vec(1,:,end),3),1)];

  %[movmean(SIE_area(1,:)./10^(12),1); movmean(SIE_area(2,:)./10^(12),1); movmean(SIE_area(3,:)./10^(12),1); movmean(pancake_area./10^(12),1); mean(fsdrad_vec(1,:,:),3)];

[hAx,hLine1,hLine2] = plotyy(t, [plot_data(:,:)],t,[movmean(nanmean(fsdrad_vec(1,:,end),3),1); -movmean(mean(mean_data_wave(1,:,end),3),30)/10]); %; movmean(mean(mean_data_newi(1,:,2),3),30); movmean(mean(mean_data_weld(1,:,2),3),30)] );

%colormap(Cmap)

legend('Unsupervised MIZ','Pancake ice','Floe size','Wave breakup in the MIZ','location','northwest','orientation','vertical');
%legend('MIZ floe size','Decrease in floe size from wave breakup')
legend boxoff  
set(gca,'box','off')
%hAx(1).YLim = [0,150];
hLine1(1).Color = Cmap(1,:);
hLine1(2).Color =  Cmap(1,:); % [254,196,79]/256; %[35,139,69]/256;
hLine1(2).LineStyle = '--';
hLine1(1).LineWidth = 3;
hLine1(2).LineWidth = 2;
%hLine1(3).Color = Cmap(3,:);
hAx(1).YColor = Cmap(1,:);

%hAx(1).YColor = [0,0,0]/256;
hAx(1).YLabel.String = 'Sea ice extent [$10^6$ km$^2$]';

hAx(2).YColor = [136,65,157]/256;
hAx(2).YLabel.String = 'Floe radius [m] \& Decrease in floe radius [10 m/day]';
hAx(2).YLim = [0,150];
hLine2(1).LineWidth = 2;
hLine2(2).LineWidth = 2;

hLine2(1).Color = [136,65,157]/256;
hLine2(2).Color = [136,65,157]/256;
hLine2(2).LineStyle = '--';
%hLine2(1).LineWidth = 3;
%ylim([0,1])
exportgraphics(f,'miz-changes.pdf', 'ContentType', 'vector','BackgroundColor','none')
addpath /Users/noahday/Documents/MATLAB/matlab2tikz/src/
matlab2tikz('miz-changes.tex', 'standalone', true);

%%
swh = ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_2019.nc','wave_sig_ht');
aice = ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_2019.nc','aice');
miz_width = ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/mizWidth_2019.nc','effective');
%%


[n_lon,n_lat,n_day] = size(swh);

for lon_lp = 1:n_lon
    for day_lp = 1:n_day
        ice_edge = find(aice(lon_lp, :, day_lp) > 0.15);
        if isempty(ice_edge)
            ice_edge = 1;
        else
            ice_edge = ice_edge(end); % + for ~0.5 lat north of ice edge
        end
        swh_north(lon_lp, day_lp) = swh(lon_lp, ice_edge, day_lp);
    end
end

%%
clc
close all
addpath /Users/noahday/GitHub/CICE-plotting-tools/functions/
addpath /Users/noahday/Documents/MATLAB/matlab2tikz/src/

plot_swh = swh_north(:,244:244+30);
plot_swh(plot_swh < 10^-1) = NaN;
plot_miz_width = miz_width(:,244:244+30);
plot_miz_width(plot_miz_width == 0) = NaN;
%plot_swh(isnan(plot_swh)) = 0;
%
conFigure(30,1)
f = figure;
scatter_kde(reshape(plot_swh(:),[],1), reshape(plot_miz_width(:)./1000,[],1), 'filled', 'MarkerSize', 100);
%xlim([0,10])
cb = colorbar();
cb.Label.String = 'Density';
set(gca,'ColorScale','log')
clim([10^-5, 10^-2])
xlabel('Significant wave height [m]')
ylabel('Effective MIZ width [km]')

exportgraphics(f,'swh_vs_miz.pdf', 'ContentType', 'vector','BackgroundColor','none')


%swh_width

plot_swh = swh_north(:,244:244+30);
plot_swh(plot_swh < 10^-1) = NaN;
plot_miz_width = swh_width(244:244+30,:,9)';
plot_miz_width(plot_miz_width == 0) = NaN;
%plot_swh(isnan(plot_swh)) = 0;
%
conFigure(30,1)
f = figure;
scatter_kde(reshape(plot_swh(:),[],1), reshape(plot_miz_width(:)./1000,[],1), 'filled', 'MarkerSize', 100);
%xlim([0,10])
cb = colorbar();
cb.Label.String = 'Density';
set(gca,'ColorScale','log')
clim([10^-5, 10^-2])
xlabel('Significant wave height [m]')
ylabel('Wave penetration distance ($H_s > \sigma$) [km]')
%hold on
%plot(mdl)

exportgraphics(f,'swh_vs_swhMIZ.pdf', 'ContentType', 'vector','BackgroundColor','none')
%%
close all


mdl = fitlm(reshape(plot_swh(:),[],1), reshape(plot_miz_width(:)./1000,[],1));

%%
% Generate data
x = normrnd(10,1,1000,1);
y = x*3 + normrnd(10,1,1000,1);
x(1) = NaN;
% Plot data using probability density estimate as function
figure(1);
scatter_kde(miz_width(1:1000)', y, 'filled', 'MarkerSize', 100);
% Add Color bar
cb = colorbar();
cb.Label.String = 'Probability density estimate';



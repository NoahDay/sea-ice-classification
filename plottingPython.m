%% Plotting tools for python formatted files
close all
clear all
cd /Users/noahday/GitHub/sea-ice-classification/
Cmap = [0.9805, 0.5000, 0.4453; 0.4416, 0.7490, 0.4322; 0.3639, 0.5755, 0.748];
    
    
%% World map
filename = "/Users/noahday/GitHub/sea-ice-classification/data/kmean_2018.nc";
ncdisp(filename)
k = ncread(filename,"k");

close all
C1 = linspecer(3);
Cmap = C1([2,3,1],:);
Cmap(1,:) = [251,128,114]/256;
lat = ncread(filename,"LAT");
lon = ncread(filename,"LON");


plot_lon = lon;
plot_lon(end+1,:) = (lon(1,:));% + lon(end,:))/2;
plot_lat = lat;
plot_lat(end+1,:) = (lat(1,:) + lat(end,:))/2;
plot_data = k(:,:,90);
plot_data(end+1,:) = plot_data(1,:);
landmask = ncread(filename,'tmask');
landmask(end+1,:) = landmask(1,:);

close all
f = figure;
w = worldmap('world');
    axesm eqaazim;
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-53]);
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
        'fontangle','italic')%,'FontSize',font_size)
exportgraphics(f,'worldmap_2018-03-01.pdf', 'ContentType', 'vector')
%% Mean variables in each class

% Read in CSV
var_list = ["aice", "hi", "hs", "fsdrad", "iage", "alvl", "longitude", "latitude", "date", "k"];
filename = "/Users/noahday/GitHub/sea-ice-classification/data/kmeans_2018.csv";
standard_table = readtable(filename);


filename = "/Users/noahday/GitHub/sea-ice-classification/data/raw_2018.csv";
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
%exportgraphics(f,'kmean_hist_sic_total.pdf','ContentType','vector')
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
%exportgraphics(f,'kmean_hist_sic_total.pdf','ContentType','vector')
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
exportgraphics(f,'kmean_hist_all_total.pdf','ContentType','vector')
matlab2tikz('kmean_hist_all_total.tex', 'standalone', true);



%% Time series
clc
[a n] = unique(datenum(raw_table.date));
D1 = sortrows([a n],2);
unique_dates = datestr(D1(:,1));
unique_dates = datetime(unique_dates);
for i = 1:length(unique_dates)
    date_lp = unique_dates(i);
    mean_data(1,i) = mean(raw_table.alvl(raw_table.date==date_lp & raw_table.k==0));
    mean_data(2,i) = mean(raw_table.alvl(raw_table.date==date_lp & raw_table.k==1));
    mean_data(3,i) = mean(raw_table.alvl(raw_table.date==date_lp & raw_table.k==2));

    std_data(1,i) = std(raw_table.alvl(raw_table.date==date_lp & raw_table.k==0));
    std_data(2,i) = std(raw_table.alvl(raw_table.date==date_lp & raw_table.k==1));
    std_data(3,i) = std(raw_table.alvl(raw_table.date==date_lp & raw_table.k==2));
end
disp('Done')



%% Plot the ts
%% Floe size
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
ylabel('Floe size [m]')
set(gca,'YScale','linear')
exportgraphics(f,'ts_floesize.pdf','ContentType','vector')

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
exportgraphics(f,'ts_age.pdf','ContentType','vector')

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
exportgraphics(f,'ts_itd.pdf','ContentType','vector')

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
exportgraphics(f,'ts_sice.pdf','ContentType','vector')



%% 
%% Mean variables in each class

% Read in CSV
var_list = ['daidtt','daidtd','Tsfc','shear','divu','strength','frazil',...
    'congel','Tair','trsig','uvel','vvel','strairx','strairy','strocnx',...
    'strocny','strintx','strinty','strcorx','strcory','wave_sig_ht',...
    'peak_period','sst','frzmlt', "longitude", "latitude", "date", "k"];
filename = "/Users/noahday/GitHub/sea-ice-classification/data/analysis_raw_2018.csv";

analysis_table = readtable(filename);

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
        mean_data(class+1,i) = mean(analysis_table.wave_sig_ht(analysis_table.date==date_lp & analysis_table.k==class));
        std_data(class+1,i) = std(analysis_table.wave_sig_ht(analysis_table.date==date_lp & analysis_table.k==class));
        %mean_data(2,i) = mean(raw_table.daidtt(raw_table.date==date_lp & raw_table.k==1));
        %mean_data(3,i) = mean(raw_table.daidtt(raw_table.date==date_lp & raw_table.k==2));
        %std_data(2,i) = std(raw_table.daidtt(raw_table.date==date_lp & raw_table.k==1));
        %std_data(3,i) = std(raw_table.daidtt(raw_table.date==date_lp & raw_table.k==2));
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
    plot(unique_dates, movmean(mean_data(i,:),30),"Color",Cmap(i,:))
    hold on
%      fill(ts_dates_CI, ts_data_CI(i,:) , 1,....
%          'facecolor',Cmap(i,:), ...
%          'edgecolor','none', ...
%          'facealpha', 0.2);
end
%ylim([0,1])
ylabel('Significant wave height [m]')
set(gca,'YScale','linear')
exportgraphics(f,'ts_swh.pdf','ContentType','vector')


%% Time series velocities
clc
[a n] = unique(datenum(analysis_table.date));
D1 = sortrows([a n],2);
unique_dates = datestr(D1(:,1));
unique_dates = datetime(unique_dates);
for i = 1:length(unique_dates)
    date_lp = unique_dates(i);
    for class = 0:2
        temp_datax = analysis_table.strairx(analysis_table.date==date_lp & analysis_table.k==class);
        temp_datay = analysis_table.strairy(analysis_table.date==date_lp & analysis_table.k==class);
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
    plot(unique_dates, mean_data(i,:),"Color",Cmap(i,:))
    hold on
    fill(ts_dates_CI, ts_data_CI(i,:) , 1,....
        'facecolor',Cmap(i,:), ...
        'edgecolor','none', ...
        'facealpha', 0.2);
end
ylabel('Ice-air stress [N/m$^2$]')
set(gca,'YScale','linear')
exportgraphics(f,'ts_strair.pdf','ContentType','vector')




%% Changes in FSD

% Read in CSV
%var_list = ['daidtt','daidtd','Tsfc','shear','divu','strength','frazil',...
%    'congel','Tair','trsig','uvel','vvel','strairx','strairy','strocnx',...
%    'strocny','strintx','strinty','strcorx','strcory','wave_sig_ht',...
%    'peak_period','sst','frzmlt', "longitude", "latitude", "date", "k"];
filename = "/Users/noahday/GitHub/sea-ice-classification/data/analysis_fsd_raw_2018.csv";

analysis_fsd_table = readtable(filename);

analysis_fsd_table.k = standard_table.k;


%% Time series
clc
[a n] = unique(datenum(analysis_fsd_table.date));
D1 = sortrows([a n],2);
unique_dates = datestr(D1(:,1));
unique_dates = datetime(unique_dates);
for i = 1:length(unique_dates)
    date_lp = unique_dates(i);
    for class = 0:2
        mean_data(class+1,i) = mean(analysis_fsd_table.dafsd_weld(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
        std_data(class+1,i) = std(analysis_fsd_table.dafsd_weld(analysis_fsd_table.date==date_lp & analysis_fsd_table.k==class));
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
    plot(unique_dates, movmean(mean_data(i,:),30),"Color",Cmap(i,:))
     hold on
%      fill(ts_dates_CI, ts_data_CI(i,:) , 1,....
%          'facecolor',Cmap(i,:), ...
%          'edgecolor','none', ...
%          'facealpha', 0.2);
end
ylabel('Change in $r_a$ from welding [m/day]')
%ylim([0,500])
set(gca,'YScale','linear')
exportgraphics(f,'ts_dafsd_weld.pdf','ContentType','vector')
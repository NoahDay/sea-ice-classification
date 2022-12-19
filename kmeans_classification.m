% Classification of CICE ice covers

% We want to experiment with different combination of variables over
% multiple years of data. 
% 1. Only ice cover "snapshot" variables
% 2. Add dynamic terms
% 3. Add thermodynamic terms
% 4. All of the terms

% Each case will have its own method:
% 1. Get filenames
% 2. Read in data
% 3. Clear NaNs
% 4. Standardise the data
% 5.Classify using kmeans

% Setup
clear all
clc
close all

%%
addpath functions
addpath /Users/noahday/GitHub/CICE-analyser/processing

historydir = '/Volumes/NoahDay5TB/WIM_on/history/fullyears/';

a = dir([historydir '/*.nc']);
n_files = numel(a);

for i = 1:n_files
   filenames(i,:) = strcat(historydir,a(i).name);
   dirdates(i,:) = a(i).name(6:end-3);
end


%% Select the relevant variables
% Columns are data
clc
sector = "SH";
close all
%var_list = {'aice','hi','hs','fsdrad','sice','iage','vlvl','vrdg', 'uvel','vvel','strength','divu','shear','daidtd','daidtt','dagedtd','dagedtt','dafsd_latg','dafsd_latm','dafsd_newi','dafsd_weld','dafsd_wave'};
var_list = {'aice','hi','hs','fsdrad','sice','iage','vlvl','vrdg'};
% Static: {'aice','hi','hs','fsdrad','sice','iage','vlvl','vrdg'}
% Dynamics: {'aice','hi','hs','fsdrad','sice','iage','vlvl','vrdg', 'uvel','vvel','strength','divu','shear','daidtd','daidtt','dagedtd','dagedtt','dafsd_latg','dafsd_latm','dafsd_newi','dafsd_weld','dafsd_wave'};

[X_raw, row_idx]= read_data_vec(filenames,sector,var_list); % [var_list, lon, lat]

%clear X_temp
label_vec = variable_dict(var_list);
size(X_raw)
%

%label_vec = variable_dict(var_list);
%
data.Xunstandard = X_raw;
data.row_idx = row_idx;
save_filename = strcat('cover_5percent_2015-19.mat');
save(save_filename,'data','-v7.3');
clear data

%%
var_list = {'aice','hi','hs','fsdrad','sice','iage','vlvl','vrdg'};

label_vec = variable_dict(var_list);
%% Load in data
clear X row_idx
load('cover_5percent_2013-14.mat')
%X.waves = data.Xunstandard;
%row_idx.waves = data.row_idx;



%% Clean the data
clear Xnan_temp 
load('cover_5percent_2015-19.mat')
SIC = 0.15;
[~,wid,dep] = size(data.Xunstandard);
X_temp = data.Xunstandard;%X_raw;
    
for i = 1:dep
    ice_mask = X_temp(:,1,i) > SIC;
    % Step 1: Apply ice mask
    for j = 1:wid
        X_temp(~ice_mask,j,i) = NaN;
    end 
end
idx = [];
Xnan = [];
row_idx = [];
for i = 1:dep
    idxtemp = isnan(X_temp(:,1,i));
    for k = 1:length(idxtemp(:,1))
        idx(k) = prod(~idxtemp(k,:));
    end
    idx = logical(idx);
    for j = 1:wid
        Xnan_temp(:,j) = X_temp(idx',j,i);
    end
    Xnan = [Xnan; Xnan_temp];
    [row, ~] = size(Xnan_temp);
    row_idx(i) = row;
    clear Xnan_temp row
end
clear Xtemp
%%
close all
conFigure(12,10)
f = figure('Position',[0, 0, 60, 10]);
for i = 1:length(label_vec)
  subplot(2,ceil(length(label_vec)/2),i)
  Xtemp(:,i) = Xnan(idx,i);
  hist(Xtemp(:,i),20)
  %hist(Xnan(idx,i),20)
  %xticks(1:length(average_stats(:,i)))
  title(label_vec{i})
end
exportgraphics(f,'distributionVariables.pdf','ContentType','vector')
%%


[Xnan,row_idx] = clearNaN(X_temp);
disp('Clear NaN done!')

dimension = size(Xnan);


% Xnan = NaN's have been removed.
%%
%Xnan.all = Xnan.waves;
%%
cleaned_data.Xnan = Xnan;
cleaned_data.label_vec = label_vec;
cleaned_data.row_idx = row_idx;
cleaned_data.dimension = dimension;
save('cleanded_cover_data2015-19.mat','cleaned_data','-v7.3');

%%
load('cleanded_cover_data2015-19.mat')
Xnan = cleaned_data.Xnan;
dimension = cleaned_data.dimension;
row_idx = cleaned_data.row_idx;
clear cleaned_data
%% Standardise the data
X_standard_all = Xnan;
[~,wid] = size(X_standard_all);

for j = 1:wid-2 % Don't standardise latitude and longitude
    if min(X_standard_all(:,j)) < 0
        % If negatives, then recentre so its all positive
        max_X(j) = max(X_standard_all(:,j));
        min_X(j) = min(X_standard_all(:,j));
        X_standard_all(:,j) = (X_standard_all(:,j) - min_X(j))/(max_X(j) - min_X(j));
    end
    %Step 1: Log transformation as almost all of the data is highly skewed
        X_standard_all(:,j) = log(X_standard_all(:,j)+1);
    %Step 2: Standardization
    % Calculate mean
        mean_X(j) = mean(X_standard_all(:,j));
        % Calculate standard deviation
        std_X(j) = std(X_standard_all(:,j));
        X_standard_all(:,j) = X_standard_all(:,j)/std_X(j);
    % Standardise the data
        max_X(j) = max(X_standard_all(:,j));
        min_X(j) = min(X_standard_all(:,j));
        X_standard_all(:,j) = (X_standard_all(:,j) - min_X(j))/(max_X(j) - min_X(j));
end
disp('Standardisation done!')
%clear X_temp
%%
eva = evalclusters(X_standard_all(:,1:end-2),'kmeans','CalinskiHarabasz','KList',1:5);
temp = eva.CriterionValues;
conFigure(11)
f = figure;
plot(eva)
ylim([0,1.2*max(eva.CriterionValues)])
exportgraphics(f,'calinskiHarabasz_cover.pdf','ContentType','vector')

%% k-means clustering
rng(2022)
idx_temp = 10000000:16339048;
aice_idx = X_standard_all(:,1)>0.15;
X = X_standard_all(:,1:end-2);
%X = X_standard_all(:,[1:3 5:end-2]);
%
num_clusters = 4;
tic
[kmeans_idx,C] = kmeans(X,num_clusters,'MaxIter',300);
%[kmeans_idx_without_floe,C] = kmeans(X,num_clusters,'MaxIter',300);
toc


% %%
% kmeans_cluster.idx = idx;
% kmeans_cluster.row_idx = row_idx;
% kmeans_cluster.label_vec = label_vec;
% kmeans_cluster.C = C;
% kmeans_cluster.X_standard_all = X_standard_all;
% kmeans_cluster.num_clusters = num_clusters;
% 
% save_filename = strcat('kmeans_3cases_',sprintf('%g',num_clusters),'_classes.mat');
% save(save_filename,'kmeans_cluster','-v7.3');
% %%
% %save_filename = strcat('kmeans_3cases_',sprintf('%g',num_clusters),'_classes.mat');
% %load(save_filename);
% 
% C = kmeans_cluster.C;
% idx = kmeans_cluster.idx;
% row_idx = kmeans_cluster.row_idx;
% X_standard_all = kmeans_cluster.X_new_unstandard;%X_standard_all;
% num_clusters = kmeans_cluster.num_clusters;
% label_vec = kmeans_cluster.label_vec;

%% Change the numbers of the kmeans classes
% Swap 1 and 3

% idx1  = kmeans_idx == 1;
% idx3 = kmeans_idx == 3;
% kmeans_idx(idx1) = 3;
% kmeans_idx(idx3) = 1;
% 
% idx2  = kmeans_idx == 2;
% idx3 = kmeans_idx == 3;
% kmeans_idx(idx2) = 3;
% kmeans_idx(idx3) = 2;

idx1  = kmeans_idx == 1;
idx3 = kmeans_idx == 3;
kmeans_idx(idx1) = 3;
kmeans_idx(idx3) = 1;

idx2  = kmeans_idx == 2;
idx3 = kmeans_idx == 3;
kmeans_idx(idx2) = 3;
kmeans_idx(idx3) = 2;


clear idx1 idx3
%%
idx1  = kmeans_idx_without_floe == 2;
idx3 = kmeans_idx_without_floe == 3;
kmeans_idx_without_floe(idx1) = 3;
kmeans_idx_without_floe(idx3) = 2;

idx2  = kmeans_idx_without_floe == 1;
idx3 = kmeans_idx_without_floe == 3;
kmeans_idx_without_floe(idx2) = 3;
kmeans_idx_without_floe(idx3) = 1;
%% Average stats
%Xnan = Xnan(aice_idx,:);
average_stats = [];
X_temp = [];
X_temp = Xnan;
X_temp = [X_temp, kmeans_idx]; % _without_floe

%X_temp = [X_temp2(:,[1:3 5:end-2]), kmeans_idx];
clear X_temp2
ylabs = {"SIC [\%]", "Ice thickness [m]", "Snow thickness [m]", "Mean floe size [m]",  "Bulk ice salinity [ppt]", "Ice age [years]", "Level ice volume", "Ridged ice volume"};

%"Mean floe radius [m]",
for i = 1:num_clusters
    index_class_data = X_temp(:,end) == i;
    class_data = X_temp(index_class_data,1:end-1);
    average_stats(i,:) = mean(class_data);
end


%%
close all
conFigure(11)
f = figure('Position',[0, 0,40, 6]);% 6, 24]);
label_vec_temp = label_vec;
%label_vec = label_vec([1:3 5:end]);
var_idx = [1,6,2,3,4];%[1,2,5];%

for i = 1:length(var_idx) %1:length(label_vec) % 
  %subplot(2,ceil(length(label_vec(var_idx))/2),i)
  subplot(1,length(var_idx),i)
  if i == 1
    b = bar(average_stats(:,var_idx(i)).*100, 'facecolor', 'flat');
    b.CData = Cmap;
  else
      b = bar(average_stats(:,var_idx(i)), 'facecolor', 'flat');
      b.CData = Cmap;
  end
  xticks(1:length(average_stats(:,i)))
  %title(label_vec{i})
  ylabel(ylabs(var_idx(i)),'FontSize',14)% ylabs
  %xlabel('Sea ice class','FontSize',14)
  if var_idx(i) == 6
      ylim([0,1])
  elseif var_idx(i) == 3
      ylim([0,0.5])
  elseif var_idx(i) == 4
      ylim([0,850])
  end
end
%sgtitle("\textbf{without floe size}",'FontSize',15)
exportgraphics(f,strcat('stat_comparison_',sprintf('%g',num_clusters),'_with_fsd_cover_15_percent_clusters.pdf'),'ContentType','vector')
label_vec = label_vec_temp;

%% Spider plot
% Initialize data points
P = average_stats(:,1:8);
%{'SIC','Ice thick.','Snow thick.','FSD','Salinity'}
% Spider plot
close all
figure
spider_plot(P,...
    'AxesLabels', label_vec,...
    'AxesInterval', 2,...
    'FillOption', {'on', 'on','on','on'},...
    'FillTransparency', 0.2*ones(1,num_clusters),...
    'AxesLimits', [zeros(1,8);1,3,1,850,15,1,1,2]);
%% Comparison of MIZ with and without floe size
ylabs = {"SIC [\%]", "Ice thickness [m]", "Snow thickness [m]", "Mean floe size [m]",  "Bulk ice salinity [ppt]", "Ice age [years]", "Level ice volume", "Ridged ice volume"};

% With floe size
average_stats = [];
X_temp = [];
X_temp = Xnan;
X_temp = [X_temp, kmeans_idx]; 
index_class_data = X_temp(:,end) == 3;
class_data = X_temp(index_class_data,1:end-1);
average_stats(1,:) = mean(class_data);
clear Xtemp class_data

% Without floe size
X_temp = [];
X_temp = Xnan;
X_temp = [X_temp, kmeans_idx_without_floe]; 
index_class_data = X_temp(:,end) == 3;
class_data = X_temp(index_class_data,1:end-1);
average_stats(2,:) = mean(class_data);
clear Xtemp

%
Cmap_comp = [Cmap(3,:); [141,211,199]./258];
close all
conFigure(11)
f = figure('Position',[0, 0,40, 6]);% 6, 24]);
label_vec_temp = label_vec;
%label_vec = label_vec([1:3 5:end]);
var_idx = [1,6,2,3,4];%[1,2,5];%

for i = 1:length(var_idx) %1:length(label_vec) % 
  %subplot(2,ceil(length(label_vec(var_idx))/2),i)
  subplot(1,length(var_idx),i)
  if i == 1
    b = bar(average_stats(:,var_idx(i)).*100, 'facecolor', 'flat');
    b.CData = Cmap_comp;
  else
      b = bar(average_stats(:,var_idx(i)), 'facecolor', 'flat');
      b.CData = Cmap_comp;
  end
  xticks(1:length(average_stats(:,i)))
  %title(label_vec{i})
  ylabel(ylabs(var_idx(i)),'FontSize',14)% ylabs
  %xlabel('Sea ice class','FontSize',14)
  if var_idx(i) == 6
      ylim([0,1])
  elseif var_idx(i) == 3
      ylim([0,0.5])
  elseif var_idx(i) == 4
      ylim([0,850])
  elseif var_idx(i) == 1
      ylim([0,100])
  elseif var_idx(i) == 2
      ylim([0,1])
  end
  xticklabels(["" ""])
  
end
%sgtitle("\textbf{without floe size}",'FontSize',15)
exportgraphics(f,strcat('FSD_comparison_',sprintf('%g',num_clusters),'_cover_15_percent_clusters.pdf'),'ContentType','vector')
label_vec = label_vec_temp;

%% Distributions between classifications
count = 1;
nBins = 11;
close all
f = figure('Position',[0, 0, 10, 10]);
for j = [1,4]
     % With floe size

    X_temp = [];
    X_temp = Xnan;
    average_stats = [];
    binEdges = linspace(min(X_temp(:,j)),max(X_temp(:,j)),nBins+1);
   
    X_temp = [X_temp, kmeans_idx];
    
    index_class_data = X_temp(:,end) == 3;
    class_data = X_temp(index_class_data,1:end-1);
    average_stats(1,:) = histcounts(class_data(:,j),binEdges,'Normalization','probability');

    % Without floe size
    X_temp = [];
    X_temp = Xnan;
    X_temp = [X_temp, kmeans_idx_without_floe];
    index_class_data = X_temp(:,end) == 3;
    class_data = X_temp(index_class_data,1:end-1);
    average_stats(2,:) = histcounts(class_data(:,j),binEdges,'Normalization','probability');

    subplot(2,1,count)
    for i = 1:2
      Xtemp(:,i) = Xnan(idx,i);
      histogram('BinEdges',binEdges,'BinCounts',average_stats(i,:),'FaceAlpha',.7,'FaceColor',Cmap_comp(i,:))
      grid on
      hold on
      %hist(Xnan(idx,i),20)
      %xticks(1:length(average_stats(:,i)))
      %title(label_vec{j})
      label_hist(i,:) = num2str(i);
      
    end
    ylabel('Probability')
    ylim([0,1])
    xlabel(ylabs{j})
    %legend(label_hist,'Location','bestoutside')
    count = count + 1;
end
exportgraphics(f,'MIZcompdistribution.pdf','ContentType','vector')
%%
line_widths = 40;
close all
f = figure('Position',[0, 0, 60, 10]);
plot(1:5,1:5,'color',Cmap_comp(1,:),'linewidth',line_widths)
hold on
plot(1:5,1:5,'color',Cmap_comp(2,:),'linewidth',line_widths)
legend("MIZ with floe size","MIZ without floe size",'Location','bestoutside','Orientation','vertical','FontSize',20)
exportgraphics(f,strcat('legend_MIZ.pdf'),'ContentType','vector')

%% Distribution of each class

close all

X_temp = [];
[~,~,dep] = size(Xnan);
X_temp = Xnan;
X_temp = [X_temp(:,1:end-2), kmeans_idx];
nBins = 10;
C1 = linspecer(num_clusters);
%Cmap = C1([3,2,1],:);
Cmap = C1([3,2,1],:);
conFigure(12,3)
clear label_hist
f = figure%('Position',[0, 0, 30, 10]);
count = 1;
for j = [1,4]
    average_stats = [];
    binEdges = linspace(min(X_temp(:,j)),max(X_temp(:,j)),nBins+1);
    for i = 1:num_clusters
        index_class_data = X_temp(:,end) == i;
        class_data = X_temp(index_class_data,1:end-1);
        average_stats(i,:) = histcounts(class_data(:,j),binEdges,'Normalization','probability');
    end
    %subplot(2,ceil(length(label_vec)/2),j)
    subplot(1,2,count)
    for i = num_clusters:-1:1
      Xtemp(:,i) = Xnan(idx,i);
      histogram('BinEdges',binEdges,'BinCounts',average_stats(i,:),'FaceAlpha',.7,'FaceColor',Cmap(i,:))
      hold on
      %hist(Xnan(idx,i),20)
      %xticks(1:length(average_stats(:,i)))
      %title(label_vec{j})
      label_hist(i,:) = num2str(i);
      
    end
    ylabel('Probability')
    ylim([0,1])
    xlabel(ylabs{j})
    %legend(label_hist,'Location','bestoutside')
    count = count + 1;
end

exportgraphics(f,'distributionClassVariables.pdf','ContentType','vector')

%%
clear Xtemp
idx = Xnan(:,1)>0.01;


close all
conFigure(12,10)
f = figure('Position',[0, 0, 60, 10]);
for i = 1:length(label_vec)
  subplot(2,ceil(length(label_vec)/2),i)
  Xtemp(:,i) = Xnan(idx,i);
  hist(Xtemp(:,i),20)
  %hist(Xnan(idx,i),20)
  %xticks(1:length(average_stats(:,i)))
  title(label_vec{i})
end


%%
X_new = Xnan;
index_lat = X_new(:,end-1) == X_new(1,end-1);
index_lon = X_new(:,end) == X_new(1,end);
index_both = index_lat.*index_lon;
sum(index_both)

index_lat = X_new(:,end-1) == X_new(200,end-1);
index_lon = X_new(:,end) == X_new(200,end);
index_both = index_lat.*index_lon;
sum(index_both)
clear index_lat index_lon index_both
%% Worldmap
close all
SIC = 0.15;
pram.ice_edge_color = 0.7*[0.4660 0.6740 0.1880];
line_width = 1;

addpath functions
[lat,lon,~,ulat,ulon] = grid_read('om2');
clear X_map MIZ_width
% Make worldmap with colour matching
hdd = "on";

if hdd == "on"
    uarea = data_format(filenames(1,:),"uarea");
end

[num_files,~] = size(filenames);
sector = "SH";
coords = sector_coords(sector);
font_size = 7;
plot_type = "kmeans";
plotting = "on";
conFigure(11)
X_new = Xnan;

[len,wid] = size(lat);
file_number = 1705;%1705%1826;%1705;%366;%609;



row_file = [0,cumsum(row_idx)];
row_vec = row_file(file_number)+1:row_file(file_number+1);
%temp_idx = kmeans_idx_without_floe;
temp_idx = kmeans_idx;
file_idx = temp_idx(row_vec);
file_idx = reshape(file_idx,length(file_idx),1);
if hdd == "on"
    [aice, sector_mask] = data_format_sector(filenames(file_number,:),'aice',sector);
end
X_map = [file_idx, X_new(row_vec,end-1:end)];
k_means = NaN.*ones(len,wid);
for i = 1:length(file_idx)
    [lon_pos,lat_pos,~] = near2(lon,lat,X_map(i,3),X_map(i,2));
    k_means(lon_pos,lat_pos) = file_idx(i); 
    k_means(~sector_mask) = NaN;
end
if hdd == "on"
    ice_mask =  aice > 0.15;
     [lat_ice_edge, lon_ice_edge, edge] = find_ice_edge(aice,0.15,sector,lat,lon);
end
if plotting == "on"
    f = figure;
    set(gcf,'Visible', 'on')
    w = worldmap('world');
        axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
        setm(w, 'Origin', [-90 0 0]);
        setm(w, 'maplatlimit', [-90,-53]);
        setm(w, 'grid', 'off');
        setm(w, 'frame', 'off');
        %setm(w, 'MLineLimit',[-60,-75]);
        %setm(w,'MLineException',[-90 0 90 180])
        setm(w, "FontColor",[0.5, 0.5, 0.5])
        setm(w, 'labelrotation', 'on')
        setm(w, 'meridianlabel', 'on','FontSize',font_size)
        setm(w, 'parallellabel', 'off','FontSize',font_size)
        setm(w, 'mlabellocation', 60);
        setm(w, 'plabellocation', 5);
        pcolorm(lat,lon,k_means,'FaceAlpha',0.99)
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(w, land, 'FaceColor', [0.5 0.5 0.5],'FaceAlpha',.5)
        %colorbar; %cmocean('deep');
        %colormap(Cmap)
        cb = colorbar;
        title(dirdates(file_number,:),'Color','black','FontSize',font_size+3)
        %plotm(lat_ice_edge(1:end-2),lon_ice_edge(1:end-2),'color','#7A7A7A','LineWidth',0.5,'LineStyle','-')
        %textm(lat_ice_edge(190)+2, lon_ice_edge(190),'15\% ice edge', 'FontSize', font_size, 'Color','#7A7A7A')
        scalebar('length',1000,...
            'units','km',...
            'color','k','location','sw',...
            'fontangle','italic','FontSize',font_size)
        if plot_type == "pc1"
            caxis([0,1])
        elseif plot_type == "kmeans"
            %cb = colorbar; cmocean('deep',num_clusters)
            %cb.TickLabels = region_label; 
            cb.Ticks = 1:num_clusters;
            %cb.Location = 'southoutside';%cb.Location = 'northoutside';
            cb.Location = 'eastoutside';
            caxis([0.5,num_clusters+0.5])
            cb.AxisLocation = 'out';
            cb.FontSize = font_size+3;
        end
end

exportgraphics(f,'3cluster.pdf','ContentType','vector')
clear X_map X_temp sb
%% MOVIES
clear area_region MIZ_width
SIC = 0.15;
pram.ice_edge_color = 0.7*[0.4660 0.6740 0.1880];
line_width = 1;

addpath functions
[lat,lon,~,ulat,ulon] = grid_read('om2');
clear X_map MIZ_width
% Make worldmap with colour matching
uarea = data_format(filenames(1,:),"uarea");

[num_files,~] = size(filenames);
sector = "SH";
coords = sector_coords(sector);
font_size = 5;
plot_type = "kmeans";
plotting = "off";
conFigure(11)
X_new = Xnan(:,:);
[len,wid] = size(lat);


file_number_vec = 1826-365:10:1826 ;
for file_number = file_number_vec 
    row_file = [0,cumsum(row_idx)];
    row_vec = row_file(file_number)+1:row_file(file_number+1);
    temp_idx = kmeans_idx;
    file_idx = temp_idx(row_vec);
    file_idx = reshape(file_idx,length(file_idx),1);

    [aice, sector_mask] = data_format_sector(filenames(file_number,:),'aice',sector);
    %[lat_ice_edge, lon_ice_edge, edge] = find_ice_edge(aice,SIC,sector,lat,lon);
    X_map = [file_idx, X_new(row_vec,end-1:end)];%[file_idx, X_new(:,end-1:end)];%
    k_means = NaN.*ones(len,wid);
    for i = 1:length(file_idx)
        [lon_pos,lat_pos,~] = near2(lon,lat,X_map(i,3),X_map(i,2));
        k_means(lon_pos,lat_pos) = file_idx(i); 
        k_means(~sector_mask) = NaN;
    end
    
    ice_mask =  aice > 0.15;
   % k_means(~ice_mask) = NaN;
    if plotting == "on"
        f = figure;
        set(gcf,'Visible', 'off')
        w = worldmap('world');
            %axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
            %setm(w, 'Origin', [-90 0 0]);
            %setm(w, 'maplatlimit', [-90,-50]);
            
            %setm(w, 'maplonlimit', [coords(1,2),coords(3,2)]);
            if sector == "EA"
                axesm eqdcylin; 
                setm(w, 'Origin', [0 28 0]); 
                setm(w, 'maplatlimit', [-75,-50]); setm(w, 'maplonlimit', [1,150]); 
                setm(w, 'meridianlabel', 'on'); setm(w, 'parallellabel', 'on'); 
                setm(w, 'mlabellocation', 30); setm(w, 'plabellocation', 10); 
                setm(w, 'mlabelparallel', 'south','FontColor','black','FontSize',3);
                setm(w, 'grid', 'off');
                setm(w, 'labelrotation', 'on')
            else
                axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
                setm(w, 'Origin', [-90 0 0]);
                setm(w, 'maplatlimit', [-90,-50]);
                setm(w, 'grid', 'off');
                setm(w, 'labelrotation', 'on')
            end
            pcolorm(lat,lon,k_means)
            land = shaperead('landareas', 'UseGeoCoords', true);
            geoshow(w, land, 'FaceColor', [0.5 0.5 0.5])
            %colorbar; %cmocean('deep');
            colormap(Cmap)
            cb = colorbar;
            title(dirdates(file_number,:),'Color','black','FontSize',font_size+5)
            %plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)

            %cb = colorbar; cmocean('deep',num_clusters)
            %cb.TickLabels = region_label; 
            cb.Ticks = 1:num_clusters;
            %cb.Location = 'southoutside';%cb.Location = 'northoutside';
            cb.Location = 'eastoutside';
            caxis([1,num_clusters])
            cb.AxisLocation = 'in';

    end
    % Calculate the width of the MIZ
    % CHECK WHAT number MIZ is !!!!!!!!
    [MIZ_width(file_number,:), miz_class, MIZ_zone] = calculate_miz_width(convertStringsToChars(filenames(file_number,:)),sector,k_means,3);
    

    % Calculate the perimeter and area of the zones
    for rgn_number = 1:num_clusters
        region_data = k_means == rgn_number;
        [len, wid] = size(k_means);
        perimeter_region(rgn_number,file_number) = find_perimeter(region_data,len,wid);

        % Area
        area_region(rgn_number, file_number) = sum(sum(uarea.*region_data));
    end
    if plotting == "on"
        if file_number == file_number_vec(1)
            gif('kmeans_SH_SA_05_stndzd_3class2014.gif','DelayTime',0.5,'resolution',500,'overwrite',true)
        else
            gif
        end
    end
    k_means_map(:,:,file_number) = k_means;
    % Clear the excess
    clear aice ice_mask k_means file_idx X_map
end
%%
clear temp
file_number_vec = 1826-365:1826;
for file_number = file_number_vec
    aice = data_format_sector(filenames(file_number,:),'aice',sector);
    mask = aice < 0.8;
    temp = aice > 0.15;
    SIC15 = temp.*mask;
    row_file = [0,cumsum(row_idx)];
    row_vec = row_file(file_number)+1:row_file(file_number+1);
    temp_idx = kmeans_idx_without_floe;
    file_idx = temp_idx(row_vec);
    file_idx = reshape(file_idx,length(file_idx),1);

    [aice, sector_mask] = data_format_sector(filenames(file_number,:),'aice',sector);
    X_map = [file_idx, X_new(row_vec,end-1:end)];
    k_means_without = NaN.*ones(len,wid);
    for i = 1:length(file_idx)
        [lon_pos,lat_pos,~] = near2(lon,lat,X_map(i,3),X_map(i,2));
        k_means_without(lon_pos,lat_pos) = file_idx(i); 
        k_means_without(~sector_mask) = NaN;
    end
     [MIZ_without_width(file_number,:), miz_class, MIZ_zone] = calculate_miz_width(convertStringsToChars(filenames(file_number,:)),sector,k_means_without,1);
     [SIC_width(file_number,:), miz_class, MIZ_zone] = calculate_miz_width(convertStringsToChars(filenames(file_number,:)),sector,SIC15,1);
     [MIZ_width(file_number,:), miz_class, MIZ_zone] = calculate_miz_width(convertStringsToChars(filenames(file_number,:)),sector,k_means_map(:,:,file_number),3);
    

    
end
%%

BrouwerData = [49.26108374384236 183.2512315270936 242.36453201970443 74.8768472906404];
BrouwerDates = datetime(['2019-02-15'; '2019-05-15'; '2019-09-15'; '2019-12-15']);
close all
   Colors = ['#78c679',"#fbb4b9",'#fecc5c','#d95f0e','#41b6c4'];
file_number_vec = 1826-365:10:1826;
%file_number_vec =  file_number_vec(1:end-1);
f = figure('Position',[0, 0, 30, 5]);
miz_data_vec = mean(MIZ_width(file_number_vec,:)');
miz_without_vec = mean(MIZ_without_width(file_number_vec,:)');

SIC_width_vec = mean(SIC_width(file_number_vec,:)');
%%
conFigure(11,2) 
%axes('NextPlot','replacechildren', 'ColorOrder',Colors());
plot(datetime(k_means_dirdates(file_number_vec,:)),miz_without_vec,'LineWidth',3,'LineStyle','-','Color','#8dd3c7') 
hold on
%hold on
plot(datetime(k_means_dirdates(file_number_vec,:)),SIC_width_vec,'LineWidth',3,'LineStyle','-','Color','#fdcdac') 
plot(datetime(k_means_dirdates(file_number_vec,:)),miz_data_vec(:,:),'LineWidth',3,'LineStyle','-','Color',Cmap(3,:)) 
plot(BrouwerDates,BrouwerData,'pentagram', 'MarkerFaceColor','yellow', 'MarkerSize',15) 
grid on
xlim([datetime(k_means_dirdates(file_number_vec(1+1),:)) datetime(k_means_dirdates(file_number_vec(end),:))])


ylim([0,500])
ylabel('MIZ width [km]')
legend('MIZ without floe size','15--80$\%$ SIC','MIZ with floe size','FIRF Exponential $^{[3]}$','Location','bestoutside')


% ylim(1.2*[-10^4,10^4])
        
       %    yline(0,'--')
        %title(varWant)
 exportgraphics(f,'LEGENDmiz_width_comp_2019.pdf','ContentType','vector')

%%
k_means_dirdates = dirdates(file_number_vec,:);
kmean_map.data = k_means_map;
kmean_map.dates = k_means_dirdates;
kmean_map.area_region = area_region;
kmean_map.MIZ_width = MIZ_width;
save(strcat('kmeans_map_',sprintf('%g',num_clusters),'_cover_clusters_2015_2019.mat'),'kmean_map','-v7.3')
%%
load('kmeans_map_3_clusters_2013_2014.mat')
%dirdates(file_number_vec,:) = k_means_dirdates;
k_means_map = kmean_map.data;
k_means_dirdates = kmean_map.dates;

%% Look at what is going on in each of these classes
% Pick a date
day = 365;
% Find what file this corresponds to
date_day = kmean_map.dates(day,:);
temp = date_day == filenames(:,end-12:end-3);
date_idx = find(sum(temp')==length(date_day));
clear temp

% Read in the data we want
varWant = 'wave_sig_ht';
data = data_format_sector(filenames(date_idx,:),varWant,sector);

idx_data = data < 0.5;
data(idx_data) = NaN;

data_temp = kmean_map.data(:,:,date_idx);
data_temp(idx_data) = NaN;
clear idx_data
close all
f = figure;
    set(gcf,'Visible', 'on')
    w = worldmap('world');
        axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
        setm(w, 'Origin', [-90 0 0]);
        setm(w, 'maplatlimit', [-90,-55]);
        setm(w, 'grid', 'on');
        setm(w, 'labelrotation', 'on')
        setm(w, 'meridianlabel', 'on','FontSize',font_size)
        setm(w, 'parallellabel', 'on','FontSize',font_size)
        setm(w, 'mlabellocation', 20);
        setm(w, 'plabellocation', 10);
        %pcolorm(lat,lon,data,'FaceAlpha',0.9)
        %[c1,h] = contourm(lat,lon,data,'LineColor','r','LevelStep',1,'LineWidth',1.2,'ShowText','on');
        %hold on
        pcolorm(lat,lon,kmean_map.data(:,:,date_idx),'FaceAlpha',0.4)
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(w, land, 'FaceColor', [0.5 0.5 0.5],'FaceAlpha',.5)
        pcolorm(lat,lon,data_temp,'FaceAlpha',1)
        
        %colorbar; %cmocean('deep');
        colormap(Cmap)
        cb = colorbar;
        title(date_day,'Color','black','FontSize',font_size+5)
        %plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
        
        cb.Ticks = 1:num_clusters;
        %cb.Location = 'southoutside';%cb.Location = 'northoutside';
        cb.Location = 'eastoutside';
        caxis([1,num_clusters])
        cb.AxisLocation = 'in';
       
%% Average statistics of variables which don't go into the classifier
close all
clear average_stats temp data ts_average_stats
% Pick a date
nBins = 21;
average_stats = zeros(num_clusters,nBins);
binEdges = linspace(-0.1,0.1,nBins+1);
normalise_stats = "off"; % on, off
time_series = "fsd"; % on, fsd, itd
sector = "SH";
time_series_type = "accum"%"accum"; % accum or mean
NFSD = ncread(filenames(1,:),"NFSD");
NCAT = ncread(filenames(1,:),"NCAT");
uarea = data_format_sector(filenames(1,:),'uarea',sector);
[floe_binwidth, floe_rad_l, floe_rad_h, floe_area_binwidth] = cice_parameters(NFSD);

% 2017
% 121 = May 1st
% 273 = Sep 30st
% 335 = Dec 1st
% 425 = March 1st
% Summer = October to March 1370:1551
% Winter = April to September 1552:1734
day_vec = 1370:10:1551; %1826-365:5:1826; %1552:10:1734;%121:273;


varVec = ['fsd_latg'; 'fsd_latm'; 'fsd_newi'; 'fsd_weld'; 'fsd_wave'];
for fsdCount = 1:5
    if time_series == "fsd" 
        average_stats = zeros(num_clusters,length(NFSD));
    elseif time_series == "itd" 
        average_stats = zeros(num_clusters,length(NCAT));
    elseif time_series == "on"
        average_stats = zeros(num_clusters,length(day_vec));
    elseif time_series == "miz_area"
        average_stats = zeros(2,length(day_vec));
    end
    varWant = varVec(fsdCount,:)
for day = day_vec
    % Find what file this corresponds to
    date_day = kmean_map.dates(day,:);
    temp = date_day == filenames(:,end-12:end-3);
    date_idx = find(sum(temp')==length(date_day));
    clear temp
 
    % Read in the data we want
    %varWant = 'pancake_welding';
    %varWant = 'pancake_formation';
    
    if length(varWant) > 4
        if varWant(1:5) == "dafsd"
            % Change in FSD
            temp = data_format_sector(filenames(date_idx,:),varWant,sector);
            [len,wid,~] = size(temp);
            for i = 1:len
                for j = 1:wid
                    data(i,j) = sum(squeeze(temp(i,j,:)).*floe_binwidth');
                end
            end
        elseif varWant == "pancake_formation"
            % Ice formation in the smallest floe size category
            temp = data_format_sector(filenames(date_idx,:),'dafsd_newi',sector);
            [len,wid,~] = size(temp);
            for i = 1:len
                for j = 1:wid
                    data(i,j) = squeeze(temp(i,j,1)).*floe_binwidth(1);
                end
            end
        elseif varWant == "nilas_formation"
            % Ice formation in the largest floe size category
            temp = data_format_sector(filenames(date_idx,:),'dafsd_newi',sector);
            [len,wid,~] = size(temp);
            for i = 1:len
                for j = 1:wid
                    data(i,j) = squeeze(temp(i,j,end)).*floe_binwidth(end);
                end
            end
        elseif varWant == "pancake_welding"
            % Ice formation in the smallest floe size category
            temp = data_format_sector(filenames(date_idx,:),'dafsd_weld',sector);
            [len,wid,~] = size(temp);
            for i = 1:len
                for j = 1:wid
                    data(i,j) = squeeze(temp(i,j,1)).*floe_binwidth(1);
                end
            end
        elseif varWant == "pancake_latm"
            % Ice formation in the smallest floe size category
            temp = data_format_sector(filenames(date_idx,:),'dafsd_latm',sector);
            [len,wid,~] = size(temp);
            for i = 1:len
                for j = 1:wid
                    data(i,j) = squeeze(temp(i,j,1)).*floe_binwidth(1);
                end
            end
        elseif varWant == "pancake"
            % SIC in FSTD(1,1)
            temp = data_format_sector(filenames(date_idx,:),'afsdn',sector);
            [len,wid,~] = size(temp);
            for i = 1:len
                for j = 1:wid
                    data(i,j) = squeeze(temp(i,j,1,1)).*floe_binwidth(1);
                end
            end
        elseif varWant == "brash"
            % Ice formation in the smallest floe size category
            temp = data_format_sector(filenames(date_idx,:),'dafsd_wave',sector);
            [len,wid,~] = size(temp);
            for i = 1:len
                for j = 1:wid
                    data(i,j) = squeeze(temp(i,j,1)).*floe_binwidth(1);
                end
            end
        elseif varWant == "afsdn1"
            % Ice formation in the smallest floe size category
            temp = data_format_sector(filenames(date_idx,:),'afsdn',sector);
            aice = data_format_sector(filenames(date_idx,:),'aice',sector);
            [len,wid,~] = size(temp);
            for i = 1:len
                for j = 1:wid
                    data(i,j) = squeeze(temp(i,j,1,1)).*floe_binwidth(1)./aice(i,j);
                   % data(i,j) = squeeze(temp(i,j,1))./aice(i,j);
                end
            end
        elseif varWant == "pancake_proportion"
            % Proportion of FSTD(1,1)
            aice = data_format_sector(filenames(date_idx,:),'aice',sector);
            temp = data_format_sector(filenames(date_idx,:),'afsdn',sector);
            [len,wid,~] = size(temp);
            for i = 1:len
                for j = 1:wid
                    data(i,j) = squeeze(temp(i,j,1,1)).*floe_binwidth(1)./aice(i,j);
                end
            end
        elseif varWant == "thick_pancake"
            % Proportion of FSTD(1,1)
            aice = data_format_sector(filenames(date_idx,:),'aice',sector);
            temp = data_format_sector(filenames(date_idx,:),'afsdn',sector);
            [len,wid,~] = size(temp);
            for i = 1:len
                for j = 1:wid
                    data(i,j) = 0;
                    for nc = 2:length(NCAT)
                        data(i,j) = data(i,j) + squeeze(temp(i,j,1,nc)).*floe_binwidth(1)./aice(i,j);
                    end
                end
            end
        elseif varWant(1:4) == "fsd_"
            data = data_format_sector(filenames(date_idx,:),strcat("dafsd",varWant(4:end)),sector);

        elseif varWant == "MIZ_width"
            class_idx = 3;
            %for class_idx = 1:3
                [MIZ_width(day,:), miz_class, MIZ_zone] = calculate_miz_width(convertStringsToChars(filenames(day,:)),sector,kmean_map.data(:,:,date_idx),class_idx);
                %ts_average_stats(class_idx,day) = mean(MIZ_width');
            %end

        else
            data = data_format_sector(filenames(date_idx,:),varWant,sector);
        end
        else
        data = data_format_sector(filenames(date_idx,:),varWant,sector);
    end
    clear temp len wid

    if time_series == "on"
        if time_series_type == "mean"
         for n = 1:num_clusters
            mask = kmean_map.data(:,:,date_idx) == n;
            ts_average_stats(n,day) = mean(reshape(data(mask),numel(data(mask)),[]));
         end
        elseif time_series_type == "accum"
            for n = 1:num_clusters
            mask = kmean_map.data(:,:,date_idx) == n;
            ts_average_stats(n,day) = sum(reshape(data(mask),numel(data(mask)),[]));
            end
        elseif varWant == "MIZ_width"
            ts_average_stats(day) = mean(MIZ_width(day,:)');
        end
    elseif time_series == "fsd"
        if time_series_type == "mean"
            for n = 1:num_clusters
                mask = kmean_map.data(:,:,date_idx) == n;
                for nf = 1:length(NFSD)
                    temp = squeeze(data(:,:,nf));
                    average_stats(n,nf) = average_stats(n,nf)+sum(reshape(temp(mask),numel(temp(mask)),[])).*floe_binwidth(nf);
                end
            end
        elseif time_series_type == "accum"
            for n = 1:num_clusters
                mask = kmean_map.data(:,:,date_idx) == n;
                for nf = 1:length(NFSD)
                    temp = squeeze(data(:,:,nf));
                    average_stats(n,nf) = average_stats(n,nf)+sum(reshape(temp(mask),numel(temp(mask)),[])).*floe_binwidth(nf);
                end
            end
        end
    elseif time_series == "itd"
        for n = 1:num_clusters
            mask = kmean_map.data(:,:,date_idx) == n;
            for nc = 1:length(NCAT)
                temp = squeeze(data(:,:,nc));
                average_stats(n,nc) = average_stats(n,nc)+sum(reshape(temp(mask),numel(temp(mask)),[]));
            end
        end
    elseif time_series == "miz_area"
            n = 3; % Pick the class referring to MIZ
            kmask = kmean_map.data(:,:,date_idx) == n;
            itd_mask = data > 0.7;
            mask = logical(kmask.*itd_mask);
            average_stats(1,day) = sum(reshape(uarea(mask),numel(uarea(mask)),[]));
            mask = logical(kmask.*~itd_mask);
            average_stats(2,day) = sum(reshape(uarea(mask),numel(uarea(mask)),[]));
            average_stats(:,day) = average_stats(:,day)./sum(average_stats(:,day));
    else % Histogram independent of time
        for n = 1:num_clusters
            mask = kmean_map.data(:,:,date_idx) == n;
            average_stats(n,:) = average_stats(n,:) + histcounts(reshape(data(mask),numel(data(mask)),[]),binEdges);
        end
    end

end % day

% Normalise average stats
if normalise_stats == "on"
    for n = 1:num_clusters
     average_stats(n,:) = average_stats(n,:)./sum(average_stats(n,:));
    end
end
dafsd_MIZ(fsdCount,:) = average_stats(3,:);

end % FSD COUNT
%

%%
close all
f = figure('Position',[0, 0, 15, 7.5]);
if time_series == "on" 
        conFigure(11,2)
        if varWant == "MIZ_width"
            axes('NextPlot','replacechildren', 'ColorOrder',Cmap(3,:));
        else
            axes('NextPlot','replacechildren', 'ColorOrder',Cmap);
        end
        
        plot(datetime(k_means_dirdates(day_vec,:)),ts_average_stats(:,day_vec),'LineWidth',3,'LineStyle','-') 
           % ylim(1.2*[-10^4,10^4])
            yline(0,'--')
        %title(varWant)
        fig_name = strcat('ts_',varWant,'_',date_day,'.pdf');
        %ylabel('Cumulative change in FSD [m$^{-1}$]')
        exportgraphics(f,fig_name,'ContentType','vector')
elseif time_series == "fsd"
    % FSD plots');
    conFigure(13,2)
    axes('NextPlot','replacechildren', 'ColorOrder',Cmap);
    set(gcf,'Visible', 'on')
    if varWant(1:4) == "fsd_"
        clear y_lab
        logdata = sign(average_stats).*log10(abs(average_stats));
        %plot(NFSD,average_stats(:,:),'LineWidth',3,'LineStyle','-','Marker','s','MarkerSize',5)
        %bar(1:12,logdata(:,:),'stacked')
        %ylabel('Change in FSD [m$^{-1}$]')
        ylabel('Cumulative change in FSD [m$^{-1}$]')
        xlabel('Floe size [m]')
        %set(gca,'YScale','log')
        grid on
        xticklabels(string(num2str(round(NFSD,0))))
        %symlog(gca,'xy',-1.7)
        %ytick_log = round(min(ytick_log),0):5:round(max(yticks),0);
        %for i = 1:length(ytick_log)
        %    y_lab(i) = string(num2str(round(sign(ytick_log(i))*exp(sign(ytick_log(i))*ytick_log(i)),1)));%{'1', '10', '1000'};
        %end
        %yticks = ytick_log;
        %yticklabels(y_lab)


        %for i = 1:num_clusters
        % bar(1:12,average_stats(i,:),'FaceAlpha',.7,'FaceColor',Cmap(i,:))
        % hold on
        %end
         %y_vec = [-10^1  10^1];
        %yticks(y_vec)
        %y_vec = sign(y_vec).*log(abs(y_vec));
        %yticks(log(NFSD))
        %clear y_label
        %count = 1;
        %for i = y_vec
        %    y_label(count) = string(num2str(i));
        %    count = count +1;
        %end
        
        %yticklabels('auto')
        %yticklabels(["10^-1","10^1"])
        %legend(region_label,'location','southwest')
        %ytick_vec = [-10^3 -10^2 -10 1 10 100 10^3];
        %ytick_log =  sign(ytick_vec).*log(abs(ytick_vec));
        %yticks(ytick_log)
        %y_lab = {"-10^{3}", "-100", "-10", "0", "10", "100", "10^{3}"};
        %for i = 1:length(ytick_log)
        %    y_lab(i) = string(num2str(round(sign(ytick_log(i))*exp(sign(ytick_log(i))*ytick_log(i)),1)));%{'1', '10', '1000'};
        %end
        %yticklabels(y_lab)
        %temp = xticks
        %xtick_vec = [1 10 100 1000];
        %xtick_log = sign(xtick_vec).*log(abs(xtick_vec));
        %xticks(xtick_log)
        %for i = 1:length(temp)
            %x_lab(i) = string(num2str(exp(xtick_log(i))));%{'1', '10', '1000'};
        %end
        %xlim(log10([5,50000]))
        %xticklabels(x_lab)
        yline(0,'--')
    else
        loglog(NFSD,average_stats(:,:),'LineWidth',3,'LineStyle','-','Marker','s','MarkerSize',5)
        ylabel('Sea ice concentration [$\%$]')
        xlabel('Floe size [m]')
    end
    
    

    %ylim([10^-5,10^2])
    fig_name = strcat(varWant,'_',date_day,'.pdf');
    exportgraphics(f,fig_name,'ContentType','vector')
elseif time_series == "itd"
    % FSD plots
    conFigure(11,1)
    axes('NextPlot','replacechildren', 'ColorOrder',Cmap);
    set(gcf,'Visible', 'on')
    plot([0;NCAT(1:end-1)],average_stats(:,:).*100,'LineWidth',3,'LineStyle','-','Marker','s','MarkerSize',5)
    %legend(region_label,'location','southwest')
    ylabel('Sea ice concentration [$\%$]')
    xlabel('Ice thickness [m]')
    ylim([10^-5,10^2])
    fig_name = strcat(varWant,'_',date_day,'.pdf');
    exportgraphics(f,fig_name,'ContentType','vector')
elseif time_series == "miz_area" 
    conFigure(11,2)
    axes('NextPlot','replacechildren', 'ColorOrder',Cmap);
    plot(datetime(k_means_dirdates(day_vec,:)),average_stats(:,day_vec),'LineWidth',3,'LineStyle','-') 
    ylabel('Area [km$^2$]')
    title(varWant)
    fig_name = strcat('ts_',varWant,'_',date_day,'.pdf');
    exportgraphics(f,fig_name,'ContentType','vector')
elseif time_series == "bar"
    conFigure(11,2)
    axes('NextPlot','replacechildren', 'ColorOrder',Cmap);
    bar([sum(average_stats'); sum(average_stats')],'stacked')

else % Histogram independent of time
    conFigure(11,2)
    for i = num_clusters:-1:1
        histogram('BinEdges',binEdges,'BinCounts',average_stats(i,:),'FaceAlpha',.7,'FaceColor',Cmap(i,:))
        hold on
        title(varWant)
        label_hist(i,:) = num2str(i);
        if normalise_stats == "on"
            ylabel('Probability')
        else
            ylabel('Frequency')
        end
    end
    legend(label_hist,'Location','northeast')
    
    if normalise_stats == "on"
        fig_name = strcat('dist_',varWant,'_',date_day,'_norm','.pdf');
    else
        fig_name = strcat('dist_',varWant,'_',date_day,'_freq','.pdf');
    end
    exportgraphics(f,fig_name,'ContentType','vector')
end
%% DAFSD MIZ

   %Colors = ["#d7191c" "#fdae61" '#ffffbf' '#a6d96a' '#1a9641'];
   Colors = ['#78c679',"#fbb4b9",'#fecc5c','#d95f0e','#41b6c4'];
  % Colors = ["#e41a1c",'#377eb8','#4daf4a','#984ea3','#ff7f00']
    close all
    f = figure('Position',[0, 0, 20, 10]);
    %axes('NextPlot','replacechildren', 'ColorOrder',Colors);
     for i = 1:width(dafsd_MIZ')
       plot(NFSD,dafsd_MIZ(i,:),'LineWidth',3,'LineStyle','-','Marker','s','MarkerSize',5,'Color',Colors(i))
       hold on
     end
   symlog() % no harm in letting symlog operate in z axis, too.
    ylim([-4.1,4.1])
   yTick = yticks;
   %x = sgn(y).C.(−1 + 10|y|)
   backTrans = sign(yTick).*(-1+10.^(abs(yTick)));
    for i = 1:length(yTick)
        yTickLab{i} = sprintf('%g0$^{%g}$',sign(backTrans(i)),log10(abs(backTrans(i)))); 
    end
    yTickLab{find(backTrans == 0)} = '0';
    %yticks(yTick)
    yticklabels(yTickLab)
    yline(0,'--')
  
    xTick = log10([1, 10, 100, 1000]);
    for i = 1:length(xTick)
        xTickLab{i} = sprintf('%g0$^{%g}$',sign(xTick(i)),log10(ceil(10^(abs(xTick(i)))))); 
    end
    xTickLab{find(xTick == 0)} = '0';
    xlim([0.4,3])
    xticks(xTick)
    xticklabels(xTickLab)
    ylabel('Change in FSD [dimless]')%%[m$^{-1}$]')
    xlabel('Floe size [m]')
    legend('Lateral growth', 'Lateral melt','New floes','Welding','Wave breakup','Location','south','Orientation','horizontal','NumColumns',3)
    exportgraphics(f,'dafsd_MIZ_summer_2018.pdf','ContentType','vector')

    %%
    % clear xTickLab  logdata
% dafsd_MIZ =[-78.0686   -5.8085   -5.8846   -2.4448   -2.8111  -0.3008   -0.0584   -0.0122   -0.0375   -0.0340  4.0217  -16.5656;...
%      0.2016    0.1211    0.0048    0.0025    0.0006 0.0003    0.0002    0.0001    0.0000    0.0000 0.0000    0.0007;...
%      339.9192   22.8999   18.9079   17.9654   18.4002 22.9754   21.4802   17.1830   13.8310   18.9008 14.2127  526.9201;...
%      -1.0477   -0.1036   -4.8148   -2.4014   -4.8859 -13.9847  -14.6517  -15.5492  -12.7272  -19.6298 -14.7161  330.2283;...
%      46.8000   -2.1910  -23.5654  -33.1271  -37.5959 -5.3921   -3.9992   -6.6194   -1.6589   -2.9203  -1.4078 -826.9496];
%    
%     
%    
%   
%   
% 
% Cols = linspecer(5);
% close all
% f = figure;
% conFigure(13,1)
%     axes('NextPlot','replacechildren', 'ColorOrder',Cols);
%     set(gcf,'Visible', 'on')
%     clear y_lab
%     logdata = sign(dafsd_MIZ).*log10(abs(dafsd_MIZ));
%     %plot(NFSD,dafsd_MIZ(:,:),'LineWidth',3,'LineStyle','-','Marker','s','MarkerSize',5)
%     plot(log10(NFSD),logdata(:,:),'LineWidth',3,'LineStyle','-','Marker','s','MarkerSize',5)
%     yline(0,'--')
%     ylabel('Change in FSD [m$^{-1}$]')
%     xlabel('Floe size [m]')
%     legend('Lateral growth', 'Lateral melt','New floes','Welding','Wave breakup','Location','Southwest')
%     xTick = log10([1,3, 10, 30, 100, 300, 1000]);
%     for i = 1:length(xTick)
%         xTickLab{i} = num2str(10.^(xTick(i)));
%     end
%     xticks(xTick)
%     xticklabels(xTickLab)
% 
% 
%     yTick = yticks;%log10([1,3, 10, 30, 100, 300, 1000]);
%     for i = 1:length(yTick)
%         yTickLab{i} = sprintf('%g0$^{%g}$',sign(yTick(i)),abs(yTick(i))); %num2str(sign(yTick(i))*10.^(abs(yTick(i))));
%     end
%     yTickLab{find(yTick == 0)} = '0';
%     yTickLab{find(yTick == 1)} = '10';
%     yTickLab{find(yTick == -1)} = '-10';
% 
%     yticks(yTick)
%     yticklabels(yTickLab)
%     %xticklabels({'0','\pi','2\pi','3\pi','4\pi','5\pi','6\pi'})
%     %symlog(gca,'xy',-1.7)
% 
%     %for i = 1:num_clusters
%     % bar(1:12,average_stats(i,:),'FaceAlpha',.7,'FaceColor',Cmap(i,:))
%     % hold on
%     %end
%      %y_vec = [-10^1  10^1];
%     %yticks(y_vec)
%     %y_vec = sign(y_vec).*log(abs(y_vec));
%     %yticks(log(NFSD))
%     %clear y_label
%     %count = 1;
%     %for i = y_vec
%     %    y_label(count) = string(num2str(i));
%     %    count = count +1;
%     %end



   %Cols = linspecer(5);
   %Cols = cmocean('phase',6);
   %Cols  = Cols(1:end,:);
    close all
 
    axes('NextPlot','replacechildren', 'ColorOrder',Cols);
    figure
    y = dafsd_MIZ;
    x = NFSD; 
    % Do log10 but keep sign
    xlog = sign(x).*log10(abs(x));
    % Just to get axis limits
    plot(xlog,y,'o')
    % Get limits
    lims = xlim;
    wdth = diff(lims);
    % Wrap negative data around to positive side
    xlog(xlog<0) = xlog(xlog<0) + wdth;

    % Do log10 but keep sign
    ylog = sign(y).*log10(abs(y));
    % Just to get axis limits
    plot(xlog,ylog,'o')
    % Get limits
    limsy = ylim;
    wdthy = diff(limsy);
    % Wrap negative data around to positive side
    ylog(ylog<0) = ylog(ylog<0) + wdthy;

    % Plot
    plot(xlog,ylog,'LineWidth',3,'LineStyle','-','Marker','s','MarkerSize',5)
    % Mess with ticks
    tck = get(gca,'XTick')';
    % Shift those that were wrapped from negative to positive (above) back 
    % to their original values
    tck(tck>lims(2)) = tck(tck>lims(2)) - wdth;
    % Convert to string, then remove any midpoint
    tcklbl = num2str(tck);
    tcklbl(tck==lims(2),:) = ' ';
    % Update tick labels
    set(gca,'XTickLabel',tcklbl)

    % Mess with ticks
    tck = get(gca,'YTick')';
    % Shift those that were wrapped from negative to positive (above) back 
    % to their original values
    tck(tck>limsy(2)) = tck(tck>limsy(2)) - wdthy;
    % Convert to string, then remove any midpoint
    tcklbl = num2str(tck);
    tcklbl(tck==limsy(2),:) = ' ';
    % Update tick labels
    set(gca,'YTickLabel',tcklbl)
    legend('Lateral growth', 'Lateral melt','New floes','Welding','Wave breakup','Location','Southwest')
  %%



    
%% SIC baseline
sector = "SH";
SIC_area = [];
file_number_vec = 1:n_files-1;
uarea = data_format_sector(filenames(1,:),'uarea',sector);
for file_number = file_number_vec 
    aice = data_format_sector(filenames(file_number,:),'aice',sector);
    mask = aice < 0.8;
    temp = aice > 0.15;
    mask = logical(temp.*mask);
    SIC_area(1,file_number) = sum(sum(uarea(mask)));

    mask = aice > 0.8;
    SIC_area(2,file_number) = sum(sum(uarea(mask)));
end
%% PLOT the areal comparison of SIC and sea ice class
% file_number_vec k_means_dirdates

% Average annually
annual_lim = [365 366 365 365 365];
idx_annual = [1,cumsum(annual_lim)];
temp_area = zeros(3,365);
for i = 1:5
    temp_area = temp_area + area_region(:,idx_annual(i):idx_annual(i)+364);
end
temp_area = temp_area./length(annual_lim);
%%
close all
conFigure(11,1.8)
f = figure;
axes('NextPlot','replacechildren', 'ColorOrder',Cmap);
plot(datetime(k_means_dirdates(1:365,:)),temp_area(:,:).*1e-6.*1e-6,'LineWidth',3,'LineStyle','-') % sum(area_region(:,file_number_vec).*1e-6.*1e-6)
ylabel('Area [$10^6$ km$^2$]')
legend('Young','Older','MIZ','Location','southoutside')
exportgraphics(f,'temp.pdf','ContentType','vector')
%%
hold on
plot(datetime(k_means_dirdates),SIC_area(1,:).*1e-6.*1e-6,'LineWidth',3,'LineStyle','-') % sum(area_region(:,file_number_vec).*1e-6.*1e-6)

ylabel('Area [million km$^2$]')
legend('MIZ','0.15$<$SIC$<$0.8','Location','northwest')
%%
exportgraphics(f,'MIZvs15-80dynamics.pdf','ContentType','vector')

f = figure;
%axes('NextPlot','replacechildren', 'ColorOrder',Cmap);
plot(datetime(k_means_dirdates),(area_region(2,file_number_vec)+area_region(3,file_number_vec)).*1e-6.*1e-6,'LineWidth',3,'LineStyle','-') % sum(area_region(:,file_number_vec).*1e-6.*1e-6)

hold on
plot(datetime(k_means_dirdates),SIC_area(2,:).*1e-6.*1e-6,'LineWidth',3,'LineStyle','-') % sum(area_region(:,file_number_vec).*1e-6.*1e-6)

ylabel('Area [million km$^2$]')
legend('Inner','SIC$>$0.8','Location','northwest')
%legend(["1","2","3","\<80","\>80"],'Location','southeast')

exportgraphics(f,'Innervs80dynamics.pdf','ContentType','vector')

%% Annual median
%year_idx = k_means_dirdates(:,end-4:end) == k_means_dirdates(1,end-4:end);


interannual_mean = [];
interannual_mean(1,:) = (area_region(1,1:365)+area_region(1,366:2*365))/2;
interannual_mean(2,:) = (SIC_area(1,1:365)+SIC_area(1,366:2*365))/2;
%
close all
f = figure;

plot(datetime(k_means_dirdates(1:365,:)),interannual_mean.*1e-6.*1e-6,'LineWidth',3,'LineStyle','-') 
ylabel('Area [million km$^2$]')
legend('MIZ','0.15$<$SIC$<$0.8','Location','northwest')
exportgraphics(f,'meanMIZvs15-80.pdf','ContentType','vector')









%%
% LINE COLORS
close all
N=6;
X = linspace(0,pi*3,1000);
Y = bsxfun(@(x,n)sin(x+2*n*pi/N), X.', 1:N);

subplot(2,2,1);
nexttile
C = linspecer(N);
axes('NextPlot','replacechildren', 'ColorOrder',C);
p = plot(date_label_vec,[perimeter_region([1,2,4],file_number_vec); sum(perimeter_region(:,file_number_vec))],'LineWidth',5);
%ylim([-1.1 1.1]);
%% Area time series
%file_number_vec = 1:361;
% load('a_p_sheet.mat');
% area_region = ap_sheet.area_region;
%file_number_vec = 1:n_files-1;
date_label = datetime(dirdates);
date_label_vec = date_label(file_number_vec);
close all
conFigure(11,2)
% f = figure;

% axes('NextPlot','replacechildren', 'ColorOrder',Cmap);
% plot(date_label_vec,[perimeter_region([1,2,3,4],file_number_vec); ],'LineWidth',3,'LineStyle','--');
% %set(p, {'color'}, {[cbar(1,:)];[cbar(2,:)];[cbar(4,:)];[0.5,0.1,0.1]});
% ylabel('Perimeter [block edges]')
%legend(region_label,'Location','southeast')

%exportgraphics(f,'wave_perimeter.pdf','ContentType','vector')

f = figure;
axes('NextPlot','replacechildren', 'ColorOrder',Cmap);
plot(date_label_vec,[area_region(:,file_number_vec).*1e-6.*1e-6],'LineWidth',3,'LineStyle','-') % sum(area_region(:,file_number_vec).*1e-6.*1e-6)
%set(p, {'color'}, {[cbar(1,:)];[cbar(2,:)];[cbar(4,:)];[0.5,0.1,0.1]});
ylabel('Area [million km$^2$]')
%legend(region_label,'Location','southeast')

exportgraphics(f,'4_SH_area_dynamics.pdf','ContentType','vector')
%%
file_number_vec = 1:364;
Cmap(1,:) = C1(2,:);
Cmap(2,:) = C1(1,:);
f = figure;
axes('NextPlot','replacechildren', 'ColorOrder',Cmap);
plot(date_label_vec,[area_region([4],file_number_vec).*1e-6.*1e-6; sum(area_region([1,2],file_number_vec)).*1e-6.*1e-6],'LineWidth',3,'LineStyle','-') % sum(area_region(:,file_number_vec).*1e-6.*1e-6)
%set(p, {'color'}, {[cbar(1,:)];[cbar(2,:)];[cbar(4,:)];[0.5,0.1,0.1]});
ylabel('Area [million km$^2$]')
%legend(region_label,'Location','southeast')

exportgraphics(f,'pack.pdf','ContentType','vector')

%%
subplot(2,1,2)
C = linspecer(N);
file_number_vec = 1:831;
date_label_vec = date_label(file_number_vec);
axes('NextPlot','replacechildren', 'ColorOrder',C);
scatter(date_label_vec,[area_region([1,2,3],file_number_vec).*1e-6.*1e-6; sum(area_region(:,file_number_vec).*1e-6.*1e-6)],'LineWidth',5)
ylabel('Area [million km$^2$]')
%legend(region_label,'Location','southeast')




% % Create a tile on the right column to get its position
% legend_label = region_label;
% legend_label{end+1} = 'Sum of regions';
% ax = subplot(3,1,1,'Visible','off');
% axPos = ax.Position;
% delete(ax)
% % Construct a Legend with the data from the sub-plots
% hL = legend(legend_label);
% % Move the legend to the position of the extra axes
% hL.Position(1:2) = axPos(1:2);




%%


NFSD = ncread(filenames(1,:),"NFSD");
NCAT = ncread(filenames(1,:),"NCAT");
f = figure;
set(gcf,'Visible', 'on')
loglog(NFSD,average_stats([1,2,4],14:end),'LineWidth',3,'LineStyle','-','Marker','s','MarkerSize',5)
legend(region_label{[1,2,4]},'location','southwest')
ylabel('Sea ice concentration [$\%$]')
xlabel('Floe size [m]')
exportgraphics(f,'floesizecompminus3.pdf','ContentType','vector')









%% Evolution of statistics in each class
% [waves; sheet; pancake]
SIC = 0.15;
pram.ice_edge_color = 0.7*[0.4660 0.6740 0.1880];
line_width = 1;
num_clusters = 4;
%dimension.pancake(1) = cleaned_data.dimension.pancake(1);
%dimension.sheet(1) = cleaned_data.dimension.sheet(1);

addpath functions
[lat,lon,~,ulat,ulon] = grid_read('om2');
clear X_map MIZ_width
% Make worldmap with colour matching
uarea = data_format(filenames.waves(1,:),"uarea");

[num_files,~] = size(filenames.waves);
sector = "EA";
coords = sector_coords(sector);
font_size = 5;
plot_type = "kmeans";
plotting = "off";
conFigure(11)
%X_standard_all(:,end-1:end) = Xnan.all(:,end-1:end);
%X_new = X_standard_all(dimension.waves(1)+dimension.sheet(1):dimension.waves(1)+dimension.sheet(1)+dimension.pancake(1),:);
X_new = Xnan.all(1:dimension.waves(1),:);
%X_new = Xnan.all(1:dimension.pancake(1),:);

[len,wid] = size(lat);
N = num_clusters;

C1 = linspecer(N);
Ctemp(1,:) = C1(4,:);
Ctemp(2,:) = C1(3,:);
Ctemp(3,:) = C1(1,:);
Ctemp(4,:) = C1(2,:);
Cmap = Ctemp([1,2,4,3],:);


file_number_vec = 1:30:365;
for file_number = file_number_vec
    if plot_type == "pc1"
        file_idx = score(row_file(file_number)+1:row_file(file_number+1),1);
        file_idx = file_idx + abs(min(file_idx));
        file_idx = file_idx./max(file_idx);
    elseif plot_type == "kmeans"
        index_lat = X_new(:,end-1) == X_new(30,end-1);
        index_lon = X_new(:,end) == X_new(30,end);
        index_both = index_lat.*index_lon;
        diff_index = diff(find(index_both)');
        %row_file = [0,cumsum(diff_index)];
        row_file = find(index_both);%[0,cumsum(row_idx)];
        row_vec = row_file(file_number)+1:row_file(file_number+1);
        %temp_idx = idx(dimension.waves(1)+dimension.sheet(1):dimension.waves(1)+dimension.sheet(1)+dimension.pancake(1)); % waves:waves+sheet
        temp_idx = idx(1:dimension.waves(1));
        file_idx = temp_idx(row_vec);
    end
    [aice, sector_mask] = data_format_sector(filenames.waves(file_number,:),'aice',sector);
    %[lat_ice_edge, lon_ice_edge, edge] = find_ice_edge(aice,SIC,sector,lat,lon);
    X_map = [file_idx, X_new(row_vec,end-1:end)];%[file_idx, X_new(:,end-1:end)];%
    k_means = NaN.*ones(len,wid);
    for i = 1:length(file_idx)
        [lon_pos,lat_pos,~] = near2(lon,lat,X_map(i,3),X_map(i,2));
        k_means(lon_pos,lat_pos) = file_idx(i); 
        k_means(~sector_mask) = NaN;
    end
    
    ice_mask =  aice > 0.15;
   % k_means(~ice_mask) = NaN;

    % Calculate the width of the MIZ
    % CHECK WHAT number MIZ is !!!!!!!!
    %[MIZ_width(file_number,:), miz_class, MIZ_zone] = calculate_miz_width(convertStringsToChars(filenames.waves(file_number,:)),sector,k_means,1);
    

    % Calculate evolution of each statistic
    for rgn_number = 1:num_clusters
        region_data = k_means == rgn_number;
        [len, wid] = size(k_means); % Size of grid
        wave_sig_ht_temp = data_format_sector(filenames.waves(file_number,:),'wave_sig_ht',sector);
        wave_sig_ht(rgn_number, file_number) = mean(mean(wave_sig_ht_temp.*region_data));
        ppd_temp = data_format_sector(filenames.waves(file_number,:),'peak_period',sector);
        ppd(rgn_number, file_number) = mean(mean(wave_sig_ht_temp.*region_data));
    end

    clear aice ice_mask k_means file_idx X_map
end

















%%
% % cbar = cmocean('deep',num_clusters);
% % close all
% % N = 6;
% % %cmocean()
% % file_number_vec = 1:351;
% % region_label = {'1','2','4'};
% % date_label = datetime(dirdates.waves);
% % date_label_vec = date_label(file_number_vec);
% % 
% % 
% % conFigure(11,0.2)
% % f = figure;
% % subplot(3,1,2)
% % C = linspecer(N);
% % axes('NextPlot','replacechildren', 'ColorOrder',C);
% % plot(date_label_vec,[perimeter_region([1,2,4],file_number_vec); sum(perimeter_region(:,file_number_vec))],'LineWidth',5);
% % %set(p, {'color'}, {[cbar(1,:)];[cbar(2,:)];[cbar(4,:)];[0.5,0.1,0.1]});
% % ylabel('Perimeter [block edges]')
% % %legend(region_label,'Location','southeast')
% % 
% % subplot(3,1,3)
% % plot(date_label_vec,[area_region([1,2,4],file_number_vec).*1e-6.*1e-6; sum(area_region(:,file_number_vec).*1e-6.*1e-6)],'LineWidth',5)
% % ylabel('Area [million km$^2$]')
% % %legend(region_label,'Location','southeast')
% % 
% % % Create a tile on the right column to get its position
% % legend_label = region_label;
% % legend_label{end+1} = 'Sum of regions';
% % ax = subplot(3,1,1,'Visible','off');
% % axPos = ax.Position;
% % delete(ax)
% % % Construct a Legend with the data from the sub-plots
% % hL = legend(legend_label);
% % % Move the legend to the position of the extra axes
% % hL.Position(1:2) = axPos(1:2);






% subplot(4,1,4)
% perim_per_area = perimeter_region([1,2,4],file_number_vec)./(area_region([1,2,4],file_number_vec).*1e-6.*1e-6);
% plot(date_label_vec,perim_per_area,'LineWidth',5)
% ylabel('Perimeter per Area [block per million km$^2$]')



%text(1,1,sprintf("$sigma =$ %g",std(perimeter_region'./(area_region./1000)')))
% ta = annotation('textarrow', [0.76 0.83], [0.71 0.71]);
% ta.String = sprintf("$\sigma =$ %g",std(perimeter_region(4,:)'./(area_region(4,:)./1000)'));
% ta.Interpreter = "latex";
%legend(region_label,'Location','southeast')



% diff_perim = perim_per_area(:,2:end)' - perim_per_area(:,1:end-1)';
% subplot(2,3,5)
% bar(std(diff_perim))
% xticklabels(region_label)
% ylabel('S.D. $\Delta$ perimeter per area')



%exportgraphics(f,'wave_area_perimeter.pdf','ContentType','vector')
















%% FUNCTIONS
function [X_total, idx] = read_data_vec(filenames,sector,var_list)
    [nfiles,~] = size(filenames);
    if nfiles == 1
        filename = filenames;
    else
        filename = filenames(1,:);
    end
    
    SIC_threshold = 0.05;
    [~, sector_mask] = data_format_sector(filename,"aice",sector);
    [lat,lon,~,~,~] = grid_read('om2');
    lat_vec = reshape(lat(sector_mask),numel(lat(sector_mask)),1);
    lon_vec = reshape(lon(sector_mask),numel(lon(sector_mask)),1);

    NFSD = ncread(filename,"NFSD");
    NCAT = ncread(filename,"NCAT");
    [floe_binwidth, ~, ~, ~] = cice_parameters(NFSD);
    
    X_total = [];
    for j = 1:nfiles % Each file
        filename = filenames(j,:);
        [aice, sector_mask] = data_format_sector(filename,"aice",sector);
        ice_mask = aice > SIC_threshold;
        X_temp = [];
        for i = 1:numel(var_list)
            if var_list{i} == "pancake"
                temp = data_format_sector(filename,"afsdn",sector);
                [len, wid] = size(aice);
                for ii = 1:len
                    for jj = 1:wid
                        if aice(ii,jj) > SIC_threshold
                            temp2(ii,jj) =  temp(ii,jj,1,1)./aice(ii,jj).*floe_binwidth(1).*100;
                        else
                            temp2(ii,jj) =  NaN;
                        end
                    end
                end
                %idx = isnan(temp2);
                %temp2(idx) = 0;
                temp_vec = reshape(temp2(sector_mask),numel(temp2(sector_mask)),1);
                X_temp = [X_temp, temp_vec];
            elseif var_list{i} == "full_fsd"
                temp = data_format_sector(filename,"afsd",sector);
                [len, wid] = size(aice);
                for kk = 1:length(NFSD)
                    for ii = 1:len
                        for jj = 1:wid
                            if aice(ii,jj) > SIC_threshold
                                temp2(ii,jj) =  temp(ii,jj,kk)./aice(ii,jj).*floe_binwidth(kk).*100;
                            else
                                temp2(ii,jj) =  NaN;
                            end
                        end
                    end
                %idx = isnan(temp2);
                %temp2(idx) = 0;
                temp_vec = reshape(temp2(sector_mask),numel(temp2(sector_mask)),1);
                X_temp = [X_temp, temp_vec];
                end
            elseif var_list{i} == "full_fstd"
                temp = data_format_sector(filename,"afsdn",sector);
                [len, wid] = size(aice);
                for nc = 1:length(NCAT)
                for kk = 1:length(NFSD)
                    for ii = 1:len
                        for jj = 1:wid
                            if aice(ii,jj) > SIC_threshold
                                temp2(ii,jj) =  temp(ii,jj,kk,nc)./aice(ii,jj).*floe_binwidth(kk).*100;
                            else
                                temp2(ii,jj) =  NaN;
                            end
                        end
                    end
                %idx = isnan(temp2);
                %temp2(idx) = 0;
                temp_vec = reshape(temp2(sector_mask),numel(temp2(sector_mask)),1);
                X_temp = [X_temp, temp_vec];
                end
                end
            elseif var_list{i} == "thick_pan"
                temp = data_format_sector(filename,"afsdn",sector);
                [len, wid] = size(aice);
                for ii = 1:len
                    for jj = 1:wid
                        if aice(ii,jj) > SIC_threshold
                            for ll = 2:5
                                temp2(ii,jj) =  temp2(ii,jj) + temp(ii,jj,1,ll)./aice(ii,jj).*floe_binwidth(1).*100;
                            end
                        else
                            temp2(ii,jj) =  NaN;
                        end
                    end
                end
                %idx = isnan(temp2);
                %temp2(idx) = 0;
                temp_vec = reshape(temp2(sector_mask),numel(temp2(sector_mask)),1);
                X_temp = [X_temp, temp_vec];
            elseif var_list{i} == "large"
                temp = data_format_sector(filename,"afsdn",sector);
                [len, wid] = size(aice);
                for ii = 1:len
                    for jj = 1:wid
                        if aice(ii,jj) > SIC_threshold
                            temp2(ii,jj) =  temp(ii,jj,numel(NFSD),1)./aice(ii,jj).*floe_binwidth(numel(NFSD)).*100;
                        else
                            temp2(ii,jj) =  NaN;
                        end
                    end
                end
                temp_vec = reshape(temp2(sector_mask),numel(temp2(sector_mask)),1);
                X_temp = [X_temp, temp_vec];

            elseif var_list{i} == "dafsd_newi"
                temp = data_format_sector(filename,"dafsd_newi",sector);
                [len, wid] = size(aice);
                for ii = 1:len
                    for jj = 1:wid
                        if aice(ii,jj) > SIC_threshold
                            temp_fsd(:) = temp(ii,jj,:);
                            temp2(ii,jj) =  sum((temp_fsd./aice(ii,jj)).*NFSD');
                        else
                            temp2(ii,jj) =  NaN;
                        end
                    end
                end
                temp_vec = reshape(temp2(sector_mask),numel(temp2(sector_mask)),1);
                X_temp = [X_temp, temp_vec];
            elseif var_list{i} == "dafsd_weld"
                temp = data_format_sector(filename,"dafsd_weld",sector);
                [len, wid] = size(aice);
                for ii = 1:len
                    for jj = 1:wid
                        if aice(ii,jj) > SIC_threshold
                            temp_fsd(:) = temp(ii,jj,:);
                            temp2(ii,jj) =  sum((temp_fsd./aice(ii,jj)).*NFSD');
                        else
                            temp2(ii,jj) =  NaN;
                        end
                    end
                end
                temp_vec = reshape(temp2(sector_mask),numel(temp2(sector_mask)),1);
                X_temp = [X_temp, temp_vec];
            elseif var_list{i} == "dafsd_wave"
                temp = data_format_sector(filename,"dafsd_wave",sector);
                [len, wid] = size(aice);
                for ii = 1:len
                    for jj = 1:wid
                        if aice(ii,jj) > SIC_threshold
                            temp_fsd(:) = temp(ii,jj,:);
                            temp2(ii,jj) =  sum((temp_fsd./aice(ii,jj)).*NFSD');
                        else
                            temp2(ii,jj) =  NaN;
                        end
                    end
                end
                temp_vec = reshape(temp2(sector_mask),numel(temp2(sector_mask)),1);
                X_temp = [X_temp, temp_vec];
            elseif var_list{i} == "dafsd_latm"
                temp = data_format_sector(filename,"dafsd_latm",sector);
                [len, wid] = size(aice);
                for ii = 1:len
                    for jj = 1:wid
                        if aice(ii,jj) > SIC_threshold
                            temp_fsd(:) = temp(ii,jj,:);
                            temp2(ii,jj) =  sum((temp_fsd./aice(ii,jj)).*NFSD');
                        else
                            temp2(ii,jj) =  NaN;
                        end
                    end
                end
                temp_vec = reshape(temp2(sector_mask),numel(temp2(sector_mask)),1);
                X_temp = [X_temp, temp_vec];
             elseif var_list{i} == "dafsd_latg"
                temp = data_format_sector(filename,"dafsd_latg",sector);
                [len, wid] = size(aice);
                for ii = 1:len
                    for jj = 1:wid
                        if aice(ii,jj) > SIC_threshold
                            temp_fsd(:) = temp(ii,jj,:);
                            temp2(ii,jj) =  sum((temp_fsd./aice(ii,jj)).*NFSD');
                        else
                            temp2(ii,jj) =  NaN;
                        end
                    end
                end
                temp_vec = reshape(temp2(sector_mask),numel(temp2(sector_mask)),1);
                X_temp = [X_temp, temp_vec];
            elseif var_list{i} == "strair"
                tempx = data_format_sector(filename,"strairx",sector);
                tempy = data_format_sector(filename,"strairy",sector);
                tempx_vec = reshape(tempx(sector_mask),numel(tempx(sector_mask)),1);
                tempy_vec = reshape(tempy(sector_mask),numel(tempy(sector_mask)),1);
                temp_vec = sqrt(tempx_vec.^2 + tempy_vec.^2);
                X_temp = [X_temp, temp_vec];
            elseif var_list{i} == "strocn"
                tempx = data_format_sector(filename,"strocnx",sector);
                tempy = data_format_sector(filename,"strocny",sector);
                tempx_vec = reshape(tempx(sector_mask),numel(tempx(sector_mask)),1);
                tempy_vec = reshape(tempy(sector_mask),numel(tempy(sector_mask)),1);
                temp_vec = sqrt(tempx_vec.^2 + tempy_vec.^2);
                X_temp = [X_temp, temp_vec];
            elseif var_list{i} == "strint"
                tempx = data_format_sector(filename,"strintx",sector);
                tempy = data_format_sector(filename,"strinty",sector);
                tempx_vec = reshape(tempx(sector_mask),numel(tempx(sector_mask)),1);
                tempy_vec = reshape(tempy(sector_mask),numel(tempy(sector_mask)),1);
                temp_vec = sqrt(tempx_vec.^2 + tempy_vec.^2);
                X_temp = [X_temp, temp_vec];
            elseif var_list{i} == "strcor"
                tempx = data_format_sector(filename,"strcorx",sector);
                tempy = data_format_sector(filename,"strcory",sector);
                tempx_vec = reshape(tempx(sector_mask),numel(tempx(sector_mask)),1);
                tempy_vec = reshape(tempy(sector_mask),numel(tempy(sector_mask)),1);
                temp_vec = sqrt(tempx_vec.^2 + tempy_vec.^2);
                X_temp = [X_temp, temp_vec];
            elseif var_list{i} == "vel"
                tempx = data_format_sector(filename,"uvel",sector);
                tempy = data_format_sector(filename,"vvel",sector);
                tempx_vec = reshape(tempx(sector_mask),numel(tempx(sector_mask)),1);
                tempy_vec = reshape(tempy(sector_mask),numel(tempy(sector_mask)),1);
                temp_vec = sqrt(tempx_vec.^2 + tempy_vec.^2);
                X_temp = [X_temp, temp_vec];  
            elseif var_list{i} == "atmspd"
                tempx = data_format_sector(filename,"uatm",sector);
                tempy = data_format_sector(filename,"vatm",sector);
                tempx_vec = reshape(tempx(sector_mask),numel(tempx(sector_mask)),1);
                tempy_vec = reshape(tempy(sector_mask),numel(tempy(sector_mask)),1);
                temp_vec = sqrt(tempx_vec.^2 + tempy_vec.^2);
                X_temp = [X_temp, temp_vec];  
            else
                temp = data_format_sector(filename,var_list{i},sector);
                [len, wid] = size(aice);
                for ii = 1:len
                    for jj = 1:wid
                        if aice(ii,jj) > SIC_threshold
                            temp2(ii,jj) =  temp(ii,jj);
                        else
                            temp2(ii,jj) =  NaN;
                        end
                    end
                end
                temp_vec = reshape(temp2(sector_mask),numel(temp2(sector_mask)),1);
                X_temp = [X_temp, temp_vec];
            end
        end 
        X_temp = [X_temp, lat_vec, lon_vec];
        [idx(j),~] = size(X_temp);
        X_total(:,:,j) = X_temp; %[X_total; X_temp];
    end
   
end


function var_out = variable_dict(var_in)
    for i = 1:numel(var_in)
        if var_in{i} == "aice"
            var_out{i} = 'SIC';
        elseif var_in{i} == "hi"
            var_out{i} = 'Ice thickness';
        elseif var_in{i} == "iage"
            var_out{i} = 'Ice age';
        elseif var_in{i} == "fsdrad"
            var_out{i} = 'Floe mean radius';
        elseif var_in{i} == "pancake"
            var_out{i} = 'Pancake ice';
        elseif var_in{i} == "large"
            var_out{i} = 'Sheet ice';
        elseif var_in{i} == "wave_sig_ht"
            var_out{i} = 'Sig. wave height';
        elseif var_in{i} == "frazil"
            var_out{i} = 'Frazil';
        elseif var_in{i} == "frzmlt"
            var_out{i} = 'Freeze/melt pot.';
        elseif var_in{i} == "congel"
            var_out{i} = 'Congelation rate';
        elseif var_in{i} == "ice_present"
            var_out{i} = 'Ice present';
        elseif var_in{i} == "peak_period"
            var_out{i} = 'Peak period';
        elseif var_in{i} == "vel"
            var_out{i} = 'Ice speed';
        elseif var_in{i} == "strint"
            var_out{i} = 'Internal stress';
        elseif var_in{i} == "Tair"
            var_out{i} = 'Air temp.';
        elseif var_in{i} == "dafsd_latg"
            var_out{i} = 'Lat growth';
        elseif var_in{i} == "dafsd_latm"
            var_out{i} = 'Lat melt';
        elseif var_in{i} == "dafsd_newi"
            var_out{i} = 'New ice';
        elseif var_in{i} == "dafsd_weld"
            var_out{i} = 'Welding';
        elseif var_in{i} == "dafsd_wave"
            var_out{i} = 'Wave breakup';
        elseif var_in{i} == "dafsd_wave"
            var_out{i} = 'Wave breakup';
        elseif var_in{i} == "sice"
            var_out{i} = 'Ice salinity';
        elseif var_in{i} == "alvl"
            var_out{i} = 'Level ice area frac.';
        elseif var_in{i} == "vlvl"
            var_out{i} = 'Level ice vol. frac.';
        elseif var_in{i} == "krdgn"
            var_out{i} = 'Mean ridge thickness per thickness of ridging ice';
        elseif var_in{i} == "ardg"
            var_out{i} = 'Fractional area of ridged ice';
        elseif var_in{i} == "vrdg"
            var_out{i} = 'Fractional volume of ridged ice';
        elseif var_in{i} == "hs"
            var_out{i} = 'Snow thickness';
        elseif var_in{i} == "FYarea"
            var_out{i} = 'First year ice area';
        elseif var_in{i} == "full_fsd"
            for j = 1:12
                var_out{i+j-1} = sprintf('FSD cat. %g',j);
            end
        elseif var_in{i} == "full_fstd"
            count = 1;
           for nc = 1:5
                for j = 1:12
                    var_out{i+count-1} = sprintf('FSTD cat. (%g %g)',j,nc);
                    count = count + 1;
                end
           end
        elseif var_in{i} == "uvel"
            var_out{i} = 'Ice u-speed';
        elseif var_in{i} == "vvel"
            var_out{i} = 'Ice v-speed';
        elseif var_in{i} == "strength"
            var_out{i} = 'Ice strength';
        elseif var_in{i} == "divu"
            var_out{i} = 'Divergence';
        elseif var_in{i} == "shear"
            var_out{i} = 'Shear';
        elseif var_in{i} == "daidtt"
            var_out{i} = 'Change in SIC (thermo.)';
        elseif var_in{i} == "daidtd"
            var_out{i} = 'Change in SIC (dynamics)';
        elseif var_in{i} == "dagedtt"
            var_out{i} = 'Change in age (thermo.)';
        elseif var_in{i} == "dagedtd"
            var_out{i} = 'Change in age (dynamics)';
        else
            var_out{i} = var_in{i};
        end
    end
end
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



function [MIZ_width, class,MIZ] = calculate_miz_width(filename,sector,k_means,miz_def)
    aice = data_format_sector(filename,'aice',sector);
    [lat,lon,~,ulat,ulon] = grid_read('om2');
    ice_mask =  aice > 0.15;
    k_means(~ice_mask) = NaN;
    idx_miz = isnan(k_means(1,:));
    [len,~] = size(aice);
    for i = 1:len
        temp = k_means(i,~idx_miz);
        [len1,wid1] = size(temp);
        if len1 == 0
            edge_class(i) = 0;
        elseif wid1 == 0
            edge_class(i) = 0;
        else
            edge_class(i) = temp(end); % Take the last cell that > 0.15 SIC
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
            temp = pathdist([ulat_vec_temp(south_point-1),ulat_vec_temp(north_point)],[ulon_vec_temp(south_point),ulon_vec_temp(north_point)],'km');
            MIZ_width(i,1:2) = temp(1:2);
        end
    end
    MIZ_width = MIZ_width(:,2);
end









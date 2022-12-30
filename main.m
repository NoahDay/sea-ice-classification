% Classification of CICE ice covers

% We want to experiment with different combination of variables over
% multiple years of data. 
% 1. Only ice cover "snapshot" variables
% 2. Add dynamic terms
% 3. Add thermodynamic terms
% 4. All of the terms

%% Setup
clear all
clc
close all

addpath /Users/noahday/Github/cice-plotting-tools/functions
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
var_list = {'aice','hi','hs','fsdrad','sice','iage','vlvl','vrdg'};
    % Static: {'aice','hi','hs','fsdrad','sice','iage','vlvl','vrdg'}
    % Dynamics: {'aice','hi','hs','fsdrad','sice','iage','vlvl','vrdg', 
    %            'uvel','vvel','strength','divu','shear','daidtd','daidtt',
    %            'dagedtd','dagedtt','dafsd_latg','dafsd_latm','dafsd_newi',
    %            'dafsd_weld','dafsd_wave'};

[X_raw, row_idx]= read_data_vec(filenames,sector,var_list); % Indices: [var_list, lon, lat]
label_vec = variable_dict(var_list);

% Save the data
save_data = 1; % On = 1
if save_data
    data.Xunstandard = X_raw;
    data.row_idx = row_idx;
    save_filename = strcat('cover_5percent_2015-19.mat');
    save(save_filename,'data','-v7.3');
    clear data
end

%% Clean the data
clear Xnan_temp 
load('data/cover_5percent_2015-19.mat')
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

[Xnan,row_idx] = clearNaN(X_temp);
disp('Clear NaN done!')
clear Xtemp

%%


load('data/cleanded_cover_data2015-19.mat')
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
toc

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









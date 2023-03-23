% Classification of CICE ice covers

% We want to experiment with different combination of variables over
% multiple years of data. 
% 1. Only ice cover "snapshot" variables
% 2. Add dynamic terms
% 3. Add thermodynamic terms
% 4. All of the terms

%% Setup
%clear all
clc
close all
addpath functions
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
var_type = "dynamics";
var_list = variable_list(var_type);

%% 
%var_list = {'aice','hi','hs','fsdrad','sice','iage','vlvl','vrdg'};
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
    save_filename = strcat(var_type,'_5percent_2015-19.mat');
    save(save_filename,'data','-v7.3');
    clear data
end

%% Clean the data
clear Xnan_temp 
load('data/dynamics_5percent_2015-19.mat')
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
% Taking just the last year
%row_idx(1826-365:end);
Xnan = Xnan(sum(row_idx(1:1826-365)):end,:);
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
num_clusters = 3;
tic
[kmeans_idx,C] = kmeans(X,num_clusters,'MaxIter',300);
toc
%%
label_vec = variable_dict(var_list);

kmeans_cluster.idx = kmeans_idx;
kmeans_cluster.row_idx = row_idx;
kmeans_cluster.label_vec = label_vec;
kmeans_cluster.C = C;
kmeans_cluster.X_standard_all = X_standard_all;
kmeans_cluster.num_clusters = num_clusters;

save_filename = strcat('kmeans_3cases_',var_type,sprintf('%g',num_clusters),'_classes.mat');
save(save_filename,'kmeans_cluster','-v7.3');

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


%% Silhouette test

X = X_standard_all(1:100000:end,1:8);%end-2);

%

close all
conFigure(8,3)
f = figure;
tiledlayout(1,4)
for i = 2:5
    nexttile
    [kmeans_idx,C] = kmeans(X,i,'MaxIter',300);
    S = silhouette(X,kmeans_idx);
    silhouette(X,kmeans_idx)
    xline(mean(S),'--')
end

exportgraphics(f,'sil_matrix.png','ContentType','image')
%% Feature selection - forward selection
X = X_standard_all(1:10000:end,1:end-2);
[len, n_total] = size(X);
X_all = X;
var_all = var_list;
label_all = label_vec;
X_opt = [];
var_opt = {};
label_opt = {};
nvar_added = 0;
silhouette_score_n = [];
for nvar = 1:n_total-1
    % Search to find the variable with the best fit
    n_search = n_total;
    [~, nvar_added] = size(X_opt);
    for i = 2:n_search - nvar_added
        % Grab a variable
        X_vec = X_all(:,i);
        % Add that to the dataset
        X_temp = [X_opt X_vec];
        [kmeans_idx,C] = kmeans(X_temp,i,'MaxIter',300);
        silhouette_score(i) = mean(silhouette(X,kmeans_idx));
    end
    maxloc = find(max(silhouette_score)==silhouette_score);
    X_opt = [X_opt  X_all(:,maxloc)];
    X_all(:,maxloc) = [];
    var_opt(end+1) = var_all(maxloc);
    var_all(maxloc) = [];
    label_opt(end+1) = label_all(maxloc);
    label_all(maxloc) = [];
    silhouette_score_n(nvar)= max(silhouette_score);
    
end

%% Feature selection - forward selection
% https://medium.com/analytics-vidhya/k-means-algorithm-in-4-parts-4-4-42bc6c781e46
clc
row_cumsum = cumsum(row_idx);
X = X_standard_all(row_cumsum(end-180):row_cumsum(end-179),1:7);
%X = X_standard_all(1:100:end,1:7); %1:end-2);
[len, n_total] = size(X);


var_opt = {};
label_opt = {};

silhouette_score_n = [];
selected_variables = {};

maxvars = 7; % Maximum number of variables retained
kmin = 2; % Minimum number of clusters
kmax = 5; % Maximum number of clusters
for k = kmin:kmax
    % Calculate the silhouette score for each variable
    X_all = X;
    var_all = var_list;
    label_all = label_vec;
    X_temp = [];
    X_opt = [];
    nvar_added = 0;
    for i = 1:maxvars
        % Search to find the variable with the best fit
        maxloc = [];
        silhouette_score = [];
        for j = 1:n_total - nvar_added
            % Grab a variable
            X_vec = X_all(:,j);
            % Add that to the dataset
            X_temp = [X_opt X_vec];
            [kmeans_idx,C] = kmeans(X_temp,k,'MaxIter',300);
            silhouette_score(j) = mean(silhouette(X_temp,kmeans_idx));
        end
        maxloc = find(max(silhouette_score)==silhouette_score);
        X_opt = [X_opt  X_all(:,maxloc)];
        X_all(:,maxloc) = [];
        var_opt(k,i) = var_all(maxloc);
        var_all(maxloc) = [];
        selected_variables(k,i) = label_all(maxloc);
        label_all(maxloc) = [];
        silhouette_score_n(k,i)= max(silhouette_score);

        nvar_added = nvar_added + 1;
    end
end
silhouette_score_n
selected_variables
%
close all
conFigure(11,1)
f = figure;
plot(1:length(silhouette_score_n(3,:)), silhouette_score_n(3,:))
set(gca,'xtick',1:length(silhouette_score_n))
set(gca,'xticklabel',selected_variables(3,:))
ylabel('Silhouette score')
ylim([0,1])

exportgraphics(f,'sil_score_summer.png','ContentType','image')

%%

close all
conFigure(8,3)
f = figure;
tiledlayout(1,4)
for k = 2:5
    nexttile
    [kmeans_idx,C] = kmeans(X(:,:),k,'MaxIter',300);
    silhouette(X,kmeans_idx)
    xline(mean(S),'--')
end

exportgraphics(f,'sil_matrix.png','ContentType','image')

%% PCA
close all
X = X_standard_all(1:100:end,1:7);
row_cumsum = cumsum(row_idx);
X = X_standard_all(row_cumsum(end-1):row_cumsum(end),1:7);

%coeff = pca(X);
[coeff,score,latent,tsquared,explained,mu] = pca(X);

conFigure(11,1)
f = figure;
b = bar(explained);
ylabel('Total variance explained [$\%$]')
xlabel 'Principal component'
matlab2tikz('pca_var_explained_summer.tex', 'standalone', true);

label_vec_tmp = label_vec(1:7);
%label_vec_tmp(8) = [];
conFigure(11,1)
f = figure;
b = bar(coeff(:,1).^2);
ylabel('Weighting for PC1 [-]')
xticks(1:length(label_vec_tmp))
 set(gca,'xticklabel',label_vec_tmp)
 [len wid] = size(X);
 yline((1/wid),'--')

exportgraphics(f,'pca1_loadings.png','ContentType','image')
addpath /Users/noahday/Documents/MATLAB/matlab2tikz/src/
matlab2tikz('pca_loadings1_summer.tex', 'standalone', true);



conFigure(11,1)
f = figure;
b = bar(coeff(:,2).^2);
ylabel('Weighting for PC2 [-]')
xticks(1:length(label_vec_tmp))
 set(gca,'xticklabel',label_vec_tmp)
 yline((1/wid),'--')

exportgraphics(f,'pca2_loadings.png','ContentType','image')
matlab2tikz('pca_loadings2_summer.tex', 'standalone', true);
%% PCA total, summer, winter
% Winter  = June, July, August
% Summer = Jan, Feb, March
close all

winter_start = 1826-213; % 1st June
winter_end = 1826-122; % 31 August

summer_start = 1826-364; % 1st Jan
summer_end = 1826-275; % 31 March



row_cumsum = cumsum(row_idx);



% TOTAL
X = X_standard_all(row_cumsum(summer_start):end,1:7);

%coeff = pca(X);
[coeff,score,latent,tsquared,explained,mu] = pca(X);

conFigure(11,1)
f = figure;
b = bar(explained);
ylabel('Total variance explained [$\%$]')
xlabel 'Principal component'
matlab2tikz('pca_var_explained_total.tex', 'standalone', true);
exportgraphics(f,'pca_var_explained_total.png','ContentType','image')

label_vec_tmp = label_vec(1:7);
%label_vec_tmp(8) = [];
conFigure(11,1)
f = figure;
b = bar(coeff(:,1).^2);
ylabel('Weighting for PC1 [-]')
xticks(1:length(label_vec_tmp))
 set(gca,'xticklabel',label_vec_tmp)
 [len wid] = size(X);
 yline((1/wid),'--')

exportgraphics(f,'pca1_loadings.png','ContentType','image')
addpath /Users/noahday/Documents/MATLAB/matlab2tikz/src/
matlab2tikz('pca_loadings1_total.tex', 'standalone', true);

% SUMMER
X = X_standard_all(row_cumsum(summer_start):row_cumsum(summer_end),1:7);

%coeff = pca(X);
[coeff,score,latent,tsquared,explained,mu] = pca(X);
close all
conFigure(11,1)
f = figure;
b = bar(explained);
ylabel('Total variance explained [$\%$]')
xlabel 'Principal component'
matlab2tikz('pca_var_explained_summer.tex', 'standalone', true);
exportgraphics(f,'pca_var_explained_summer.png','ContentType','image')

label_vec_tmp = label_vec(1:7);
%label_vec_tmp(8) = [];
conFigure(11,1)
f = figure;
b = bar(coeff(:,1).^2);
ylabel('Weighting for PC1 [-]')
xticks(1:length(label_vec_tmp))
 set(gca,'xticklabel',label_vec_tmp)
 [len wid] = size(X);
 yline((1/wid),'--')

exportgraphics(f,'pca1_loadings_summer.png','ContentType','image')
addpath /Users/noahday/Documents/MATLAB/matlab2tikz/src/
matlab2tikz('pca_loadings1_summer.tex', 'standalone', true);


% WINTER
X = X_standard_all(row_cumsum(winter_start):row_cumsum(winter_end),1:7);

%coeff = pca(X);
[coeff,score,latent,tsquared,explained,mu] = pca(X);

conFigure(11,1)
f = figure;
b = bar(explained);
ylabel('Total variance explained [$\%$]')
xlabel 'Principal component'
matlab2tikz('pca_var_explained_winter.tex', 'standalone', true);
exportgraphics(f,'pca_var_explained_winter.png','ContentType','image')

label_vec_tmp = label_vec(1:7);
%label_vec_tmp(8) = [];
conFigure(11,1)
f = figure;
b = bar(coeff(:,1).^2);
ylabel('Weighting for PC1 [-]')
xticks(1:length(label_vec_tmp))
 set(gca,'xticklabel',label_vec_tmp)
 [len wid] = size(X);
 yline((1/wid),'--')

exportgraphics(f,'pca1_loadings_winter.png','ContentType','image')
addpath /Users/noahday/Documents/MATLAB/matlab2tikz/src/
matlab2tikz('pca_loadings1_winter.tex', 'standalone', true);




%%
rng(0,'twister'); % For reproducibility
N = 100;
X = rand(N,20);
y = -ones(N,1);
y(X(:,3).*X(:,9)./X(:,15) < 0.4) = 1;

mdl = fscnca(X,y,'Solver','sgd','Verbose',1);


%%

for i = 1:wid
    X_temp = X(:,i);
    [kmeans_idx,C] = kmeans(X_temp,i,'MaxIter',300);
    silhouette_score(i) = mean(silhouette(X,kmeans_idx));
end
silhouette_score
% Add hs
for i = [1 3 4 5 6 7 8]
    if i < 2
        X_temp = X(:,[i 2]);
        [kmeans_idx,C] = kmeans(X_temp,i,'MaxIter',300);
        silhouette_score(i) = mean(silhouette(X,kmeans_idx));
    else
        X_temp = X(:,[2 i]);
        [kmeans_idx,C] = kmeans(X_temp,i,'MaxIter',300);
        silhouette_score(i) = mean(silhouette(X,kmeans_idx));
    end
end
silhouette_score
% Add fsdrad
for i = [1 3 5 6 7 8]
    if i < 2
        X_temp = X(:,[i 2 4]);
        [kmeans_idx,C] = kmeans(X_temp,i,'MaxIter',300);
        silhouette_score(i) = mean(silhouette(X,kmeans_idx));
    elseif i == 3 
        X_temp = X(:,[2 i 4]);
        [kmeans_idx,C] = kmeans(X_temp,i,'MaxIter',300);
        silhouette_score(i) = mean(silhouette(X,kmeans_idx));
    else
        X_temp = X(:,[2 4 i]);
        [kmeans_idx,C] = kmeans(X_temp,i,'MaxIter',300);
        silhouette_score(i) = mean(silhouette(X,kmeans_idx));
    end
end
var_list
silhouette_score

%%
label_vec

close all
figure
plot(X(:,2),X(:,8),'o')
%corr(X(:,6),X(:,16))

%%
rng('default')  % For reproducibility
X = [randn(10,2)+3;randn(10,2)-3];
scatter(X(:,1),X(:,2));
title('Randomly Generated Data');
clust = kmeans(X,2);
silhouette(X,clust)


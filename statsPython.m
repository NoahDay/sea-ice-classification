% Read in CSV
var_list = ["aice", "hi", "hs", "fsdrad", "iage", "alvl", "longitude", "latitude", "date", "k"];
filename = "/Volumes/NoahDay1TB/sea-ice-classification/data/kmeans_2011.csv";
standard_table = readtable(filename);

for year = 2012:2019
    filename = strcat("/Volumes/NoahDay1TB/sea-ice-classification/data/kmeans_",num2str(year),".csv");
    new_table = readtable(filename);
    standard_table = [ standard_table; new_table];
end

%% Read in the data
year = 2019;
filename = strcat("/Volumes/NoahDay1TB/sea-ice-classification/data/kmeans_",num2str(year),".csv");
standard_table = readtable(filename);

summary(standard_table)


%% PCA
close all
% X = X_standard_all(1:100:end,1:7);
% row_cumsum = cumsum(row_idx);
% X = X_standard_all(row_cumsum(end-1):row_cumsum(end),1:7);

X = [standard_table.aice, standard_table.hi, standard_table.hs, standard_table.fsdrad, standard_table.iage, standard_table.alvl];
X = X - mean(X);

%coeff = pca(X);
[coeff,score,latent,tsquared,explained,mu] = pca(X);

conFigure(11,1)
f = figure;
b = bar(explained);
ylabel('Total variance explained [$\%$]')
xlabel 'Principal component'
matlab2tikz('pca_var_explained.tex', 'standalone', true);

label_vec_tmp = ["SIC", "Ice thickness", "Snow thickness", "Floe size", "Age", "Level ice"];
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

%%

summer_idx = standard_table.date(:) < datetime(2019,03,01) | standard_table.date(:) > datetime(2019,11,30);
winter_idx = standard_table.date(:) >= datetime(2019,06,01) & standard_table.date(:) <= datetime(2019,08,31);

[coeff_total,score,latent,tsquared,explained_total,mu] = pca(X);

[coeff_summer,score,latent,tsquared,explained_summer,mu] = pca(X(summer_idx,:));

[coeff_winter,score,latent,tsquared,explained_winter,mu] = pca(X(winter_idx,:));

plot_data = [explained_total'; explained_summer'; explained_winter']';
plot_coeff_all = [coeff_total(:,1), coeff_summer(:,1), coeff_winter(:,1)].^2;

%plot_data = explained;
%plot_coeff_all = coeff(:,1); 

cmap_bar = ["#1c9099"; "#fdbb84"; "#bcbddc"];
cmap_bar = ["#1c9099"; "#fdbb84"; "#bcbddc"];

nIDs = 4;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); % {'(a)','(b)','(c)','(d)'}
%text(0.025,0.95,charlbl{1},'Units','normalized')

conFigure(30)
close all
%conFigure(15,2)
f = figure;
set(gcf, 'Position',  [100, 100, 100, 100])
tiledlayout(1,2)
nexttile
b = bar(plot_data);
ylabel('Variance explained [$\%$]')
xlabel 'Principal component'

b(1).FaceColor = cmap_bar(1,:);
b(2).FaceColor = cmap_bar(2,:);
b(3).FaceColor = cmap_bar(3,:);
text(0.025,0.95,charlbl{1},'Units','normalized')
 legend('Full year','Summer','Winter','')

nexttile

b = bar(plot_coeff_all);

b(1).FaceColor = cmap_bar(1,:);
b(2).FaceColor = cmap_bar(2,:);
b(3).FaceColor = cmap_bar(3,:);
ylabel('PC1 weighting [-]')
xticks(1:length(label_vec_tmp))
set(gca,'xticklabel',label_vec_tmp)
[len wid] = size(X);
 

%yline((1/wid),'--','HandleVisibility','off')
%yregion(0, 1/wid);
ybars = [0, 1/wid];
p1 = patch([min(xlim) max(xlim) max(xlim) min(xlim)], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8]);
p1.FaceAlpha = 0.2;
p1.LineStyle = '--';
%alpha(0.3)
text(0.025,0.95,charlbl{2},'Units','normalized')

 exportgraphics(f,'pca_analysis.png','ContentType','image')
matlab2tikz('pca_analysis.tex', 'standalone', true);
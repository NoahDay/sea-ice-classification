% Read in data
%clear all
effective = ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/mizWidth_2019.nc','effective')./1000;
absolute = ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/mizWidth_2019.nc','absolute')./1000;
sic = ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/mizWidth_2019.nc','sic')./1000;
lon = ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/mizWidth_2019.nc','LON');

effective_nofsd = ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/mizWidth_nofsd_2019.nc','effective')./1000;

t1 = datetime(2019,1,1);
t2 = datetime(2019,12,31);
ts_dates = t1:t2;


Cmap = [0.9805, 0.5000, 0.4453; 0.4416, 0.7490, 0.4322; 0.3639, 0.5755, 0.748];
    
  
%% Calculate SIC width from the ice edge
miz_idx = [0 1 1 0 0 0 1 1 1 0 0];
clear swh_width k_width sic_width
count = 0; 
for year_lp = 2011:2019
    count = count + 1;
    filename = strcat("/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_",num2str(year_lp),'.nc');
    aice_data = ncread(filename,'aice');
    swh_data = ncread(filename,'wave_sig_ht');
    k_data = ncread(filename,'k');
    %ppd_data = ncread(filename,'peak_period');
    [n_lon, n_lat, n_days] = size(aice_data);
    hte = ncread(filename,'HTE');

    for day_lp = 1:n_days
        for lon_lp = 1:n_lon
            sic_vec = aice_data(lon_lp,:,day_lp);
            miz_idx = (sic_vec > 0.15) & (sic_vec < 0.8);
            if sum(miz_idx) == 0
                sic_width(day_lp,lon_lp) = 0;
                
            else
                diff_idx = diff(miz_idx);
                temp = find(diff_idx == 1);
                inner_miz = temp(end) +1;
                temp = find(miz_idx == 1);
                outer_miz = temp(end);
                sic_width(day_lp,lon_lp) = sum(hte(lon_lp,inner_miz:outer_miz).*aice_data(lon_lp,inner_miz:outer_miz,day_lp));
               
            end
            
            % SWH
            swh_vec = swh_data(lon_lp,:,day_lp);
            swh_idx = (swh_vec > 0.3805);
            sic_idx = (sic_vec > 0.15);
            if sum(sic_idx) == 0
                swh_width(day_lp,lon_lp,count) = NaN;
            elseif sum(swh_idx) == 0
                swh_width(day_lp,lon_lp,count) = NaN;
            else
                swh_width(day_lp,lon_lp,count) = sum(hte(lon_lp,swh_idx).*aice_data(lon_lp,swh_idx,day_lp));
            end
    
    
            % Peak period
    %         ppd_vec = ppd_data(lon_lp,:,day_lp);
    %         ppd_idx = (ppd_vec < 15);
    %         ppd_idx = (sic_vec > 0.15);
    %         if sum(sic_idx) == 0
    %             ppd_width(day_lp,lon_lp) = NaN;
    %         elseif sum(swh_idx) == 0
    %             ppd_width(day_lp,lon_lp) = 0;
    %         else
    %             ppd_width(day_lp,lon_lp) = sum(hte(lon_lp,ppd_idx).*aice_data(lon_lp,ppd_idx,day_lp));
    %         end
    
            
    
            % K-means
            k_vec = k_data(lon_lp,:,day_lp);
            k_idx = k_vec == 0;
            if  sum(sic_vec > 0.15) == 0
                k_width(day_lp,lon_lp,count) = NaN;
                
            else
     
                k_width(day_lp,lon_lp,count) = sum(hte(lon_lp,k_idx).*aice_data(lon_lp,k_idx,day_lp));
            end
    
            %if day_lp == 180 & lon_lp== 180
            %    return
            %end
        end
    
    end
end
sprintf('SIC error rate is: %g',sum(sum(isnan(sic_width())))./numel(sic_width()))
sprintf('SWH error rate is: %g',sum(sum(isnan(swh_width())))./numel(swh_width()))
sprintf('k error rate is: %g',sum(sum(isnan(k_width())))./numel(k_width()))
%%

addpath /Users/noahday/GitHub/CICE-plotting-tools/functions/
addpath /Users/noahday/Documents/MATLAB/matlab2tikz/src/
effective = [];
sic_MIZ_width_mean = [];
swh_MIZ_width_mean = [];
for i = 1:9
    %MIZ_width_mean =nanmean(k_width(1:365,:,:),2)./1000;
    effective = [effective; k_width(1:365,:,i)']; 
    %sic_MIZ_width_mean = [sic_MIZ_width_mean; sic_width(1:365,:,i)'];
    swh_MIZ_width_mean = [swh_MIZ_width_mean; swh_width(1:365,:,i)'];
end

sic_MIZ_width_mean = nanmean(sic_width(1:365,:),2);
effective = effective./1000;
%effective = %mean(k_width,[2,3])./1000;

effective_feb = effective(:,31:31+28);
effective_may = effective(:,31+28+31+30:31+28+31+30+31);
effective_sep = effective(:,31+28+31+30+31+30+30+31:31+28+31+30+31+30+30+31+30);
effective_dec = effective(:,end-31:end);
effective_month_centre = [nanmedian(effective_feb,'all') nanmedian(effective_may,'all') nanmedian(effective_sep,'all') nanmedian(effective_dec,'all')];

effective_month_IQR = [quantile(effective_feb(:),0.25), quantile(effective_may(:),0.25), quantile(effective_sep(:),0.25), quantile(effective_dec(:),0.25);...
                       quantile(effective_feb(:),0.75), quantile(effective_may(:),0.75), quantile(effective_sep(:),0.75), quantile(effective_dec(:),0.75)];




MIZ_width_mean = mean(absolute);
MIZ_width_sd = 1.96*std(absolute);

ts_dates_CI  = [ts_dates, fliplr(ts_dates)];
MIZ_width_CI = [MIZ_width_mean-MIZ_width_sd, fliplr(MIZ_width_mean+MIZ_width_sd)];

corrected_MIZ_width_sd = nanstd(effective);
corrected_MIZ_width_mean = nanmean(effective);
MIZ_width_corrected_CI = [corrected_MIZ_width_mean-corrected_MIZ_width_sd, fliplr(corrected_MIZ_width_mean+corrected_MIZ_width_sd)];
corrected_MIZ_width_IQR = iqr(effective);
MIZ_IQR = [quantile(effective,0.25), fliplr(quantile(effective,0.75))];


% No FSD
corrected_MIZ_width_sd_nofsd = nanstd(effective_nofsd);
corrected_MIZ_width_mean_nofsd = nanmean(effective_nofsd);
MIZ_width_corrected_CI_nofsd = [corrected_MIZ_width_mean_nofsd-corrected_MIZ_width_sd_nofsd, fliplr(corrected_MIZ_width_mean_nofsd+corrected_MIZ_width_sd_nofsd)];
corrected_MIZ_width_IQR_nofsd = iqr(effective_nofsd);
MIZ_IQR_nofsd = [quantile(effective_nofsd,0.25), fliplr(quantile(effective_nofsd,0.75))];



corrected_MIZ_width_mean_nofsd = nanmean(effective_nofsd);
sic_MIZ_width_mean = nanmean(sic_width./1000,2); %nanmean(sic);
swh_MIZ_width_mean = nanmean(swh_width(1:365,:,1:end)./1000,2);
k_MIZ_width_mean = movmean(nanmean(effective'./1000),1);
MIZ_IQR = [quantile(effective'./1000,0.25), fliplr(quantile(effective'./1000,0.75))];

%brouwer_dates = [31,31+28+31+30,31+28+31+30+31+30+31+31,31+28+31+30+31+30+31+31+30+31+30];
brouwer_linear_centre = [36.187, 143.457, 192.569, 60.743];
brouwer_linear_lower = [20.679, 100.808, 125.363, 33.603];
brouwer_linear_upper = [60.743, 186.107, 275.283, 148.627];
brouwer_dates = [datetime(2019,2,15), datetime(2019,5,15), datetime(2019,9,15), datetime(2019,12,15)];%[6,15,27,36];

font_size = 12;
legend_labels = {'Unsupervised MIZ \textbf{with} floe size','Wave propagation observations$^{ 1}$','Unsupervised MIZ \textbf{without} floe size','Sea ice concentration (15--80$\%$)','Significant wave height $ > \sigma$'};

close all
f = figure('Position',[0, 0, 40, 10]);
conFigure(11,4)
set(gca,'FontSize',20) 
brouwer_color = [0.4660 0.6740 0.1880];
MIZ_color = Cmap(1,:); 
hold on


%plot(datenum(ts_dates),MIZ_width_mean,'linewidth',2)
plot(datenum(ts_dates),movmean(corrected_MIZ_width_mean,2),'color',MIZ_color,'linewidth',2)
%plot(datenum(brouwer_dates),brouwer_linear_centre,'s','color',[0.4660 0.6740 0.1880],'markersize',10,"MarkerFaceColor",brouwer_color)
errorbar(datenum(brouwer_dates),brouwer_linear_centre,...
    brouwer_linear_centre-brouwer_linear_lower,...
    brouwer_linear_upper-brouwer_linear_centre,...
    "o","MarkerSize",10,...
    "MarkerEdgeColor",brouwer_color,"MarkerFaceColor",brouwer_color,'linewidth',2,'Color',brouwer_color)


plot(datenum(ts_dates),movmean(corrected_MIZ_width_mean_nofsd,1),'linewidth',2,'color',[254,196,79]/256)

plot(datenum(ts_dates),sic_MIZ_width_mean(1:365),'linewidth',2,'color',[189,201,225]/256)


plot(datenum(ts_dates),movmean(swh_MIZ_width_mean(1:365),10),'linewidth',2,'color',[66,146,198]/256)


h = fill(datenum(ts_dates_CI), MIZ_width_corrected_CI, 1,....
       'facecolor',MIZ_color, ...%
       'edgecolor',MIZ_color, ...
       'facealpha', 0.2);



h = fill(datenum(ts_dates_CI), MIZ_width_corrected_CI_nofsd, 1,....
       'facecolor',[254,196,79]/256, ...%
       'edgecolor',[254,196,79]/256, ...
       'facealpha', 0.05);



 ax = ancestor(h, 'axes');
 ax.XAxis.Exponent = 0;
 xtickformat('%.0f');


legend(legend_labels,'Location','northwest','Orientation','vertical')
legend boxoff 
ylabel('Corrected MIZ width [km]')

ylim([0,300])
datetick('x','mmm')
%Fix ticks
%xlim([ts_dates(1),ts_dates(end)])
%xlabel('Months of 2019')
%ax=gca;
%ts_dates_months = cumsum([ts_dates(1), 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]);
%ax.XTick = ts_dates_months;
%ax.XTickLabels = {'J','F','M','A','M','J','J','A','S','O','N','D'};%cellstr(dirdates(cumsum([1826-365, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]),:));
%ax.XTickLabels = cellstr(datetime(dirdates(cumsum([1826-365, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]),:),"Format","MMM")+1);
exportgraphics(f,'ts_mizwidth1.pdf','ContentType','vector','BackgroundColor','none')
addpath /Users/noahday/Documents/MATLAB/matlab2tikz/src/
matlab2tikz('miz_comp_brouwer2.tex', 'standalone', true);


%ax.XTick=ax.XTick(1:12:end)
%Set xtickformat -- read more about possible formats: 
%https://se.mathworks.com/help/matlab/ref/datestr.html
%datetick('x','mm','keepticks');

%%






















%% MIZ WIDTH Produce tikz for EGU 2023
addpath /Users/noahday/GitHub/CICE-plotting-tools/functions/
addpath /Users/noahday/Documents/MATLAB/matlab2tikz/src/

MIZ_width_mean = mean(absolute);
MIZ_width_sd = std(absolute);

ts_dates_CI  = [ts_dates, fliplr(ts_dates)];
MIZ_width_CI = [MIZ_width_mean-MIZ_width_sd, fliplr(MIZ_width_mean+MIZ_width_sd)];

MIZ_width_mean =nanmean(k_width(1:365,:,:),2)./1000;
effective = squeeze(MIZ_width_mean)'; 
corrected_MIZ_width_sd = nanstd(effective);

corrected_MIZ_width_mean = nanmean(effective);
MIZ_width_corrected_CI = [corrected_MIZ_width_mean-corrected_MIZ_width_sd, fliplr(corrected_MIZ_width_mean+corrected_MIZ_width_sd)];
corrected_MIZ_width_IQR = iqr(effective);
MIZ_IQR = [quantile(effective,0.25), fliplr(quantile(effective,0.75))];

corrected_MIZ_width_mean_nofsd = nanmean(effective_nofsd);
sic_MIZ_width_mean = nanmean(sic_width'./1000); %nanmean(sic);

brouwer_linear_centre = [36.187, 143.457, 192.569, 60.743];
brouwer_linear_lower = [20.679, 100.808, 125.363, 33.603];
brouwer_linear_upper = [60.743, 186.107, 275.283, 148.627];
brouwer_dates = [datetime(2019,2,15), datetime(2019,5,15), datetime(2019,9,15), datetime(2019,12,15)];%[6,15,27,36];
font_size = 12;
legend_labels = {'Unsupervised MIZ \textbf{with} floe size','Wave penetration$^{ 1}$','Unsupervised MIZ \textbf{without} floe size','Sea ice concentration (15--80$\%$)'};



close all
% 1st line
f = figure('Position',[0, 0, 30, 10]);
conFigure(12,3)
set(gca,'FontSize',18) 
brouwer_color = [0.4660 0.6740 0.1880];
MIZ_color = Cmap(1,:); %[0.9290 0.6940 0.1250];
hold on
plot(datenum(ts_dates),movmean(corrected_MIZ_width_mean,1),'color',MIZ_color,'linewidth',2)
errorbar(datenum(brouwer_dates),brouwer_linear_centre,brouwer_linear_centre-brouwer_linear_lower,brouwer_linear_upper-brouwer_linear_centre,"s","MarkerSize",10,...
    "MarkerEdgeColor",brouwer_color,"MarkerFaceColor",brouwer_color,'linewidth',2,'Color',brouwer_color)


plot(datenum(ts_dates),movmean(corrected_MIZ_width_mean_nofsd,1),'linewidth',2,'color',[201,148,199]/256)
%plot(datenum(ts_dates),sic_MIZ_width_mean,'linewidth',2,'color',[189,201,225]/256)

h = fill(datenum(ts_dates_CI), MIZ_width_corrected_CI, 1,....
       'facecolor',MIZ_color, ...%
       'edgecolor',MIZ_color, ...
       'facealpha', 0.2);

% ax = ancestor(h, 'axes');
% ax.XAxis.Exponent = 0;
% xtickformat('%.0f');
 
legend(legend_labels(1:3),'Location','northwest','Orientation','vertical')
ylabel('Corrected MIZ width [km]')

ylim([0,500])

%Fix ticks
datetick('x','mmm')

matlab2tikz('miz_comp_brouwer3.tex', 'standalone', true);
%exportgraphics(f,'ts_mizwidth.pdf','ContentType','vector','BackgroundColor','none')
%addpath /Users/noahday/Documents/MATLAB/matlab2tikz/src/



%%
%clear
close all;
figure;

theta = linspace( 0.0, 2.0*pi, 20 );
%h = fill( cos(theta), sin(theta), 'Red' );
fill(datenum(ts_dates_CI), MIZ_width_corrected_CI, 'Red')

matlab2tikz( 'test.tikz', 'standalone',true );
system( 'pdflatex test.tikz' );

%%
upper_loc = ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/mizWidth_2019.nc','upper_location');
lower_loc = ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/mizWidth_2019.nc','lower_location');

upper_ts = mean(upper_loc);
lower_ts = mean(lower_loc);

close all
conFigure(11)
figure
plot(ts_dates,lower_ts)
hold on
plot(ts_dates,upper_ts)
ylabel("MIZ location [$^\circ$ S]")
legend("Southern boundary","Northern boundary")

upper_ts_dt = diff(upper_loc,1);
lower_ts_dt = diff(lower_loc,1);

conFigure(11)
figure
%plot(ts_dates,movmean(lower_ts_dt(1:30:end,:),1),'color','b')
%hold on
plot(ts_dates,movmean(upper_ts_dt(1:30:end,:),1))
ylabel("Change in MIZ location [$^\circ$ S]")
%legend("Southern boundary","Northern boundary")
legend
%% Change in MIZ width
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
%yticks = 1:10:350;
%yticklabels = lon(1:10:350);
%% Change in Upper location over the year
close all
plot_data = (diff(lower_loc(:,:),1));
plot_data = [plot_data(281:end,:); plot_data(1:280,:)];

conFigure(30)
figure
h = pcolor(ts_dates(1:end),plot_lon(1:end-1), movmean(movmean(plot_data,10,1),10,2));
 %h = pcolor();
 set(h, 'EdgeColor', 'none');
cmocean('balance','pivot');
%clim([-200,200])
cb = colorbar;
cb.Label.String = 'Change in MIZ lower [$^\circ$/day]';
xlabel 'Days of the year'
ylabel 'Longitude'

plot_data = (diff(upper_loc(:,:),1));
plot_data = [plot_data(281:end,:); plot_data(1:280,:)];

conFigure(30)
figure
h = pcolor(ts_dates(1:end),plot_lon(1:end-1), movmean(movmean(plot_data,10,1,'omitnan'),10,2));
 %h = pcolor();
 set(h, 'EdgeColor', 'none');
cmocean('balance','pivot');
%clim([-200,200])
cb = colorbar;
cb.Label.String = 'Change in MIZ Upper [$^\circ$/day]';
xlabel 'Days of the year'
ylabel 'Longitude'

%% Location of northern boundary
close all
plot_data = upper_loc(:,:);
plot_data = [plot_data(281:end,:); plot_data(1:280,:)];

conFigure(30)
figure
h = pcolor(ts_dates(1:end),plot_lon(1:end), movmean(plot_data,80,1));
 %h = pcolor();
 set(h, 'EdgeColor', 'none');
cmocean('thermal');
%clim([-200,200])
cb = colorbar;
cb.Label.String = 'MIZ Upper boundary [$^\circ$]';
xlabel 'Days of the year'
ylabel 'Longitude'
%% Change in MIZ width
lower_loc = ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/mizWidth_2019.nc','lower_location');

close all
plot_data = upper_loc(:,:);
plot_data = [plot_data(281:end,:); plot_data(1:280,:)];

conFigure(30)
figure
h = pcolor(ts_dates(1:end),plot_lon(1:end), movmean(plot_data,20,1));
 %h = pcolor();
 set(h, 'EdgeColor', 'none');
cmocean('balance');
%clim([-200,200])
cb = colorbar;
cb.Label.String = 'Change in MIZ Upper [$^\circ$/day]';
xlabel 'Days of the year'
ylabel 'Longitude'
%%
filename = '/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_2019.nc';
fsd_data = ncread(filename, 'fsdrad');
swh_data = ncread(filename, 'wave_sig_ht');
k_data = ncread(filename, 'k');
lon = ncread(filename,'LON');
[n_lon, n_lat, n_days] = size(fsd_data);

clear plot_data
for day_lp = 1:n_days
    for lon_lp = 1:n_lon
        k_mask = k_data(lon_lp,:,day_lp) == 0;
        plot_data_fsd(lon_lp,day_lp) = nanmean(fsd_data(lon_lp,k_mask,day_lp));
        plot_data_swh(lon_lp,day_lp) = nanmean(swh_data(lon_lp,k_mask,day_lp));
    end
end
%%
close all
conFigure(30)
figure
%h = pcolor(ts_dates(1:end),plot_lon(1:end-1), movmean(plot_data,80,1));
h = pcolor(movmean(plot_data_fsd,1,1))
 %h = pcolor();
 set(h, 'EdgeColor', 'none');
cmocean('thermal');
%clim([-200,200])
cb = colorbar;
cb.Label.String = 'MIZ floe size [m]';
xlabel 'Days of the year'
ylabel 'Longitude'


conFigure(30)
figure
%h = pcolor(ts_dates(1:end),plot_lon(1:end-1), movmean(plot_data,80,1));
h = pcolor(movmean(movmean(diff(plot_data_fsd,1,2),20,1,"omitnan"),20,2,"omitnan"))
 %h = pcolor();
 set(h, 'EdgeColor', 'none');
cmocean('balance','pivot');
clim([-10,10])
cb = colorbar;
cb.Label.String = 'Change MIZ floe size [m]';
xlabel 'Days of the year'
ylabel 'Longitude'

%%
close all
conFigure(30)
figure
%h = pcolor(ts_dates(1:end),plot_lon(1:end-1), movmean(plot_data,80,1));
h = pcolor(movmean(squeeze(nanmean(fsd_data,2)),10,2,"omitnan"))
 %h = pcolor();
 set(h, 'EdgeColor', 'none');
cmocean('thermal');
%clim([-50,50])
cb = colorbar;
cb.Label.String = 'Longitudinal average floe size [m]';
xlabel 'Days of the year'
ylabel 'Longitude'


conFigure(30)
figure
%h = pcolor(ts_dates(1:end),plot_lon(1:end-1), movmean(plot_data,80,1));
h = pcolor(movmean(diff(squeeze(nanmean(fsd_data,2)),1,2),10,2,"omitnan"))
 %h = pcolor();
 set(h, 'EdgeColor', 'none');
cmocean('balance');
clim([-50,50])
cb = colorbar;
cb.Label.String = 'Change in lon average floe size [m]';
xlabel 'Days of the year'
ylabel 'Longitude'


%% SWH
close all
conFigure(30)
figure
%h = pcolor(ts_dates(1:end),plot_lon(1:end-1), movmean(plot_data,80,1));
h = pcolor(movmean(movmean(diff(squeeze(plot_data_swh),1,2),5,1,"omitnan"),20,2,"omitnan"))
 %h = pcolor();
 set(h, 'EdgeColor', 'none');
cmocean('balance');
clim([-0.5,0.5])
cb = colorbar;
cb.Label.String = 'Change in SWH [m]';
xlabel 'Days of the year'
ylabel 'Longitude'

figure
%h = pcolor(ts_dates(1:end),plot_lon(1:end-1), movmean(plot_data,80,1));
h = pcolor(movmean(movmean((squeeze(plot_data_swh)),20,1,"omitnan"),20,2,"omitnan"))
 %h = pcolor();
 set(h, 'EdgeColor', 'none');
cmocean('balance','pivot');
%clim([0,5])
cb = colorbar;
cb.Label.String = 'SWH [m]';
xlabel 'Days of the year'
ylabel 'Longitude'

%%
winter_MIZ = effective(:,90:180);
figure
hist(winter_MIZ(:),101)

miz_lon = mean(effective');
figure
plot(lon,miz_lon)


%%
nofsd_k = ncread("/Volumes/NoahDay1TB/sea-ice-classification/datakmean_nofsd_2019.nc","k");
lat = ncread("/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_2019.nc","LAT");
lon = ncread("/Volumes/NoahDay1TB/sea-ice-classification/data/kmean_2019.nc","LON");

close all
figure
w = worldmap('world');
    axesm eqaazim;
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-53]);
    setm(w, 'maplonlimit', [0,360]);
    setm(w, 'grid', 'off');
    setm(w, 'frame', 'off')
    pcolorm(lat,lon,nofsd_k(:,:,180));
    colormap(Cmap)
    bedmap2('patchgl');
    %outlineashelf('all','color','k');
    %scalebar('color','k','units','km','location','se')
    scalebar('length',1000,...
        'units','km',...
        'color','k','location','sw',...
        'fontangle','italic')%,'FontSize',font_size)
    colorbar
%exportgraphics(f,'worldmap_2018-03-01.pdf', 'ContentType', 'vector','BackgroundColor','none')


%% MIZ width over years
%effective = ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/mizWidth_2018.nc','effective')./1000;
effective_2010 = nanmean(ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/mizWidth_2010.nc','sic')./1000);

effective_2011 = nanmean(ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/mizWidth_2011.nc','effective')./1000);

effective_2012 = nanmean(ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/mizWidth_2012.nc','effective')./1000);

effective_2013 = nanmean(ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/mizWidth_2013.nc','effective')./1000);

effective_2014 = nanmean(ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/mizWidth_2014.nc','effective')./1000);

effective_2015 = nanmean(ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/mizWidth_2015.nc','effective')./1000);

effective_2016 = nanmean(ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/mizWidth_2016.nc','effective')./1000);

effective_2017 = nanmean(ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/mizWidth_2017.nc','effective')./1000);

%effective_2018 = nanmean(ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/mizWidth_2018.nc','effective')./1000);

effective_2019 = nanmean(ncread('/Volumes/NoahDay1TB/sea-ice-classification/data/mizWidth_2019.nc','effective')./1000);




close all
addpath /Users/noahday/GitHub/CICE-plotting-tools/functions/
conFigure(30)
%MIZ_width_mean =nanmean(k_width(1:365,:,:),2)./1000;
plot_data = [effective_2010];
plot_data = plot_data(plot_data>0);
figure
plot(plot_data,'o')
%ylim([0,300])

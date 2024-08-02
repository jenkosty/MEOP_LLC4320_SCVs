%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EXPLANATION OF LILLY DETECTION ALGORITHM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% LLC4320 snapshot dates
snapshot_dates = {'01-Oct-2011','01-Nov-2011', '01-Dec-2011', '01-Jan-2012', '01-Feb-2012', '01-Mar-2012', '01-Apr-2012', '01-May-2012', ...
    '01-Jun-2012', '01-Jul-2012', '01-Aug-2012', '01-Sep-2012'};
date = snapshot_dates{2};

input_path = '/Volumes/Elements/LLCsealdata/Snapshot_';
load(string(input_path) + string(date) + '/LLCsealdata'); %%% LLCsealdata

LLC = cell(4,1);
sectors = {'LLC_1', 'LLC_2', 'LLC_4', 'LLC_5'};
for i = 1:4
    load(string(input_path) + string(date) + '/' + string(sectors{i}) + '/lats.mat');
    LLC{i}.lats = lats;
    load(string(input_path) + string(date) + '/' + string(sectors{i}) + '/lons.mat');
    LLC{i}.lons = lons;
    load(string(input_path) + string(date) + '/' + string(sectors{i}) + '/OW.mat');
    LLC{i}.OW = OW;
    load(string(input_path) + string(date) + '/' + string(sectors{i}) + '/vort.mat');
    LLC{i}.vort = vort;
    load(string(input_path) + string(date) + '/' + string(sectors{i}) + '/salt.mat');
    LLC{i}.salt = salt;
    load(string(input_path) + string(date) + '/' + string(sectors{i}) + '/temp.mat');
    LLC{i}.temp = temp;
    load(string(input_path) + string(date) + '/' + string(sectors{i}) + '/depth.mat');
    LLC{i}.depth = depth;
    load(string(input_path) + string(date) + '/' + string(sectors{i}) + '/polygon.mat')
    LLC{i}.polygon = polygon;
    LLC{i}.date = date;
    disp('Snapshot ' + string(i))
    clear lats lons OW vort depth polygon salt temp
end
LLC_1 = LLC{1};
LLC_2 = LLC{2};
LLC_4 = LLC{3};
LLC_5 = LLC{4};
clear LLC i sectors

load(string(input_path) + string(date) + '/lilly_data_final'); %%% Lilly Data
for tag_no = 1:length(LLCsealdata)
    LLCsealdata(tag_no).contourdata = lilly_data(tag_no).contourdata;
end


%%

tag_no = 8;
i = 112;

%%% Note: cyclone from the November 1, 2011 snapshot
casts = [LLCsealdata(tag_no).ref_ind{1,i} LLCsealdata(tag_no).ref_ind{2,i} i];
ind = min(casts):max(casts);
LLC_ind = LLCsealdata(tag_no).sector(i);
if LLC_ind == 2
    LLC = LLC_2;
end

figure()
load('AntarcticCoastline_rtopo2.mat')
load('rtopo_1080x310.mat')
axesm('stereo', 'Origin', [-90 0], 'MapLatLimit', [-90 -57]);
axis off; framem on; hold on;
plotm(cntrs_sub{1}(2,:),cntrs_sub{1}(1,:),'Color','k','LineWidth',2,'LineStyle','--');
contourm(YC, XC, coastline, [0 0], 'Fill', 'off', 'Color', 'k', 'LineWidth', 2)
plotm(LLCsealdata(tag_no).lat, LLCsealdata(tag_no).lon, 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'b');

%%

fs = 15;
no_colors = 20;
depth_ind = 33;
disp(LLC.depth(depth_ind))

f = figure('Position', [100 100 1000 900]);
tiledlayout(3,2,'TileSpacing','Compact')
%sgtitle({'LLC4320 "Seal" ' + string(LLCsealdata(tag_no).tag) + ', ' + string(date) + ' Snapshot', "  "}, 'FontWeight', 'bold', 'FontSize', 20)

%%%
lon_min = min(LLCsealdata(tag_no).lon(ind));
lon_max = max(LLCsealdata(tag_no).lon(ind));
delta_lon = lon_max - lon_min;
lon_min = lon_min - 0.05*delta_lon;
lon_max = lon_max + 0.05*delta_lon;
lat_min = min(LLCsealdata(tag_no).lat(ind));
lat_max = max(LLCsealdata(tag_no).lat(ind));
delta_lat = lat_max - lat_min;
lat_min = lat_min - 0.05*delta_lat;
lat_max = lat_max + 0.05*delta_lat;

%%% Okubo-Weiss subplot
ax1 = nexttile;
pcolor(LLC.lons, LLC.lats, LLC.OW);
shading flat
colormap(ax1, cmocean('balance', no_colors)); 
cmin = -1e-9;
cmax = 1e-9;
clim([cmin cmax])
cticks = linspace(cmin, cmax, (no_colors+4) / 4);
h = colorbar; 
h.Ticks = cticks;
h.Label.String = "Okubo-Weiss parameter [s^{-2}]"; h.Label.Rotation = 270; h.Label.VerticalAlignment = "bottom"; h.Label.FontSize = fs;
hold on
plot(LLCsealdata(tag_no).lon, LLCsealdata(tag_no).lat, 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'color', 'k')
plot(LLCsealdata(tag_no).lon(i), LLCsealdata(tag_no).lat(i), 'Marker', 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', 'color', 'k')
plot([lon_min lon_max lon_max lon_min lon_min], [lat_min lat_min lat_max lat_max lat_min], 'color', 'r', 'LineWidth', 2)
drawnow
xlabel('Longitude', 'FontSize', fs)
xmin = 0.99*min(LLCsealdata(tag_no).lon(:));
xmax = 1.01*max(LLCsealdata(tag_no).lon(:));
xlim([xmin xmax])
ylabel('Latitude', 'FontSize', fs)
ymin = 1.01*min(LLCsealdata(tag_no).lat(:));
ymax = 0.99*max(LLCsealdata(tag_no).lat(:));
ylim([ymin ymax])
pos1 = get(gca, 'Position'); % Position of subplot 1

pos1_x_arrow = interp1([xmin xmax], [pos1(1) pos1(1)+pos1(3)], [lon_min + 0.25*(lon_max-lon_min)]);
pos1_y_arrow = interp1([ymin ymax], [pos1(2) pos1(2)+pos1(4)], [lat_min+0.4]);
ax = gca; ax.FontSize = fs;

%%% Temperature subplot
ax2 = nexttile;
pcolor(LLC.lons, LLC.lats, LLC.temp(:,:,depth_ind));
colormap(ax2, cmocean('thermal', no_colors));
cmin = -1.75;
cmax = 1.75;
clim([cmin cmax])
cticks = linspace(cmin, cmax, (no_colors+4) / 4);
h = colorbar; 
h.Ticks = cticks;
h.TickLabels = round(cticks, 1);
shading flat
h = colorbar; 
h.Label.String = "Temperature [" + char(176) + "C]"; h.Label.Rotation = 270; h.Label.VerticalAlignment = "bottom"; h.Label.FontSize = fs;
hold on
plot(LLCsealdata(tag_no).lon, LLCsealdata(tag_no).lat, 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'color', 'k')
plot(LLCsealdata(tag_no).lon(i), LLCsealdata(tag_no).lat(i), 'Marker', 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', 'color', 'k')
plot([lon_min lon_max lon_max lon_min lon_min], [lat_min lat_min lat_max lat_max lat_min], 'color', 'r', 'LineWidth', 2)
drawnow
xlabel('Longitude', 'FontSize', fs)
xmin = 0.99*min(LLCsealdata(tag_no).lon(:));
xmax = 1.01*max(LLCsealdata(tag_no).lon(:));
xlim([xmin xmax])
ylabel('Latitude', 'FontSize', fs)
ymin = 1.01*min(LLCsealdata(tag_no).lat(:));
ymax = 0.99*max(LLCsealdata(tag_no).lat(:));
ylim([ymin ymax])
ax = gca; ax.FontSize = fs;

%%%
% lon_min = min(LLCsealdata(tag_no).contourdata(i).lon(:));
% lon_max = max(LLCsealdata(tag_no).contourdata(i).lon(:));
% lat_min = min(LLCsealdata(tag_no).contourdata(i).lat(:));
% lat_max = max(LLCsealdata(tag_no).contourdata(i).lat(:));

%%% Okubo-Weiss subplot (zoomed)
ax3 = nexttile;
pos2 = get(gca, 'Position'); % Position of subplot 1
pcolor(LLC.lons, LLC.lats, LLC.OW);
colormap(ax3, cmocean('balance', no_colors));
shading flat
cmin = -1e-9;
cmax = 1e-9;
clim([cmin cmax])
cticks = linspace(cmin, cmax, (no_colors+4) / 4);
h = colorbar; 
h.Ticks = cticks;
h.Label.String = "Okubo-Weiss parameter [s^{-2}]"; h.Label.Rotation = 270; h.Label.VerticalAlignment = "bottom"; h.Label.FontSize = fs;
hold on
plot(LLCsealdata(tag_no).lon(ind), LLCsealdata(tag_no).lat(ind), 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'color', 'k')
plot(LLCsealdata(tag_no).lon(i), LLCsealdata(tag_no).lat(i), 'Marker', 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', 'color', 'k')
plot(LLCsealdata(tag_no).contourdata(i).ell_lon, LLCsealdata(tag_no).contourdata(i).ell_lat, '-', 'Color', 'g', 'LineWidth', 2)
plot([lon_min lon_max lon_max lon_min lon_min], [lat_min lat_min lat_max lat_max lat_min], '-', 'color', 'r', 'LineWidth', 2)
drawnow
xlabel('Longitude', 'FontSize', fs)
xmin = 0.995*min(LLCsealdata(tag_no).lon(ind));
xmax = 1.005*max(LLCsealdata(tag_no).lon(ind));
xlim([xmin xmax])
ylabel('Latitude', 'FontSize', fs)
ymin = 1.002*min(LLCsealdata(tag_no).lat(ind));
ymax = 0.998*max(LLCsealdata(tag_no).lat(ind));
ylim([ymin ymax])
pos2_x_arrow = interp1([xmin xmax], [pos2(1) pos2(1)+pos2(3)], [lon_min + 0.7*(lon_max-lon_min)]);
pos2_y_arrow = interp1([ymin ymax], [pos2(2) pos2(2)+pos2(4)], [lat_max]);
ax = gca; ax.FontSize = fs;

%%% Temperature subplot (zoomed)
ax4 = nexttile;
pcolor(LLC.lons, LLC.lats, LLC.temp(:,:,depth_ind));
colormap(ax4, cmocean('thermal', no_colors));
shading flat
cmin = -1.75;
cmax = 1.75;
clim([cmin cmax])
cticks = linspace(cmin, cmax, (no_colors+4) / 4);
h = colorbar; 
h.Ticks = cticks;
h.TickLabels = round(cticks, 1);
h.Label.String = "Temperature [" + char(176) + "C]"; h.Label.Rotation = 270; h.Label.VerticalAlignment = "bottom"; h.Label.FontSize = fs;
hold on
plot(LLCsealdata(tag_no).lon(ind), LLCsealdata(tag_no).lat(ind), 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'color', 'k')
plot(LLCsealdata(tag_no).lon(i), LLCsealdata(tag_no).lat(i), 'Marker', 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', 'color', 'k')
plot(LLCsealdata(tag_no).contourdata(i).ell_lon, LLCsealdata(tag_no).contourdata(i).ell_lat, '-', 'Color', 'g', 'LineWidth', 2)
plot([lon_min lon_max lon_max lon_min lon_min], [lat_min lat_min lat_max lat_max lat_min], '-', 'color', 'r', 'LineWidth', 2)
drawnow
xlabel('Longitude', 'FontSize', fs)
xmin = 0.995*min(LLCsealdata(tag_no).lon(ind));
xmax = 1.005*max(LLCsealdata(tag_no).lon(ind));
xlim([xmin xmax])
ylabel('Latitude', 'FontSize', fs)
ymin = 1.002*min(LLCsealdata(tag_no).lat(ind));
ymax = 0.998*max(LLCsealdata(tag_no).lat(ind));
ylim([ymin ymax])
ax = gca; ax.FontSize = fs;

%%% Temperature timeseries
ax5 = nexttile([1 2]);
wgs84 = wgs84Ellipsoid("km");
dist_array(1) = 0;
for ii = 1:length(ind)-1
    dist_array(ii+1) = dist_array(ii) + distance(LLCsealdata(tag_no).lat(ind(ii)), LLCsealdata(tag_no).lon(ind(ii)), LLCsealdata(tag_no).lat(ind(ii+1)), LLCsealdata(tag_no).lon(ind(ii+1)), wgs84);
end
contourf(dist_array, LLCsealdata(tag_no).ps.pres(:,1), LLCsealdata(tag_no).ps.temp(:,ind), 15)
set(gca, 'YDir', 'reverse'); ylim([0 500])
colormap(ax5, cmocean('thermal', no_colors)); shading flat
cmin = -1.75;
cmax = 1.75;
clim([cmin cmax])
cticks = linspace(cmin, cmax, (no_colors+4) / 4);
h = colorbar; 
h.Ticks = cticks;
h.TickLabels = round(cticks, 1);
h.Label.String = "Temperature [" + char(176) + "C]"; h.Label.Rotation = 270; h.Label.VerticalAlignment = "bottom"; h.Label.FontSize = fs;
hold on
xline(dist_array(13), 'g-', 'LineWidth', 3)
xlabel('Distance [km]', 'FontSize', fs);
ylabel('Pressure [dbar]', 'FontSize', fs)
ax = gca; ax.FontSize = fs;

annotation('arrow',  [pos1_x_arrow(1), pos2_x_arrow(1)], [pos1_y_arrow(1), pos2_y_arrow(1)], 'Color', 'r', 'LineWidth', 10, 'HeadWidth', 30, 'HeadLength', 20);

exportgraphics(f,'Paper Figures/LLCsealtracks.png','Resolution',600)


% hold on

% annotation('line',  [pos1_x_fig(2), pos2_x_fig(2)], [pos1_y_fig(2), pos2_y_fig(2)], 'Color', 'r');
% annotation('line',  [pos1_x_fig(2), pos2_x_fig(1)], [pos1_y_fig(2), pos2_y_fig(1)], 'Color', 'r');


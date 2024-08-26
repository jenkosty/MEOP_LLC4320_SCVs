%%% Loading MEOP data
load('qc_ts.mat')

%%% Summarizing Geographic Location of MEOP Profiles
meop_lats = [qc_ts.lat];
meop_lons = [qc_ts.lon];
lats_grid = -90:1:-60; % latitude grid
lons_grid = -180:1:180; % longitude grid
density_meop_profs = NaN(length(lats_grid), length(lons_grid));
for i = 1:length(lats_grid)-1
    for j = 1:length(lons_grid)-1

        %%% Extracting lat/lon data for each grid cell
        cell_lats = [lats_grid(i) lats_grid(i) lats_grid(i+1) lats_grid(i+1) lats_grid(i)];
        cell_lons = [lons_grid(j) lons_grid(j+1) lons_grid(j+1) lons_grid(j) lons_grid(j)];

        %%% Calculating area of grid
        wgs84 = wgs84Ellipsoid("km");
        area = areaint(cell_lats, cell_lons, wgs84);

        %%% Iddentifying MEOP profiles in each cell
        density_meop_profs(i,j) = 100 * length(find(inpolygon(meop_lats, meop_lons, cell_lats, cell_lons))) / area;

    end
end
density_meop_profs(density_meop_profs == 0) = NaN;

%%
MEOPdatetime = [];
for tag_no = 1:length(qc_ts)
    MEOPdatetime = vertcat(MEOPdatetime, datetime(qc_ts(tag_no).time));
end

%%
%%% Coastline data
load('AntarcticCoastline_rtopo2.mat')
load('rtopo_1080x310.mat')

fs = 15;
no_colors = 18;
clrmap = flipud(cmocean('solar', no_colors));

%%% Creating heatmap of meop data distribution
f = figure('Position', [500 100 1200 500]);
tiledlayout(1,2,'TileSpacing','tight')

nexttile()
inds = datetime(2004,01,01):calmonths(3):datetime(2019,01,31);
histogram(MEOPdatetime,inds, 'FaceColor', [0.5 0.5 0.5])
xline(inds(1:4:end), 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
xlim([inds(1) inds(end)])
ylabel('Profiles Per Quarter', 'FontSize', fs)
xlabel('Year', 'FontSize', fs)
ax = gca; ax.FontSize = fs;

nexttile()
axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -60],'MLineLocation', 30, 'PLineLocation', 10, 'FontSize',fs);
axis off; framem on; gridm on; mlabel on; plabel on;
hold on
pcolorm(lats_grid, lons_grid, density_meop_profs)
colormap(clrmap);
cmin = min(density_meop_profs(:));
cmax = max(density_meop_profs(:));
clim([cmin cmax])
%cticks = linspace(cmin, cmax, (no_colors+3) / 3);
set(gca,'ColorScale','log')
h = colorbar; 
%h.Ticks = cticks;
%h.TickLabels = round(cticks, 0);
h.Label.String = "Profile Density (# profiles per 100 km^2)"; h.Label.Rotation = 270; h.Label.VerticalAlignment = "bottom"; h.Label.FontSize = fs;
contourm(YC, XC, coastline, [0 0], 'Fill', 'off', 'Color', 'k', 'LineWidth', 2)
plotm(cntrs_sub{1}(2,:),cntrs_sub{1}(1,:),'Color','k','LineWidth',2,'LineStyle','--');
ax = gca; ax.FontSize = fs;

exportgraphics(f,'Paper Figures/MEOPdistribution.png','Resolution',600)

%% Getting average distance between adjacent profiles

wgs84 = wgs84Ellipsoid("km");
u = 1;
for tag_no = 1:length(qc_ts)
    for i = 1:length(qc_ts(tag_no).cast)-1
        dist_array(u) = distance(qc_ts(tag_no).lat(i), qc_ts(tag_no).lon(i),qc_ts(tag_no).lat(i+1), qc_ts(tag_no).lon(i+1), wgs84);
        u = u + 1;
    end
end

disp('Median: ' + string(median(dist_array)) + ' km');
disp('Mean: ' + string(mean(dist_array)) + ' km')
disp('Min: ' + string(min(dist_array)) + ' km');
disp('Max: ' + string(max(dist_array)) + ' km');

%% Getting distribution of original un-QCed dataset

%%% Loading MEOP data
load('/Volumes/Elements/Raw Data/SealData_All')

%%% Summarizing Geographic Location of MEOP Profiles
meop_lats = [sealdata_all.LAT];
meop_lons = [sealdata_all.LON];
lats_grid = -90:1:-60; % latitude grid
lons_grid = -180:1:180; % longitude grid
density_meop_profs_og = NaN(length(lats_grid), length(lons_grid));
for i = 1:length(lats_grid)-1
    for j = 1:length(lons_grid)-1

        %%% Extracting lat/lon data for each grid cell
        cell_lats = [lats_grid(i) lats_grid(i) lats_grid(i+1) lats_grid(i+1) lats_grid(i)];
        cell_lons = [lons_grid(j) lons_grid(j+1) lons_grid(j+1) lons_grid(j) lons_grid(j)];

        %%% Calculating area of grid
        wgs84 = wgs84Ellipsoid("km");
        area = areaint(cell_lats, cell_lons, wgs84);

        %%% Iddentifying MEOP profiles in each cell
        density_meop_profs_og(i,j) = 100 * length(find(inpolygon(meop_lats, meop_lons, cell_lats, cell_lons))) / area;

    end
end
density_meop_profs_og(density_meop_profs_og == 0) = NaN;


%% Figure to compare distributions of original and QCed datasets
load('AntarcticCoastline_rtopo2.mat')
load('rtopo_1080x310.mat')

fs = 15;
no_colors = 18;
clrmap = flipud(cmocean('solar', no_colors));

tiledlayout(1,2)

%%% Original
nexttile()
axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -60],'MLineLocation', 30, 'PLineLocation', 10, 'FontSize',fs);
axis off; framem on; gridm on; mlabel on; plabel on;
hold on
pcolorm(lats_grid, lons_grid, density_meop_profs_og)
colormap(clrmap);
cmin = min(density_meop_profs_og(:));
cmax = max(density_meop_profs_og(:));
clim([cmin cmax])
set(gca,'ColorScale','log')
h = colorbar; 
h.Label.String = "Profile Density (# profiles per 100 km^2)"; h.Label.Rotation = 270; h.Label.VerticalAlignment = "bottom"; h.Label.FontSize = fs;
contourm(YC, XC, coastline, [0 0], 'Fill', 'off', 'Color', 'k', 'LineWidth', 2)
plotm(cntrs_sub{1}(2,:),cntrs_sub{1}(1,:),'Color','k','LineWidth',2,'LineStyle','--');
title('Original Dataset', 'FontSize', fs)
ax = gca; ax.FontSize = fs;

%%% QCed
nexttile()
axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -60],'MLineLocation', 30, 'PLineLocation', 10, 'FontSize',fs);
axis off; framem on; gridm on; mlabel on; plabel on;
hold on
pcolorm(lats_grid, lons_grid, density_meop_profs)
colormap(clrmap);
cmin = min(density_meop_profs(:));
cmax = max(density_meop_profs(:));
clim([cmin cmax])
set(gca,'ColorScale','log')
h = colorbar; 
h.Label.String = "Profile Density (# profiles per 100 km^2)"; h.Label.Rotation = 270; h.Label.VerticalAlignment = "bottom"; h.Label.FontSize = fs;
contourm(YC, XC, coastline, [0 0], 'Fill', 'off', 'Color', 'k', 'LineWidth', 2)
plotm(cntrs_sub{1}(2,:),cntrs_sub{1}(1,:),'Color','k','LineWidth',2,'LineStyle','--');
title('Post-QC', 'FontSize', fs)
ax = gca; ax.FontSize = fs;

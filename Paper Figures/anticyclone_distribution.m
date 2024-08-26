%% Loading Data
save_fig = 1;

%%% Coastline data
load('AntarcticCoastline_rtopo2.mat')
load('rtopo_1080x310.mat')
load('ACCfronts.mat')
front_lats = {'LatSACCF'};
front_lons = {'LonSACCF'};
input_path = '/Volumes/Elements/MEOPdata';
load(string(input_path) + '/optimized_anticyclones.mat')

LLC = LLCanticyclones_optimized;
for i = 1:length(LLC)
    if isempty(LLC(i).spice_A)
        LLC(i).spice_A = NaN;
    end
end
MEOP = MEOPanticyclones_optimized;

%% Creating Figure

no_colors = 18;
clrmap = turbo(no_colors);
MEOPshape = 's';
MEOPsize = 100;
LLCshape = 'o';
LLCsize = 50;

figure('Position', [100 100 1000 1000])
f = tiledlayout(2,2,'TileSpacing','Compact');
fs = 15;

%%% Spice Anomaly
param = 'spice_A';
nexttile()
axesm('stereo', 'Origin', [-90 0], 'MapLatLimit', [-90 -57]);
axis off; framem on; hold on;
contourm(YC, XC, coastline, [0 0], 'Fill', 'off', 'Color', 'k', 'LineWidth', 2)
for u = 1:length(front_lons)
    plotm(ACCfronts.(front_lats{u}), ACCfronts.(front_lons{u}), 'b', 'LineWidth', 2)
end
plotm(cntrs_sub{1}(2,:),cntrs_sub{1}(1,:),'Color','k','LineWidth',2,'LineStyle','--');
a = scatterm([LLC.lat], [LLC.lon], LLCsize, 'k', LLCshape, 'filled', 'DisplayName', 'LLC4320'); % Just for legend
b = scatterm([MEOP.lat], [MEOP.lon], MEOPsize, 'k', MEOPshape, 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'MEOP'); % Just for legend
scatterm([LLC.lat], [LLC.lon], LLCsize, [LLC.(param)], LLCshape, 'filled', 'DisplayName', 'LLC4320');
scatterm([MEOP.lat], [MEOP.lon], MEOPsize, [MEOP.(param)], MEOPshape, 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'MEOP');
title('Spice Anomaly', 'FontSize', fs)
colormap(clrmap);
cmin = prctile(horzcat([LLC.(param)], [MEOP.(param)]), 0);
cmax = prctile(horzcat([LLC.(param)], [MEOP.(param)]), 100);
clim([cmin cmax])
cticks = linspace(cmin, cmax, (no_colors+3) / 3);
h = colorbar; 
h.Ticks = cticks;
h.TickLabels = round(cticks, 2);
h.Label.String = "[kg/m^3]"; h.Label.Rotation = 270; h.Label.VerticalAlignment = "bottom"; h.Label.FontSize = fs;
legend([b a], 'FontSize', fs, 'Location', 'southwest' );
ax = gca; ax.FontSize = fs;

%%% DHA Anomaly
param = 'dha_A';
nexttile
axesm('stereo', 'Origin', [-90 0], 'MapLatLimit', [-90 -57]);
axis off; framem on; hold on;
contourm(YC, XC, coastline, [0 0], 'Fill', 'off', 'Color', 'k', 'LineWidth', 2)
for u = 1:length(front_lons)
    plotm(ACCfronts.(front_lats{u}), ACCfronts.(front_lons{u}), 'b', 'LineWidth', 2)
end
plotm(cntrs_sub{1}(2,:),cntrs_sub{1}(1,:),'Color','k','LineWidth',2,'LineStyle','--');
scatterm([LLC.lat], [LLC.lon], LLCsize, [LLC.(param)], LLCshape, 'filled')
scatterm([MEOP.lat], [MEOP.lon], MEOPsize, [MEOP.(param)], MEOPshape, 'filled', 'MarkerEdgeColor', 'k')
title('Dynamic Height Anomaly', 'FontSize', fs)
colormap(clrmap);
cmin = prctile(horzcat([LLC.(param)], [MEOP.(param)]), 0);
cmax = prctile(horzcat([LLC.(param)], [MEOP.(param)]), 90);
clim([cmin cmax])
cticks = linspace(cmin, cmax, (no_colors+3) / 3);
h = colorbar; 
h.Ticks = cticks;
h.TickLabels = round(cticks, 2);
h.Label.String = "[m^2/s^2]"; h.Label.Rotation = 270; h.Label.VerticalAlignment = "bottom"; h.Label.FontSize = fs;
ax = gca; ax.FontSize = fs;

%%% Core Depth
param = 'dha_P';
nexttile
axesm('stereo', 'Origin', [-90 0], 'MapLatLimit', [-90 -57]);
axis off; framem on; hold on;
contourm(YC, XC, coastline, [0 0], 'Fill', 'off', 'Color', 'k', 'LineWidth', 2)
for u = 1:length(front_lons)
    plotm(ACCfronts.(front_lats{u}), ACCfronts.(front_lons{u}), 'b', 'LineWidth', 2)
end
plotm(cntrs_sub{1}(2,:),cntrs_sub{1}(1,:),'Color','k','LineWidth',2,'LineStyle','--');
scatterm([LLC.lat], [LLC.lon], LLCsize, [LLC.(param)], LLCshape, 'filled')
scatterm([MEOP.lat], [MEOP.lon], MEOPsize, [MEOP.(param)], MEOPshape, 'filled', 'MarkerEdgeColor', 'k')
title('Core Depth', 'FontSize', fs)
colormap(clrmap);
cmin = prctile(horzcat([LLC.(param)], [MEOP.(param)]), 0);
cmax = prctile(horzcat([LLC.(param)], [MEOP.(param)]), 90);
clim([cmin cmax])
cticks = linspace(cmin, cmax, (no_colors+3) / 3);
h = colorbar; 
h.Ticks = cticks;
h.TickLabels = round(cticks, 0);
h.Label.String = "[dbar]"; h.Label.Rotation = 270; h.Label.VerticalAlignment = "bottom"; h.Label.FontSize = fs;
ax = gca; ax.FontSize = fs;

%%% Core Width
param = 'dha_Hcore';
nexttile
axesm('stereo', 'Origin', [-90 0], 'MapLatLimit', [-90 -57]);
axis off; framem on; hold on;
contourm(YC, XC, coastline, [0 0], 'Fill', 'off', 'Color', 'k', 'LineWidth', 2)
for u = 1:length(front_lons)
    plotm(ACCfronts.(front_lats{u}), ACCfronts.(front_lons{u}), 'b', 'LineWidth', 2)
end
plotm(cntrs_sub{1}(2,:),cntrs_sub{1}(1,:),'Color','k','LineWidth',2,'LineStyle','--');
scatterm([LLC.lat], [LLC.lon], LLCsize, [LLC.(param)], LLCshape, 'filled')
scatterm([MEOP.lat], [MEOP.lon], MEOPsize, [MEOP.(param)], MEOPshape, 'filled', 'MarkerEdgeColor', 'k')
title('Core Width', 'FontSize', fs)
colormap(clrmap);
cmin = prctile(horzcat([LLC.(param)], [MEOP.(param)]), 0);
cmax = prctile(horzcat([LLC.(param)], [MEOP.(param)]), 90);
clim([cmin cmax])
cticks = linspace(cmin, cmax, (no_colors+3) / 3);
h = colorbar; 
h.Ticks = cticks;
h.TickLabels = round(cticks, 0);
h.Label.String = "[dbar]"; h.Label.Rotation = 270; h.Label.VerticalAlignment = "bottom"; h.Label.FontSize = fs;
ax = gca; ax.FontSize = fs;

if save_fig == 1
    exportgraphics(f,'Paper Figures/anticyclone_map.png','Resolution',600)
end

%%
lons = [];
for tag_no = 1:467
     lons = [lons qc_ts(tag_no).lon];
end
meop_dist = histcounts(lons, -180:1:180);

lons = [];
for u = 1:length(MEOPanticyclones_optimized)
    lons = [lons MEOPanticyclones_optimized(u).lon];
end
eddy_dist = histcounts(lons, -180:1:180);

lons = [];
for u = 1:length(LLCanticyclones_optimized)
    lons = [lons LLCanticyclones_optimized(u).lon];
end
llc_eddy_dist = histcounts(lons, -180:1:180);

figure('Position', [100 100 1000 500])
fs = 15;
lw = 2;
colororder({'k','r'})
hold on
yyaxis left
plot(-179:1:180, meop_dist, 'k-', 'LineWidth', lw)
ylabel('MEOP Profiles', 'FontSize', fs)
yyaxis right
plot(-179:1:180, llc_eddy_dist, 'r-', 'LineWidth', lw)
ylabel('MEOP Cyclones', 'FontSize', fs, 'Rotation', -90)
xlim([-180 180])
xlabel('Longitude', 'FontSize', fs)
grid on
ax = gca; ax.FontSize = fs;


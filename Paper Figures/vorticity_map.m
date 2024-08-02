
%%% Coastline data
load('AntarcticCoastline_rtopo2.mat')
load('rtopo_1080x310.mat')

load('/Volumes/Elements/MEOP Optimization/cyclones/optimized.mat');
load('/Volumes/Elements/MEOP Optimization/anticyclones/optimized.mat');

ind = [LLCcyclones_optimized.True_SCV] == 1;
LLCcyclones_optimized = LLCcyclones_optimized(ind);

ind = [LLCanticyclones_optimized.True_SCV] == 1;
LLCanticyclones_optimized = LLCanticyclones_optimized(ind);


%% Creating Figure

figure('Position', [100 100 1000 1000])
tiledlayout(1,1,'TileSpacing','Compact')
%sgtitle('Minty Cyclones', 'FontSize', 15, 'FontWeight', 'bold')
fs = 15;

%%% Vorticity
nexttile
axesm('stereo', 'Origin', [-90 0], 'MapLatLimit', [-90 -57]);
axis off; framem on; hold on;
contourm(YC, XC, coastline, [0 0], 'Fill', 'off', 'Color', 'k', 'LineWidth', 2)
plotm(cntrs_sub{1}(2,:),cntrs_sub{1}(1,:),'Color','k','LineWidth',2,'LineStyle','--');
b = scatterm([LLCanticyclones_optimized.lat], [LLCanticyclones_optimized.lon], 50, abs([LLCanticyclones_optimized.vort_in_contour]), 'o', 'filled', 'DisplayName', 'Anticyclones');
a = scatterm([LLCcyclones_optimized.lat], [LLCcyclones_optimized.lon], 50, abs([LLCcyclones_optimized.vort_in_contour]), 's', 'filled','MarkerEdgeColor', 'k', 'DisplayName', 'Cyclones');
title('Depth-Averaged Absolute Vorticity', 'FontSize', fs)
h = colorbar; clim([0 0.4])
h.Label.String = "[s^{-1}]"; h.Label.Rotation = 270; h.Label.VerticalAlignment = "bottom"; h.Label.FontSize = fs;
ax = gca; ax.FontSize = fs;
legend([b a], 'FontSize', fs, 'Location', 'southwest')
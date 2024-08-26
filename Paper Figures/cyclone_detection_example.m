%%% Loading data
save_fig = 1;
input_path = '/Volumes/Elements/MEOPdata';
load(string(input_path) + '/qc_ts.mat');
load(string(input_path) + '/optimized_cyclones_deep.mat')
ts_data = qc_ts;
u = 41; %59; %58; %23;
tag_no = MEOPcyclones_optimized(u).tag_no;
i = MEOPcyclones_optimized(u).cast;
% tag_no = 182; %406; %386;
% i = 41; %109; %76;
%%
clc;
run('LLCseals_algorithm_settings.m');
ts_data = preppingTimeSeriesForDetectionAlgorithm(ts_data, depth_grid, tag_no, prms, 0);
ts_data(tag_no).spice_gauss(i) = gaussian_fit_spice(ts_data, tag_no, i, prms);
ts_data(tag_no).dha_gauss_cyclones(i) = gaussian_fit_dha_modes(ts_data, tag_no, i, 1, prms);

disp('MEOP Seal ' + string(ts_data(tag_no).tag) + ', Date: ' + string(datetime(ts_data(tag_no).time(i,:))));
casts = [ts_data(tag_no).ref_ind{1,i} ts_data(tag_no).ref_ind{2,i} i];
ind = min(casts):max(casts);
x_datenum = datenum(qc_ts(tag_no).time(ind,:));

no_colors = 20;
lw = 2;
lw_box = 5;
fs = 15;
letter_fs = 25;
main_clr = 'b';
ymin = 0;
ymax = 400;
label_placement = [0.03, 0.95];
label_placement_ts = [0.02, 0.9];
rectangle_width = 0.75;

use_subplot = 1;

%%% Boxes for reference profiles
xmin_1 = x_datenum(1); xmax_1 = x_datenum(11);
xmin_2 = x_datenum(15); xmax_2 = x_datenum(end);
ymin_box = 3; ymax_box = ymax-3;

%%% Cyclone Detection
f = figure('Position', [100 100 1200 800]);
pos = get(f, 'Position');
if use_subplot == 0
    t = tiledlayout(4,3, 'TileSpacing', 'tight');
end

%%% Sigma0 array
u = ind(1);
sigma0_arr = ts_data(tag_no).ps.sigma0(1:50:end,u);
sigma0_arr = sigma0_arr(~isnan(sigma0_arr));
sigma0_arr = round(sigma0_arr,2);
sigma0_arr = unique(sigma0_arr);
clear u

%%% pcolor settings
panel_width = 0.05;

%%% Temperature subplot
if use_subplot == 0
    ax1 = nexttile([1 2]);
else
    ax1 = subplot(4,3,[1 2]);
    pos1 = get(ax1, 'Position');
    pos1(1) = 0.1;
    pos1(3) = pos1(3) + panel_width;
    set(ax1, 'Position', pos1)
end
contourf(x_datenum, ts_data(tag_no).ps.pres(:,1), ts_data(tag_no).ps.temp(:,ind), no_colors, 'LineColor', 'none')
hold on
[C, hh] = contour(x_datenum, ts_data(tag_no).ps.pres(:,1), ts_data(tag_no).ps.sigma0(:,ind), sigma0_arr, 'k');
xline(datenum(qc_ts(tag_no).time(i,:)), 'b', 'LineWidth', lw)
letter = char('a');
text(label_placement_ts(1), label_placement_ts(2), ['(' letter ')'], 'Units', 'normalized', 'FontSize', letter_fs);
rectangle('Position', [xmin_1, 0, rectangle_width, 100], 'FaceColor', 'white')
clabel(C,hh, sigma0_arr(2:2:end-2),'FontSize', fs, 'Color', 'k', 'LabelSpacing', 700);
set(gca, 'YDir', 'reverse'); ylim([ymin ymax])
colormap(ax1, cmocean('thermal', no_colors)); shading flat
cmin = minmin(ts_data(tag_no).ps.temp(:,ind));
cmax = maxmax(ts_data(tag_no).ps.temp(:,ind));
clim([cmin cmax])
cticks = linspace(cmin, cmax, (no_colors+4) / 4);
h = colorbar; 
h.Ticks = cticks;
h.TickLabels = round(cticks, 1);
h.Label.String = "Temperature [" + char(176) + "C]"; h.Label.Rotation = 270; h.Label.VerticalAlignment = "bottom"; h.Label.FontSize = fs;
hold on
plot([xmin_1 xmax_1 xmax_1 xmin_1 xmin_1],[ymin_box ymin_box ymax_box ymax_box ymin_box], 'LineWidth', lw_box, 'Color', 'r')
plot([xmin_2 xmax_2 xmax_2 xmin_2 xmin_2],[ymin_box ymin_box ymax_box ymax_box ymin_box], 'LineWidth', lw_box, 'Color', 'r')
datetick('x', 'mm/dd'); xlim([x_datenum(1) x_datenum(end)])
ylabel('Pressure [dbar]', 'FontSize', fs)
ax = gca; ax.FontSize = fs;

%%% Salinity subplot
if use_subplot == 0
    ax2 = nexttile([1 2]);
else
    ax2 = subplot(4,3,[4 5]);
    pos2 = get(ax2, 'Position');
    pos2(1) = 0.1;
    pos2(2) = pos1(2) - pos2(4);
    pos2(3) = pos2(3) + panel_width;
    set(ax2, 'Position',pos2);
end
contourf(x_datenum, ts_data(tag_no).ps.pres(:,1), ts_data(tag_no).ps.salt(:,ind), no_colors, 'LineColor', 'none')
hold on
[C, hh] = contour(x_datenum, ts_data(tag_no).ps.pres(:,1), ts_data(tag_no).ps.sigma0(:,ind), sigma0_arr, 'k');
xline(datenum(qc_ts(tag_no).time(i,:)), 'b', 'LineWidth', lw)
letter = char('b');
text(label_placement_ts(1), label_placement_ts(2), ['(' letter ')'], 'Units', 'normalized', 'FontSize', letter_fs);
rectangle('Position', [xmin_1, 0, rectangle_width, 100], 'FaceColor', 'white')
clabel(C,hh, sigma0_arr(2:2:end-2),'FontSize', fs, 'Color', 'k', 'LabelSpacing', 700);
set(gca, 'YDir', 'reverse'); ylim([ymin ymax])
colormap(ax2, cmocean('haline', no_colors)); shading flat
cmin = minmin(ts_data(tag_no).ps.salt(:,ind));
cmax = maxmax(ts_data(tag_no).ps.salt(:,ind));
clim([cmin cmax])
cticks = linspace(cmin, cmax, (no_colors+4) / 4);
cticks = cticks(1:end-1);
h = colorbar; 
h.Ticks = cticks;
h.TickLabels = round(cticks, 1);
h.Label.String = "Salinity [psu]"; h.Label.Rotation = 270; h.Label.VerticalAlignment = "bottom"; h.Label.FontSize = fs;
hold on
plot([xmin_1 xmax_1 xmax_1 xmin_1 xmin_1],[ymin_box ymin_box ymax_box ymax_box ymin_box], 'LineWidth', lw_box, 'Color', 'r')
plot([xmin_2 xmax_2 xmax_2 xmin_2 xmin_2],[ymin_box ymin_box ymax_box ymax_box ymin_box], 'LineWidth', lw_box, 'Color', 'r')
datetick('x', 'mm/dd'); xlim([x_datenum(1) x_datenum(end)])
xlabel('Date [mm/dd]', 'FontSize', fs)
ylabel('Pressure [dbar]', 'FontSize', fs)
ax = gca; ax.FontSize = fs;

%%% Map subplot
if use_subplot == 0
    ax0 = nexttile(3, [2 1]);
else
    ax0 = subplot(4,3,[3 6]);
    pos0 = get(ax0, 'Position');
    pos0(2) = pos0(2) + 0.02;
    set(ax0, 'Position', pos0);
end
load('AntarcticCoastline_rtopo2.mat')
load('rtopo_1080x310.mat') 
axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -57]);
axis off; framem on; hold on;
contourm(YC,XC,coastline,[0 0], 'Fill', 'off', 'Color','k', 'LineWidth', 2);
plotm(cntrs_sub{1}(2,:),cntrs_sub{1}(1,:),'Color','k','LineWidth',2,'LineStyle','--'); %position of th 1000m contour (shelf break)
scatterm(ts_data(tag_no).lat(i), ts_data(tag_no).lon(i), 100, 'ko', 'MarkerFaceColor', main_clr)
ax = gca; ax.FontSize = fs;

%%% Settings for bottom 3 panels
panel_width = 0.85/3;
panel_height = 0.05;

%%% IQR Check
if use_subplot == 0
    nexttile([2 1]);
else
    ax3 = subplot(4,3,[7 10]);
    pos3 = get(ax3, 'Position');
    pos3(1) = 0.1;
    pos3(3) = panel_width;
    pos3(2) = pos3(2) + panel_height;
    set(ax3, 'Position', pos3);
end
hold on
for ii = 1:length(ind)
    a = plot(ts_data(tag_no).ds.anoms.spice(:,ind(ii)), ts_data(tag_no).ds.pres(:,i), 'Color', [0.5 0.5 0.5], 'DisplayName', 'Background Profiles');
end
letter = char('c');
text(label_placement(1), label_placement(2), ['(' letter ')'], 'Units', 'normalized', 'FontSize', letter_fs);
b = plot(ts_data(tag_no).ds.iqrs.spice_anom_lim_hi(:,i), ts_data(tag_no).ds.pres(:,i), 'k', 'LineWidth', lw, 'DisplayName', 'IQR-Based Threshold');
plot(ts_data(tag_no).ds.iqrs.spice_anom_lim_lo(:,i), ts_data(tag_no).ds.pres(:,i), 'k', 'LineWidth', lw)
c = plot(ts_data(tag_no).ds.anoms.spice(:,i), ts_data(tag_no).ds.pres(:,i), main_clr, 'LineWidth', lw, 'DisplayName','Cyclone Profile');
xline(0, '--k', 'LineWidth', 1);
set(gca, 'YDir', 'reverse')
grid on
xlabel('Spice Anomaly [kg/m^3]', 'FontSize', fs)
xticks(-0.2:0.1:0.2)
ylabel('Pressure [dbar]', 'FontSize', fs)
ylim([ymin ymax])
legend([c b a], 'Location', 'southeast', 'FontSize', fs)
ax = gca; ax.FontSize = fs;

%%% Spice Gaussian Fit
if use_subplot == 0
    ax4 = nexttile([2 1]);
else
    ax4 = subplot(4,3,[8 11]);
    pos4 = get(ax4, 'Position');
    pos4(1) = pos3(1) + pos3(3);
    pos4(3) = panel_width;
    pos4(2) = pos4(2) + panel_height;
    set(ax4, 'Position', pos4)
end
hold on
b = plot(ts_data(tag_no).ds.anoms.spice(:,i), ts_data(tag_no).ds.pres(:,i), main_clr, 'DisplayName', 'Cyclone Profile', 'LineWidth', lw);
if ~isempty(ts_data(tag_no).spice_gauss) && (length(ts_data(tag_no).spice_gauss) >= i) && ~isempty(ts_data(tag_no).spice_gauss(i).rejected)
    a = plot(ts_data(tag_no).spice_gauss(i).X,ts_data(tag_no).spice_gauss(i).Y,'--','Color', 'm', 'LineWidth',lw, 'DisplayName','Gaussian Fit');
end
letter = char('d');
text(label_placement(1), label_placement(2), ['(' letter ')'], 'Units', 'normalized', 'FontSize', letter_fs);
xline(0, '--k', 'LineWidth', 1);
hold off
xlabel('Spice Anomaly [kg/m^3]', 'FontSize', fs);
xlim([min(ts_data(tag_no).ds.anoms.spice(:,i)) max(ts_data(tag_no).ds.anoms.spice(:,i))])
xticks(-0.15:0.05:0.05)
set(gca, 'YDir', 'reverse');
set(gca, 'yticklabel', []);
ylim([ymin ymax]);
grid on
legend([b a], 'Location', 'southeast', 'FontSize', fs)
ax = gca; ax.FontSize = fs;

%%% DHA Gaussian Fit
if use_subplot == 0
    nexttile([2 1])
else
    ax5 = subplot(4,3,[9 12]);
    pos5 = get(ax5, 'Position');
    pos5(1) = pos4(1) + pos4(3);
    pos5(3) = panel_width;
    pos5(2) = pos5(2) + panel_height;
    set(ax5, 'Position', pos5)
end
hold on
xline(0, '--k', 'LineWidth', 1);
a = plot(ts_data(tag_no).ps.anoms.dyn_height_anom(:,i), ts_data(tag_no).ps.pres(:,i),main_clr, 'DisplayName', 'Cyclone Profile','LineWidth',lw);
if ~isempty(ts_data(tag_no).dha_gauss_cyclones) && (length(ts_data(tag_no).dha_gauss_cyclones) >= i) && ~isempty(ts_data(tag_no).dha_gauss_cyclones(i).rejected)
    b = plot(ts_data(tag_no).dha_gauss_cyclones(i).dataX_orig,ts_data(tag_no).dha_gauss_cyclones(i).dataY_orig,'Color', [0.5 0.5 0.5],'LineWidth',lw, 'DisplayName', '1st BC Mode Removed');
    c = plot(ts_data(tag_no).dha_gauss_cyclones(i).X,ts_data(tag_no).dha_gauss_cyclones(i).Y,'--','Color', 'm', 'LineWidth',lw, 'DisplayName', 'Gaussian Fit');
    legend([a b c], 'Location', 'southeast', 'FontSize', fs)
end
letter = char('e');
text(label_placement(1), label_placement(2), ['(' letter ')'], 'Units', 'normalized', 'FontSize', letter_fs);
hold off
xlabel('Dynamic Height Anomaly [m^2/s^2]', 'FontSize', fs);
xticks(-0.02:0.02:0.06)
set(gca, 'YDir', 'reverse');
set(gca, 'yticklabel', []);
ylim([ymin ymax]);
grid on
%title('Dynamic Height Anomaly Gaussian Fit', 'FontSize', fs, 'FontWeight', 'Bold')
ax = gca; ax.FontSize = fs;

if save_fig == 1
    exportgraphics(f,'Paper Figures/cyclone_detection_example.png','Resolution',600)
end

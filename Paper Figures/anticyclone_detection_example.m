%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ANTICYCLONE DETECTION EXAMPLE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Loading data
save_figs = 0;
input_path = '/Volumes/Elements/MEOPdata';
load(string(input_path) + '/qc_ts.mat');
ts_data = qc_ts;
load(string(input_path) + '/optimized_anticyclones_shallow.mat')
% tag_no = 115;
% i = 154;
u = 6; %8 18 26
tag_no = MEOPanticyclones_optimized(u).tag_no;
i = MEOPanticyclones_optimized(u).cast;

run('LLCseals_algorithm_settings.m');
ts_data = preppingTimeSeriesForDetectionAlgorithm(ts_data, depth_grid, tag_no, prms, 0);
ts_data(tag_no).isa_gauss(i) = gaussian_fit_isa(ts_data, tag_no, i, prms);
ts_data(tag_no).spice_gauss(i) = gaussian_fit_spice(ts_data, tag_no, i, prms);
ts_data(tag_no).dha_gauss_anticyclones(i) = gaussian_fit_dha_modes(ts_data, tag_no, i, 0, prms);

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
ymax = 500;
label_placement = [0.85, 0.95];
label_placement_ts = [0.02, 0.9];
rectangle_width = 0.75;

use_subplot = 1;

%%% Boxes for reference profiles
xmin_1 = x_datenum(1); xmax_1 = x_datenum(11);
xmin_2 = x_datenum(15); xmax_2 = x_datenum(end);
ymin_box = 3; ymax_box = ymax-3;

%%% Anticyclone Detection
f = figure('Position', [100 100 1200 1000]);
pos = get(f, 'Position');
if use_subplot == 0
    tiledlayout(4,3, 'TileSpacing', 'tight');
end

%%% Sigma0 array
u = ind(1);
sigma0_arr = ts_data(tag_no).ps.sigma0(1:25:end,u);
sigma0_arr = sigma0_arr(~isnan(sigma0_arr));
sigma0_arr = round(sigma0_arr,2);
sigma0_arr = unique(sigma0_arr);
clear u

%%% pcolor settings
panel_width = 0.05;

%%% Temperature subplot
if use_subplot == 0
    ax1 = nexttile([1 3]);
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
clabel(C,hh, sigma0_arr(2:2:end-2),'FontSize', fs, 'Color', 'k', 'LabelSpacing', 700);
letter = char('a');
text(label_placement_ts(1), label_placement_ts(2), ['(' letter ')'], 'Units', 'normalized', 'FontSize', letter_fs);
rectangle('Position', [xmin_1, 0, rectangle_width, 100], 'FaceColor', 'white')
xline(datenum(qc_ts(tag_no).time(i,:)), 'b', 'LineWidth', lw)
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
xlabel('Date (mm/dd)', 'FontSize', fs)
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
letter = char('b');
text(label_placement_ts(1), label_placement_ts(2), ['(' letter ')'], 'Units', 'normalized', 'FontSize', letter_fs);
rectangle('Position', [xmin_1, 0, rectangle_width, 100], 'FaceColor', 'white')
clabel(C,hh, sigma0_arr(2:2:end-2),'FontSize', fs, 'Color', 'k', 'LabelSpacing', 700);
xline(datenum(qc_ts(tag_no).time(i,:)), 'b', 'LineWidth', lw)
set(gca, 'YDir', 'reverse'); ylim([ymin ymax])
colormap(ax2, cmocean('haline', no_colors)); colorbar; shading flat
cmin = minmin(ts_data(tag_no).ps.salt(:,ind));
cmax = maxmax(ts_data(tag_no).ps.salt(:,ind));
clim([cmin cmax])
cticks = linspace(cmin, cmax, (no_colors+4) / 4);
cticks = cticks(1:end-1);
h = colorbar; 
h.Ticks = cticks;
h.TickLabels = round(cticks, 1);
h.Label.String = "Salinity [psu]"; h.Label.Rotation = 270; h.Label.VerticalAlignment = "bottom"; h.Label.FontSize = fs;
plot([xmin_1 xmax_1 xmax_1 xmin_1 xmin_1],[ymin_box ymin_box ymax_box ymax_box ymin_box], 'LineWidth', lw_box, 'Color', 'r')
plot([xmin_2 xmax_2 xmax_2 xmin_2 xmin_2],[ymin_box ymin_box ymax_box ymax_box ymin_box], 'LineWidth', lw_box, 'Color', 'r')
datetick('x', 'mm/dd'); xlim([x_datenum(1) x_datenum(end)])
ylabel('Pressure [dbar]', 'FontSize', fs)
xlabel('Date [mm/dd]', 'FontSize', fs)
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
    ax3 = subplot(4,4,[9 13]);
    pos3 = get(ax3, 'Position');
    pos3(1) = 0.1;
    pos3(3) = panel_width;
    pos3(2) = pos3(2) + panel_height;
    set(ax3, 'Position', pos3);
end
hold on
for ii = 1:length(ind)
    if ind(ii) ~= i
        a = plot(ts_data(tag_no).ds.anoms.isopycnal_separation_normalized(:,ind(ii)), ts_data(tag_no).ds.pres(:,i), 'Color', [0.5 0.5 0.5], 'DisplayName', 'Background Profiles');
    end
end
letter = char('c');
text(label_placement(1), label_placement(2), ['(' letter ')'], 'Units', 'normalized', 'FontSize', letter_fs);
b = plot(ts_data(tag_no).ds.iqrs.isopycnal_separation_anom_normalized_lim_hi(:,i), ts_data(tag_no).ds.pres(:,i), 'k', 'LineWidth', 2, 'DisplayName', 'IQR-Based Threshold');
plot(ts_data(tag_no).ds.iqrs.isopycnal_separation_anom_normalized_lim_lo(:,i), ts_data(tag_no).ds.pres(:,i), 'k', 'LineWidth', 2)
c = plot(ts_data(tag_no).ds.anoms.isopycnal_separation_normalized(:,i), ts_data(tag_no).ds.pres(:,i), main_clr, 'LineWidth', 2, 'DisplayName', 'Anticyclonic SCV');
set(gca, 'YDir', 'reverse')
grid on
xlabel('NISA [-]', 'FontSize', fs)
ylabel('Pressure [dbar]', 'FontSize', fs)
ylim([ymin ymax])
xlim([min(ts_data(tag_no).ds.iqrs.isopycnal_separation_anom_normalized_lim_lo(:,i))-1 max(ts_data(tag_no).ds.iqrs.isopycnal_separation_anom_normalized_lim_hi(:,i))+20])
legend([c b a], 'Location', 'southeast', 'FontSize', fs);
ax = gca; ax.FontSize = fs;

%%% ISA Gaussian Fit
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
letter = char('d');
text(label_placement(1), label_placement(2), ['(' letter ')'], 'Units', 'normalized', 'FontSize', letter_fs);
plot(ts_data(tag_no).ds.anoms.isopycnal_separation_normalized(:,i), ts_data(tag_no).ds.pres(:,i), main_clr, 'DisplayName', 'Profile', 'LineWidth', lw)
if ~isempty(ts_data(tag_no).isa_gauss) && (length(ts_data(tag_no).isa_gauss) >= i) && ~isempty(ts_data(tag_no).isa_gauss(i).rejected)
    a = plot(ts_data(tag_no).isa_gauss(i).X,ts_data(tag_no).isa_gauss(i).Y,'--','Color', 'm', 'LineWidth',lw, 'DisplayName','Gaussian Fit');
end
xline(0, '--k', 'LineWidth', 1);
hold off
xlabel('NISA [-]', 'FontSize', fs);
set(gca, 'yticklabel', []);
set(gca, 'YDir', 'reverse');
ylim([ymin ymax]);
grid on
legend([a], 'Location', 'southeast', 'FontSize', fs)
ax = gca; ax.FontSize = fs;

% %%% Spice Gaussian Fit
% if use_subplot == 0
%     nexttile([2 1])
% else
%     ax5 = subplot(4,4,[11 15]);
%     pos5 = get(ax5, 'Position');
%     pos5(1) = pos4(1) + pos4(3);
%     pos5(3) = panel_width;
%     pos5(2) = pos5(2) + panel_height;
%     set(ax5, 'Position', pos5)
% end
% hold on
% letter = char('c');
% text(label_placement(1), label_placement(2), ['(' letter ')'], 'Units', 'normalized', 'FontSize', letter_fs);
% plot(ts_data(tag_no).ds.anoms.spice(:,i), ts_data(tag_no).ds.pres(:,i), main_clr, 'DisplayName', 'Profile', 'LineWidth', lw)
% if ~isempty(ts_data(tag_no).spice_gauss) && (length(ts_data(tag_no).spice_gauss) >= i) && ~isempty(ts_data(tag_no).isa_gauss(i).rejected)
%     a = plot(ts_data(tag_no).spice_gauss(i).X,ts_data(tag_no).spice_gauss(i).Y,'--','Color', 'm', 'LineWidth',lw, 'DisplayName','Gaussian Fit');
% end
% xline(0, '--k', 'LineWidth', 1);
% hold off
% xlabel('Spice Anomaly [kg/m^3]', 'FontSize', fs);
% set(gca, 'yticklabel', []);
% set(gca, 'YDir', 'reverse');
% ylim([ymin ymax]);
% grid on
% legend([a], 'Location', 'southeast', 'FontSize', fs)
% ax = gca; ax.FontSize = fs;


%%% DHA Gaussian Fit
if use_subplot == 0
    nexttile([2 1])
else
    ax6 = subplot(4,3,[9 12]);
    pos6 = get(ax6, 'Position');
    pos6(1) = pos4(1) + pos4(3);
    pos6(3) = panel_width;
    pos6(2) = pos6(2) + panel_height;
    set(ax6, 'Position', pos6)
end
hold on
letter = char('e');
text(label_placement(1), label_placement(2), ['(' letter ')'], 'Units', 'normalized', 'FontSize', letter_fs);
xline(0, '--k', 'LineWidth', 1);
a = plot(ts_data(tag_no).ps.anoms.dyn_height_anom(:,i), ts_data(tag_no).ps.pres(:,i),main_clr, 'DisplayName', 'Anticyclone Profile','LineWidth',lw);
if ~isempty(ts_data(tag_no).dha_gauss_anticyclones) && (length(ts_data(tag_no).dha_gauss_anticyclones) >= i) && ~isempty(ts_data(tag_no).dha_gauss_anticyclones(i).rejected)
    b = plot(ts_data(tag_no).dha_gauss_anticyclones(i).dataX_orig,ts_data(tag_no).dha_gauss_anticyclones(i).dataY_orig,'Color', [0.5 0.5 0.5],'LineWidth',lw, 'DisplayName', '1st BC Mode Removed');
    c = plot(ts_data(tag_no).dha_gauss_anticyclones(i).X,ts_data(tag_no).dha_gauss_anticyclones(i).Y,'--','Color', 'm', 'LineWidth',lw, 'DisplayName', 'Gaussian Fit');
    legend([a b c], 'Location', 'southeast', 'FontSize', fs)
end
hold off
xlabel('Dynamic Height Anomaly [m^2/s^2]', 'FontSize', fs);
set(gca, 'yticklabel', []);
set(gca, 'YDir', 'reverse');
ylim([ymin ymax]);
grid on
ax = gca; ax.FontSize = fs;

if save_figs == 1
    exportgraphics(f,'Paper Figures/anticyclone_detection_example.png','Resolution',600)
end

%% Loading data
save_fig = 1;
input_path = '/Volumes/Elements/MEOPdata';
load(string(input_path) + '/optimized_anticyclones.mat')
LLC = LLCanticyclones_optimized;
MEOP = MEOPanticyclones_optimized;

load('qc_ts.mat')

%%

MEOPclr = [208,28,139] ./ 255;
LLCclr = [77,172,38] ./ 255;
alpha = 0.2;
lw = 4;

figure('Position', [100 100 1000 500])
f = tiledlayout(1,3, 'TileSpacing', 'compact');
fs = 15;

%%%%%%%%%%%%%%%%%%%%%
%%% Spice Anomaly %%%
%%%%%%%%%%%%%%%%%%%%%

nexttile;
grid on
hold on

%%% MEOP
MEOPdata = NaN(length(depth_grid), length(MEOP));
for ii = 1:length(MEOP)
    %plot(MEOP(ii).spice_anom_profile_fitted, MEOP(ii).spice_pres_profile_fitted, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    ind = ~isnan(MEOP(ii).spice_anom_profile_fitted) & ~isnan(MEOP(ii).spice_pres_profile_fitted);
    MEOPdata(:,ii) = interp1(MEOP(ii).spice_pres_profile_fitted(ind), MEOP(ii).spice_anom_profile_fitted(ind), depth_grid);
end
mean_profile = mean(MEOPdata, 2, 'omitnan');
rmse_profile = rmse(MEOPdata, mean_profile, 2, 'omitnan');
%plot(mean_profile - rmse_profile, depth_grid, 'r--', 'LineWidth', 3, 'DisplayName', 'MEOP');
%plot(mean_profile + rmse_profile, depth_grid, 'r--', 'LineWidth', 3, 'DisplayName', 'MEOP');
x_fill = vertcat(mean_profile+rmse_profile, flipud(mean_profile-rmse_profile)); 
y_fill = vertcat(depth_grid, flipud(depth_grid)); 
fill(x_fill, y_fill, MEOPclr, 'FaceAlpha', alpha, 'LineStyle', 'none');
a = plot(mean(MEOPdata(:,strcmp([MEOP.type], "deep")), 2, 'omitnan'), depth_grid, 'Color', MEOPclr, 'LineWidth', lw, 'DisplayName', 'MEOP');
plot(mean(MEOPdata(:,strcmp([MEOP.type], "shallow")), 2, 'omitnan'), depth_grid, '--', 'Color', MEOPclr, 'LineWidth', lw, 'DisplayName', 'MEOP');

%%% LLC
LLCdata = NaN(length(depth_grid), length(LLC));
for ii = 1:length(LLC)
    ind = ~isnan(LLC(ii).spice_anom_profile_fitted) & ~isnan(LLC(ii).spice_pres_profile_fitted);
    if isempty(ind)
        continue
    end
    LLCdata(:,ii) = interp1(LLC(ii).spice_pres_profile_fitted(ind), LLC(ii).spice_anom_profile_fitted(ind), depth_grid);
end
mean_profile = mean(LLCdata, 2, 'omitnan');
rmse_profile = rmse(LLCdata, mean_profile, 2, 'omitnan');
%plot(mean_profile - rmse_profile, depth_grid, 'b--', 'LineWidth', 3, 'DisplayName', 'MEOP');
%plot(mean_profile + rmse_profile, depth_grid, 'b--', 'LineWidth', 3, 'DisplayName', 'MEOP');
x_fill = vertcat(mean_profile+rmse_profile, flipud(mean_profile-rmse_profile)); 
y_fill = vertcat(depth_grid, flipud(depth_grid)); 
fill(x_fill, y_fill, LLCclr, 'FaceAlpha', alpha, 'LineStyle', 'none');
b = plot(mean(LLCdata(:,strcmp([LLC.type], "deep")), 2, 'omitnan'), depth_grid, 'Color', LLCclr, 'LineWidth', lw, 'DisplayName', 'LLC4320');
plot(mean(LLCdata(:,strcmp([LLC.type], "shallow")), 2, 'omitnan'), depth_grid, '--', 'Color', LLCclr, 'LineWidth', lw, 'DisplayName', 'LLC4320');

xline(0, '--k', 'LineWidth', 2)
set(gca, 'YDir', 'reverse')
ylim([0 600]);
ylabel('Pressure [dbar]', 'FontSize', fs)
xlabel('Spice Anomaly [kg/m^3]', 'FontSize', fs)
ax = gca; ax.FontSize = fs;
legend([a b], 'Location', 'southwest', 'FontSize', fs)
clear LLCdata MEOPdata mean_profile rmse_profile a b ind

%%%%%%%%%%%
%%% DHA %%%
%%%%%%%%%%%
nexttile;
grid on
hold on

%%% MEOP
MEOPdata = NaN(length(depth_grid), length(MEOP));
for ii = 1:length(MEOP)
    %plot(MEOP(ii).dha_anom_profile_fitted, MEOP(ii).dha_pres_profile_fitted, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.1);
    ind = ~isnan(MEOP(ii).dha_anom_profile_fitted) & ~isnan(MEOP(ii).dha_pres_profile_fitted);
    MEOPdata(:,ii) = interp1(MEOP(ii).dha_pres_profile_fitted(ind), MEOP(ii).dha_anom_profile_fitted(ind), depth_grid);
end
mean_profile = mean(MEOPdata, 2, 'omitnan');
rmse_profile = rmse(MEOPdata, mean_profile, 2, 'omitnan');
%plot(mean_profile - rmse_profile, depth_grid, 'r--', 'LineWidth', 3, 'DisplayName', 'MEOP');
%plot(mean_profile + rmse_profile, depth_grid, 'r--', 'LineWidth', 3, 'DisplayName', 'MEOP');
x_fill = vertcat(mean_profile+rmse_profile, flipud(mean_profile-rmse_profile)); 
y_fill = vertcat(depth_grid, flipud(depth_grid)); 
fill(x_fill, y_fill, MEOPclr, 'FaceAlpha', alpha, 'LineStyle', 'none'); 
a = plot(mean(MEOPdata(:,strcmp([MEOP.type], "deep")), 2, 'omitnan'), depth_grid, 'Color', MEOPclr, 'LineWidth', lw, 'DisplayName', 'MEOP');
plot(mean(MEOPdata(:,strcmp([MEOP.type], "shallow")), 2, 'omitnan'), depth_grid, '--', 'Color', MEOPclr, 'LineWidth', lw, 'DisplayName', 'MEOP');

%%% LLC
LLCdata = NaN(length(depth_grid), length(LLC));
for ii = 1:length(LLC)
    %plot(LLC(ii).dha_anom_profile_fitted, LLC(ii).dha_pres_profile_fitted, 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
    ind = ~isnan(LLC(ii).dha_anom_profile_fitted) & ~isnan(LLC(ii).dha_pres_profile_fitted);
    LLCdata(:,ii) = interp1(LLC(ii).dha_pres_profile_fitted(ind), LLC(ii).dha_anom_profile_fitted(ind), depth_grid);
end
mean_profile = mean(LLCdata, 2, 'omitnan');
rmse_profile = rmse(LLCdata, mean_profile, 2, 'omitnan');
%plot(mean_profile - rmse_profile, depth_grid, 'b--', 'LineWidth', 3, 'DisplayName', 'MEOP');
%plot(mean_profile + rmse_profile, depth_grid, 'b--', 'LineWidth', 3, 'DisplayName', 'MEOP');
x_fill = vertcat(mean_profile+rmse_profile, flipud(mean_profile-rmse_profile)); 
y_fill = vertcat(depth_grid, flipud(depth_grid)); 
fill(x_fill, y_fill, LLCclr, 'FaceAlpha', alpha, 'LineStyle', 'none'); 
b = plot(mean(LLCdata(:,strcmp([LLC.type], "deep")), 2, 'omitnan'), depth_grid, 'Color', LLCclr, 'LineWidth', lw, 'DisplayName', 'LLC4320');
plot(mean(LLCdata(:,strcmp([LLC.type], "shallow")), 2, 'omitnan'), depth_grid, '--', 'Color', LLCclr,'LineWidth', lw, 'DisplayName', 'LLC4320');

xline(0, '--k', 'LineWidth', 2)
set(gca, 'YDir', 'reverse')
ylim([0 600]);
xlabel('Dynamic Height Anomaly [m^2/s^2]', 'FontSize', fs)
ax = gca; ax.FontSize = fs;

%%%%%%%%%%%
%%% ISA %%%
%%%%%%%%%%%

nexttile;
grid on
hold on

%%% MEOP
MEOPdata = NaN(length(depth_grid), length(MEOP));
for ii = 1:length(MEOP)
    %plot(MEOP(ii).isa_anom_profile_fitted, MEOP(ii).isa_pres_profile_fitted, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    ind = ~isnan(MEOP(ii).isa_anom_profile_fitted) & ~isnan(MEOP(ii).isa_pres_profile_fitted);
    MEOPdata(:,ii) = interp1(MEOP(ii).isa_pres_profile_fitted(ind), MEOP(ii).isa_anom_profile_fitted(ind), depth_grid);
end
mean_profile = mean(MEOPdata, 2, 'omitnan');
rmse_profile = rmse(MEOPdata, mean_profile, 2, 'omitnan');
%plot(mean_profile - rmse_profile, depth_grid, 'r--', 'LineWidth', 3, 'DisplayName', 'MEOP');
%plot(mean_profile + rmse_profile, depth_grid, 'r--', 'LineWidth', 3, 'DisplayName', 'MEOP');
x_fill = vertcat(mean_profile+rmse_profile, flipud(mean_profile-rmse_profile)); 
y_fill = vertcat(depth_grid, flipud(depth_grid)); 
fill(x_fill, y_fill, MEOPclr, 'FaceAlpha', alpha, 'LineStyle', 'none'); 
a = plot(mean(MEOPdata(:,strcmp([MEOP.type], "deep")), 2, 'omitnan'), depth_grid, 'Color', MEOPclr, 'LineWidth', lw, 'DisplayName', 'MEOP');
plot(mean(MEOPdata(:,strcmp([MEOP.type], "shallow")), 2, 'omitnan'), depth_grid, '--', 'Color', MEOPclr,'LineWidth', lw, 'DisplayName', 'MEOP');

%%% LLC
LLCdata = NaN(length(depth_grid), length(LLC));
for ii = 1:length(LLC)
        %plot(LLC(ii).isa_anom_profile_fitted, LLC(ii).isa_pres_profile_fitted, 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
        ind = ~isnan(LLC(ii).isa_anom_profile_fitted) & ~isnan(LLC(ii).isa_pres_profile_fitted);
        LLCdata(:,ii) = interp1(LLC(ii).isa_pres_profile_fitted(ind), LLC(ii).isa_anom_profile_fitted(ind), depth_grid);
end
mean_profile = mean(LLCdata, 2, 'omitnan');
rmse_profile = rmse(LLCdata, mean_profile, 2, 'omitnan');
%plot(mean_profile - rmse_profile, depth_grid, 'b--', 'LineWidth', 3, 'DisplayName', 'MEOP');
%plot(mean_profile + rmse_profile, depth_grid, 'b--', 'LineWidth', 3, 'DisplayName', 'MEOP');
x_fill = vertcat(mean_profile+rmse_profile, flipud(mean_profile-rmse_profile)); 
y_fill = vertcat(depth_grid, flipud(depth_grid)); 
fill(x_fill, y_fill, LLCclr, 'FaceAlpha', alpha, 'LineStyle', 'none'); 
b = plot(mean(LLCdata(:,strcmp([LLC.type], "deep")), 2, 'omitnan'), depth_grid, 'Color', LLCclr, 'LineWidth', lw, 'DisplayName', 'LLC4320');
plot(mean(LLCdata(:,strcmp([LLC.type], "shallow")), 2, 'omitnan'), depth_grid, '--', 'Color', LLCclr, 'LineWidth', lw, 'DisplayName', 'LLC4320');
xlim([0 20])

xline(0, '--k', 'LineWidth', 2)
set(gca, 'YDir', 'reverse')
ylim([0 600]);
xlabel('Norm. Isopycnal Separation Anomaly [-]', 'FontSize', fs)
ax = gca; ax.FontSize = fs;
clear LLCisa ind ii

if save_fig == 1
    exportgraphics(f,'Paper Figures/anticyclone_profiles.png','Resolution',600)
end


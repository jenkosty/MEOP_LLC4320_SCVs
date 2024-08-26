%%% Loading MEOP data
input_path = '/Volumes/Elements/MEOPdata';
load(string(input_path) + '/qc_ts_full.mat')
depth_threshold = 1000;

%%
%%% Getting N2 data from MEOP
u = 1;
for tag_no = 1:length(qc_ts)
    for i = 1:length(qc_ts(tag_no).cast)
        N2(:,u) = qc_ts(tag_no).ps.N2(:,i);
        N2ref(:,u) = qc_ts(tag_no).ps.ref_N2(:,i);
        IS(:,u) = qc_ts(tag_no).ds.isopycnal_separation(:,i);
        NISA(:,u) = qc_ts(tag_no).ds.anoms.isopycnal_separation_normalized(:,i);
        DH(:,u) = qc_ts(tag_no).ps.dyn_height_anom(:,i);
        DHA(:,u) = qc_ts(tag_no).ps.anoms.dyn_height_anom(:,i);
        bathymetry(u) = qc_ts(tag_no).bathymetry(i);
        u = u + 1;
    end
end

%%% N2
MEOP_N2{1,1} = 'MEOP';
MEOP_N2{2,1} = mean(N2, 2, 'omitnan');
MEOP_N2{3,1} = mean(N2(:,abs(bathymetry) < depth_threshold), 2, 'omitnan');
MEOP_N2{4,1} = mean(N2(:,abs(bathymetry) >= depth_threshold), 2, 'omitnan');

%%% Isopycnal Separation
MEOP_IS{1,1} = 'MEOP';
MEOP_IS{2,1} = mean(IS, 2, 'omitnan');
MEOP_IS{3,1} = mean(IS(:,abs(bathymetry) < depth_threshold), 2, 'omitnan');
MEOP_IS{4,1} = mean(IS(:,abs(bathymetry) >= depth_threshold), 2, 'omitnan');

%%% Normalized isopycnal separation anomaly
MEOP_NISA{1,1} = 'MEOP';
MEOP_NISA{2,1} = mean(NISA, 2, 'omitnan');
MEOP_NISA{3,1} = mean(NISA(:,abs(bathymetry) < depth_threshold), 2, 'omitnan');
MEOP_NISA{4,1} = mean(NISA(:,abs(bathymetry) >= depth_threshold), 2, 'omitnan');

%%% Dynamic height
MEOP_DH{1,1} = 'MEOP';
MEOP_DH{2,1} = mean(DH, 2, 'omitnan');
MEOP_DH{3,1} = mean(DH(:,abs(bathymetry) < depth_threshold), 2, 'omitnan');
MEOP_DH{4,1} = mean(DH(:,abs(bathymetry) >= depth_threshold), 2, 'omitnan');

%%% Dynamic height anomaly
MEOP_DHA{1,1} = 'MEOP';
MEOP_DHA{2,1} = mean(DHA, 2, 'omitnan');
MEOP_DHA{3,1} = mean(DHA(:,abs(bathymetry) < depth_threshold), 2, 'omitnan');
MEOP_DHA{4,1} = mean(DHA(:,abs(bathymetry) >= depth_threshold), 2, 'omitnan');

clear bathymetry N2 N2ref IS NISA DH DHA

%%
%%% Loading LLC data
snapshot_dates = {'01-Oct-2011','01-Nov-2011', '01-Dec-2011', '01-Jan-2012', '01-Feb-2012', '01-Mar-2012', '01-Apr-2012', '01-May-2012', ...
    '01-Jun-2012', '01-Jul-2012', '01-Aug-2012', '01-Sep-2012'};
input_path = '/Volumes/Elements/LLCsealdata'

for uu = 1:length(snapshot_dates)
    date = snapshot_dates{uu};
    disp(date);

    %%% Loading LLC seal track data
    load(string(input_path) + '/Snapshot_' + string(date) + '/LLCsealdata_full.mat')

    %%% Extracting N2 data
    u = 1;
    for tag_no = 1:length(LLCsealdata)
        for i = 1:length(LLCsealdata(tag_no).cast)
            N2(:,u) = LLCsealdata(tag_no).ps.N2(:,i);
            N2ref(:,u) = LLCsealdata(tag_no).ps.ref_N2(:,i);
            IS(:,u) = LLCsealdata(tag_no).ds.isopycnal_separation(:,i);
            NISA(:,u) = LLCsealdata(tag_no).ds.anoms.isopycnal_separation_normalized(:,i);
            DH(:,u) = LLCsealdata(tag_no).ps.dyn_height_anom(:,i);
            DHA(:,u) = LLCsealdata(tag_no).ps.anoms.dyn_height_anom(:,i);
            bathymetry(u) = LLCsealdata(tag_no).bathymetry(i);
            u = u + 1;
        end
    end

    LLC_N2{1,uu} = date;
    LLC_N2{2,uu} = mean(N2, 2, 'omitnan');
    LLC_N2{3,uu} = mean(N2(:,abs(bathymetry) < 1000), 2, 'omitnan');
    LLC_N2{4,uu} = mean(N2(:,abs(bathymetry) >= 1000), 2, 'omitnan');

    %%% Isopycnal Separation
    LLC_IS{1,uu} = date;
    LLC_IS{2,uu} = mean(IS, 2, 'omitnan');
    LLC_IS{3,uu} = mean(IS(:,abs(bathymetry) < depth_threshold), 2, 'omitnan');
    LLC_IS{4,uu} = mean(IS(:,abs(bathymetry) >= depth_threshold), 2, 'omitnan');

    %%% Normalized isopycnal separation anomaly
    LLC_NISA{1,uu} = date;
    LLC_NISA{2,uu} = mean(NISA, 2, 'omitnan');
    LLC_NISA{3,uu} = mean(NISA(:,abs(bathymetry) < depth_threshold), 2, 'omitnan');
    LLC_NISA{4,uu} = mean(NISA(:,abs(bathymetry) >= depth_threshold), 2, 'omitnan');

    %%% Dynamic height
    LLC_DH{1,uu} = date;
    LLC_DH{2,uu} = mean(DH, 2, 'omitnan');
    LLC_DH{3,uu} = mean(DH(:,abs(bathymetry) < depth_threshold), 2, 'omitnan');
    LLC_DH{4,uu} = mean(DH(:,abs(bathymetry) >= depth_threshold), 2, 'omitnan');

    %%% Dynamic height anomaly
    LLC_DHA{1,uu} = date;
    LLC_DHA{2,uu} = mean(DHA, 2, 'omitnan');
    LLC_DHA{3,uu} = mean(DHA(:,abs(bathymetry) < depth_threshold), 2, 'omitnan');
    LLC_DHA{4,uu} = mean(DHA(:,abs(bathymetry) >= depth_threshold), 2, 'omitnan');

    clear LLCsealdata N2 IS NISA DH DHA

end

save('Paper Figures/background_profiles.mat', 'MEOP_N2', 'MEOP_IS', 'MEOP_NISA', 'MEOP_DH', 'MEOP_DHA', 'LLC_N2', 'LLC_IS', 'LLC_NISA', 'LLC_DH', 'LLC_DHA', '-v7.3')
%%
%load('backgrounds_new.mat')
%load('qc_ts.mat')
snapshot_dates = {'01-Oct-2011','01-Nov-2011', '01-Dec-2011', '01-Jan-2012', '01-Feb-2012', '01-Mar-2012', '01-Apr-2012', '01-May-2012', ...
    '01-Jun-2012', '01-Jul-2012', '01-Aug-2012', '01-Sep-2012'};
date = snapshot_dates{1};
LLC = cell(4,1);
sectors = {'LLC_1'};
for i = 1
    load('/Volumes/Elements/LLCsealdata/Snapshot_' + string(date) + '/' + string(sectors{i}) + '/depth.mat');
    LLC{i}.depth = depth;
end
LLC_1 = LLC{1};

%%
save_fig = 0;
f = figure('Position', [100 100 1000 800]);
tiledlayout(1,2)
fs = 15;

%%% LLC colors
clrs = distinguishable_colors(length(LLC_N2));
density_grid = qc_ts(1).ds.sigma0(:,1);

LLC_x = LLC_DH;
MEOP_x = MEOP_DH;
y = depth_grid;
xlabel_text = 'DHA';
ylabel_text = 'z';

nexttile
hold on
for i = 1:length(LLC_x)
    a(i) = plot(LLC_x{3,i}, y, 'Color', clrs(i,:), 'LineWidth', 2, 'DisplayName', LLC_N2{1,i});
end
b = plot(MEOP_x{3,1}, y, 'k', 'LineWidth', 6, 'DisplayName', 'MEOP');
set(gca, 'YDir', 'reverse')
ylim([0 500])
%yline(-LLC_1.depth, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5)
grid on
xlabel(string(xlabel_text), 'FontSize', fs)
ylabel(string(ylabel_text), 'FontSize', fs)
legend([a,b], 'Location', 'SouthEast','FontSize', fs)
title('Shelf')
ax = gca; ax.FontSize = fs;

nexttile
hold on
for i = 1:length(LLC_x)
    a(i) = plot(LLC_x{4,i},y, 'Color', clrs(i,:), 'LineWidth', 2, 'DisplayName', LLC_N2{1,i});
end
b = plot(MEOP_x{4,1}, y, 'k', 'LineWidth', 6, 'DisplayName', 'MEOP');
set(gca, 'YDir', 'reverse')
ylim([0 500])
%yline(-LLC_1.depth, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5)
grid on
xlabel(string(xlabel_text), 'FontSize', fs)
ylabel(string(ylabel_text), 'FontSize', fs)
%legend([a,b], 'Location', 'SouthEast','FontSize', fs)
title('Off-Shelf')
ax = gca; ax.FontSize = fs;

if save_fig == 1
    exportgraphics(f,'Paper Figures/background_stratification_w_bathymetry.png','Resolution',600)
end

% nexttile
% hold on
% for i = 1:length(LLC_IS)
%     plot(LLC_IS{2,i}, 1:length(LLC_IS{2,1}), 'Color', clrs(i,:), 'LineWidth', 2, 'DisplayName', LLC_N2ref{1,i})
% end
% set(gca, 'YDir', 'reverse')
% xlabel('Isopycnal Separation', 'FontSize', fs)
% ylabel('Pressure [dbar]', 'FontSize', fs)
% grid on
% legend('Location', 'northeast', 'FontSize', fs)

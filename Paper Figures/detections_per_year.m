save_fig = 0;
input_path = '/Volumes/Elements/MEOPdata';

%%% Getting temporal distribution of all MEOP profiles 
load(string(input_path) + "/qc_ts.mat")
MEOPdates = [];
for tag_no = 1:length(qc_ts)
    MEOPdates = vertcat(MEOPdates, qc_ts(tag_no).time);
end
MEOPmonths_all = MEOPdates(:,2);
clear MEOPdates
for i = 1:12
    m(i) = sum(MEOPmonths_all == i);
end

%%% LLC4320 snapshot dates
snapshot_dates = {'01-Oct-2011','01-Nov-2011', '01-Dec-2011', '01-Jan-2012', '01-Feb-2012', '01-Mar-2012', '01-Apr-2012', '01-May-2012', ...
    '01-Jun-2012', '01-Jul-2012', '01-Aug-2012', '01-Sep-2012'};
snapshot_months = {'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep'};

%%% Creating figure
fs = 15;  lw = 2;
f = figure('Position', [100 100 800 600]);
tiledlayout(3,1, 'TileSpacing', 'tight')

%%%%%%%%%%%%%%%
%%% LLC4320 %%%
%%%%%%%%%%%%%%%

nexttile()
title('LLC4320', 'FontSize', fs)
hold on
%%% CYCLONES

%%% Loading data
load(string(input_path) + '/optimized_cyclones_shallow.mat')
LLCcyclones_shallow = LLCcyclones_optimized;
load(string(input_path) + '/optimized_cyclones_deep.mat')
LLCcyclones_deep = LLCcyclones_optimized;
LLCcyclones_optimized = [LLCcyclones_deep LLCcyclones_shallow];

%%% Getting detections per month
for i = 1:length(snapshot_dates)
    n(i) = sum(strcmp(string(vertcat(LLCcyclones_optimized.date)), snapshot_dates{i}));
end
plot(1:12, n, 'r-o', 'MarkerFaceColor', 'r','LineWidth', lw, 'DisplayName', 'Cyclones')

%%% ANTICYCLONES

%%% Loading data
load(string(input_path) + '/optimized_anticyclones_shallow.mat')
LLCanticyclones_shallow = LLCanticyclones_optimized;
load(string(input_path) + '/optimized_anticyclones_deep.mat')
LLCanticyclones_deep = LLCanticyclones_optimized;
LLCanticyclones_optimized = [LLCanticyclones_deep LLCanticyclones_shallow];

%%% Getting detections per month
for i = 1:length(snapshot_dates)
    n(i) = sum(strcmp(string(vertcat(LLCanticyclones_optimized.date)), snapshot_dates{i}));
end
plot(1:12, n, 'b-o', 'MarkerFaceColor', 'b', 'LineWidth', lw,'DisplayName', 'Anticyclones')

grid on
xlim([1 12])
xticks(1:12)
xticklabels([])
ylabel('Number of Detections', 'FontSize', fs)
legend('Location', 'northeast', 'FontSize', fs)
ax = gca; ax.FontSize = fs;

%%%%%%%%%%%%
%%% MEOP %%%
%%%%%%%%%%%%

nexttile()
title('MEOP', 'FontSize', fs)
hold on

%%% CYCLONES

%%% Loading data
load(string(input_path) + '/optimized_cyclones_shallow.mat')
MEOPcyclones_shallow = MEOPcyclones_optimized;
load(string(input_path) + '/optimized_cyclones_deep.mat')
MEOPcyclones_deep = MEOPcyclones_optimized;
if ~isempty(MEOPcyclones_shallow)
    MEOPcyclones_optimized = [MEOPcyclones_deep MEOPcyclones_shallow];
else
    MEOPcyclones_optimized = MEOPcyclones_deep;
end

%%% Getting detections per month
MEOPmonths = vertcat(MEOPcyclones_optimized.date);
MEOPmonths = MEOPmonths(:,2);
for i = 1:12
    n(i) = sum(MEOPmonths == i);
end
plot(1:12, n, 'r-o', 'MarkerFaceColor', 'r', 'LineWidth', lw,'DisplayName', 'Cyclones')

%%% ANTICYCLONES

%%% Loading data
load(string(input_path) + '/optimized_anticyclones_shallow.mat')
MEOPanticyclones_shallow = MEOPanticyclones_optimized;
load(string(input_path) + '/optimized_anticyclones_deep.mat')
MEOPanticyclones_deep = MEOPanticyclones_optimized;
if ~isempty(MEOPanticyclones_shallow)
    MEOPanticyclones_optimized = [MEOPanticyclones_deep MEOPanticyclones_shallow];
else
    MEOPanticyclones_optimized = MEOPanticyclones_deep;
end

%%% Getting detections per month
MEOPmonths = vertcat(MEOPanticyclones_optimized.date);
MEOPmonths = MEOPmonths(:,2);
for i = 1:12
    n(i) = sum(MEOPmonths == i);
end
plot(1:12, n, 'b-o', 'MarkerFaceColor', 'b', 'LineWidth', lw,'DisplayName', 'Anticyclones')
grid on
xlim([1 12])
xticks(1:12)
xticklabels([]);
ylabel('Number of Detections', 'FontSize', fs)
ax = gca; ax.FontSize = fs;

%%%%%%%%%%%%%%%%%%%%%%%
%%% MEOP Normalized %%%
%%%%%%%%%%%%%%%%%%%%%%%

nexttile()
title('MEOP Normalized by Data Coverage', 'FontSize', fs)
hold on

%%% CYCLONES

%%% Detections normalized
MEOPmonths = vertcat(MEOPcyclones_optimized.date);
MEOPmonths = MEOPmonths(:,2);
for i = 1:12
    n(i) = sum(MEOPmonths == i);
end
plot(1:12, (n ./ m) * 100, 'r-o', 'MarkerFaceColor', 'r', 'LineWidth', lw,'DisplayName', 'Cyclones')

%%% ANTICYCLONES

%%% Detections normalized
MEOPmonths = vertcat(MEOPanticyclones_optimized.date);
MEOPmonths = MEOPmonths(:,2);
for i = 1:12
    n(i) = sum(MEOPmonths == i);
end
plot(1:12, (n ./ m) * 100, 'b-o', 'MarkerFaceColor', 'b', 'LineWidth', lw,'DisplayName', 'Anticyclones')
grid on
xlim([1 12])
xticks(1:12)
ylabel({'Perecentage of Profiles'; 'Identified as SCVs [%]'}, 'FontSize', fs)
ax = gca; ax.FontSize = fs;

xticklabels(snapshot_months)

if save_fig == 1
    exportgraphics(f,'Paper Figures/detections_per_snapshot.png','Resolution',600)
end

load('/Volumes/Elements/MEOP Optimization/Anticyclones/final_stats_anticyclones_deep.mat');
anticyclone_deep_precision = precision;
anticyclone_deep_f = f_all;
anticyclone_deep_no_iteration = 1:length(anticyclone_deep_f);

load('/Volumes/Elements/MEOP Optimization/Anticyclones/final_stats_anticyclones_shallow.mat');
anticyclone_shallow_precision = precision;
anticyclone_shallow_f = f_all;
anticyclone_shallow_no_iteration = 1:length(anticyclone_shallow_f);

load('/Volumes/Elements/MEOP Optimization/Cyclones/final_stats_cyclones_deep.mat')
cyclone_deep_precision = precision;
cyclone_deep_f = f_all;
cyclone_deep_no_iteration = 1:length(cyclone_deep_f);

load('/Volumes/Elements/MEOP Optimization/Cyclones/final_stats_cyclones_shallow.mat')
cyclone_shallow_precision = precision;
cyclone_shallow_f = f_all;
cyclone_shallow_no_iteration = 1:length(cyclone_shallow_f);


%%

fs = 15;
lw = 2;

f = figure();
tiledlayout(2,1,"TileSpacing","tight")

nexttile()
hold on
a = plot(anticyclone_deep_no_iteration, anticyclone_deep_precision*100,'k',  'LineWidth', lw, 'DisplayName', 'Anticyclones');
plot(anticyclone_shallow_no_iteration, anticyclone_shallow_precision*100, 'k--','LineWidth', lw, 'DisplayName', 'Anticyclones');
b = plot(cyclone_deep_no_iteration, cyclone_deep_precision*100,'r', 'LineWidth', lw,  'DisplayName', 'Cyclones');
plot(cyclone_shallow_no_iteration, cyclone_shallow_precision*100,'r--', 'LineWidth', lw,  'DisplayName', 'Cyclones');
hold off
ylabel('Precision [%]', 'FontSize', fs)
legend([a b], 'Location','southeast','FontSize', fs)
ax = gca; ax.FontSize = fs;
grid on
xlim([0 50])

nexttile()
hold on
a = plot(anticyclone_deep_no_iteration, anticyclone_deep_f, 'k', 'LineWidth', lw, 'DisplayName', 'Anticyclones');
plot(anticyclone_shallow_no_iteration, anticyclone_shallow_f, 'k--', 'LineWidth', lw, 'DisplayName', 'Anticyclones');
b = plot(cyclone_deep_no_iteration, cyclone_deep_f,'r', 'LineWidth', lw,  'DisplayName', 'Cyclones');
plot(cyclone_shallow_no_iteration, cyclone_shallow_f,'r--', 'LineWidth', lw,  'DisplayName', 'Cyclones');
hold off
ylabel('F-Score', 'FontSize', fs)
xlabel('Iteration Number', 'FontSize', fs)
ax = gca; ax.FontSize = fs;
grid on
xlim([0 50])

exportgraphics(f,'Paper Figures/optimization_results.png','Resolution',600)
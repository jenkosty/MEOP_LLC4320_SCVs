clear
input_path = '/Volumes/Elements/MEOPdata';
load(string(input_path) + '/optimized_cyclones_shallow.mat')
load(string(input_path) + '/final_stats_cyclones_shallow.mat')
clc
disp('CYCLONES')

%%% LLC cyclones
disp('---LLC---')
disp('Average number of detections: ' + string(length(LLCcyclones_optimized)/12))
disp('Average precision: ' + string(precision(end)));
disp('Number of minty eddies: ' + string(length(LLCcyclones_optimized([LLCcyclones_optimized.anom_A] < 0))));
disp('Average spice anomaly: ' + string(round(mean([LLCcyclones_optimized.anom_A]), 2, 'significant')))
disp('Average DHA anomaly: ' + string(round(mean([LLCcyclones_optimized.dha_A]), 2, 'significant')))
disp('Average width: ' + string(round(mean([LLCcyclones_optimized.dha_Hcore]), 2, 'significant')))
disp('Average depth: ' + string(round(mean([LLCcyclones_optimized.dha_P]), 2, 'significant')))

disp('   ')

%%% MEOP cyclones
disp('---MEOP---')
disp('Number of detections: ' + string(length(MEOPcyclones_optimized)))
disp('Number of minty eddies: ' + string(length(MEOPcyclones_optimized([MEOPcyclones_optimized.anom_A] < 0))));
disp('Average spice anomaly: ' + string(round(mean([MEOPcyclones_optimized.anom_A]), 2, 'significant')))
disp('Average DHA anomaly: ' + string(round(mean([MEOPcyclones_optimized.dha_A]), 2, 'significant')))
disp('Average width: ' + string(round(mean([MEOPcyclones_optimized.dha_Hcore]), 2, 'significant')))
disp('Average depth: ' + string(round(mean([MEOPcyclones_optimized.dha_P]), 2, 'significant')))

disp('   ')

%%% Min Values
disp('---Minimum Values---')
for i = 1:length(param_names_min)
    disp(string(param_names_min(i)) + ': ' + string(round(final_stats_min.(param_names_min(i)), 2, 'significant')))
end

disp('   ')

%%% Max Values
disp('---Maximum Values---')
for i = 1:length(param_names_max)
    disp(string(param_names_max(i)) + ': ' + string(round(final_stats_max.(param_names_max(i)), 2, 'significant')))
end

%%
clear
input_path = '/Volumes/Elements/MEOPdata';
load(string(input_path) + '/optimized_anticyclones_shallow.mat')
load(string(input_path) + '/final_stats_anticyclones_shallow.mat')
clc
disp('ANTICYCLONES')

%%% LLC anticyclones
disp('---LLC---')
disp('Average number of detections: ' + string(length(LLCanticyclones_optimized)/12))
disp('Average precision: ' + string(precision(end)));
disp('Number of minty eddies: ' + string(length(LLCanticyclones_optimized([LLCanticyclones_optimized.spice_A] < 0))));
disp('Average spice anomaly: ' + string(round(mean([LLCanticyclones_optimized.spice_A]), 2, 'significant')))
disp('Average ISA anomaly: ' + string(round(mean([LLCanticyclones_optimized.anom_A]), 2, 'significant')))
disp('Average DHA anomaly: ' + string(round(mean([LLCanticyclones_optimized.dha_A]), 2, 'significant')))
disp('Average width: ' + string(round(mean([LLCanticyclones_optimized.dha_Hcore]), 2, 'significant')))
disp('Average depth: ' + string(round(mean([LLCanticyclones_optimized.dha_P]), 2, 'significant')))

disp(' ')

%%% MEOP anticyclones
disp('---MEOP---')
disp('Number of detections: ' + string(length(MEOPanticyclones_optimized)))
disp('Number of minty eddies: ' + string(length(MEOPanticyclones_optimized([MEOPanticyclones_optimized.spice_A] < 0))));
disp('Average spice anomaly: ' + string(round(mean([MEOPanticyclones_optimized.spice_A]), 2, 'significant')))
disp('Average ISA anomaly: ' + string(round(mean([MEOPanticyclones_optimized.anom_A]), 2, 'significant')))
disp('Average DHA anomaly: ' + string(round(mean([MEOPanticyclones_optimized.dha_A]), 2, 'significant')))
disp('Average width: ' + string(round(mean([MEOPanticyclones_optimized.dha_Hcore]), 2, 'significant')))
disp('Average depth: ' + string(round(mean([MEOPanticyclones_optimized.dha_P]), 2, 'significant')))

disp(' ')

%%% Min Values
disp('---Minimum Values---')
for i = 1:length(param_names_min)
    disp(string(param_names_min(i)) + ': ' + string(round(final_stats_min.(param_names_min(i)), 2, 'significant')))
end

disp('   ')

%%% Max Values
disp('---Maximum Values---')
for i = 1:length(param_names_max)
    disp(string(param_names_max(i)) + ': ' + string(round(final_stats_max.(param_names_max(i)), 2, 'significant')))
end


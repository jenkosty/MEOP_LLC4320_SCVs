clear
input_path = '/Volumes/Elements/MEOPdata';
load(string(input_path) + '/optimized_cyclones_deep.mat')
load(string(input_path) + '/final_stats_cyclones_deep.mat')
clc
disp('CYCLONES')

%%% LLC cyclones
disp('---LLC---')
disp('Average number of detections: ' + string(round(length(LLCcyclones_optimized)/12, 2, 'significant')))
disp('Average precision: ' + string(round(precision(end), 2, 'significant')));
disp('Number of minty eddies: ' + string(length(LLCcyclones_optimized([LLCcyclones_optimized.anom_A] < 0))));
disp('Average spice anomaly: ' + string(round(mean([LLCcyclones_optimized.anom_A]), 2, 'significant')))
disp('Average DHA anomaly: ' + string(round(mean([LLCcyclones_optimized.dha_A]), 2, 'significant')))
disp('Average width: ' + string(round(mean([LLCcyclones_optimized.dha_Hcore]), 2, 'significant')))
disp('Average depth: ' + string(round(mean([LLCcyclones_optimized.dha_P]), 2, 'significant')))

disp('   ')

%%% MEOP cyclones
disp('---MEOP---')
disp('Number of detections: ' + string(round(length(MEOPcyclones_optimized), 2, 'significant')))
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

clear; clc;
input_path = '/Volumes/Elements/MEOPdata';

%%% Shallow detections
load(string(input_path) + '/optimized_cyclones_shallow.mat')
load(string(input_path) + '/final_stats_cyclones_shallow.mat')
LLCcyclones_optimized_shallow = LLCcyclones_optimized;
Pshallow = precision(end);
MEOPcyclones_optimized_shallow = MEOPcyclones_optimized;

%%% Deep detections
load(string(input_path) + '/optimized_cyclones_deep.mat')
load(string(input_path) + '/final_stats_cyclones_deep.mat')
LLCcyclones_optimized_deep = LLCcyclones_optimized;
Pdeep = precision(end);
MEOPcyclones_optimized_deep = MEOPcyclones_optimized;

%%% All detecctions
load(string(input_path) + '/optimized_cyclones.mat')

disp('--- LLC ---')
Nshallow = length(LLCcyclones_optimized_shallow);
Ndeep = length(LLCcyclones_optimized_deep);
disp('Total number of detections: ' + string(Nshallow+Ndeep))
disp('Detections per snapshot: ' + string((Nshallow+Ndeep) / 12))
disp('Detections off-shelf: ' + string((Ndeep)))
disp('Detections on-shelf: ' + string((Nshallow)))
disp('Percentage off-shelf: ' + string(Ndeep/(Nshallow+Ndeep) * 100))
disp('Precision off-shelf: ' + string(Pdeep*100))
disp('Precision on-shelf: ' + string(Pshallow*100))
disp('Average precision: ' + string(((Ndeep*Pdeep) + (Nshallow*Pshallow)) / (Ndeep+Nshallow)))
ind = [LLCcyclones_optimized.anom_A] < 0;
disp('Minty core: ' + string(length(LLCcyclones_optimized(ind)) / length(LLCcyclones_optimized) * 100))

disp('--- MEOP ---')
Nshallow = length(MEOPcyclones_optimized_shallow);
Ndeep = length(MEOPcyclones_optimized_deep);
disp('Total number of detections: ' + string(Nshallow+Ndeep))
disp('Percentage off-shelf: ' + string(Ndeep/(Nshallow+Ndeep) * 100))
ind = [MEOPcyclones_optimized.anom_A] < 0;
disp('Minty core: ' + string(length(MEOPcyclones_optimized(ind)) / length(MEOPcyclones_optimized) * 100))


%%
clear
input_path = '/Volumes/Elements/MEOPdata';
load(string(input_path) + '/optimized_anticyclones_deep.mat')
load(string(input_path) + '/final_stats_anticyclones_deep.mat')
clc
disp('ANTICYCLONES')

%%% LLC anticyclones
disp('---LLC---')
disp('Average number of detections: ' + string(round(length(LLCanticyclones_optimized)/12, 2, 'significant')))
disp('Average precision: ' + string(round(precision(end), 2, 'significant')));
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

%%
clear; clc;
input_path = '/Volumes/Elements/MEOPdata';

%%% Shallow detections
load(string(input_path) + '/optimized_anticyclones_shallow.mat')
load(string(input_path) + '/final_stats_anticyclones_shallow.mat')
LLCanticyclones_optimized_shallow = LLCanticyclones_optimized;
Pshallow = precision(end);
MEOPanticyclones_optimized_shallow = MEOPanticyclones_optimized;

%%% Deep detections
load(string(input_path) + '/optimized_anticyclones_deep.mat')
load(string(input_path) + '/final_stats_anticyclones_deep.mat')
LLCanticyclones_optimized_deep = LLCanticyclones_optimized;
Pdeep = precision(end);
MEOPanticyclones_optimized_deep = MEOPanticyclones_optimized;

%%% All detecctions
load(string(input_path) + '/optimized_anticyclones.mat')

disp('--- LLC ---')
Nshallow = length(LLCanticyclones_optimized_shallow);
Ndeep = length(LLCanticyclones_optimized_deep);
disp('Total number of detections: ' + string(Nshallow+Ndeep))
disp('Detections per snapshot: ' + string((Nshallow+Ndeep) / 12))
disp('Detections off-shelf: ' + string((Ndeep)))
disp('Detections on-shelf: ' + string((Nshallow)))
disp('Percentage off-shelf: ' + string(Ndeep/(Nshallow+Ndeep) * 100))
disp('Precision off-shelf: ' + string(Pdeep*100))
disp('Precision on-shelf: ' + string(Pshallow*100))
disp('Average precision: ' + string(((Ndeep*Pdeep) + (Nshallow*Pshallow)) / (Ndeep+Nshallow)))
ind = [LLCanticyclones_optimized.anom_A] < 0;
disp('Minty core: ' + string(length(LLCanticyclones_optimized(ind)) / length(LLCanticyclones_optimized) * 100))

disp('--- MEOP ---')
Nshallow = length(MEOPanticyclones_optimized_shallow);
Ndeep = length(MEOPanticyclones_optimized_deep);
disp('Total number of detections: ' + string(Nshallow+Ndeep))
disp('Percentage off-shelf: ' + string(Ndeep/(Nshallow+Ndeep) * 100))
ind = [MEOPanticyclones_optimized.anom_A] < 0;
disp('Minty core: ' + string(length(MEOPanticyclones_optimized(ind)) / length(MEOPanticyclones_optimized) * 100))

spice_anoms = [MEOPanticyclones_optimized_deep.spice_A MEOPanticyclones_optimized_shallow.spice_A];
[h,p] = ttest(spice_anoms);
disp('p-value: ' + string(p));


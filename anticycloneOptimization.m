%% Loading characteristics from anticyclone detections (all snapshots)

%%% LLC4320 snapshot dates
snapshot_dates = {'01-Oct-2011','01-Nov-2011', '01-Dec-2011', '01-Jan-2012', '01-Feb-2012', '01-Mar-2012', '01-Apr-2012', '01-May-2012', ...
    '01-Jun-2012', '01-Jul-2012', '01-Aug-2012', '01-Sep-2012'};

%%% Path for LLCsealdata and algorithm output
input_path = '/Volumes/Elements/LLCsealdata/Snapshot_';
output_path = '/Volumes/Elements/MEOPdata';

%%% Loading LLC seal track data (only need 1 year)
date = snapshot_dates{1};
load(string(input_path) + string(date) + '/LLCsealdata_full.mat')

%%% Looping through snapshots
u = 1;
uu = 1;
for ii = 1:length(snapshot_dates)
    disp(string(snapshot_dates{ii}))
    load(string(input_path) + string(snapshot_dates{ii}) + '/anticyclone_data.mat')
    load(string(input_path) + string(snapshot_dates{ii}) + '/lilly_data_final');
    load(string(input_path) + string(snapshot_dates{ii}) + '/background_data')
    date = snapshot_dates{ii};

    %%% Looping through tracks
    for tag_no = 1:length(LLCsealdata)
        for i = anticyclones(tag_no).dha
            
            isa_gauss_tag = isa_gauss{tag_no};
            dha_gauss_
            tag = dha_gauss{tag_no};
            spice_gauss_tag = spice_gauss{tag_no};

            % %%% Excluding detections with a bad dha fit
            % if median(dha_gauss_tag(i).dataX) < 0
            %  continue
            % end

            %%% Basic info
            anticyclone_data(u).date = date;
            anticyclone_data(u).tag_no = tag_no;
            anticyclone_data(u).cast = i;
            anticyclone_data(u).bathymetry = LLCsealdata(tag_no).bathymetry(i);

            %%% ISA fit info
            anticyclone_data(u).anom_A = isa_gauss_tag(i).A;
            anticyclone_data(u).anom_R2 = isa_gauss_tag(i).R2;
            anticyclone_data(u).anom_nrmse = isa_gauss_tag(i).nrmse;
            anticyclone_data(u).anom_Hcore = isa_gauss_tag(i).Hcore;
            anticyclone_data(u).anom_P = isa_gauss_tag(i).P;
            
            %%% DHA fit info
            anticyclone_data(u).dha_A = dha_gauss_tag(i).A;
            anticyclone_data(u).dha_R2 = dha_gauss_tag(i).R2;
            anticyclone_data(u).dha_nrmse = dha_gauss_tag(i).nrmse;
            anticyclone_data(u).dha_Hcore = dha_gauss_tag(i).Hcore;
            anticyclone_data(u).dha_P = dha_gauss_tag(i).P;

            %%% Spice fit info
            if i <= length(spice_gauss_tag)
                anticyclone_data(u).spice_A = spice_gauss_tag(i).A;
                anticyclone_data(u).spice_R2 = spice_gauss_tag(i).R2;
                anticyclone_data(u).spice_nrmse = spice_gauss_tag(i).nrmse;
                anticyclone_data(u).spice_Hcore = spice_gauss_tag(i).Hcore;
                anticyclone_data(u).spice_P = spice_gauss_tag(i).P;
            end

            %%% Background data
            anticyclone_data(u).bathymetric_var = background(tag_no).bathymetric_var(i);
            anticyclone_data(u).shelf_break_ratio = background(tag_no).shelf_break_ratio(i);
            anticyclone_data(u).isopycnal_stability = background(tag_no).isopycnal_stablity(i);
            anticyclone_data(u).MLD = background(tag_no).MLD(i);
            anticyclone_data(u).spice_std = background(tag_no).spice_std(i);
            anticyclone_data(u).max_pres = background(tag_no).max_pres(i);

            %%% Lilly SCV data
            anticyclone_data(u).True_SCV = lilly_data(tag_no).scv_deep(i);
            anticyclone_data(u).scv_reason = lilly_data(tag_no).scv_reason_deep(i);
            if anticyclone_data(u).True_SCV == 1
                anticyclone_data(u).OW_in_contour = lilly_data(tag_no).contourdata(i).OW_in_contour;
                anticyclone_data(u).area = lilly_data(tag_no).contourdata(i).area;
                anticyclone_data(u).ecc = lilly_data(tag_no).contourdata(i).ecc;
                anticyclone_data(u).area_ratio = lilly_data(tag_no).contourdata(i).area_ratio;
                anticyclone_data(u).vort_in_contour = lilly_data(tag_no).contourdata(i).vort_in_contour;
                if (lilly_data(tag_no).contourdata(i).vort_in_contour < 0)
                    anticyclone_data(u).True_anticyclone = 1;
                else
                    anticyclone_data(u).True_anticyclone = 0;
                end
            else
                anticyclone_data(u).mean_OW = NaN;
                anticyclone_data(u).area = NaN;
                anticyclone_data(u).True_anticyclone = 0;
            end

            u = u + 1;

        end
    end

    %%% Missed Cyclones
    for tag_no = 1:length(lilly_data)
        detected = [anticyclones(tag_no).dha];
        OW_flagged = find([lilly_data(tag_no).contourdata.vort_in_contour] < 0);
        ind = setdiff(OW_flagged, detected);
        for i = 1:length(ind)
            missed_anticyclones(uu).date = date;
            missed_anticyclones(uu).tag_no = tag_no;
            missed_anticyclones(uu).cast = ind(i);

            %%% Background data
            missed_anticyclones(uu).bathymetric_var = background(tag_no).bathymetric_var(ind(i));
            missed_anticyclones(uu).shelf_break_ratio = background(tag_no).shelf_break_ratio(ind(i));
            missed_anticyclones(uu).isopycnal_stability = background(tag_no).isopycnal_stablity(ind(i));
            missed_anticyclones(uu).MLD = background(tag_no).MLD(ind(i));
            missed_anticyclones(uu).spice_std = background(tag_no).spice_std(i);
            missed_anticyclones(uu).max_pres = background(tag_no).max_pres(i);

            missed_anticyclones(uu).mean_OW = lilly_data(tag_no).contourdata(ind(i)).OW_in_contour;
            missed_anticyclones(uu).area = lilly_data(tag_no).contourdata(ind(i)).area;
            uu = uu + 1;
        end
     end

    clear anticyclones isa_gauss dha_gauss N2_gauss N2_gauss_tag isa_gauss_tag dha_gauss_tag spice_gauss_tag lilly_data background
end

detected_anticyclones = anticyclone_data;
save(string(output_path) + '/LLCanticyclones', 'detected_anticyclones', 'missed_anticyclones', '-v7.3')

%% Applying a minimum OW threshold 

input_path = '/Volumes/Elements/MEOPdata';
output_path = '/Volumes/Elements/MEOPdata';
load(string(input_path) + '/LLCanticyclones.mat')

OW_requirement = 0.2e-9;

for i = 1:length(detected_anticyclones)
    if abs(detected_anticyclones(i).OW_in_contour) < OW_requirement
        detected_anticyclones(i).True_anticyclone = 0;
    end
end
clear i

%%% Shallow or Deep?
deep = 1;
depth_threshold = 1000;
if deep == 1
    detected_anticyclones = detected_anticyclones(abs([detected_anticyclones.bathymetry]) >= depth_threshold);
else
    detected_anticyclones = detected_anticyclones(abs([detected_anticyclones.bathymetry]) < depth_threshold);
end

%%% Loading seal data
load(string(input_path) + '/qc_ts.mat');

%%% Removing edge cases
for u = 1:length(detected_anticyclones)
    tag_no = detected_anticyclones(u).tag_no;
    i = detected_anticyclones(u).cast;

    first_cast = 1;
    last_cast = length(qc_ts(tag_no).cast);

    if i <= (first_cast + 9) | i >= (last_cast - 9)
        ind(u) = 1;
    else
        ind(u) = 0;
    end
end
detected_anticyclones = detected_anticyclones(~ind);

% %%% Excluding detections based on pressure difference of fits
% clear ind
% for u = 1:length(detected_anticyclones)
%     pres_diff = abs(detected_anticyclones(u).anom_P - detected_anticyclones(u).dha_P);
% 
%     if pres_diff <= 150
%         ind(u) = 1;
%     else
%         ind(u) = 0;
%     end
% 
% end
% detected_anticyclones = detected_anticyclones(~ind);

%%% Formatting data for optimization
detected_scvs = detected_anticyclones;
missed_scvs = missed_anticyclones;

%%% Additional requirements (making sure spice fits are reasonable)
detected_scvs = detected_scvs(find(~cellfun(@isempty,{detected_scvs.spice_A})));
detected_scvs = detected_scvs([detected_scvs.spice_P] >= 0);
detected_scvs = detected_scvs([detected_scvs.max_pres] >= 400);

%%% Getting original results
detected_scvs_og = detected_scvs;
missed_scvs_og = missed_scvs;
OG_TP = length(detected_scvs_og([detected_scvs_og.True_anticyclone] == 1));
OG_FP = length(detected_scvs_og([detected_scvs_og.True_anticyclone] == 0));
OG_FN = length(missed_scvs_og);

fn = fieldnames(detected_scvs);
fn = fn([5:14, 20:24]);

clear u tag_no last_cast first_cast i ind

%% Getting initial guesses

%%% Calculating f-score for each parameter
clear guess_min_og guess_all_min guess_max_og guess_all_max

betas_all = [0.2, 0.3, 0.4, 0.5];

%%% one loop for maxes, one loop for mins
for cntr = 1:2
    for uu = 1:numel(fn)
        disp(uu)

        %%% Finding values to test (excluding extreme values)
        param = abs([detected_scvs.(fn{uu})]);
        idx = (prctile(param, 90)-prctile(param, 10))/200;
        arr = prctile(param, 10):idx:prctile(param, 90);

        u = 1;
        clear stats
        for i = arr

            if cntr == 1
                ind = (param >= i);
            else
                ind = (param <= i);
            end

            detected_scvs_opt = detected_scvs(ind);

            TP = length(detected_scvs_opt([detected_scvs_opt.True_anticyclone] == 1));
            FP = length(detected_scvs_opt([detected_scvs_opt.True_anticyclone] == 0));
            FN = OG_FN + (OG_TP - TP);

            stats(u).precision = TP / (TP + FP);
            stats(u).recall = TP / (TP + FN);
            ii = 1;
            beta = betas_all(ii);
            stats(u).f_2 = (1+beta^2) * (stats(u).precision * stats(u).recall) / ((beta^2 * stats(u).precision) + stats(u).recall);
            ii = ii + 1;
            stats(u).f_3 = (1+beta^2) * (stats(u).precision * stats(u).recall) / ((beta^2 * stats(u).precision) + stats(u).recall);
            ii = ii + 1;
            beta = betas_all(ii);
            stats(u).f_4 = (1+beta^2) * (stats(u).precision * stats(u).recall) / ((beta^2 * stats(u).precision) + stats(u).recall);
            ii = ii + 1;
            stats(u).f_5 = (1+beta^2) * (stats(u).precision * stats(u).recall) / ((beta^2 * stats(u).precision) + stats(u).recall);
            stats(u).no_detections = length(detected_scvs_opt);
            stats(u).val = i;
            stats(u).param = fn{uu};
            u = u + 1;
        end

        n = 1;
        f = [stats.f_2];
        [m,~] = max(f, [], "all", "omitnan");
        ind = find([stats.f_2] == m);
        initial_stats_all{n,uu} = stats(ind);

        n = n + 1;
        f = [stats.f_3];
        [m,~] = max(f, [], "all", "omitnan");
        ind = find([stats.f_3] == m);
        initial_stats_all{n,uu} = stats(ind);

        n = n + 1;
        f = [stats.f_4];
        [m,~] = max(f, [], "all", "omitnan");
        ind = find([stats.f_4] == m);
        initial_stats_all{n,uu} = stats(ind);

        n = n + 1;
        f = [stats.f_5];
        [m,~] = max(f, [], "all", "omitnan");
        ind = find([stats.f_5] == m);
        initial_stats_all{n,uu} = stats(ind);

    end

    %%% Saving initial guesses for all f-scores
    for uu = 1:size(initial_stats_all,2)
        initial_stats_tmp(1).(initial_stats_all{1,uu}(1).param) = initial_stats_all{1,uu}(1).val;
        initial_stats_tmp(2).(initial_stats_all{2,uu}(1).param) = initial_stats_all{2,uu}(1).val;
        initial_stats_tmp(3).(initial_stats_all{3,uu}(1).param) = initial_stats_all{3,uu}(1).val;
        initial_stats_tmp(4).(initial_stats_all{4,uu}(1).param) = initial_stats_all{4,uu}(1).val;
    end
    initial_stats_all_tmp(1,:) = initial_stats_all(1,:);
    initial_stats_all_tmp(2,:) = initial_stats_all(2,:);
    initial_stats_all_tmp(3,:) = initial_stats_all(3,:);
    initial_stats_all_tmp(4,:) = initial_stats_all(4,:);

    if cntr == 1
        guess_all_min = initial_stats_all_tmp;
    else
        guess_all_max = initial_stats_all_tmp;
    end

end

%%% Savining initial guess for chosen f-score
for cntr = 1:2
    if cntr == 1
        initial_stats_all = guess_all_min;
    else
        initial_stats_all = guess_all_max;
    end
    for uu = 1:size(initial_stats_all,2)
        initial_stats_tmp(1).(initial_stats_all{1,uu}(1).param) = initial_stats_all{1,uu}(1).val;
        initial_stats_tmp(2).(initial_stats_all{2,uu}(1).param) = initial_stats_all{2,uu}(1).val;
        initial_stats_tmp(3).(initial_stats_all{3,uu}(1).param) = initial_stats_all{3,uu}(1).val;
        initial_stats_tmp(4).(initial_stats_all{4,uu}(1).param) = initial_stats_all{4,uu}(1).val;
    end

    ii = 4;
    names = fieldnames(initial_stats_tmp);
    if cntr == 1
        for uu = 1:length(names)
            guess_min_og.(names{uu}) = initial_stats_tmp(ii).(names{uu});
            guess_min_og.beta = betas_all(ii);
        end
    else
        for uu = 1:length(names)
            guess_max_og.(names{uu}) = initial_stats_tmp(ii).(names{uu});
            guess_max_og.beta = betas_all(ii);
        end  
    end

end

clear initial_stats_tmp initial_stats_all_tmp f FN FP i idx ind m missed_scvs missed_scvs_og param stats TP u uu...
    u initial_stats_all detected_scvs_og detected_scvs_opt beta arr n ii

%% Getting initial delta values

%%% List of parameters to optimize
% param_names_min = ["anom_A",  "anom_R2", "anom_P", "anom_Hcore", "dha_A", "dha_R2", "dha_P", "dha_Hcore","spice_A", "spice_R2", "spice_P", "spice_Hcore"];
% param_names_max = ["bathymetric_var", "isopycnal_stability", "spice_std", "MLD", "anom_nrmse", "anom_P", "anom_Hcore", "dha_nrmse", "dha_P", "dha_Hcore", "spice_nrmse", "spice_Hcore", "spice_P"];

param_names_min = ["anom_A",  "anom_R2", "anom_P", "anom_Hcore", "dha_A", "dha_R2", "dha_P", "dha_Hcore"];
param_names_max = ["bathymetric_var", "isopycnal_stability", "spice_std", "MLD", "anom_nrmse", "anom_P", "anom_Hcore", "dha_nrmse", "dha_P", "dha_Hcore"];

param_names = horzcat(param_names_min, param_names_max);

clear deltas_min deltas_max deltas deltas_og
prct = 0.2;
for uu = 1:length(param_names)
    if uu <= length(param_names_min)
        idx = find(strcmp(param_names(uu), fn));
        param = abs([detected_scvs.(fn{idx})]);
        deltas_min(uu) = prct*(prctile(param, 90)-prctile(param, 10));
    else
        idx = find(strcmp(param_names(uu), fn));
        param = abs([detected_scvs.(fn{idx})]);
        deltas_max(uu - length(param_names_min)) = prct*(prctile(param, 90)-prctile(param, 10));
    end
end
deltas = horzcat(deltas_min, deltas_max);

disp(deltas)
deltas_og = deltas;

%% Optimization loop

guess_min = guess_min_og;
guess_max = guess_max_og;
deltas = deltas_og;

% Define hyperparameters for gradient descent
max_iterations = 100;...; % Maximum number of iterations

if deep == 1
    beta0 = 0.03;
else
    beta0 = 0.02;
end
guess_min.beta = beta0;
guess_max.beta = beta0;
last_iteration = 5;

% Gradient descent optimization
clear tested_values inds f_all precision recall no_detections deltas_all param_values_all
for i = 1:max_iterations
    disp(i)

    clear tested_values inds stats
    for u = 1:length(param_names)

        %%% Getting indices for each parameter tested
        if u <= length(param_names_min)
            [tested_values{u,1}, inds{u,1}] = indexing_each_parameter(detected_scvs, param_names(u), guess_min.(param_names(u)), deltas(u), "min");
        else
            [tested_values{u,1}, inds{u,1}] = indexing_each_parameter(detected_scvs, param_names(u), guess_max.(param_names(u)), deltas(u), "max");
        end
    end

    %%% Calculating f-score and getting parameter value for next iteration
    for u = 1:length(param_names)

        all_other_params = setdiff(1:length(param_names), u);
        tmp = vertcat(inds{all_other_params,1});
        tmp = double(tmp(5:9:size(tmp,1),:));
        ind_all_other_params = prod(tmp,1);

        for uu = 1:size(inds{u,1})

            ind = ind_all_other_params .* inds{u,1}(uu,:); 
            ind = (ind == 1);

            detected_scvs_opt = detected_scvs(ind);

            TP = length(detected_scvs_opt([detected_scvs_opt.True_anticyclone] == 1));
            FP = length(detected_scvs_opt([detected_scvs_opt.True_anticyclone] == 0));
            FN = OG_FN + (OG_TP - TP);

            stats.precision(uu) = TP / (TP + FP);
            stats.recall(uu) = TP / (TP + FN);
            beta = beta0;
            stats.f(uu) = (1+beta^2) * (stats.precision(uu) * stats.recall(uu)) / ((beta^2 * stats.precision(uu)) + stats.recall(uu));

        end

        %%% Getting max f-score
        max_f = max(stats.f, [], 'omitnan');
        uu = find(stats.f == max_f);
       
        if length(uu) > 1
            uu = last_iteration;
        end

        %%% Updating parameter value
        if u <= length(param_names_min)
            guess_min.(param_names(u)) = tested_values{u,1}(uu);

        else
            guess_max.(param_names(u)) = tested_values{u,1}(uu);
        end

        %%% Setting min values for R^2
        if strcmp(param_names(u), "anom_R2") || strcmp(param_names(u), "spice_R2") || strcmp(param_names(u), "dha_R2")
            if guess_min.(param_names(u)) < 0.5
                guess_min.(param_names(u)) = 0.5;
            end
        end

        %%% Setting max value for NRMSE
        if strcmp(param_names(u), "anom_nrmse") || strcmp(param_names(u), "spice_nrmse") || strcmp(param_names(u), "dha_nrmse")
            if guess_max.(param_names(u)) > 0.5
                guess_max.(param_names(u)) = 0.5;
            end
        end

        %%% Updating delta value
        if abs(stats.f(uu) - stats.f(last_iteration)) > 0
            gamma = 1-((abs(stats.f(last_iteration) - stats.f(uu)) / stats.f(last_iteration)));
            deltas(u) = gamma*deltas(u);
        end

    end

    tmp = vertcat(inds{1:size(inds),1});
    tmp = double(tmp(5:9:size(tmp,1),:));
    ind_all_params = prod(tmp,1);
    ind = (ind_all_params == 1);

    detected_scvs_opt = detected_scvs(ind);
    TP = length(detected_scvs_opt([detected_scvs_opt.True_anticyclone] == 1));
    FP = length(detected_scvs_opt([detected_scvs_opt.True_anticyclone] == 0));
    FN = OG_FN + (OG_TP - TP);

    %%% Saving values for each iteration
    precision(i) = TP / (TP + FP);
    recall(i) = TP / (TP + FN);
    no_detections(i) = length(detected_scvs_opt);
    disp(string(precision(i)) + ', ' + string(no_detections(i)))
    beta = beta0;
    f_all(i) = (1+beta^2) * (precision(i) * recall(i)) / ((beta^2 * precision(i)) + recall(i));
    final_stats_min = guess_min;
    final_stats_max = guess_max;

    %%% Stopping when convergence is achieved
    if i > 10
        if abs(mean(precision(end-5:end)) - precision(end)) / precision(end) < 1e-10
            break
        end
    end

end

%%% Saving optimization results
if deep == 1
    save(string(output_path) + '/final_stats_anticyclones_deep', 'final_stats_max', 'final_stats_min', 'guess_max_og', 'guess_min_og', 'param_names_min', 'param_names_max', 'detected_scvs_opt', 'precision', 'f_all', 'OW_requirement', 'prct')
else
    save(string(output_path) + '/final_stats_anticyclones_shallow', 'final_stats_max', 'final_stats_min', 'guess_max_og', 'guess_min_og', 'param_names_min', 'param_names_max', 'detected_scvs_opt', 'precision', 'f_all', 'OW_requirement', 'prct')
end

%%
%%% Getting "reasonable" threshold values
final_stats_max_old = final_stats_max;
param_names_min = ["anom_A",  "anom_R2", "dha_A", "dha_R2"];
param_names_max = ["bathymetric_var", "isopycnal_stability", "spice_std", "MLD", "anom_nrmse", "dha_nrmse"];
param_names = horzcat(param_names_min, param_names_max);

clear final_stats_min final_stats_max
final_stats_min.anom_A = 5;
final_stats_min.anom_R2 = 0.6;
final_stats_min.dha_A = 0.02;
final_stats_min.dha_R2 = 0.6;

final_stats_max.bathymetric_var = final_stats_max_old.bathymetric_var;
final_stats_max.isopycnal_stability = final_stats_max_old.isopycnal_stability;
final_stats_max.spice_std = final_stats_max_old.spice_std;
final_stats_max.MLD = final_stats_max_old.MLD;
final_stats_max.anom_nrmse = 0.3;
final_stats_max.dha_nrmse = 0.3;

%%% Saving reasonable data
if deep == 1
    save(string(output_path) + '/final_stats_anticyclones_deep_reasonable', 'final_stats_max', 'final_stats_min', 'param_names_min', 'param_names_max');
else
    save(string(output_path) + '/final_stats_anticyclones_shallow_reasonable', 'final_stats_max', 'final_stats_min', 'param_names_min', 'param_names_max');
end

%% Figure to show convergence

figure()
yyaxis left
plot(1:length(precision), precision .*100)
ylabel('Precision')

yyaxis right
plot(1:length(no_detections), no_detections)
ylabel('Number of Detections')
xlabel('Iteration')


%%

function [tested_vals, inds] = indexing_each_parameter(detected_scvs, param_name, param_val, delta, extremum)
i = 1;
for u = -2:0.5:2
    tested_vals(i) = param_val + (u*delta);
    if strcmp(extremum, "min")
        inds(i,:) = abs([detected_scvs.(string(param_name))]) >= tested_vals(i);
        i = i + 1;
    else
        inds(i,:) = abs([detected_scvs.(string(param_name))]) <= tested_vals(i);
        i = i + 1;
    end
end
end



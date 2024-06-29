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
    load(string(input_path) + string(snapshot_dates{ii}) + '/cyclone_data.mat')
    load(string(input_path) + string(snapshot_dates{ii}) + '/lilly_data_final');
    load(string(input_path) + string(snapshot_dates{ii}) + '/background_data')
    date = snapshot_dates{ii};

     %%% Looping through tracks
    for tag_no = 1:467
        for i = cyclones(tag_no).dha
            
            spice_gauss_tag = spice_gauss{tag_no};
            dha_gauss_tag = dha_gauss{tag_no};

            if median(dha_gauss_tag(i).dataX) > 0
                continue
            end

            cyclone_data(u).date = date;
            cyclone_data(u).tag_no = tag_no;
            cyclone_data(u).cast = i;
            cyclone_data(u).bathymetry = LLCsealdata(tag_no).bathymetry(i);

            %%% Spice gauss
            cyclone_data(u).anom_A = spice_gauss_tag(i).A;
            cyclone_data(u).anom_R2 = spice_gauss_tag(i).R2;
            cyclone_data(u).anom_nrmse = spice_gauss_tag(i).nrmse;
            cyclone_data(u).anom_H = spice_gauss_tag(i).H;
            cyclone_data(u).anom_Hcore = spice_gauss_tag(i).Hcore;
            cyclone_data(u).anom_P = spice_gauss_tag(i).P;
            cyclone_data(u).anom_Plow = spice_gauss_tag(i).Plow;
            cyclone_data(u).anom_Phih = spice_gauss_tag(i).Phih;

            %%% DHA gauss
            cyclone_data(u).dha_A = dha_gauss_tag(i).A;
            cyclone_data(u).dha_R2 = dha_gauss_tag(i).R2;
            cyclone_data(u).dha_nrmse = dha_gauss_tag(i).nrmse;
            cyclone_data(u).dha_H = dha_gauss_tag(i).H;
            cyclone_data(u).dha_Hcore = dha_gauss_tag(i).Hcore;
            cyclone_data(u).dha_P = dha_gauss_tag(i).P;
            cyclone_data(u).dha_Plow = dha_gauss_tag(i).Plow;
            cyclone_data(u).dha_Phih = dha_gauss_tag(i).Phih;

            %%% Background data
            cyclone_data(u).bathymetric_var = background(tag_no).bathymetric_var(i);
            cyclone_data(u).shelf_break_ratio = background(tag_no).shelf_break_ratio(i);
            cyclone_data(u).isopycnal_stability = background(tag_no).isopycnal_stablity(i);
            cyclone_data(u).MLD = background(tag_no).MLD(i);
            cyclone_data(u).spice_std = background(tag_no).spice_std(i);
            cyclone_data(u).max_pres = background(tag_no).max_pres(i);

            %%% Lilly SCV data
            cyclone_data(u).True_SCV = lilly_data(tag_no).scv_deep(i);
            cyclone_data(u).scv_reason = lilly_data(tag_no).scv_reason_deep(i);
            if cyclone_data(u).True_SCV == 1
                cyclone_data(u).OW_in_contour = lilly_data(tag_no).contourdata(i).OW_in_contour;
                cyclone_data(u).area = lilly_data(tag_no).contourdata(i).area;
                cyclone_data(u).ecc = lilly_data(tag_no).contourdata(i).ecc;
                cyclone_data(u).area_ratio = lilly_data(tag_no).contourdata(i).area_ratio;
                cyclone_data(u).vort_in_contour = lilly_data(tag_no).contourdata(i).vort_in_contour;
                if (cyclone_data(u).vort_in_contour > 0)
                    cyclone_data(u).True_cyclone = 1;
                else
                    cyclone_data(u).True_cyclone = 0;
                end
            else
                cyclone_data(u).OW_in_contour = NaN;
                cyclone_data(u).area = NaN;
                cyclone_data(u).True_cyclone = 0;
            end
            u = u + 1;

        end
    end

    %%% Missed Cyclones
    for tag_no = 1:length(lilly_data)
        detected = [cyclones(tag_no).dha];
        OW_flagged = find([lilly_data(tag_no).contourdata.vort_in_contour] > 0);
        ind = setdiff(OW_flagged, detected);
        for i = 1:length(ind)
            missed_cyclones(uu).date = date;
            missed_cyclones(uu).tag_no = tag_no;
            missed_cyclones(uu).cast = ind(i);

            %%% Background data
            missed_cyclones(uu).bathymetric_var = background(tag_no).bathymetric_var(ind(i));
            missed_cyclones(uu).shelf_break_ratio = background(tag_no).shelf_break_ratio(ind(i));
            missed_cyclones(uu).isopycnal_stability = background(tag_no).isopycnal_stablity(ind(i));
            missed_cyclones(uu).MLD = background(tag_no).MLD(ind(i));
            missed_cyclones(uu).spice_std = background(tag_no).spice_std(i);
            missed_cyclones(uu).max_pres = background(tag_no).max_pres(i);

            %%% Lilly data
            missed_cyclones(uu).mean_OW = lilly_data(tag_no).contourdata(ind(i)).OW_in_contour;
            missed_cyclones(uu).area = lilly_data(tag_no).contourdata(ind(i)).area;
            uu = uu + 1;
        end
     end

    clear cyclones spice_gauss dha_gauss spice_gauss_tag dha_gauss_tag lilly_data background
end

detected_cyclones = cyclone_data;
save(string(output_path) + '/LLCcyclones.mat', 'detected_cyclones', 'missed_cyclones', '-v7.3')

%% Extra requirements for SCV classification (optional)!!

output_path = '/Volumes/Elements/MEOPdata';
load(string(output_path) + '/LLCcyclones.mat')

OW_requirement = 0.1e-9;
for i = 1:length(detected_cyclones)
    if abs(detected_cyclones(i).OW_in_contour) < OW_requirement
        detected_cyclones(i).True_cyclone = 0;
    end
end
clear i

%%% Shallow or Deep?
deep = 1;
depth_threshold = 1000;

if deep == 1
    detected_cyclones = detected_cyclones(abs([detected_cyclones.bathymetry]) >= depth_threshold);
else
    detected_cyclones = detected_cyclones(abs([detected_cyclones.bathymetry]) < depth_threshold);
end

%%% Removing edge cases

%%% Loading seal data
load('qc_ts.mat');

for u = 1:length(detected_cyclones)
    tag_no = detected_cyclones(u).tag_no;
    i = detected_cyclones(u).cast;

    first_cast = qc_ts(tag_no).cast(1);
    last_cast = qc_ts(tag_no).cast(end);

    if i <= (first_cast + 9) | i >= (last_cast - 9)
        ind(u) = 1;
    else
        ind(u) = 0;
    end
end
detected_cyclones = detected_cyclones(~ind);

%%% Formatting data for optimization
detected_scvs = detected_cyclones;
missed_scvs = missed_cyclones;

%%% Additional requirements
detected_scvs = detected_scvs([detected_scvs.anom_A] < 0); %%% Minty Cyclones
detected_scvs = detected_scvs([detected_scvs.max_pres] >= 400);

%%% Getting original results
detected_scvs_og = detected_scvs;
missed_scvs_og = missed_scvs;
OG_TP = length(detected_scvs_og([detected_scvs_og.True_cyclone] == 1));
OG_FP = length(detected_scvs_og([detected_scvs_og.True_cyclone] == 0));
OG_FN = length(missed_scvs_og);

fn = fieldnames(detected_scvs);
fn = fn(4:26);

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

            TP = length(detected_scvs_opt([detected_scvs_opt.True_cyclone] == 1));
            FP = length(detected_scvs_opt([detected_scvs_opt.True_cyclone] == 0));
            FN = OG_FN + (OG_TP - TP);

            stats(u).precision = TP / (TP + FP);
            stats(u).recall = TP / (TP + FN);
            ii = 1;
            beta = betas_all(ii);
            stats(u).f_2 = (1+beta^2) * (stats(u).precision * stats(u).recall) / ((beta^2 * stats(u).precision) + stats(u).recall);
            ii = ii + 1;
            beta = betas_all(ii);
            stats(u).f_3 = (1+beta^2) * (stats(u).precision * stats(u).recall) / ((beta^2 * stats(u).precision) + stats(u).recall);
            ii = ii + 1;
            beta = betas_all(ii);
            stats(u).f_4 = (1+beta^2) * (stats(u).precision * stats(u).recall) / ((beta^2 * stats(u).precision) + stats(u).recall);
            ii = ii + 1;
            beta = betas_all(ii);
            stats(u).f_5 = (1+beta^2) * (stats(u).precision * stats(u).recall) / ((beta^2 * stats(u).precision) + stats(u).recall);
            stats(u).no_detections = length(detected_scvs_opt);
            stats(u).val = i;
            stats(u).param = fn{uu};
            u = u + 1;
        end
        clear ii

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

%save('initial_guesses_cyclones.mat', 'guess_all_max', 'guess_all_min');

%%
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

%save('final_stats_cyclones.mat', 'final_stats_max', 'final_stats_min');

clear initial_stats_tmp initial_stats_all_tmp f FN FP i idx ind m missed_scvs missed_scvs_og param stats TP u uu...
    u initial_stats_all detected_scvs_og detected_scvs_opt beta arr n cntr ii names

%%% Getting initial delta values
param_names_min = ["anom_A", "anom_R2",  "anom_P", "anom_Hcore",  "dha_A", "dha_R2", "dha_P", "dha_Hcore"];
param_names_max = [ "spice_std", "bathymetric_var", "isopycnal_stability", "MLD", "anom_nrmse", "anom_P", "anom_Hcore", "dha_nrmse", "dha_P", "dha_Hcore"];

param_names = horzcat(param_names_min, param_names_max);

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

clear deltas_min deltas_max idx param uu 

%%

guess_min = guess_min_og;
guess_max = guess_max_og;
deltas = deltas_og;

% Define hyperparameters for gradient descent
max_iterations = 100; % Maximum number of iterations
if deep == 1
    beta0 = 0.035;
else
    beta0 = 0.025;  
end
guess_min.beta = beta0;
guess_max.beta = beta0;
last_iteration = 5;

% Gradient descent optimization
clear tested_values inds f_all precision recall no_detections deltas_all param_values_all

for i = 1:max_iterations
    disp(i)

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

        % if u <= length(param_names_min)
        %     param_values_all(i,u) = guess_min.(param_names(u));
        % else
        %     param_values_all(i,u) = guess_max.(param_names(u));
        % end
        % deltas_all(i,:) = deltas;

        all_other_params = setdiff(1:length(param_names), u);
        tmp = vertcat(inds{all_other_params,1});
        tmp = double(tmp(5:9:size(tmp,1),:));
        ind_all_other_params = prod(tmp,1);

        for uu = 1:size(inds{u,1})

            ind = ind_all_other_params .* inds{u,1}(uu,:);
            ind = (ind == 1);

            detected_scvs_opt = detected_scvs(ind);

            TP = length(detected_scvs_opt([detected_scvs_opt.True_cyclone] == 1));
            FP = length(detected_scvs_opt([detected_scvs_opt.True_cyclone] == 0));
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

        %%% Updating delta value
        if abs(stats.f(uu) - stats.f(last_iteration)) > 0
            gamma = 1-(abs(stats.f(last_iteration) - stats.f(uu)) / (stats.f(last_iteration)));
            deltas(u) = gamma*deltas(u);
        end

    end

    clear detected_scvs_opt ind
    tmp = vertcat(inds{1:size(inds),1});
    tmp = double(tmp(5:9:size(tmp,1),:));
    ind_all_params = prod(tmp,1);
    ind = (ind_all_params == 1);

    detected_scvs_opt = detected_scvs(ind);
    TP = length(detected_scvs_opt([detected_scvs_opt.True_cyclone] == 1));
    FP = length(detected_scvs_opt([detected_scvs_opt.True_cyclone] == 0));
    FN = OG_FN + (OG_TP - TP);

    precision(i) = TP / (TP + FP);
    recall(i) = TP / (TP + FN);
    no_detections(i) = length(detected_scvs_opt);
    disp(string(precision(i)) + ', ' + string(no_detections(i)))
    beta = beta0;
    f_all(i) = (1+beta^2) * (precision(i) * recall(i)) / ((beta^2 * precision(i)) + recall(i));
    final_stats_min = guess_min;
    final_stats_max = guess_max;

end

if deep == 1
    save(string(output_path) + '/final_stats_cyclones_deep', 'final_stats_max', 'final_stats_min', 'guess_max_og', 'guess_min_og', 'param_names_min', 'param_names_max', 'detected_scvs_opt', 'precision', 'f_all', 'OW_requirement', 'prct')
else
    save(string(output_path) + '/final_stats_cyclones_shallow', 'final_stats_max', 'final_stats_min', 'guess_max_og', 'guess_min_og', 'param_names_min', 'param_names_max', 'detected_scvs_opt', 'precision', 'f_all', 'OW_requirement', 'prct')
end

%%

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





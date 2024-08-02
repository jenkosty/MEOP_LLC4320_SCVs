%% Loading data

deep = 1;
depth_threshold = 1000;
input_path = '/Volumes/Elements/MEOPdata';

%%% MEOP anticyclone Data
load(string(input_path) + '/background_data.mat')
load(string(input_path) + '/anticyclone_data.mat')
load('/Users/jenkosty/Documents/Research/SCV_Project/Seal Data/qc_ts.mat')

%%% LLC anticyclone Data
load(string(input_path) + '/LLCanticyclones.mat')
detected_anticyclones = detected_anticyclones([detected_anticyclones.max_pres] >= 400);
detected_anticyclones = detected_anticyclones(find(~cellfun(@isempty,{detected_anticyclones.spice_A})));
detected_anticyclones = detected_anticyclones([detected_anticyclones.spice_P] >= 0);

if deep == 1
    detected_anticyclones = detected_anticyclones(abs([detected_anticyclones.bathymetry]) >= depth_threshold);
elseif deep == 0
    detected_anticyclones = detected_anticyclones(abs([detected_anticyclones.bathymetry]) < depth_threshold);
end
LLCanticyclones_all = detected_anticyclones;

%%% Adding lat/lon data to LLC anticyclones
for uu = 1:length(LLCanticyclones_all)
    tag_no = LLCanticyclones_all(uu).tag_no;
    i = LLCanticyclones_all(uu).cast;
    LLCanticyclones_all(uu).lat = qc_ts(tag_no).lat(i);
    LLCanticyclones_all(uu).lon = qc_ts(tag_no).lon(i);
end

%%% Removing edge cases
clear ind
for u = 1:length(LLCanticyclones_all)
    tag_no = LLCanticyclones_all(u).tag_no;
    i = LLCanticyclones_all(u).cast;

    first_cast = 1;
    last_cast = length(qc_ts(tag_no).cast);

    if i <= (first_cast + 9) | i >= (last_cast - 9)
        ind(u) = 1;
    else
        ind(u) = 0;
    end
end
LLCanticyclones_all = LLCanticyclones_all(~ind);

%%% Formating MEOP Data
u = 1;
for tag_no = 1:length(qc_ts)
    for i = anticyclones(tag_no).dha

        isa_gauss_tag = isa_gauss{tag_no};
        spice_gauss_tag = spice_gauss{tag_no};
        dha_gauss_tag = dha_gauss{tag_no};

        % %%% Excluding detections based on pressure difference of fits
        % pres_diff = abs(isa_gauss_tag(i).P - dha_gauss_tag(i).P);
        % if pres_diff > 100
        %     continue
        % end

        anticyclone_data(u).date = qc_ts(tag_no).time(i,:);
        anticyclone_data(u).tag_no = tag_no;
        anticyclone_data(u).cast = i;
        anticyclone_data(u).lat = qc_ts(tag_no).lat(i);
        anticyclone_data(u).lon = qc_ts(tag_no).lon(i);
        anticyclone_data(u).bathymetry = qc_ts(tag_no).bathymetry(i);

        %%% ISA gauss
        anticyclone_data(u).anom_A = isa_gauss_tag(i).A;
        anticyclone_data(u).anom_R2 = isa_gauss_tag(i).R2;
        anticyclone_data(u).anom_nrmse = isa_gauss_tag(i).nrmse;
        anticyclone_data(u).anom_H = isa_gauss_tag(i).H;
        anticyclone_data(u).anom_Hcore = isa_gauss_tag(i).Hcore;
        anticyclone_data(u).anom_P = isa_gauss_tag(i).P;
        anticyclone_data(u).anom_Plow = isa_gauss_tag(i).Plow;
        anticyclone_data(u).anom_Phih = isa_gauss_tag(i).Phih;

        %%% Spice gauss
        if i > length(spice_gauss_tag)
            anticyclone_data(u).spice_A = NaN;
            anticyclone_data(u).spice_R2 = NaN;
            anticyclone_data(u).spice_nrmse = NaN;
            anticyclone_data(u).spice_H = NaN;
            anticyclone_data(u).spice_Hcore = NaN;
            anticyclone_data(u).spice_P = NaN;
            anticyclone_data(u).spice_Plow = NaN;
            anticyclone_data(u).spice_Phih = NaN;
        else
            anticyclone_data(u).spice_A = spice_gauss_tag(i).A;
            anticyclone_data(u).spice_R2 = spice_gauss_tag(i).R2;
            anticyclone_data(u).spice_nrmse = spice_gauss_tag(i).nrmse;
            anticyclone_data(u).spice_H = spice_gauss_tag(i).H;
            anticyclone_data(u).spice_Hcore = spice_gauss_tag(i).Hcore;
            anticyclone_data(u).spice_P = spice_gauss_tag(i).P;
            anticyclone_data(u).spice_Plow = spice_gauss_tag(i).Plow;
            anticyclone_data(u).spice_Phih = spice_gauss_tag(i).Phih;
        end

        %%% DHA gauss
        anticyclone_data(u).dha_A = dha_gauss_tag(i).A;
        anticyclone_data(u).dha_R2 = dha_gauss_tag(i).R2;
        anticyclone_data(u).dha_nrmse = dha_gauss_tag(i).nrmse;
        anticyclone_data(u).dha_H = dha_gauss_tag(i).H;
        anticyclone_data(u).dha_Hcore = dha_gauss_tag(i).Hcore;
        anticyclone_data(u).dha_P = dha_gauss_tag(i).P;
        anticyclone_data(u).dha_Plow = dha_gauss_tag(i).Plow;
        anticyclone_data(u).dha_Phih = dha_gauss_tag(i).Phih;

        %%% Background data
        anticyclone_data(u).bathymetric_var = background(tag_no).bathymetric_var(i);
        anticyclone_data(u).shelf_break_ratio = background(tag_no).shelf_break_ratio(i);
        anticyclone_data(u).isopycnal_stability = background(tag_no).isopycnal_stablity(i);
        anticyclone_data(u).MLD = background(tag_no).MLD(i);
        anticyclone_data(u).spice_std = background(tag_no).spice_std(i);
        anticyclone_data(u).max_pres = background(tag_no).max_pres(i);

        u = u + 1;
    end
end
MEOPanticyclones_all = anticyclone_data;

%%% Removing edge cases
clear ind
for u = 1:length(MEOPanticyclones_all)
    tag_no = MEOPanticyclones_all(u).tag_no;
    i = MEOPanticyclones_all(u).cast;

    first_cast = 1;
    last_cast = length(qc_ts(tag_no).cast);

    if i <= (first_cast + 9) | i >= (last_cast - 9)
        ind(u) = 1;
    else
        ind(u) = 0;
    end
end
MEOPanticyclones_all = MEOPanticyclones_all(~ind);

%%% Extra requirements (spice anomaly and bathymetry)
if deep == 1
    MEOPanticyclones_all = MEOPanticyclones_all(abs([MEOPanticyclones_all.bathymetry]) >= depth_threshold);
elseif deep == 0
    MEOPanticyclones_all = MEOPanticyclones_all(abs([MEOPanticyclones_all.bathymetry]) < depth_threshold);
end
MEOPanticyclones_all = MEOPanticyclones_all(find(~cellfun(@isempty,{MEOPanticyclones_all.spice_A})));
MEOPanticyclones_all = MEOPanticyclones_all([MEOPanticyclones_all.spice_P] > 0);

clear anticyclone_data u spice_gauss_tag dha_gauss_tag dha_gauss anticyclones background spice_gauss...
    tag_no i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Applying optimization %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% SEPATATE FOR LLC AND MEOP DATA

%%% Anticyclone Optimization Parameters
if deep == 1
    load(string(input_path) + '/final_stats_anticyclones_deep.mat')
else
    load(string(input_path) + '/final_stats_anticyclones_shallow.mat')
end
param_names = horzcat(param_names_min, param_names_max);

%%% LLC Snapshots
ind = ones(1,length(LLCanticyclones_all));
for uu = 1:length(param_names)
    if uu <= length(param_names_min)
        ind_tmp = abs([LLCanticyclones_all.(param_names(uu))]) >= final_stats_min.(param_names(uu));
        ind = ind .* double(ind_tmp);
    else
        ind_tmp = abs([LLCanticyclones_all.(param_names(uu))]) <= final_stats_max.(param_names(uu));
        ind = ind .* double(ind_tmp);
    end
end
ind = (ind == 1);
LLCanticyclones_optimized = LLCanticyclones_all(ind);

%%% Anticyclone Optimization Parameters
if deep == 1
    load(string(input_path) + '/final_stats_anticyclones_deep_reasonable.mat')
else
    load(string(input_path) + '/final_stats_anticyclones_shallow_reasonable.mat')
end
param_names = horzcat(param_names_min, param_names_max);

%%% MEOP Snapshots 
ind = ones(1,length(MEOPanticyclones_all));
for uu = 1:length(param_names)
    if uu <= length(param_names_min)
        ind_tmp = abs([MEOPanticyclones_all.(param_names(uu))]) >= final_stats_min.(param_names(uu));
        ind = ind .* double(ind_tmp);
        disp(sum(ind))
    else
        ind_tmp = abs([MEOPanticyclones_all.(param_names(uu))]) <= final_stats_max.(param_names(uu));
        ind = ind .* double(ind_tmp);
        disp(sum(ind))
    end
end
ind = (ind == 1);
MEOPanticyclones_optimized = MEOPanticyclones_all(ind);

%clear ind LLCanticyclones_all MEOPanticyclones_all

%%% Getting Vertical Profiles of LLC anticyclones (Takes a little bit of time...)

%%% Formating LLC Data
snapshot_dates = {'01-Oct-2011','01-Nov-2011', '01-Dec-2011', '01-Jan-2012', '01-Feb-2012', '01-Mar-2012', '01-Apr-2012', '01-May-2012', ...
    '01-Jun-2012', '01-Jul-2012', '01-Aug-2012', '01-Sep-2012'};

u = 1;
for ii = 1:length(snapshot_dates)
    disp(snapshot_dates{ii})
    ind = strcmp(snapshot_dates{ii}, {LLCanticyclones_optimized.date});
    snapshot_data = LLCanticyclones_optimized(ind);
    load('/Volumes/Elements/LLCsealdata/Snapshot_' + string(snapshot_dates{ii}) + '/anticyclone_data.mat')
    
    for uu = 1:length(snapshot_data)
        tag_no = snapshot_data(uu).tag_no;
        i = snapshot_data(uu).cast;

        %%% Fitted
        LLCanticyclones_optimized(u).isa_anom_profile_fitted = isa_gauss{1,tag_no}(i).X;
        LLCanticyclones_optimized(u).dha_anom_profile_fitted = dha_gauss{1,tag_no}(i).X;
        LLCanticyclones_optimized(u).spice_anom_profile_fitted = spice_gauss{1,tag_no}(i).X;
        LLCanticyclones_optimized(u).isa_pres_profile_fitted = isa_gauss{1,tag_no}(i).Y;
        LLCanticyclones_optimized(u).dha_pres_profile_fitted =  dha_gauss{1,tag_no}(i).Y;
        LLCanticyclones_optimized(u).spice_pres_profile_fitted =  spice_gauss{1,tag_no}(i).Y;

        %%% Raw
        LLCanticyclones_optimized(u).isa_anom_profile = isa_gauss{1,tag_no}(i).dataX;
        LLCanticyclones_optimized(u).dha_anom_profile = dha_gauss{1,tag_no}(i).dataX;
        LLCanticyclones_optimized(u).spice_anom_profile = spice_gauss{1,tag_no}(i).dataX;
        LLCanticyclones_optimized(u).isa_pres_profile = isa_gauss{1,tag_no}(i).dataY;
        LLCanticyclones_optimized(u).dha_pres_profile =  dha_gauss{1,tag_no}(i).dataY;
        LLCanticyclones_optimized(u).spice_pres_profile =  spice_gauss{1,tag_no}(i).dataY;

        if deep == 1
            LLCanticyclones_optimized(u).type = "deep";
        elseif deep == 0
            LLCanticyclones_optimized(u).type = "shallow";
        end

        u = u + 1;
    end

    clear dha_gauss anticyclones spice_gauss uu tag_no i

end

clear u ind snapshot_data snapshot_dates date ii

%%% Getting vertical profiles of MEOP anticyclones

load(string(input_path) + '/anticyclone_data.mat')
load(string(input_path) + '/qc_ts_full.mat')

%%% Formatting MEOP anticyclones
for u = 1:length(MEOPanticyclones_optimized)
    tag_no = MEOPanticyclones_optimized(u).tag_no;
    i = MEOPanticyclones_optimized(u).cast;

    %%% Fitted
    MEOPanticyclones_optimized(u).isa_anom_profile_fitted = isa_gauss{1,tag_no}(i).X;
    MEOPanticyclones_optimized(u).dha_anom_profile_fitted = dha_gauss{1,tag_no}(i).X;
    MEOPanticyclones_optimized(u).spice_anom_profile_fitted = spice_gauss{1,tag_no}(i).X;
    MEOPanticyclones_optimized(u).isa_pres_profile_fitted = isa_gauss{1,tag_no}(i).Y;
    MEOPanticyclones_optimized(u).dha_pres_profile_fitted =  dha_gauss{1,tag_no}(i).Y;
    MEOPanticyclones_optimized(u).spice_pres_profile_fitted =  spice_gauss{1,tag_no}(i).Y;
    
    %%% Raw
    MEOPanticyclones_optimized(u).isa_anom_profile = isa_gauss{1,tag_no}(i).dataX;
    MEOPanticyclones_optimized(u).dha_anom_profile = dha_gauss{1,tag_no}(i).dataX;
    MEOPanticyclones_optimized(u).spice_anom_profile = spice_gauss{1,tag_no}(i).dataX;
    MEOPanticyclones_optimized(u).isa_pres_profile = isa_gauss{1,tag_no}(i).dataY;
    MEOPanticyclones_optimized(u).dha_pres_profile =  dha_gauss{1,tag_no}(i).dataY;
    MEOPanticyclones_optimized(u).spice_pres_profile =  spice_gauss{1,tag_no}(i).dataY;

    MEOPanticyclones_optimized(u).N2_anom_profile = qc_ts(tag_no).ds.anoms.N2(:,i);
    MEOPanticyclones_optimized(u).salt_anom_profile = qc_ts(tag_no).ds.anoms.salt(:,i);
    MEOPanticyclones_optimized(u).temp_anom_profile = qc_ts(tag_no).ds.anoms.temp(:,i);
    MEOPanticyclones_optimized(u).N2_pres_profile = qc_ts(tag_no).ds.pres(:,i);

    if deep == 1
        MEOPanticyclones_optimized(u).type = "deep";
    elseif deep == 0
        MEOPanticyclones_optimized(u).type = "shallow";
    end

end

clear dha_gauss anticyclones spice_gauss uu tag_no i u ii


LLCanticyclones_optimized = LLCanticyclones_optimized([LLCanticyclones_optimized.dha_A]>=0);
MEOPanticyclones_optimized = MEOPanticyclones_optimized([MEOPanticyclones_optimized.dha_A]>=0);

if deep == 1
   save(string(input_path) + '/optimized_anticyclones_deep.mat', 'MEOPanticyclones_optimized', 'LLCanticyclones_optimized')
else
   save(string(input_path) + '/optimized_anticyclones_shallow.mat', 'MEOPanticyclones_optimized', 'LLCanticyclones_optimized')
end

%% Combining shallow and deep detections (must run both versions first!)

load(string(input_path) + '/optimized_anticyclones_deep.mat')
LLCanticyclones_optimized_tmp1 = LLCanticyclones_optimized;
MEOPanticyclones_optimized_tmp1 = MEOPanticyclones_optimized;

load(string(input_path) + '/optimized_anticyclones_shallow.mat')
LLCanticyclones_optimized_tmp2 = LLCanticyclones_optimized;
MEOPanticyclones_optimized_tmp2 = MEOPanticyclones_optimized;

LLCanticyclones_optimized = [LLCanticyclones_optimized_tmp1 LLCanticyclones_optimized_tmp2];
MEOPanticyclones_optimized = [MEOPanticyclones_optimized_tmp1 MEOPanticyclones_optimized_tmp2];
save(string(input_path) + '/optimized_anticyclones.mat', 'MEOPanticyclones_optimized', 'LLCanticyclones_optimized')


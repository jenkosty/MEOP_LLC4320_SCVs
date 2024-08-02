%% Loading data

deep = 0;
depth_threshold = 1000;
input_path = '/Volumes/Elements/MEOPdata';

%%% Cyclone Optimization Parameters
if deep == 1
    load(string(input_path) + '/final_stats_cyclones_deep.mat')
else
    load(string(input_path) + '/final_stats_cyclones_shallow.mat')
end

%%% MEOP Cyclone Data
load(string(input_path) + '/background_data.mat')
load(string(input_path) + '/cyclone_data.mat')
load(string(input_path) + '/qc_ts.mat')

%%% LLC Cyclone Data
load(string(input_path) + '/LLCcyclones.mat')
detected_cyclones = detected_cyclones([detected_cyclones.max_pres] >= 400);
if deep == 1
    detected_cyclones = detected_cyclones(abs([detected_cyclones.bathymetry]) >= depth_threshold);
else
    detected_cyclones = detected_cyclones(abs([detected_cyclones.bathymetry]) < depth_threshold);
end
LLCcyclones_all = detected_cyclones;

%%% Adding lat/lon data to LLC Cyclones
for uu = 1:length(LLCcyclones_all)
    tag_no = LLCcyclones_all(uu).tag_no;
    i = LLCcyclones_all(uu).cast;
    LLCcyclones_all(uu).lat = qc_ts(tag_no).lat(i);
    LLCcyclones_all(uu).lon = qc_ts(tag_no).lon(i);
end

%%% Removing edge cases
clear ind
for u = 1:length(LLCcyclones_all)
    tag_no = LLCcyclones_all(u).tag_no;
    i = LLCcyclones_all(u).cast;

    first_cast = 1;
    last_cast = length(qc_ts(tag_no).cast);

    if i <= (first_cast + 9) | i >= (last_cast - 9)
        ind(u) = 1;
    else
        ind(u) = 0;
    end
end
LLCcyclones_all = LLCcyclones_all(~ind);

clear detected_cyclones uu i tag_no

%%% Formating MEOP Data
u = 1;
for tag_no = 1:length(qc_ts)
    for i = cyclones(tag_no).dha

        spice_gauss_tag = spice_gauss{tag_no};
        dha_gauss_tag = dha_gauss{tag_no};

        cyclone_data(u).date = qc_ts(tag_no).time(i,:);
        cyclone_data(u).tag_no = tag_no;
        cyclone_data(u).cast = i;
        cyclone_data(u).lat = qc_ts(tag_no).lat(i);
        cyclone_data(u).lon = qc_ts(tag_no).lon(i);
        cyclone_data(u).bathymetry = qc_ts(tag_no).bathymetry(i);

        %%% Spice gauss
        cyclone_data(u).anom_A = spice_gauss_tag(i).A;
        cyclone_data(u).anom_R2 = spice_gauss_tag(i).R2;
        cyclone_data(u).anom_nrmse = spice_gauss_tag(i).nrmse;
        cyclone_data(u).anom_Hcore = spice_gauss_tag(i).Hcore;
        cyclone_data(u).anom_P = spice_gauss_tag(i).P;
        cyclone_data(u).anom_Plow = spice_gauss_tag(i).Plow;
        cyclone_data(u).anom_Phih = spice_gauss_tag(i).Phih;

        %%% DHA gauss
        cyclone_data(u).dha_A = dha_gauss_tag(i).A;
        cyclone_data(u).dha_R2 = dha_gauss_tag(i).R2;
        cyclone_data(u).dha_nrmse = dha_gauss_tag(i).nrmse;
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

        u = u + 1;
    end
end

MEOPcyclones_all = cyclone_data;

%%% Removing edge cases
clear ind
for u = 1:length(MEOPcyclones_all)
    tag_no = MEOPcyclones_all(u).tag_no;
    i = MEOPcyclones_all(u).cast;

    first_cast = 1;
    last_cast = length(qc_ts(tag_no).cast);

    if i <= (first_cast + 9) | i >= (last_cast - 9)
        ind(u) = 1;
    else
        ind(u) = 0;
    end
end
MEOPcyclones_all = MEOPcyclones_all(~ind);

%%% Extra requirements (spice anomaly and bathymetry)
if deep == 1
    MEOPcyclones_all = MEOPcyclones_all(abs([MEOPcyclones_all.bathymetry]) >= depth_threshold);
else
    MEOPcyclones_all = MEOPcyclones_all(abs([MEOPcyclones_all.bathymetry]) < depth_threshold);
end

clear cyclone_data u spice_gauss_tag dha_gauss_tag dha_gauss cyclones background spice_gauss...
    tag_no i

%%% Applying optimization

param_names = horzcat(param_names_min, param_names_max);

%%% LLC Snapshots
ind = ones(1,length(LLCcyclones_all));
for uu = 1:length(param_names)
    if uu <= length(param_names_min)
        ind_tmp = abs([LLCcyclones_all.(param_names(uu))]) >= final_stats_min.(param_names(uu));
        ind = ind .* double(ind_tmp);
    else
        ind_tmp = abs([LLCcyclones_all.(param_names(uu))]) <= final_stats_max.(param_names(uu));
        ind = ind .* double(ind_tmp);
    end
end
ind = (ind == 1);
LLCcyclones_optimized = LLCcyclones_all(ind);

%%% MEOP Snapshots
ind = ones(1,length(MEOPcyclones_all));
for uu = 1:length(param_names)
    if uu <= length(param_names_min)
        ind_tmp = abs([MEOPcyclones_all.(param_names(uu))]) >= final_stats_min.(param_names(uu));
        ind = ind .* double(ind_tmp);
    else
        ind_tmp = abs([MEOPcyclones_all.(param_names(uu))]) <= final_stats_max.(param_names(uu));
        ind = ind .* double(ind_tmp);
    end
end
ind = (ind == 1);
MEOPcyclones_optimized = MEOPcyclones_all(ind);


%%% Getting Vertical Profiles of LLC cyclones (Takes a little bit of time...)

snapshot_dates = {'01-Oct-2011','01-Nov-2011', '01-Dec-2011', '01-Jan-2012', '01-Feb-2012', '01-Mar-2012', '01-Apr-2012', '01-May-2012', ...
    '01-Jun-2012', '01-Jul-2012', '01-Aug-2012', '01-Sep-2012'};

u = 1;
for ii = 1:length(snapshot_dates)
    disp(snapshot_dates{ii})
    ind = strcmp(snapshot_dates{ii}, {LLCcyclones_optimized.date});
    snapshot_data = LLCcyclones_optimized(ind);
    load('/Volumes/Elements/LLCsealdata/Snapshot_' + string(snapshot_dates{ii}) + '/cyclone_data.mat')
    load('/Volumes/Elements/LLCsealdata/Snapshot_' + string(snapshot_dates{ii}) + '/LLCsealdata_full.mat')
    
    for uu = 1:length(snapshot_data)

        tag_no = snapshot_data(uu).tag_no;
        i = snapshot_data(uu).cast;

        %%% Fitted
        LLCcyclones_optimized(u).dha_anom_profile_fitted = dha_gauss{1,tag_no}(i).X;
        LLCcyclones_optimized(u).spice_anom_profile_fitted = spice_gauss{1,tag_no}(i).X;
        LLCcyclones_optimized(u).dha_pres_profile_fitted =  dha_gauss{1,tag_no}(i).Y;
        LLCcyclones_optimized(u).spice_pres_profile_fitted =  spice_gauss{1,tag_no}(i).Y;

        %%% Raw
        LLCcyclones_optimized(u).dha_anom_profile = dha_gauss{1,tag_no}(i).dataX;
        LLCcyclones_optimized(u).spice_anom_profile = spice_gauss{1,tag_no}(i).dataX;
        LLCcyclones_optimized(u).dha_pres_profile =  dha_gauss{1,tag_no}(i).dataY;
        LLCcyclones_optimized(u).spice_pres_profile =  spice_gauss{1,tag_no}(i).dataY;

        %%% 
        LLCcyclones_optimized(u).N2_anom_profile = LLCsealdata(tag_no).ds.anoms.N2(:,i);
        LLCcyclones_optimized(u).salt_anom_profile = LLCsealdata(tag_no).ds.anoms.salt(:,i);
        LLCcyclones_optimized(u).temp_anom_profile = LLCsealdata(tag_no).ds.anoms.temp(:,i);
        LLCcyclones_optimized(u).N2_pres_profile = LLCsealdata(tag_no).ds.pres(:,i);
        LLCcyclones_optimized(u).vort_profile = LLCsealdata(tag_no).vort(:,i);

        if deep == 1
            LLCcyclones_optimized(u).type = "deep";
        else
            LLCcyclones_optimized(u).type = "shallow";
        end
        
        u = u + 1;
    end

    clear dha_gauss cyclones spice_gauss uu tag_no i

end

clear u ind snapshot_data snapshot_dates date ii

%%% Getting vertical profiles of MEOP cyclones

%%% Formatting MEOP data
load(string(input_path) + '/cyclone_data.mat')
load(string(input_path) + '/qc_ts_full.mat')

for u = 1:length(MEOPcyclones_optimized)
    tag_no = MEOPcyclones_optimized(u).tag_no;
    i = MEOPcyclones_optimized(u).cast;

    %%% Fitted
    MEOPcyclones_optimized(u).dha_anom_profile_fitted = dha_gauss{1,tag_no}(i).X;
    MEOPcyclones_optimized(u).spice_anom_profile_fitted = spice_gauss{1,tag_no}(i).X;
    MEOPcyclones_optimized(u).dha_pres_profile_fitted =  dha_gauss{1,tag_no}(i).Y;
    MEOPcyclones_optimized(u).spice_pres_profile_fitted =  spice_gauss{1,tag_no}(i).Y;

    %%% Raw
    MEOPcyclones_optimized(u).dha_anom_profile = dha_gauss{1,tag_no}(i).dataX;
    MEOPcyclones_optimized(u).spice_anom_profile = spice_gauss{1,tag_no}(i).dataX;
    MEOPcyclones_optimized(u).dha_pres_profile =  dha_gauss{1,tag_no}(i).dataY;
    MEOPcyclones_optimized(u).spice_pres_profile =  spice_gauss{1,tag_no}(i).dataY;

    MEOPcyclones_optimized(u).N2_anom_profile = qc_ts(tag_no).ds.anoms.N2(:,i);
    MEOPcyclones_optimized(u).salt_anom_profile = qc_ts(tag_no).ds.anoms.salt(:,i);
    MEOPcyclones_optimized(u).temp_anom_profile = qc_ts(tag_no).ds.anoms.temp(:,i);
    MEOPcyclones_optimized(u).N2_pres_profile = qc_ts(tag_no).ds.pres(:,i);

    if deep == 1
       MEOPcyclones_optimized(u).type = "deep";
    else
       MEOPcyclones_optimized(u).type = "shallow";
    end
end

%%% Saving data

clear dha_gauss cyclones spice_gauss uu tag_no i u ii

LLCcyclones_optimized = LLCcyclones_optimized([LLCcyclones_optimized.dha_A]<=0);
MEOPcyclones_optimized = MEOPcyclones_optimized([MEOPcyclones_optimized.dha_A]<=0);

if deep == 1
   save(string(input_path) + '/optimized_cyclones_deep.mat', 'MEOPcyclones_optimized', 'LLCcyclones_optimized')
else
   save(string(input_path) + '/optimized_cyclones_shallow.mat', 'MEOPcyclones_optimized', 'LLCcyclones_optimized')
end

%% Combining shallow and deep detections (must run both versions first!)

input_path = '/Volumes/Elements/MEOPdata';

load(string(input_path) + '/optimized_cyclones_deep.mat')
LLCcyclones_optimized_tmp1 = LLCcyclones_optimized;
MEOPcyclones_optimized_tmp1 = MEOPcyclones_optimized;

load(string(input_path) + '/optimized_cyclones_shallow.mat')
LLCcyclones_optimized_tmp2 = LLCcyclones_optimized;
MEOPcyclones_optimized_tmp2 = MEOPcyclones_optimized;

LLCcyclones_optimized = [LLCcyclones_optimized_tmp1 LLCcyclones_optimized_tmp2];
MEOPcyclones_optimized = [MEOPcyclones_optimized_tmp1 ];
save(string(input_path) + '/optimized_cyclones.mat', 'MEOPcyclones_optimized', 'LLCcyclones_optimized')

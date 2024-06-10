%%% LLC4320 snapshot dates
snapshot_dates = {'01-Oct-2011','01-Nov-2011', '01-Dec-2011', '01-Jan-2012', '01-Feb-2012', '01-Mar-2012', '01-Apr-2012', '01-May-2012', ...
    '01-Jun-2012', '01-Jul-2012', '01-Aug-2012', '01-Sep-2012'};
date = snapshot_dates{1};

input_path = '/Volumes/Elements/LLCsealdata/Snapshot_';
load(string(input_path) + string(date) + '/LLCsealdata');

output_path = '/Volumes/Elements/LLCsealdata/Snapshot_';

ts_data = LLCsealdata;
modeldata = 1;
test_prof = 1:length(LLCsealdata);

%%% Loading algorithm settings
run('LLCseals_algorithm_settings.m')

%%
%%% Prepping time series for detection algorithm(s) (calculating variables
%%% etc.)
ts_data = preppingTimeSeriesForDetectionAlgorithm(ts_data, depth_grid, test_prof, prms, modeldata);

%%% Saving data
if modeldata == 1
    LLCsealdata = ts_data;
    % save(output_path + string(date) + '/LLCsealdata_full', 'LLCsealdata', 'depth_grid', '-v7.3')
else
    qc_ts = ts_data;
    % save(output_path'/qc_ts_full', 'qc_ts', 'depth_grid', '-v7.3')
end

%%
%%% Anticyclonic-specific algorithm checks
ts_data = anticyclonic_checks(LLCsealdata, test_prof, prms);

%%% Extracting anticyclone-specific data
isa_gauss = cell(1,length(test_prof));
dha_gauss = cell(1,length(test_prof));
spice_gauss = cell(1,length(test_prof));
for tag_no = 1:max(test_prof)
    isa_gauss{tag_no} = ts_data(tag_no).isa_gauss;
    dha_gauss{tag_no} = ts_data(tag_no).dha_gauss_anticyclones;
    spice_gauss{tag_no} = ts_data(tag_no).spice_gauss_anticyclones;
    anticyclones(tag_no) = ts_data(tag_no).anticyclones;
end

%%% Saving data
% if modeldata == 1
%     date = ts_data(1).date;
%     save(string(output_data) + string(date) + '/anticyclone_data.mat', 'isa_gauss','dha_gauss', 'spice_gauss', 'anticyclones', 'date', '-v7.3')
% end

clear isa_gauss dha_gauss spice_gauss anticyclones date

%%
%%% Cyclonic-specific algorithm checks
ts_data = cyclonic_checks(LLCsealdata, test_prof, prms);

%%% Extracting cyclonic-specific data
spice_gauss = cell(1,length(test_prof));
dha_gauss = cell(1,length(test_prof));
for tag_no = 1:max(test_prof)
    spice_gauss{tag_no} = ts_data(tag_no).spice_gauss;
    dha_gauss{tag_no} = ts_data(tag_no).dha_gauss_cyclones;
    cyclones(tag_no) = ts_data(tag_no).cyclones;
end

%%% Saving data
% if modeldata == 1
%     date = ts_data(1).date;
%     save(string(output_data) + string(date) + '/cyclone_data.mat', 'spice_gauss', 'dha_gauss', 'cyclones', 'date', '-v7.3')
% end

clear spice_gauss dha_gauss cyclones date


%%
%%% Background-specific algorithm checks
ts_data = background_checks(LLCsealdata, test_prof);

%%% Extracting background-specific data
for tag_no = 1:max(test_prof)
    background(tag_no).bathymetric_var = ts_data(tag_no).bathymetric_var;
    background(tag_no).shelf_break_ratio = ts_data(tag_no).shelf_break_ratio;
    background(tag_no).jump = ts_data(tag_no).jump;
    background(tag_no).isopycnal_stablity = ts_data(tag_no).isopycnal_stability;
    background(tag_no).MLD = ts_data(tag_no).MLD;
    background(tag_no).spice_std = ts_data(tag_no).spice_std;
    background(tag_no).max_pres = ts_data(tag_no).max_pres;
end

%%% Saving data
% if modeldata == 1
%     date = ts_data(1).date;
%     save(string(output_data) + string(date) + '/background_data.mat', 'background', 'date', '-v7.3')
% end

%clear background date




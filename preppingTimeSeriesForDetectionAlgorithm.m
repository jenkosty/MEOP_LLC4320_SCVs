function ts_data = preppingTimeSeriesForDetectionAlgorithm(ts_data, depth_grid, test_prof, prms, modeldata)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating Variables %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('---------------------')
disp('Calculating Variables')
disp('---------------------')

for tag_no = test_prof

    %%% Creating arrays to hold rejection and justification data
    rejected.anticyclones = zeros(size(ts_data(tag_no).cast)); % Anticyclones
    reason.anticyclones = strings(size(ts_data(tag_no).cast));
    rejected.cyclones = zeros(size(ts_data(tag_no).cast));     % Cyclones
    reason.cyclones = strings(size(ts_data(tag_no).cast));
    ts_data(tag_no).rejected = rejected;                       
    ts_data(tag_no).reason = reason;     

    %%% Creating array to hold anticyclonic detections
    anticyclones.isa_iqr = [];
    anticyclones.isa_gaussian = [];
    anticyclones.dha = [];
    ts_data(tag_no).anticyclones = anticyclones;

    %%% Creating array to hold cyclonic detections
    cyclones.spice_iqr = [];
    cyclones.spice_gaussian = [];
    cyclones.dha = [];
    ts_data(tag_no).cyclones = cyclones;

    %%% Calculating BedMachineAntarctica bathymetry
    [x,y] = ll2xy(ts_data(tag_no).lat, ts_data(tag_no).lon, -1);
    ts_data(tag_no).bathymetry = interpBedmachineAntarctica(x, y, 'bed');

    %%% Calculating pressure space variables
    ts_data(tag_no).ps = calc_pres_space_vars(ts_data, tag_no, depth_grid, modeldata);

end

clear rejected reason x y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Checking for Density Inversions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-------------------------------')
disp('Checking for Density Inversions')
disp('-------------------------------')

for tag_no = test_prof

    %%% Checking for density inversions. Small inversions are corrected.
    %%% Profiles with large inversions will be excluded.
    [ts_data(tag_no).ps.sigma0, ts_data(tag_no).ps.max_sigma0_inversion,...
        ts_data(tag_no).rejected, ts_data(tag_no).reason]...
        = checking_for_density_inversions(ts_data, tag_no);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Interpolating to Density Grid %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('------------------------------')
disp('Interpolating to Density Space')
disp('------------------------------')

%%% Creating density grid
density_grid = (26.0:0.001:28.5)';

for tag_no = test_prof
    
    %%% Interpolating pressure space variables to density grid
    ts_data(tag_no).ds = interp_to_sigma0_space(ts_data, tag_no, density_grid);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating Isopycnal Separation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--------------------------------')
disp('Calculating Isopycnal Separation')
disp('--------------------------------')

for tag_no = test_prof
    
    %%% Calculating isopycnal separation
    ts_data(tag_no).ds.isopycnal_separation = calc_isopycnal_separation(ts_data, tag_no);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Finding Indices to Build Reference Profiles %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--------------------------------------')
disp('Getting Indices for Reference Profiles')
disp('--------------------------------------')
 
for tag_no = test_prof
        
    %%% Finding indices to build reference profiles
    ts_data(tag_no).ref_ind = indices_for_ref_profiles(ts_data,tag_no, prms.refprof);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Building Reference Profiles %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('---------------------------')
disp('Building Reference Profiles')
disp('---------------------------')

for tag_no = test_prof

    %%% Building density-space reference profiles
    [ts_data(tag_no).ds.ref_salt, ts_data(tag_no).ds.ref_temp,...
        ts_data(tag_no).ds.ref_N2, ts_data(tag_no).ds.ref_spice,...
        ts_data(tag_no).ds.ref_isopycnal_separation] = build_density_space_ref_profiles(ts_data, tag_no);

    %%% Building pressure-space reference profiles
    [ts_data(tag_no).ps.ref_dyn_height_anom, ts_data(tag_no).ps.ref_N2]...
        = build_pres_space_ref_profiles(ts_data, tag_no);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating Anomalies %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('---------------------')
disp('Calculating Anomalies')
disp('---------------------')

for tag_no = test_prof

    %%% Calculating density-space anomalies
    ts_data(tag_no).ds.anoms = calc_density_space_anomalies(ts_data, tag_no);

    %%% Calculating pressure-space anomalies
    ts_data(tag_no).ps.anoms = calc_pres_space_anomalies(ts_data, tag_no);

end

%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating IQR %%%
%%%%%%%%%%%%%%%%%%%%%%%

disp('---------------')
disp('Calculating IQR')
disp('---------------')

for tag_no = test_prof

    %%% Calculating iqr, lower-threshold, and upper threshold
    ts_data(tag_no).ds.iqrs = calc_iqr(ts_data, tag_no);
    
end

end

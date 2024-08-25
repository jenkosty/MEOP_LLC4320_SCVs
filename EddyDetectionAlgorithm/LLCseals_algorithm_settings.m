
clear prms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME SERIES QUALITY CONTROL THRESHOLDS
qc.min_depth    = 350; % minimum depth to be considered a 'good' profile 
qc.min_profiles = 50;   % minimum number of 'good' profiles for a time series 
qc.max_time_gap = 5;   % max gap (days) between 'good' profiles before rejection
prms.qc = qc;
clear qc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REFERENCE PROFILE SETTINGS
refprof.inner_window = 2;
refprof.outer_window = 12; %%% number of days away from profile of interest
refprof.no_profiles = 10;
prms.refprof = refprof;
clear refprof

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ISOPYCNAL SEPARATION ANOMALY IQR CHECK SETTINGS
isa_iqr_check.min_pres = 50;
isa_iqr_check.max_pres = 600;
isa_iqr_check.min_thickness = 10;
isa_iqr_check.min_density_levels = 2;
isa_iqr_check.magnitude = 1;
prms.isa_iqr_check = isa_iqr_check;
clear isa_iqr_check

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ISOPYCNAL SEPARATION ANOMALY GAUSSIAN FIT CHECK SETTINGS (anticyclones)
isa_gauss.amplitude = 1;
isa_gauss.height = 10;
isa_gauss.min_pres = 50;
isa_gauss.max_pres = 600;
isa_gauss.R2 = 0;
isa_gauss.nrmse = 1;
isa_gauss.min_MLD = 100;
prms.isa_gauss = isa_gauss;
clear isa_gauss

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPICE IQR CHECK SETTINGS
spice_iqr_check.min_pres = 50;
spice_iqr_check.max_pres = 600;
spice_iqr_check.min_thickness = 10;
spice_iqr_check.min_density_levels = 2;
%spice_iqr_check.no_profiles = 15;
spice_iqr_check.magnitude = 0.02;
prms.spice_iqr_check = spice_iqr_check;
clear spice_iqr_check

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPICE ANOMALY GAUSSIAN FIT CHECK SETTINGS 
spice_gauss.amplitude = 0.02;
spice_gauss.height = 10;
spice_gauss.min_pres = 50;
spice_gauss.max_pres = 600;
spice_gauss.R2 = 0;
spice_gauss.nrmse = 1;
spice_gauss.min_MLD = 50;
prms.spice_gauss = spice_gauss;
clear spice_gauss

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DYNAMIC HEIGHT ANOMALY CHECK SETTINGS 
dha_check.amplitude = 0.001;
dha_check.height = 10;
dha_check.min_pres = 50;
dha_check.max_pres = 600;
dha_check.R2 = 0;
dha_check.nrmse = 1;
dha_check.min_MLD = 100;
prms.dha_gauss = dha_check;
clear dha_check

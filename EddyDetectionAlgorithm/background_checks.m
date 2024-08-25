function ts_data = background_checks(ts_data, test_prof)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Isopycnal Stability Check (anticyclones AND cyclones) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-------------------------')
disp('Isopycnal Stability Check')
disp('-------------------------')

for tag_no = test_prof

    %%% Creating array to hold isopycnal separation variance data
    ts_data(tag_no).isopycnal_stability = NaN(size(ts_data(tag_no).cast));

    %%% Isopycnal stability check
    ts_data(tag_no).isopycnal_stability = isopycnal_stability_check(ts_data, tag_no);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bathymetric Stability Check (anticyclones AND cyclones) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('---------------------------');
disp('Bathymetric Stability Check');
disp('---------------------------');

for tag_no = test_prof

    %%% Creating array to hold bathymetric variance data
    ts_data(tag_no).bathymetric_var = NaN(size(ts_data(tag_no).cast));

    %%% Bathymetric stability check
    [ts_data(tag_no).bathymetric_var, ts_data(tag_no).shelf_break_ratio, ts_data(tag_no).jump] = bathymetry_check(ts_data, tag_no);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mixed Layer Depth Check %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-----------------------')
disp('Mixed Layer Depth Check')
disp('-----------------------')

for tag_no = test_prof

    %%% Creating array to hold MLD data
    ts_data(tag_no).MLD = NaN(size(ts_data(tag_no).cast));

    %%% MLD check
    ts_data(tag_no).MLD = MLD_check(ts_data, tag_no);

end


%%%%%%%%%%%%%%%%%%%%%%
%%% Spice STD Test %%%
%%%%%%%%%%%%%%%%%%%%%%

disp('--------------')
disp('Spice STD Test')
disp('--------------')

for tag_no = test_prof

    %%% Creating array to hold spice std data
    ts_data(tag_no).spice_std = NaN(size(ts_data(tag_no).cast));

    %%% Spice STD check
    ts_data(tag_no).spice_std = spice_std_check(ts_data, tag_no);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Checking Max Depth %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('------------------')
disp('Checking Max Depth')
disp('------------------')

for tag_no = test_prof

    for i = 1:length(ts_data(tag_no).cast)

        tmp = max(ts_data(tag_no).ps.pres(~isnan(ts_data(tag_no).ps.temp(:,i)),i));
        if isempty(tmp)
            tmp = 0;
        end
        ts_data(tag_no).max_pres(i) = tmp;
    end
end

end
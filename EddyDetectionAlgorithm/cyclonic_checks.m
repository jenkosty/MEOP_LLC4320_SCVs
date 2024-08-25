function ts_data = cyclonic_checks(ts_data,test_prof, prms)

    %%%%%%%%%%%%%%%%%%%%%%%
    %%% Spice IQR Check %%%
    %%%%%%%%%%%%%%%%%%%%%%%

    disp('---------------');
    disp('Spice IQR Check');
    disp('---------------');

    for tag_no = test_prof

        %%% IQR check for cyclone detection
        [ts_data(tag_no).cyclones, ts_data(tag_no).rejected, ...
            ts_data(tag_no).reason] = spice_iqr_check(ts_data, tag_no, prms);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Spice Gaussian Check %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('--------------------')
    disp('Spice Gaussian Check')
    disp('--------------------')

    for tag_no = test_prof

        disp('Tag ' + string(tag_no))

        %%% Spice gaussian check
        [ts_data(tag_no).cyclones, ts_data(tag_no).rejected, ...
            ts_data(tag_no).reason, ts_data(tag_no).spice_gauss] = spice_gaussian_check(ts_data, tag_no, prms);

    end
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Dynamic Height Anomaly Gaussian Check %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('---------')
    disp('DHA Check')
    disp('---------')

    for tag_no = test_prof

         disp('Tag ' + string(tag_no))

        %%% DHA check
        [ts_data(tag_no).cyclones, ts_data(tag_no).rejected, ...
            ts_data(tag_no).reason, ts_data(tag_no).dha_gauss_cyclones] = dha_check_cyclones(ts_data, tag_no, prms);

    end

end
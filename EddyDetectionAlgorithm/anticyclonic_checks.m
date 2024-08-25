function ts_data = anticyclonic_checks(ts_data, test_prof, prms)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Isopycnal Separation Anomaly IQR Check %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('-------------');
    disp('ISA IQR Check');
    disp('-------------');

    for tag_no = test_prof

        %%% IQR check for anticyclone detection
        [ts_data(tag_no).anticyclones, ts_data(tag_no).rejected, ...
            ts_data(tag_no).reason] = isa_iqr_check(ts_data, tag_no, prms);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Isopycnal Separation Anomaly Gaussian Check %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('------------------')
    disp('ISA Gaussian Check')
    disp('------------------')

    for tag_no = test_prof

         disp('Tag ' + string(tag_no))

        %%% Isopycnal separation anomaly gaussian check
        [ts_data(tag_no).anticyclones, ts_data(tag_no).rejected,...
            ts_data(tag_no).reason, ts_data(tag_no).isa_gauss] = isa_gaussian_check(ts_data, tag_no, prms);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Dynamic Height Anomaly Check %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('---------')
    disp('DHA Check')
    disp('---------')

    for tag_no = test_prof

         disp('Tag ' + string(tag_no))

        %%% DHA check
        [ts_data(tag_no).anticyclones, ts_data(tag_no).rejected,...
            ts_data(tag_no).reason, ts_data(tag_no).dha_gauss_anticyclones] = dha_check_anticyclones(ts_data, tag_no, prms);

    end

    

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Spice Anomaly Fit %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%

    disp('---------')
    disp('Spice Fit')
    disp('---------')

    for tag_no = test_prof

         disp('Tag ' + string(tag_no))

        clear gauss_all
        u = 0;
        for i = ts_data(tag_no).anticyclones.dha
            clear gauss
            gauss = gaussian_fit_spice(ts_data, tag_no, i, prms);
            if gauss.rejected == 0
                gauss_all(i) = gauss;
                u = u + 1;
            end
        end
        if isempty(ts_data(tag_no).anticyclones.dha) | (u == 0)
            gauss_all = struct([]);
        end
        ts_data(tag_no).spice_gauss_anticyclones = gauss_all;
    end

end
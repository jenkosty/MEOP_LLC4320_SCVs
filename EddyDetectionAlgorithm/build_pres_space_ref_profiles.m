function [ref_dyn_height_anom, ref_N2] = build_pres_space_ref_profiles(ts_data, tag_no)

    %%% Creating reference profiles for each of the time series profiles
    ref_dyn_height_anom = NaN(size(ts_data(tag_no).ps.dyn_height_anom));
    ref_N2 = NaN(size(ts_data(tag_no).ps.N2));

    for i = 1:length(ts_data(tag_no).cast)
    
        %%% Extracting the profiles to build the reference profile
        tmp_dyn_height_anom_ps = ts_data(tag_no).ps.dyn_height_anom(:,[ts_data(tag_no).ref_ind{1,i} ts_data(tag_no).ref_ind{2,i}]);
        tmp_N2_ps = ts_data(tag_no).ps.N2(:,[ts_data(tag_no).ref_ind{1,i} ts_data(tag_no).ref_ind{2,i}]);

        %%% Calculating a climatological value for each pressure level if at least
        %%% 75% of the potential reference profiles have data at said
        %%% level
        for j = 1:size(tmp_dyn_height_anom_ps, 1)
            tmp_dyn_height_anom_ps_level = tmp_dyn_height_anom_ps(j,~isnan(tmp_dyn_height_anom_ps(j,:)));
            if length(tmp_dyn_height_anom_ps_level) > 0.75*size(tmp_dyn_height_anom_ps,2)
                ref_dyn_height_anom(j,i) = mean(tmp_dyn_height_anom_ps_level);  
            end
        end

        for j = 1:size(tmp_N2_ps, 1)
            tmp_N2_ps_level = tmp_N2_ps(j,~isnan(tmp_N2_ps(j,:)));
            if length(tmp_N2_ps_level) > 0.75*size(tmp_N2_ps,2)
                ref_N2(j,i) = median(tmp_N2_ps_level);
            end
        end
        
    end

end


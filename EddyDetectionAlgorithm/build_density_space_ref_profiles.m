function [ref_salt, ref_temp, ref_N2, ref_spice, ref_isopycnal_separation] = build_density_space_ref_profiles(ts_data,tag_no)
 
    %%% Creating reference profiles for each of the time series profiles
    ref_salt = NaN(size(ts_data(tag_no).ds.salt));
    ref_temp = NaN(size(ts_data(tag_no).ds.temp));
    ref_spice = NaN(size(ts_data(tag_no).ds.spice));
    ref_N2 = NaN(size(ts_data(tag_no).ds.N2));
    ref_isopycnal_separation = NaN(size(ts_data(tag_no).ds.isopycnal_separation));

    for i = 1:length(ts_data(tag_no).cast)
    
        %%% Extracting the profiles to build the reference profile
        tmp_salt_ds = ts_data(tag_no).ds.salt(:,[ts_data(tag_no).ref_ind{1,i} ts_data(tag_no).ref_ind{2,i}]);
        tmp_temp_ds = ts_data(tag_no).ds.temp(:,[ts_data(tag_no).ref_ind{1,i} ts_data(tag_no).ref_ind{2,i}]);
        tmp_N2_ds = ts_data(tag_no).ds.N2(:,[ts_data(tag_no).ref_ind{1,i} ts_data(tag_no).ref_ind{2,i}]);
        tmp_spice_ds = ts_data(tag_no).ds.spice(:,[ts_data(tag_no).ref_ind{1,i} ts_data(tag_no).ref_ind{2,i}]);
        tmp_isopycnal_separation_ds = ts_data(tag_no).ds.isopycnal_separation(:,[ts_data(tag_no).ref_ind{1,i} ts_data(tag_no).ref_ind{2,i}]);

        %%% Calculating a climatological value for each density level if at least
        %%% 75% of the potential reference profiles have data at said
        %%% level
        for j = 1:size(tmp_salt_ds, 1)
            tmp_salt_ds_level = tmp_salt_ds(j,~isnan(tmp_salt_ds(j,:)));
            tmp_temp_ds_level = tmp_temp_ds(j,~isnan(tmp_temp_ds(j,:)));
            tmp_spice_ds_level = tmp_spice_ds(j,~isnan(tmp_spice_ds(j,:)));
            if length(tmp_salt_ds_level) > 0.75*size(tmp_salt_ds,2)
                ref_salt(j,i) = median(tmp_salt_ds_level);
                ref_temp(j,i) = median(tmp_temp_ds_level);
                ref_spice(j,i) = median(tmp_spice_ds_level); 
            end
            
            tmp_N2_ds_level = tmp_N2_ds(j,~isnan(tmp_N2_ds(j,:)));
            if length(tmp_N2_ds_level) > 0.75*size(tmp_N2_ds,2)
                ref_N2(j,i) = median(tmp_N2_ds_level);  
            end
            
            tmp_isopycnal_separation_ds_level = tmp_isopycnal_separation_ds(j,~isnan(tmp_isopycnal_separation_ds(j,:)));
            if length(tmp_isopycnal_separation_ds_level) > 0.75*size(tmp_isopycnal_separation_ds,2)
                ref_isopycnal_separation(j,i) = mean(tmp_isopycnal_separation_ds_level);  
            end
        end
        
    end

end


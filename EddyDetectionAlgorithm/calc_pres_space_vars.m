function tmp_pres_space = calc_pres_space_vars(ts_data, tag_no, depth_grid, model)
    
    %%% Creating pressure space structure
    tmp_pres_space.pres = depth_grid .* ones(size(ts_data(tag_no).salt));
    
    %%% Assigning temperature and salinity data to new pressure space
    %%% structure
    tmp_pres_space.salt = ts_data(tag_no).salt;
    tmp_pres_space.temp = ts_data(tag_no).temp;
    if model == 1
        tmp_pres_space.vort = ts_data(tag_no).vort;
    end
    
    %%% Calculating absolute salinity and conservative temperature
    tmp_pres_space.salt_absolute = gsw_SA_from_SP(tmp_pres_space.salt, depth_grid, ts_data(tag_no).lon, ts_data(tag_no).lat);
    tmp_pres_space.temp_conservative = gsw_CT_from_t(tmp_pres_space.salt_absolute, tmp_pres_space.temp, tmp_pres_space.pres);
    
    %%% Calculating potential density anomaly
    tmp_pres_space.sigma0 = gsw_sigma0(tmp_pres_space.salt_absolute, tmp_pres_space.temp_conservative);
    
    %%% Calculating N^2
    [N2, mid_pres] = gsw_Nsquared(tmp_pres_space.salt_absolute, tmp_pres_space.temp_conservative, tmp_pres_space.pres, ts_data(tag_no).lat .* ones(size(tmp_pres_space.salt)));
    for i = 1:length(ts_data(tag_no).cast)
        tmp_pres_space.N2(:,i) = interp1(mid_pres(:,i), N2(:,i), tmp_pres_space.pres(:,i));
    end
    
    %%% Calculating dynamic height anomaly
    tmp_pres_space.dyn_height_anom = gsw_geo_strf_dyn_height(tmp_pres_space.salt_absolute, tmp_pres_space.temp_conservative, tmp_pres_space.pres, 0);
    
    %%% Calculating spice
    tmp_pres_space.spice = gsw_spiciness0(tmp_pres_space.salt_absolute, tmp_pres_space.temp_conservative);
end
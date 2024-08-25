function tmp_sigma0_space = interp_to_sigma0_space(ts_data, tag_no, density_grid)

    %%% Creating matrices for density-interpolated data
    interp_salt = NaN(length(density_grid), length(ts_data(tag_no).cast));
    interp_temp = NaN(length(density_grid), length(ts_data(tag_no).cast));
    interp_N2 = NaN(length(density_grid), length(ts_data(tag_no).cast));
    interp_spice = NaN(length(density_grid), length(ts_data(tag_no).cast));
    interp_pres = NaN(length(density_grid), length(ts_data(tag_no).cast));
    
    for i = 1:length(ts_data(tag_no).cast)
        
        %%% Removing NaNs
        tmp_sigma0 = ts_data(tag_no).ps.sigma0(~isnan(ts_data(tag_no).ps.sigma0(:,i)),i);
        tmp_salt = ts_data(tag_no).ps.salt(~isnan(ts_data(tag_no).ps.salt(:,i)) & ~isnan(ts_data(tag_no).ps.sigma0(:,i)),i);
        tmp_temp = ts_data(tag_no).ps.temp(~isnan(ts_data(tag_no).ps.temp(:,i)) & ~isnan(ts_data(tag_no).ps.sigma0(:,i)),i);
        tmp_N2 = ts_data(tag_no).ps.N2(~isnan(ts_data(tag_no).ps.N2(:,i)) & ~isnan(ts_data(tag_no).ps.sigma0(:,i)),i);
        tmp_sigma0_N2 = ts_data(tag_no).ps.sigma0(~isnan(ts_data(tag_no).ps.sigma0(:,i)) & ~isnan(ts_data(tag_no).ps.N2(:,i)),i);
        tmp_spice = ts_data(tag_no).ps.spice(~isnan(ts_data(tag_no).ps.spice(:,i)) & ~isnan(ts_data(tag_no).ps.sigma0(:,i)),i);
        tmp_pres = ts_data(tag_no).ps.pres(~isnan(ts_data(tag_no).ps.sigma0(:,i)) & ~isnan(ts_data(tag_no).ps.sigma0(:,i)),i);
        
        if isempty(tmp_salt) || length(tmp_salt) < 2
            interp_salt(:,i) = NaN(length(density_grid), 1);
            interp_temp(:,i) = NaN(length(density_grid), 1);
            interp_spice(:,i) = NaN(length(density_grid), 1);
            interp_pres(:,i) = NaN(length(density_grid), 1);
        else
            %%% Interpolating data
            interp_salt(:,i) = interp1(tmp_sigma0, tmp_salt, density_grid);
            interp_temp(:,i) = interp1(tmp_sigma0, tmp_temp, density_grid);
            interp_spice(:,i) = interp1(tmp_sigma0, tmp_spice, density_grid);
            interp_pres(:,i) = interp1(tmp_sigma0, tmp_pres, density_grid);
        end

        if length(tmp_sigma0_N2) < 2 || length(tmp_N2) < 2
            interp_N2(:,i) = NaN(length(density_grid), 1);
        else
            interp_N2(:,i) = interp1(tmp_sigma0_N2, tmp_N2, density_grid);
        end
    end
    
    %%% Saving density-interpolated data
    tmp_sigma0_space.salt = interp_salt;
    tmp_sigma0_space.temp = interp_temp;
    tmp_sigma0_space.N2 = interp_N2;
    tmp_sigma0_space.spice = interp_spice;
    tmp_sigma0_space.pres = interp_pres;
    tmp_sigma0_space.sigma0 = density_grid .* ones(size(interp_salt));
end
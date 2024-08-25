function anoms = calc_pres_space_anomalies(ts_data,tag_no)
 
    %%% Pre-allocating storage for anomaly matrices
    anoms.dyn_height_anom = NaN(size(ts_data(tag_no).ps.dyn_height_anom));

    %%% Calculating anomalies
    for i = 1:length(ts_data(tag_no).cast)
        anoms.dyn_height_anom(:,i) = ts_data(tag_no).ps.dyn_height_anom(:,i) - ts_data(tag_no).ps.ref_dyn_height_anom(:,i);
    end

end


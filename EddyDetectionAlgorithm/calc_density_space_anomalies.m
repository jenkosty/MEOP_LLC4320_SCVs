function anoms = calc_density_space_anomalies(ts_data,tag_no)

    %%% Pre-allocating storage for anomaly matrices
    anoms.salt = NaN(size(ts_data(tag_no).ds.salt));
    anoms.temp = NaN(size(ts_data(tag_no).ds.temp));
    anoms.N2 = NaN(size(ts_data(tag_no).ds.N2));
    anoms.spice = NaN(size(ts_data(tag_no).ds.spice));
    anoms.isopycnal_separation = NaN(size(ts_data(tag_no).ds.isopycnal_separation));
    anoms.isopycnal_separation_normalized = NaN(size(ts_data(tag_no).ds.isopycnal_separation));

    for i = 1:length(ts_data(tag_no).cast)

        %%% Calculating anomalies
        anoms.salt(:,i) = ts_data(tag_no).ds.salt(:,i) - ts_data(tag_no).ds.ref_salt(:,i);
        anoms.temp(:,i) = ts_data(tag_no).ds.temp(:,i) - ts_data(tag_no).ds.ref_temp(:,i);
        anoms.N2(:,i) = ts_data(tag_no).ds.N2(:,i) - ts_data(tag_no).ds.ref_N2(:,i);
        anoms.spice(:,i) = ts_data(tag_no).ds.spice(:,i) - ts_data(tag_no).ds.ref_spice(:,i);
        anoms.isopycnal_separation(:,i) = ts_data(tag_no).ds.isopycnal_separation(:,i) - ts_data(tag_no).ds.ref_isopycnal_separation(:,i);

        %%% Normalizing isopycnal separation anomaly
        anoms.isopycnal_separation_normalized(:,i) = anoms.isopycnal_separation(:,i) ./ ts_data(tag_no).ds.ref_isopycnal_separation(:,i);
    end

end


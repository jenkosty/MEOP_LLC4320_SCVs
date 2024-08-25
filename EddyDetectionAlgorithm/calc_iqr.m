function iqr_data = calc_iqr(ts_data,tag_no)

    %%% Pre-allocating storage for iqr data
    iqr_data.spice_anom_iqr = NaN(size(ts_data(tag_no).ds.spice));
    iqr_data.N2_anom_iqr = NaN(size(ts_data(tag_no).ds.N2));
    iqr_data.isopycnal_separation_anom_iqr = NaN(size(ts_data(tag_no).ds.isopycnal_separation));
    iqr_data.isopycnal_separation_anom_normalized_iqr = NaN(size(ts_data(tag_no).ds.isopycnal_separation));

    %%% Pre-allocating storage for lower threshold data
    iqr_data.spice_anom_lim_lo = NaN(size(ts_data(tag_no).ds.spice));
    iqr_data.N2_anom_lim_lo = NaN(size(ts_data(tag_no).ds.N2));
    iqr_data.isopycnal_separation_anom_normalized_lim_lo = NaN(size(ts_data(tag_no).ds.isopycnal_separation));

    %%% Pre-allocating storage for upper threshold data
    iqr_data.spice_anom_lim_hi = NaN(size(ts_data(tag_no).ds.spice));
    iqr_data.N2_anom_lim_hi = NaN(size(ts_data(tag_no).ds.N2));
    iqr_data.isopycnal_separation_anom_normalized_lim_hi = NaN(size(ts_data(tag_no).ds.isopycnal_separation));

    for i = 1:length(ts_data(tag_no).cast)

        %%% Extracting the assigned anomaly profiles
        tmp_spice_anom_ds = ts_data(tag_no).ds.anoms.spice(:,[ts_data(tag_no).ref_ind{1,i} ts_data(tag_no).ref_ind{2,i}]);
        tmp_spice_anom_ds(sum(double(~isnan(tmp_spice_anom_ds)), 2) < 0.75*size(tmp_spice_anom_ds,2),:) = NaN;
        tmp_N2_anom_ds = ts_data(tag_no).ds.anoms.N2(:,[ts_data(tag_no).ref_ind{1,i} ts_data(tag_no).ref_ind{2,i}]);
        tmp_N2_anom_ds(sum(double(~isnan(tmp_N2_anom_ds)), 2) < 0.75*size(tmp_N2_anom_ds,2),:) = NaN;
        tmp_isopycnal_separation_anom_ds = ts_data(tag_no).ds.anoms.isopycnal_separation(:,[ts_data(tag_no).ref_ind{1,i} ts_data(tag_no).ref_ind{2,i}]);
        tmp_isopycnal_separation_anom_ds(sum(double(~isnan(tmp_isopycnal_separation_anom_ds)), 2) < 0.75*size(tmp_isopycnal_separation_anom_ds,2),:) = NaN;
        tmp_isopycnal_separation_anom_normalized_ds = ts_data(tag_no).ds.anoms.isopycnal_separation_normalized(:,[ts_data(tag_no).ref_ind{1,i} ts_data(tag_no).ref_ind{2,i}]);
        tmp_isopycnal_separation_anom_normalized_ds(sum(double(~isnan(tmp_isopycnal_separation_anom_normalized_ds)), 2) < 0.75*size(tmp_isopycnal_separation_anom_normalized_ds,2),:) = NaN;
              
        %%% Calculating IQR
        iqr_data.spice_anom_iqr(:,i) = prctile(tmp_spice_anom_ds, 75, 2) - prctile(tmp_spice_anom_ds, 25, 2);
        iqr_data.N2_anom_iqr(:,i) = prctile(tmp_N2_anom_ds, 75, 2) - prctile(tmp_N2_anom_ds, 25, 2);
        iqr_data.isopycnal_separation_anom_iqr(:,i) = prctile(tmp_isopycnal_separation_anom_ds, 75, 2) - prctile(tmp_isopycnal_separation_anom_ds, 25, 2);
        iqr_data.isopycnal_separation_anom_normalized_iqr(:,i) = prctile(tmp_isopycnal_separation_anom_normalized_ds, 75, 2) - prctile(tmp_isopycnal_separation_anom_normalized_ds, 25, 2);
        
        %%% Threshold value (to be multiplied by the IQR)
        iqr_threshold = 1.5;

        %%% Calculating upper thresholds
        iqr_data.spice_anom_lim_hi(:,i) = prctile(tmp_spice_anom_ds, 75, 2) + iqr_threshold*iqr_data.spice_anom_iqr(:,i);
        iqr_data.N2_anom_lim_hi(:,i) = prctile(tmp_N2_anom_ds, 75, 2) + iqr_threshold*iqr_data.N2_anom_iqr(:,i);
        iqr_data.isopycnal_separation_anom_lim_hi(:,i) = prctile(tmp_isopycnal_separation_anom_ds, 75, 2) + iqr_threshold*iqr_data.isopycnal_separation_anom_iqr(:,i);
        iqr_data.isopycnal_separation_anom_normalized_lim_hi(:,i) = prctile(tmp_isopycnal_separation_anom_normalized_ds, 75, 2) + iqr_threshold*iqr_data.isopycnal_separation_anom_normalized_iqr(:,i);
        
        %%% Calculating lower thresholds
        iqr_data.spice_anom_lim_lo(:,i) = prctile(tmp_spice_anom_ds, 25, 2) - iqr_threshold*iqr_data.spice_anom_iqr(:,i);
        iqr_data.N2_anom_lim_lo(:,i) = prctile(tmp_N2_anom_ds, 25, 2) - iqr_threshold*iqr_data.N2_anom_iqr(:,i);
        iqr_data.isopycnal_separation_anom_lim_lo(:,i) = prctile(tmp_isopycnal_separation_anom_ds, 25, 2) - iqr_threshold*iqr_data.isopycnal_separation_anom_iqr(:,i);
        iqr_data.isopycnal_separation_anom_normalized_lim_lo(:,i) = prctile(tmp_isopycnal_separation_anom_normalized_ds, 25, 2) - iqr_threshold*iqr_data.isopycnal_separation_anom_normalized_iqr(:,i);
    end

end


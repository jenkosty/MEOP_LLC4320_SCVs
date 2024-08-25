function isopycnal_stability = isopycnal_stability_check(ts_data,tag_no)

    %%% Creating array to hold profiles that pass isopycnal stability check
    isopycnal_stability = NaN(size(ts_data(tag_no).cast));

    %%% Calculating isopycnal stability for all profiles
    for i = 1:length(ts_data(tag_no).cast)
        ind = [ts_data(tag_no).ref_ind{1,i} ts_data(tag_no).ref_ind{2,i}]; % indices of surrounding profiles
        ref_isopycnal_separation = ts_data(tag_no).ds.ref_isopycnal_separation(:,i);
        pres = ts_data(tag_no).ds.pres(:,ind);
        pres(pres < 100) = NaN;
        mean_pres = mean(pres, 2, 'omitnan');
        rmse = std(pres, 0, 2, 'omitnan');
        %rmse = sqrt(sum((pres - mean_pres).^2, 2, 'omitnan') / length(ind));
        idx = sum(~isnan(pres), 2);
        rmse(idx < (0.75*length(ind))) = NaN;
        rmse = rmse ./ ref_isopycnal_separation;
        mean_rmse = mean(rmse, 'omitnan');
        isopycnal_stability(i) = mean_rmse; 
    end

end


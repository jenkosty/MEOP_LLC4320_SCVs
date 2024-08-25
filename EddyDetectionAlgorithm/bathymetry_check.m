function [bathymetric_var, shelf_break_ratio, jump] = bathymetry_check(ts_data,tag_no)

    %%% Creating array to hold profiles that pass bathymetry check
    bathymetric_var = NaN(size(ts_data(tag_no).cast));
    shelf_break_ratio = NaN(size(ts_data(tag_no).cast));

    %%% Calculating bathymetric variance for all profiles
    for i = 1:length(ts_data(tag_no).cast)
        ind = [ts_data(tag_no).ref_ind{1,i} ts_data(tag_no).ref_ind{2,i}];
        bathymetry = ts_data(tag_no).bathymetry(:, min(ind):max(ind));
        %bathymetric_var_prof = var(bathymetry, 'omitnan') / (mean(bathymetry, 'omitnan')^2);
        bathymetric_var_prof = std(bathymetry, 'omitnan') / mean(bathymetry, 'omitnan');
        bathymetric_var(i) = bathymetric_var_prof;
    end

    %%% Checking for large bathmetric jumps in individual profile
    for i = 1:length(ts_data(tag_no).cast)
        ind = [ts_data(tag_no).ref_ind{1,i} ts_data(tag_no).ref_ind{2,i}];
        bathymetry = ts_data(tag_no).bathymetry(:,[i ind]);
        bathymetric_zscore = zscore(bathymetry);
        if abs(bathymetric_zscore(1)) > 2
            jump(i) = 1;
        else
            jump(i) = 0;
        end
    end

    %%% Checking for percentage of profiles above / below 1000 dbar
    for i = 1:length(ts_data(tag_no).cast)
        ind = [ts_data(tag_no).ref_ind{1,i} ts_data(tag_no).ref_ind{2,i}];
        bathymetry = abs(ts_data(tag_no).bathymetry(:,[i ind]));
        shallow = sum(abs(bathymetry) <= 1000);
        deep = sum(abs(bathymetry) > 1000);
        if (shallow ~= 0) && (deep ~= 0)
            shelf_break_ratio(i) = min([shallow deep]) / length(bathymetry);
        else
            shelf_break_ratio(i) = 0;
        end
    end

end


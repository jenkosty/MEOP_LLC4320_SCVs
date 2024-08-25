function isopycnal_separation = calc_isopycnal_separation(ts_data, tag_no)

    %%% Creating matrix for isopycnal separation data
    isopycnal_separation = NaN(size(ts_data(tag_no).ds.sigma0));
    
    %%% Calculating isopycnal separation
    for i = 1:length(ts_data(tag_no).cast)
        for j = 2:length(ts_data(tag_no).ds.sigma0(:,i))-1
            isopycnal_separation(j,i) = ts_data(tag_no).ds.pres(j+1, i) - ts_data(tag_no).ds.pres(j-1, i);
        end
    end

    %%% Note: Normalization will take place in the anomaly calculation
    %%% step. 

end


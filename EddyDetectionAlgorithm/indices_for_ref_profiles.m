function  mean_ind = indices_for_ref_profiles(ts_data,tag_no, refprof)

    mean_ind = cell(2,length(ts_data(tag_no).cast));
    
    %%% Finding indices to use for background profile calculation
    for i = 1:length(ts_data(tag_no).cast)
        
        mean_ind(1,i) = {(i-refprof.outer_window):(i-refprof.inner_window)};
        mean_ind{1,i}(mean_ind{1,i} < 1) = []; %%% Making sure indices remain within ts boundaries
        mean_ind(2,i) = {(i+refprof.inner_window):(i+refprof.outer_window)};
        mean_ind{2,i}(mean_ind{2,i} > length(ts_data(tag_no).cast)) = []; %%% Making sure indices remain within ts boundaries
        
        %%% Creating special reference profiles for the start and end of
        %%% the time series
        if (length(mean_ind{1,i}) < refprof.outer_window - refprof.inner_window+1)
            mean_ind{2,i} = [mean_ind{2,i} max(mean_ind{2,i})+1:max(mean_ind{2,i})+(refprof.outer_window)-(refprof.inner_window-1)-length(mean_ind{1,i})];
        end
        
        if (length(mean_ind{2,i}) < refprof.outer_window - refprof.inner_window+1)
            mean_ind{1,i} = sort([mean_ind{1,i} min(mean_ind{1,i})-1:-1:min(mean_ind{1,i})-(refprof.outer_window)+(refprof.inner_window-1)+length(mean_ind{2,i})]);
        end
        
    end
    
end


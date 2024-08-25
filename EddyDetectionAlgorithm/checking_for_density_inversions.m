function [sigma0, max_sigma0_inversion, rejected, reason] = checking_for_density_inversions(ts_data, tag_no)

    %%% Extracting arrays that hold profile rejection and justification data
    rejected = ts_data(tag_no).rejected;
    reason = ts_data(tag_no).reason;

    %%% Creating arrays to hold inversion and corrected density data
    max_sigma0_inversion = NaN(size(ts_data(tag_no).cast));
    sigma0 = NaN(size(ts_data(tag_no).ps.sigma0));
    
    for i =  1:length(ts_data(tag_no).cast)
        
        %%% Checking to see if sigma0 profile is ascending with depth
        if issorted(ts_data(tag_no).ps.sigma0(~isnan(ts_data(tag_no).ps.sigma0(:,i)),i), 'ascend') == 1

            max_sigma0_inversion(i) = 0;
            sigma0(:,i) = ts_data(tag_no).ps.sigma0(:,i);

        else

            %%% If not, sorting the profile
            sigma0_orig = ts_data(tag_no).ps.sigma0(:,i);
            idx = ~isnan(sigma0_orig);
            sigma0_sort_initial = sort(sigma0_orig, 'ascend');
            sigma0_sort_initial(isnan(sigma0_sort_initial)) = [];
            sigma0_sort = sigma0_orig;
            sigma0_sort(idx) = sigma0_sort_initial;

            %%% Checking max difference between sorted and original density
            %%% profile. If the max difference it too large, the profile
            %%% will be rejected. If the max difference is small, the
            %%% profile will be replaced with the sorted version. 
         
            max_sigma0_inversion(:,i) = max(abs(sigma0_orig-sigma0_sort));
       
            if max(abs(sigma0_orig-sigma0_sort)) > 0.1 % kg/m^3
                sigma0(:,i) = NaN;
                rejected.anticyclones(i) = 1;
                reason.anticyclones(i) = "Density Inversion";
                rejected.cyclones(i) = 1;
                reason.cyclones(i) = "Density Inversion";
            else
                sigma0(:,i) = sigma0_sort;
            end
           
        end

    end

end


function  avg_MLD = MLD_check(ts_data,tag_no)

%%% Creating array to hold profiles that pass MLD check
avg_MLD = NaN(size(ts_data(tag_no).cast));

for i = 1:length(ts_data(tag_no).cast)

    %%% Indices of surrounding profiles
    ind = [ts_data(tag_no).ref_ind{1,i} ts_data(tag_no).ref_ind{2,i}];

    %%% Creating array to hold MLD of surrounding profiles
    MLD = NaN(size(ind));

    for ii = 1:length(ind)

        %%% Extracting potential density profile
        sigma0 = ts_data(tag_no).ps.sigma0(:,ind(ii));
        sigma0_nonan = sigma0(~isnan(sigma0));
        pres = ts_data(tag_no).ps.pres(:,ind(ii));

        %%% Extracting surface density
        if ~isempty(sigma0_nonan)
            surface_density = sigma0_nonan(1);
        else
            MLD(ii) = NaN;
            continue
        end
        
        %%%
        density_diff = sigma0 - surface_density;

        %%% MLD Calculation
        MLD_ind = find(density_diff > 0.03, 1);
        if isempty(MLD_ind)
            MLD(ii) = max(pres(~isnan(sigma0)));
        else
            MLD(ii) = pres(MLD_ind);
        end

    end

    %%% Averaging the surrounding MLD
    avg_MLD(i) = mean(MLD, 'omitnan');

end

end
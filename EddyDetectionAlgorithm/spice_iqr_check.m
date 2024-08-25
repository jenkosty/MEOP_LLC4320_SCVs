function [cyclones, rejected, reason] = spice_iqr_check(ts_data, tag_no, prms)

%%% Extracting arrays that hold profile classification, rejection, and justification data
cyclones = ts_data(tag_no).cyclones;
rejected = ts_data(tag_no).rejected;
reason = ts_data(tag_no).reason;

%%% Creating array to hold profiles that pass iqr check
cyclones.spice_iqr = [];

u = 1;

for i = 1:length(ts_data(tag_no).cast)


    %%% Checking spice
    if mean(ts_data(tag_no).ds.anoms.spice(:,i), 'omitnan') > 0
        spice_check = (ts_data(tag_no).ds.anoms.spice(:,i) > ts_data(tag_no).ds.iqrs.spice_anom_lim_hi(:,i));
    else
        spice_check = (ts_data(tag_no).ds.anoms.spice(:,i) < ts_data(tag_no).ds.iqrs.spice_anom_lim_lo(:,i));
    end

    if sum(double(spice_check)) <= prms.spice_iqr_check.min_density_levels
        rejected.cyclones(i) = 1;
        reason.cyclones(i) = "No IQR Anomaly";
    else

        %%% Extracting number of continuous anomalies
        y = diff(find([0 double(spice_check(ts_data(tag_no).ds.pres(:,i) >= prms.spice_iqr_check.min_pres))' 0]==0))-1;
        y(y==0) = [];

        %%% Indices at which isopyncal separation check passes
        A = find(spice_check .* ts_data(tag_no).ds.pres(:,i) >= prms.spice_iqr_check.min_pres);

        %%% Creating cells for continuous blocks of indices
        B = cell(1, length(y));
        for j = 1:length(y)
            B{1,j} = NaN(y(j),1);
        end

        %%% Filling in cells with indices
        q = 1;
        while q <= length(B)
            for j = 1:length(y)
                for k = 1:length(B{1,j})
                    B{1,j}(k,1) = A(q);
                    q = q+1;
                end
            end
        end

        if isempty(B)
            rejected.cyclones(i) = 1;
            reason.cyclones(i) = "No IQR Anomaly";
            continue
        else
            for j = 1:length(B)

                %%% Extracting pressure levels
                pres_levels = ts_data(tag_no).ds.pres(B{1,j},i);

                %%% Extracting associated isopycnal separation anomaly values
                spice_anom = ts_data(tag_no).ds.anoms.spice(B{1,j},i);

                %%% Rejecting anomaly for various reasons
                %%% 1. Anomaly width
                if length(B{1,j}) < prms.spice_iqr_check.min_density_levels
                    rejected.cyclones(i) = 1;
                    reason.cyclones(i) = "Too Few Isopycnals (IQR)";
                    continue
                elseif (max(pres_levels) - min(pres_levels(pres_levels > prms.spice_iqr_check.min_pres))) < prms.spice_iqr_check.min_thickness
                    rejected.cyclones(i) = 1;
                    reason.cyclones(i) = "Too Short (IQR)";
                    continue
                %%% 2. Location of anomaly in water column
                elseif max(pres_levels) < prms.spice_iqr_check.min_pres 
                    rejected.cyclones(i) = 1;
                    reason.cyclones(i) = "Too Shallow (IQR)";
                    continue
                elseif min(pres_levels) > prms.spice_iqr_check.max_pres 
                    rejected.cyclones(i) = 1;
                    reason.cyclones(i) = "Too Shallow (IQR)";
                    continue
                % elseif pres_levels(abs(spice_anom) == max(abs(spice_anom))) < prms.spice_iqr_check.min_pres
                %     rejected.cyclones(i) = 1;
                %     reason.cyclones(i) = "Peak Too Shallow (IQR)";
                %     continue
                % elseif pres_levels(abs(spice_anom) == max(abs(spice_anom))) > prms.spice_iqr_check.max_pres
                %     rejected.cyclones(i) = 1;
                %     reason.cyclones(i) = "Peak Too Deep (IQR)";
                %     continue
                %%% 3. Strength of anomaly
                elseif max(abs(spice_anom)) < prms.spice_iqr_check.magnitude
                    rejected.cyclones(i) = 1;
                    reason.cyclones(i) = "Spice Anomaly Too Small (IQR)";
                    continue

                %%% Profile passes all checks so flag it for further
                %%% testing
                else
                    rejected.cyclones(i) = 0;
                    reason.cyclones(i) = strings;
                    cyclones.spice_iqr(u) = i;
                    u = u + 1;
                    break
                end

            end

        end

    end

end

end


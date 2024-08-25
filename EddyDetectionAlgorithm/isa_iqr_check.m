function [anticyclones, rejected, reason] = isa_iqr_check(ts_data,tag_no, prms)

%%% Extracting arrays that hold profile classification, rejection, and justification data
anticyclones = ts_data(tag_no).anticyclones;
rejected = ts_data(tag_no).rejected;
reason = ts_data(tag_no).reason;

%%% Creating array to hold profiles that pass iqr check
anticyclones.isa_iqr = [];

u = 1;

for i = 1:length(ts_data(tag_no).cast)

    %%% Checking isopycnal separation
    isopycnal_separation_check = (ts_data(tag_no).ds.anoms.isopycnal_separation_normalized(:,i) > ts_data(tag_no).ds.iqrs.isopycnal_separation_anom_normalized_lim_hi(:,i));

    if sum(double(isopycnal_separation_check)) <= prms.isa_iqr_check.min_density_levels
        rejected.anticyclones(i) = 1;
        reason.anticyclones(i) = "No IQR Anomaly";
    else

        %%% Extracting number of continuous anomalies
        y = diff(find([0 double(isopycnal_separation_check(ts_data(tag_no).ds.pres(:,i) >= prms.isa_iqr_check.min_pres))' 0]==0))-1;
        y(y==0) = [];

        %%% Indices at which isopyncal separation check passes
        A = find(isopycnal_separation_check .* ts_data(tag_no).ds.pres(:,i) >= prms.isa_iqr_check.min_pres);

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
            rejected.anticyclones(i) = 1;
            reason.anticyclones(i) = "No IQR Anomaly";
            continue
        else

            for j = 1:length(B)

                %%% Extracting pressure levels
                pres_levels = ts_data(tag_no).ds.pres(B{1,j},i);

                %%% Extracting associated isopycnal separation anomaly values
                isopycnal_sep_anom = ts_data(tag_no).ds.anoms.isopycnal_separation_normalized(B{1,j},i);

                %%% Rejecting anomaly for various reasons
                %%% 1. Anomaly width
                if length(B{1,j}) < prms.isa_iqr_check.min_density_levels
                    rejected.anticyclones(i) = 1;
                    reason.anticyclones(i) = "Too Few Isopycnals (IQR)";
                    continue
                elseif (max(pres_levels) - min(pres_levels(pres_levels > prms.isa_iqr_check.min_pres))) < prms.isa_iqr_check.min_thickness
                    rejected.anticyclones(i) = 1;
                    reason.anticyclones(i) = "Too Short (IQR)";
                    continue

                %%% 2. Location of anomaly in water column
                elseif max(pres_levels) < prms.isa_iqr_check.min_pres
                    rejected.anticyclones(i) = 1;
                    reason.anticyclones(i) = "Too Shallow (IQR)";
                    continue
                elseif min(pres_levels) > prms.isa_iqr_check.max_pres  
                    rejected.anticyclones(i) = 1;
                    reason.anticyclones(i) = "Too Deep (IQR)";
                    continue

                %%% 3. Strength of anomaly
                elseif max(ts_data(tag_no).ds.anoms.isopycnal_separation_normalized(B{1,j},i)) < prms.isa_iqr_check.magnitude
                    rejected.anticyclones(i) = 1;
                    reason.anticyclones(i) = "Isopycnal Separation Anomaly Too Small (IQR)";
                    continue

                %%% Profile passes all checks so flag it for further
                %%% testing
                else

                    rejected.anticyclones(i) = 0;
                    reason.anticyclones(i) = strings;
                    anticyclones.isa_iqr(u) = i;
                    u = u + 1;

                    break
                end

            end

        end

    end

end

end


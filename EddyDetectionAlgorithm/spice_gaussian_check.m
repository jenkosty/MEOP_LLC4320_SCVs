function [cyclones, rejected, reason, gauss_all] = spice_gaussian_check(ts_data,tag_no, prms)

    %%% Extracting arrays hold profile classification, rejection, and justification data
    cyclones = ts_data(tag_no).cyclones;
    rejected = ts_data(tag_no).rejected;
    reason = ts_data(tag_no).reason;

    %%% Creating array to hold profiles that pass spice gaussian check
    cyclones.spice_gaussian = [];
    u = 1;

    for i = cyclones.spice_iqr

        %%% Gaussian Fit 
        gauss = gaussian_fit_spice(ts_data, tag_no, i, prms);

        %%% Rejecting if gaussian fit is poor or amplitude is small
        if gauss.rejected == 1
            rejected.cyclones(i) = 1;
            reason.cyclones(i) = 'Poor Spice Fit';
            continue
        elseif (abs(gauss.A) < prms.spice_gauss.amplitude)
            rejected.cyclones(i) = 1;
            reason.cyclones(i) = 'Small Spice Amplitude';
            continue
        elseif (gauss.H < prms.spice_gauss.height)
            rejected.cyclones(i) = 1;
            reason.cyclones(i) = 'Width of Spice Anomaly';
            continue
        elseif (gauss.P < prms.spice_gauss.min_pres) || (gauss.P > prms.spice_gauss.max_pres)
            rejected.cyclones(i) = 1;
            reason.cyclones(i) = 'Depth of Spice Anomaly';
            continue
        % elseif isnan(mean(ts_data(tag_no).ds.anoms.spice((ts_data(tag_no).ds.pres(:,i) > gauss.Phih),i), 'omitnan'))
        %     rejected.cyclones(i) = 1;
        %     reason.cyclones(i) = 'Spice Gauss Not Closed';
        %     continue
        % elseif isnan(mean(ts_data(tag_no).ds.anoms.spice((ts_data(tag_no).ds.pres(:,i) < gauss.Plow),i), 'omitnan'))
        %     rejected.cyclones(i) = 1;
        %     reason.cyclones(i) = 'Spice Gauss Not Closed';
        %     continue
        else
            cyclones.spice_gaussian(u) = i;
            u = u + 1;
            gauss_all(i) = gauss;
        end
        
    end

    if u == 1
        gauss_all = struct([]);
    end

end
function [cyclones, rejected, reason, gauss_all] = dha_check_cyclones(ts_data, tag_no, prms)

    %%% Extracting arrays hold profile classification, rejection, and justification data
    cyclones = ts_data(tag_no).cyclones;
    rejected = ts_data(tag_no).rejected;
    reason = ts_data(tag_no).reason;

    %%% Creating array to hold profiles that pass dha check
    cyclones.dha = [];
    u = 1;

    for i = cyclones.spice_gaussian

        %%% Dynamic Height Anomaly Gaussian Fit
        gauss = gaussian_fit_dha_modes(ts_data, tag_no, i, 1, prms);

        %%% Rejecting if gaussian fit is poor or amplitude is small
        if gauss.rejected == 1
            rejected.cyclones(i) = 1;
            reason.cyclones(i) = 'Poor DHA Fit';
            continue
        elseif (gauss.A > 0)
            rejected.cyclones(i) = 1;
            reason.cyclones(i) = 'Positive DHA Amplitude';
            continue
        elseif (abs(gauss.A) < prms.dha_gauss.amplitude)
            rejected.cyclones(i) = 1;
            reason.cyclones(i) = 'Small DHA Amplitude';
            continue
        elseif(gauss.H < prms.dha_gauss.height)
            rejected.cyclones(i) = 1;
            reason.cyclones(i) = 'Width of DHA Amplitude';
            continue
        elseif (gauss.P > prms.dha_gauss.max_pres) || (gauss.P < prms.dha_gauss.min_pres)
            rejected.cyclones(i) = 1;
            reason.cyclones(i) = 'Depth of DHA Anomaly';
            continue
        % elseif isnan(mean(gauss.dataX(gauss.dataY > gauss.Phih), 'omitnan'))
        %     rejected.cyclones(i) = 1;
        %     reason.cyclones(i) = 'DHA Gauss Not Closed';
        %     continue
        % elseif isnan(mean(gauss.dataX(gauss.dataY < gauss.Plow), 'omitnan'))
        %     rejected.cyclones(i) = 1;
        %     reason.cyclones(i) = 'DHA Gauss Not Closed';
        %     continue
        else
            cyclones.dha(u) = i;
            u = u + 1;
            gauss_all(i) = gauss;
        end

    end

    if u == 1
        gauss_all = struct([]);
    end

end
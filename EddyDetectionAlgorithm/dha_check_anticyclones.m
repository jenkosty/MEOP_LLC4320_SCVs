function [anticyclones, rejected, reason, gauss_all] = dha_check_anticyclones(ts_data, tag_no, prms)

    %%% Extracting arrays hold profile classification, rejection, and justification data
    anticyclones = ts_data(tag_no).anticyclones;
    rejected = ts_data(tag_no).rejected;
    reason = ts_data(tag_no).reason;

    %%% Creating array to hold profiles that pass dha check
    anticyclones.dha = [];
    u = 1;

    for i = anticyclones.isa_gaussian

        %%% Dynamic Height Anomaly Gaussian Fit
        gauss = gaussian_fit_dha_modes(ts_data, tag_no, i, 0, prms);

        %%% Rejecting if gaussian fit is poor or amplitude is small
        if gauss.rejected == 1
            rejected.anticyclones(i) = 1;
            reason.anticyclones(i) = 'Poor DHA Fit';
            continue
        elseif (gauss.A < prms.dha_gauss.amplitude)
            rejected.anticyclones(i) = 1;
            reason.anticyclones(i) = 'Small DHA Amplitude';
            continue
        elseif(gauss.H < prms.dha_gauss.height)
            rejected.anticyclones(i) = 1;
            reason.anticyclones(i) = 'Width of DHA Amplitude';
            continue
        elseif (gauss.P > prms.dha_gauss.max_pres) || (gauss.P < prms.dha_gauss.min_pres)
            rejected.anticyclones(i) = 1;
            reason.anticyclones(i) = 'Depth of DHA Anomaly';
            continue
        % elseif isnan(mean(gauss.dataX(gauss.dataY > gauss.Phih), 'omitnan'))
        %     rejected.anticyclones(i) = 1;
        %     reason.anticyclones(i) = 'DHA Gauss Not Closed';
        %     continue
        % elseif isnan(mean(gauss.dataX(gauss.dataY < gauss.Plow), 'omitnan'))
        %     rejected.anticyclones(i) = 1;
        %     reason.anticyclones(i) = 'DHA Gauss Not Closed';
        %     continue
        % elseif isnan(mean(gauss.prof.dyn_height_anom_BC1(gauss.prof.dyn_height_pres_BC1 > gauss.Phih), 'omitnan'))
        %     rejected.anticyclones(i) = 1;
        %     reason.anticyclones(i) = 'DHA Gauss Not Closed';
        %     continue
        % elseif isnan(mean(gauss.prof.dyn_height_anom_BC1(gauss.prof.dyn_height_pres_BC1 < gauss.Plow), 'omitnan'))
        %     rejected.anticyclones(i) = 1;
        %     reason.anticyclones(i) = 'DHA Gauss Not Closed';
        %     continue
        else
            anticyclones.dha(u) = i;
            u = u + 1;
            gauss_all(i) = gauss;
        end

    end
    
    if u == 1
        gauss_all = struct([]);
    end

end
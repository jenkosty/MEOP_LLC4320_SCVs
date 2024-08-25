function [anticyclones, rejected, reason, gauss_all] = isa_gaussian_check(ts_data,tag_no, prms)

    %%% Extracting arrays hold profile classification, rejection, and justification data
    anticyclones = ts_data(tag_no).anticyclones;
    rejected = ts_data(tag_no).rejected;
    reason = ts_data(tag_no).reason;

    %%% Creating array to hold profiles that pass isa gaussian check
    anticyclones.isa_gaussian = [];
    u = 1;

    for i = anticyclones.isa_iqr

        %%% Gaussian Fit 
        gauss = gaussian_fit_isa(ts_data, tag_no, i, prms);
        isa_profile = ts_data(tag_no).ds.anoms.isopycnal_separation_normalized(:,i);
        pres_profile = ts_data(tag_no).ds.pres(:,i);

        %%% Rejecting if gaussian fit is poor or amplitude is small
        if gauss.rejected == 1
            rejected.anticyclones(i) = 1;
            reason.anticyclones(i) = 'Poor ISA Fit (Gauss)';
            continue
        elseif (gauss.A < prms.isa_gauss.amplitude)
            rejected.anticyclones(i) = 1;
            reason.anticyclones(i) = 'Small ISA Amplitude (Gauss)';
            continue
        elseif (gauss.H < prms.isa_gauss.height)
            rejected.anticyclones(i) = 1;
            reason.anticyclones(i) = 'Width of ISA (Gauss)';
            continue
        elseif (gauss.P < prms.isa_gauss.min_pres) || (gauss.P > prms.isa_gauss.max_pres)
            rejected.anticyclones(i) = 1;
            reason.anticyclones(i) = 'Depth of ISA (Gauss)';
            continue
        % elseif isnan(mean(isa_profile(pres_profile > gauss.Phih), 'omitnan'))
        %     rejected.anticyclones(i) = 1;
        %     reason.anticyclones(i) = 'ISA Gauss Not Closed';
        %     continue
        % elseif isnan(mean(isa_profile(pres_profile < gauss.Plow), 'omitnan'))
        %     rejected.anticyclones(i) = 1;
        %     reason.anticyclones(i) = 'ISA Gauss Not Closed';
        %     continue
        else
            anticyclones.isa_gaussian(u) = i;
            u = u + 1;
            gauss_all(i) = gauss;
        end
        
    end

    if u == 1
        gauss_all = struct([]);
    end

end
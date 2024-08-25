function results = gaussian_fit_dha_modes(ts_data,tag_no,i,cyclonic, prms)

prof = dha_modes_decomp(ts_data, tag_no, i);
if prof.rejected == 1
    results.rejected = 1;
    results.reason = "DHA Modes";
    return
end
prof_final.dyn_pres = prof.dyn_pres;
prof_final.dyn_height_anom = prof.dyn_height_anom;
prof_final.dyn_height_anom_BC1 = prof.dyn_height_anom_BC1;
prof_final.dyn_height_pres_BC1 = prof.dyn_height_pres_BC1;

%%% Extracting pressure and anomaly profiles for fitting
pres = prof.dyn_height_pres_BC1;
anom = prof.dyn_height_anom_BC1;
pres_orig = pres;
anom_orig = anom;

%%% Calculating MLD
sigma0 = ts_data(tag_no).ps.sigma0(:,i);
pres_ps = ts_data(tag_no).ps.pres(:,i);
sigma0_nonan = sigma0(~isnan(sigma0));
surface_density = sigma0_nonan(1);
density_diff = sigma0 - surface_density;
MLD_ind = find(density_diff > 0.03, 1);
if isempty(MLD_ind)
    results.rejected = 1;
    results.reason = "No Mixed Layer";
    return
end
MLD = pres_ps(MLD_ind);
if MLD < prms.dha_gauss.min_MLD
    MLD = prms.dha_gauss.min_MLD;
end

%%% Extracting pressure and anomaly profiles below MLD && removing
%%% NaNs
idx = (pres > MLD);
pres = pres(idx);
anom = anom(idx);

idx = (pres < 0.9*max(pres));
pres = pres(idx);
anom = anom(idx);

idy = ~isnan(pres) & ~isnan(anom);
pres = pres(idy);
anom = anom(idy);

if isempty(pres) || isempty(anom)
    results.rejected = 1;
    results.reason = 'No Data Below Mixed Layer';
    return
end

%%% flipping anomaly direction if cyclonic
if cyclonic == 1
    anom = -anom;
end

% Grab amplitude and depth of max anomaly
spike.A  = max(anom);
% if spike.A < 0
%     results.rejected = 1;
%     results.reason = "Negative Peak";
%     return
% end
spike.P = pres(find(anom == spike.A));

% Get range of allowable parameters
prng = [-0.3:0.05:0.3];  % allow pressure peak to vary between +- 20% of height
arng  = [0.7:0.05:1.3];  % allow amplitude range of +- 20% of anomaly peak
hrng  = [10:10:300];     % allow height to vary between 50 and 500 meters

% Goodness-of-fit thresholds
R2_thresh    = prms.dha_gauss.R2;
NRMSE_thresh = prms.dha_gauss.nrmse;

% Set up matrix for least-squared error calculation
lse = nan(length(prng),length(arng),length(hrng));
R2    = nan(size(lse));
NRMSE = nan(size(lse));

% Go through all possible combinations
hcnt = 0; % reset h counter
for h = hrng
    hcnt = hcnt + 1; % increase 'h' counter
    acnt = 0;        % reset 'a' counter
    for a = arng
        acnt = acnt + 1; % increase 'a' counter
        pcnt = 0;        % reset 'p' counter
        for p = prng
            pcnt = pcnt + 1; % increase 'p'

            % Center Gaussian model around spike.P + p*h
            zo = [];
            zo = double(pres - [spike.P + p*(4)*sqrt(h^2/2)]);
            sa = double(anom);

            % Reduce to where data exists
            dat      = sa + zo;
            sa       = sa(~isnan(dat));
            zo       = zo(~isnan(dat));

            % Generate gaussian model using updated amplitude, center, and height
            gauss = (spike.A*a)*exp((-(zo.^2))/(h.^2));

            % Grab results (conducting lse calc over entire
            % profile!!)
            dataX = anom;
            dataY = pres;
            dataY  = round(dataY, 6);
            zp     = [zo + spike.P + p*(4)*sqrt(h^2/2)];
            modelX = gauss;
            modelY = zp;
            modelY = round(modelY, 6);

            % Check that depths of model and data intersect
            if length(dataX) < length(modelX) | length(modelX) < length(dataX)
                [c,ia,ib] = intersect(dataY,modelY);
                if ~isempty(c)
                    ind     = find((min(c) <= dataY) & (dataY <= max(c)));
                    dataX   = dataX(ia);   dataY = dataY(ia);
                    ind     = find((min(c) <= modelY) & (modelY <= max(c)));
                    modelX  = modelX(ib); modelY = modelY(ib);
                else
                    R2(pcnt, acnt, hcnt) = NaN;
                    lse(pcnt, acnt, hcnt) = NaN;
                    continue
                end
            end

            % Calculate R^2
            R2(pcnt,acnt,hcnt) = corr2(dataX,modelX).^2;

            % Calculate NRMSE (must be < 0.5);
            RMSE                  = sqrt(sum((dataX - modelX).^2)/length(dataX));
            NRMSE(pcnt,acnt,hcnt) = RMSE/(max(dataX) - min(dataX));

            % Save least-squared error results (ignore if bad R2 value (i.e. < 0.5))
            if (R2(pcnt,acnt,hcnt) < R2_thresh) || isnan(R2(pcnt,acnt,hcnt))
                lse(pcnt,acnt,hcnt) = NaN;
            elseif NRMSE(pcnt,acnt,hcnt) > NRMSE_thresh
                lse(pcnt,acnt,hcnt) = NaN;
            else
                lse(pcnt,acnt,hcnt) = sum([dataX-modelX].^2);
            end

        end
    end
end

% Find best zo,A,H combo according to lse
[minlse,idxlse] = min(lse(:));
if isnan(minlse)
    results.rejected = 1;
    results.reason = "Bad R^2 Value (Gaussian)";
    return
end
[a,b,c] = ind2sub(size(lse),idxlse);

% Update parameters
results.A    = spike.A*arng(b);
results.H    = hrng(c);
results.P    = spike.P + prng(a)*(4)*sqrt(results.H^2/2);
results.Hcore = sqrt(2)*results.H;
results.Plow = results.P - 2*sqrt((results.H^2)/2);
results.Phih = results.P + 2*sqrt((results.H^2)/2);
results.Plow = round(results.Plow/10)*10;
results.Phih = round(results.Phih/10)*10;
results.minlse = minlse;
results.R2 = R2(idxlse);
results.nrmse = NRMSE(idxlse);
results.rejected = 0;
results.dataX = anom;
results.dataY = pres;
results.dataX_orig = anom_orig;
results.dataY_orig = pres_orig;
results.prof = prof_final;
results.MLD = MLD;

% Update zo,zp,gauss for final model
zo    = [0:1:800] - [results.P];
zp    = zo + results.P;
gauss = results.A*exp((-(zo.^2))/(results.H.^2));

% Save final model
results.X = gauss;
results.Y = zp;

% Finally, fix orientation of model if cyclonic
if cyclonic == 1
    results.X = -results.X;
    results.A = -results.A;
    results.dataX = -results.dataX;
    anom = -anom;
end

% figure()
% plot(ts_data(tag_no).ps.anoms.dyn_height_anom(:,i),ts_data(tag_no).ps.pres(:,i),'k','linewidth',3)
% hold on; grid on; set(gca,'YDir','Reverse')
% yline(MLD, 'b--', 'LineWidth', 3)
% plot(anom, pres, 'k--', 'linewidth', 3)
% plot(anom_orig, pres_orig, 'k--', 'LineWidth', 3);
% plot(results.X,results.Y,'Color','r','LineWidth',3,'LineStyle','-.')
% xline(0, 'k', 'LineWidth', 3)
% xlabel('Dynamic Height Anomaly (m^2/s^2)', 'FontSize', 12)
% ylabel('Pressure (dbar)', 'FontSize', 12);
% %title('Dynamic Height Anomaly Gaussian Fit', 'FontSize', 15, 'FontWeight', 'bold')
% title('Tag #: ' + string(tag_no) + ', Cast #: ' + string(i) + ', R^2: ' + string(results.R2) + ', nrmlse: ' + string(results.nrmse) + ', ' + string(median(results.dataX)))
% ylim([0 700])
% grid on

end

function results = gaussian_fit_spice(ts_data,tag_no,i,prms)

clear results

%%% Extracting pressure and anomaly profiles for fitting
pres = ts_data(tag_no).ds.pres(:,i);
anom = ts_data(tag_no).ds.anoms.spice(:,i);
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
if MLD < prms.spice_gauss.min_MLD
    MLD = prms.spice_gauss.min_MLD;
end

%%% Extracting pressure and anomaly profiles below MLD && removing
%%% NaNs
idx = pres > MLD;
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

%%% flipping anomaly direction if minty
if mean(anom, 'omitnan') < 0
    minty = 1;
    anom = -anom;
else
    minty = 0;
end

% Get range of allowable parameters
prng = [-0.3:0.05:0.3];  % allow pressure peak to vary between +- 20% of height
arng  = [0.7:0.05:1.3];  % allow amplitude range of +- 20% of anomaly peak
hrng  = [10:10:300];     % allow height to vary between 50 and 500 meters

% Goodness-of-fit thrsholds
R2_thresh    = prms.spice_gauss.R2;
NRMSE_thresh = prms.spice_gauss.nrmse;

% Grab amplitude and depth of max anomaly
ind = islocalmax(anom, 'MaxNumExtrema', 5);
possible_maxes = anom(ind);
possible_pres = pres(ind);

%%% Gaussian fit for each possible max
nrmse_vals = NaN(length(possible_maxes),1);
R2_vals = NaN(length(possible_maxes),1);
for run_no = 1:length(possible_maxes)

    clear results

    spike.A  = possible_maxes(run_no);
    spike.P = possible_pres(run_no);

    % % Grab amplitude and depth of max anomaly
    % spike.A  = max(anom);
    % spike.P = pres(find(anom == spike.A));
    % spike.P = spike.P(1);

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

                % % Get gaussian limits for testing
                % pl = [spike.P + p*(4)*sqrt(h^2/2)] - 2*sqrt((h^2)/2); pl  = round(pl/10)*10;
                % ph = [spike.P + p*(4)*sqrt(h^2/2)] + 2*sqrt((h^2)/2); ph  = round(ph/10)*10;
                %
                % % Grab results
                % zp     = [zo + spike.P + p*(4)*sqrt(h^2/2)];
                % dataX  = anom(pl <= pres & pres <= ph);
                % dataY  = pres(pl <= pres & pres <= ph);
                % dataY  = round(dataY, 6);
                % modelX = gauss(pl <= zp & zp <= ph);
                % modelY = zp(pl <= zp & zp <= ph);
                % modelY = round(modelY, 6);

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
    [minlse,idxlse] = min(lse(:), [], 'omitnan');
    if isnan(minlse)
        results.rejected = 1;
        results.reason = "Bad R^2 Value (Gaussian)";
        continue
    end
    [a,b,c] = ind2sub(size(lse),idxlse);

    % Update parameters
    results.A    = spike.A*arng(b); % amplitude of gaussian
    results.H    = hrng(c); % height of gaussian
    results.P    = spike.P + prng(a)*(4)*sqrt(results.H^2/2); % pressure peak of gaussian
    results.Hcore = sqrt(2)*results.H;
    results.Plow = results.P - 2*sqrt((results.H^2)/2);
    results.Phih = results.P + 2*sqrt((results.H^2)/2);
    results.Plow = round(results.Plow/10)*10;
    results.Phih = round(results.Phih/10)*10;
    results.minlse = minlse;
    results.R2 = R2(a, b, c);
    results.nrmse = NRMSE(a, b, c);
    results.rejected = 0;
    results.dataX = dataX;
    results.dataY = dataY;
    results.dataX_orig = anom_orig;
    results.dataY_orig = pres_orig;
    results.MLD = MLD;

    results_all{run_no} = results;
    R2_vals(run_no) = results.R2;
    nrmse_vals(run_no) = results.nrmse;

end

if isempty(possible_maxes) || sum(isnan(nrmse_vals)) == length(nrmse_vals)
    results.rejected = 1;
    results.reason = "Bad R^2 Value (All Gaussians)";
    return
end

%%% Saving best fit (based on nrmse value)
[~,ind] = min(nrmse_vals);
results = results_all{ind};

% Update zo,zp,gauss for final model
zo    = double([0:1:800] - [results.P]);
zp    = zo + results.P;
gauss = results.A*exp((-(zo.^2))/(results.H.^2));

% Save final model
results.X = gauss;
results.Y = zp;

% Finally, fix orientation of model if minty
if minty == 1
    results.X = -results.X;
    results.A = -results.A;
    results.dataX = -results.dataX;
    possible_maxes = -possible_maxes;
    anom = -anom;
end

% figure()
% plot(anom,pres,'k','linewidth',2)
% hold on; grid on; set(gca,'YDir','Reverse')
% plot(results.X,results.Y,'Color','r','LineWidth',3,'LineStyle','-.')
% scatter(possible_maxes, possible_pres, 70, 'o', 'filled', 'MarkerFaceColor', 'k');
% % yline(results.Plow, '--g');
% % yline(results.Phih, '--r');
% xlabel('Spice Anomaly (kg/m^3)', 'FontSize', 12)
% ylabel('Pressure (dbar)', 'FontSize', 12);
% %title('Tag #: ' + string(tag_no) + ', Cast #: ' + string(i) + ', R^2: ' + string(results.R2) + ', lse: ' + string(results.minlse))
% title('Spice Anomaly Gaussian Fit', 'FontSize', 15, 'FontWeight', 'bold')
% grid on

end
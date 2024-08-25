function results = gaussian_fit_isa_ds(ts_data,tag_no,i,prms)

clear results

%%% Extracting pressure and anomaly profiles for fitting
pres = ts_data(tag_no).ds.pres(:,i);
sigma0 = ts_data(tag_no).ds.sigma0(:,i);
anom = ts_data(tag_no).ds.anoms.isopycnal_separation_normalized(:,i);

%%% Calculating MLD
sigma0_ps = ts_data(tag_no).ps.sigma0(:,i);
pres_ps = ts_data(tag_no).ps.pres(:,i);
sigma0_nonan = sigma0_ps(~isnan(sigma0_ps));
surface_density = sigma0_nonan(1);
density_diff = sigma0_ps - surface_density;
MLD_ind = find(density_diff > 0.03, 1);
if isempty(MLD_ind)
    results.rejected = 1;
    results.reason = "No Mixed Layer";
    return
end
MLD = pres_ps(MLD_ind);
if MLD < 100
    MLD = 100;
end

%%% Extracting pressure and anomaly profiles below MLD && removing
%%% NaNs
idx = (pres > MLD) & (pres < 0.9*max(pres));
sigma0 = sigma0(idx);
anom = anom(idx);

idy = ~isnan(sigma0) & ~isnan(anom);
sigma0 = sigma0(idy);
anom = anom(idy);

if isempty(sigma0) || isempty(anom)
    results.rejected = 1;
    results.reason = 'No Data Below Mixed Layer';
    return
end

% Get range of allowable parameters
prng = [-0.2:0.05:0.2];  % allow pressure peak to vary between +- 20% of height
arng  = [0.8:0.05:1.2];  % allow amplitude range of +- 20% of anomaly peak
hrng  = [0.001:0.001:0.3];     % allow height to vary between 50 and 500 meters

% Goodness-of-fit thrsholds
R2_thresh    = 0; %prms.isa_gauss.R2;
NRMSE_thresh = 1; %prms.isa_gauss.nrmse;

% Grab amplitude and depth of max anomaly
ind = islocalmax(anom, 'MaxNumExtrema', 5) & anom > 0;
possible_maxes = anom(ind);
possible_sigma0 = sigma0(ind);

for run_no = 1:length(possible_maxes)

spike.A  = possible_maxes(run_no);
spike.P = possible_sigma0(run_no);

% Set up matrices for least-squared error calculation
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
            zo = double(sigma0 - [spike.P + p*(4)*sqrt(h^2/2)]);
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
            dataY = sigma0;
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

            % Calculate NRMSE;
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
    return
end
[a,b,c] = ind2sub(size(lse),idxlse);

% Update parameters
results.A    = spike.A*arng(b); % amplitude of gaussian
results.H    = hrng(c); % height of gaussian
results.P    = spike.P + prng(a)*(4)*sqrt(results.H^2/2); % pressure peak of gaussian
results.a = a;
results.b = b;
results.c = c;
results.Plow = results.P - 1*sqrt((results.H^2)/2);
results.Phih = results.P + 1*sqrt((results.H^2)/2);
results.Plow = round(results.Plow/10)*10;
results.Phih = round(results.Phih/10)*10;
results.minlse = minlse;
results.R2 = R2(a, b, c);
results.nrmse = NRMSE(a, b, c);
results.rejected = 0;
results.lse = lse;
results.dataX = dataX;
results.dataY = dataY;
results.MLD = MLD;

results_all{run_no} = results;
R2_vals(run_no) = results.R2;
nrmse_vals(run_no) = results.nrmse;

end

if isempty(possible_maxes)
    results.rejected = 1;
    return
end

[~,ind] = min(nrmse_vals);
results = results_all{ind};

% Update zo,zp,gauss for final model
zo    = double(pres - [results.P]);
zp    = zo + results.P;
gauss = results.A*exp((-(zo.^2))/(results.H.^2));

% Interpolating to smooth pressure grid
zo    = [27:0.001:28] - [results.P];
zp_interp    = zo + results.P;
gauss_interp = results.A*exp((-(zo.^2))/(results.H.^2));

% Save final model
results.X = gauss_interp;
results.Y = zp_interp;
% 
% figure()
% plot(ts_data(tag_no).ds.anoms.isopycnal_separation_normalized(:,i),ts_data(tag_no).ds.sigma0(:,i),'k','linewidth',2);
% hold on; grid on; set(gca,'YDir','Reverse');
% %plot(anom, pres, 'r')
% scatter(possible_maxes, possible_sigma0, 70, 'o', 'filled');
% %yline(results.MLD, 'Color', 'b', 'LineWidth', 3, 'LineStyle', '-.');
% plot(results.X,results.Y,'Color','r','LineWidth',3,'LineStyle','-.');
% %plot(anom, pres, '--c')
% %plot(spike.A, spike.P, 'o', 'MarkerFaceColor', 'c')
% xlabel('Normalized Isopycnal Separation Anomaly', 'FontSize', 12);
% ylabel('Pressure (dbar)', 'FontSize', 12);
% title('Tag #: ' + string(tag_no) + ', Cast #: ' + string(i) + ', R^2: ' + string(results.R2) + ', nrmlse: ' + string(results.nrmse))
% %title('Isopycnal Separation Anomaly Gaussian Fit', 'FontSize', 15, 'FontWeight', 'bold')
% grid on

end
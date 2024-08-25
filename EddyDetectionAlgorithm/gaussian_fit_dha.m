function results = gaussian_fit_dha(ts_data,tag_no,i,cyclonic)

        %%% Extracting pressure and anomaly profiles for fitting
        pres = ts_data(tag_no).ps.pres(:,i);
        anom = ts_data(tag_no).ps.anoms.dyn_height_anom(:,i);
        idx = pres > 75;
        pres = pres(idx);
        anom = anom(idx);
        
        %%% flipping anomaly direction if cyclonic
        if cyclonic == 1
            anom = -anom;
        end
    
        % Grab amplitude and depth of max anomaly
        spike.A  = max(anom);
        spike.P = pres(find(anom == spike.A));
        
        % Get range of allowable parameters
        prng = [-0.2:0.05:0.2];  % allow pressure peak to vary between +- 20% of height
        arng  = [0.8:0.05:1.2];  % allow amplitude range of +- 20% of anomaly peak
        hrng  = [50:10:500];     % allow height to vary between 50 and 500 meters

        % Goodness-of-fit thresholds
        R2_thresh    = 0.5;
        NRMSE_thresh = 0.5;
        
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
        results.Plow = spike.P - 2*sqrt((results.H^2)/2);
        results.Phih = spike.P + 2*sqrt((results.H^2)/2);
        results.Plow = round(results.Plow/10)*10;
        results.Phih = round(results.Phih/10)*10;
        results.minlse = minlse;
        results.R2 = R2(idxlse);
        results.nrmse = NRMSE(idxlse);
        results.rejected = 0;
        results.dataX = dataX;
        results.modelX = modelX;
        results.dataY = dataY;
        results.modelY = modelY;
        
        % Update zo,zp,gauss for final model
        zo    = double(pres - [results.P]);
        zp    = zo + results.P;
        gauss = results.A*exp((-(zo.^2))/(results.H.^2));
        
        % Save final model
        results.X = gauss;
        results.Y = zp;
        
        %// Finally, fix orientation of model if cyclonic
        if cyclonic == 1
            results.X = -results.X;
            results.A = -results.A;
        end
        
        figure()
        plot(ts_data(tag_no).ps.anoms.dyn_height_anom(:,i),ts_data(tag_no).ps.pres(:,i),'k','linewidth',2)
        hold on; grid on; set(gca,'YDir','Reverse')
        plot(results.X,results.Y,'Color','r','LineWidth',3,'LineStyle','-.')
        xlabel('Dynamic Height Anomaly (m^2/s^2)', 'FontSize', 12)
        ylabel('Pressure (dbar)', 'FontSize', 12);
        title('Dynamic Height Anomaly Gaussian Fit', 'FontSize', 15, 'FontWeight', 'bold')
        ylim([0 700])
        grid on

end


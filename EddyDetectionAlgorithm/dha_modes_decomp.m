function prof = dha_modes_decomp(ts_data, tag_no, i)

        opts1 = optimset('display','off','UseParallel',false);

        %%% Extract vertical velocity and horizontal structure modes of climatology
        BA = ~isnan(ts_data(tag_no).ps.ref_N2(:,i));
        [meop_profile(i).ref.wmodes, meop_profile(i).ref.pmodes, ~, ~] = dynmodes(ts_data(tag_no).ps.ref_N2(BA,i), ts_data(tag_no).ps.pres(BA,i),1);

        %%% Grab pressure levels of mode decomposition, add zero level (surface)
        meop_profile(i).ref.mode_pres = ts_data(tag_no).ps.pres(BA,i);
        if meop_profile(i).ref.mode_pres(1) ~= 0
            meop_profile(i).ref.mode_pres = [0;meop_profile(i).ref.mode_pres];
        end

        % Interpolate 1st baroclinic mode to pressure of SCV cast
        BA2 = ~isnan(ts_data(tag_no).ps.anoms.dyn_height_anom(:,i));
        meop_profile(i).dyn_pres = ts_data(tag_no).ps.pres(BA2,i);
        meop_profile(i).dyn_height_anom = ts_data(tag_no).ps.anoms.dyn_height_anom(BA2,i);
        meop_profile(i).ref.BC1_data = meop_profile(i).ref.pmodes(:,1);
        meop_profile(i).ref.BC1_pres = meop_profile(i).ref.mode_pres;
        meop_profile(i).ref.BC1_data = interp1(meop_profile(i).ref.mode_pres(~isnan(meop_profile(i).ref.BC1_data)),meop_profile(i).ref.BC1_data(~isnan(meop_profile(i).ref.BC1_data)),meop_profile(i).dyn_pres);
        meop_profile(i).ref.BC1_pres = meop_profile(i).dyn_pres;

        % Create function that describes residuals between projected BC1 and dyn_height_anom
        % Exclude data inbetween SCV limits for better fit to first mode
        dat = [];
        dat = [meop_profile(i).ref.BC1_data + meop_profile(i).dyn_height_anom];
        x_o = [];
        x_o = meop_profile(i).ref.BC1_data(~isnan(dat));
        x_p = [];
        x_p = meop_profile(i).ref.BC1_pres(~isnan(dat));
        x_f = [];
        x_f = meop_profile(i).dyn_height_anom(~isnan(dat));

        % % Get limits
        % pl = 100; %ts_data(tag_no).gauss_fit{1,i}.Plow;
        % ph = 300; %ts_data(tag_no).gauss_fit{1,i}.Phih;
        % 
        % %Remove values between upper/lower limits of SCV to avoid bad fit
        % ind = [];
        % ind = find(pl < x_p & x_p < ph);
        % if isempty(ind)
        %     meop_profile(i).rejected = 1;
        %     meop_profile(i).reason = "DHA Preprocessing";
        %     prof = meop_profile(i);
        %     return;
        % end
        % if ind(end) == length(x_p)
        %     ind = ind(1:end-1);
        % end
        % x_o(ind) = [];
        % x_f(ind) = [];
        % x_p(ind) = [];

        % Remove mixed layer depths (Lynne Talley method, first density greater than 0.03 from sfc value
        ind      = [];
        mld_dens = ts_data(tag_no).ps.sigma0(~isnan(ts_data(tag_no).ps.sigma0(:,i)),i);
        if isempty(mld_dens)
            meop_profile(i).rejected = 1;
            meop_profile(i).reason = "DHA Preprocessing";
            prof = meop_profile(i);
            return;
        end
        mld_pres = ts_data(tag_no).ps.pres(~isnan(ts_data(tag_no).ps.sigma0(:,i)),i);
        ind      = find(mld_dens > mld_dens(1)+0.03);
        if isempty(ind)
            meop_profile(i).rejected = 1;
            meop_profile(i).reason = "DHA Preprocessing";
            prof = meop_profile(i);
            return;
        end
        mld_pres = mld_pres(ind(1));
        ind      = find(x_p < mld_pres);
        if length(ind) >= length(x_o)
            meop_profile(i).rejected = 1;
            meop_profile(i).reason = "DHA Preprocessing";
            prof = meop_profile(i);
            return;
        end
        x_o(ind) = [];
        x_f(ind) = [];

        % f simply evaluates a given alpha (modal amplitude) and returns the
        % difference between the input DHanom profile and the projected 1st mode
        % We want to restrict our solutions such that the bottom of the projected
        % profile is equal to the bottom of the DHanom profile
        % SO let alpha2 = DHanom(end) - alpha*BT1(end)
        f = [];
        f = @(alpha) (alpha*x_o - x_f + (x_f(1) - alpha*x_o(1))); %%% CHANGED FROM END TO 1
        x0  = 0.05; % First guess

        % Solve for best modal amplitude
        alpha = [];
        alpha = lsqnonlin(f,x0,[-1],[1],opts1);

        % Redefine x_o and x_f with full profile
        x_o = meop_profile(i).ref.BC1_data(~isnan(dat));
        x_p = meop_profile(i).ref.BC1_pres(~isnan(dat)); %%% NOTE: Changed from meop_profile(i).ref.mode_pres - need to check with Danny
        x_f = meop_profile(i).dyn_height_anom(~isnan(dat));

        % Fix dynamic height anomaly by removing projected 1st mode, add back in barotopic mode
        meop_profile(i).dyn_height_anom_BC1 = [x_f] - [x_o*alpha + (x_f(1) - alpha*x_o(1))]; %%% CHANGED FROM END TO 1
        meop_profile(i).dyn_height_pres_BC1 = meop_profile(i).dyn_pres(~isnan(dat));

        % Save VMD results
        meop_profile(i).ref.VMD.x_f      = x_f;
        meop_profile(i).ref.VMD.x_o      = x_o;
        meop_profile(i).ref.VMD.x_p      = x_p;
        meop_profile(i).ref.VMD.alpha    = alpha;

        % Get mode decomposition results
        BC1 = meop_profile(i).ref.VMD.x_o*meop_profile(i).ref.VMD.alpha;
        BC1 = BC1 - BC1(end); %// Set bottom to zero
        
        %%% Check for 'peak' in dynamic height around core of SCV
        % plidx  = find(abs([pl-meop_profile(i).dyn_height_pres_BC1]) == min(abs([pl-meop_profile(i).dyn_height_pres_BC1])));
        % phidx  = find(abs([ph-meop_profile(i).dyn_height_pres_BC1]) == min(abs([ph-meop_profile(i).dyn_height_pres_BC1])));
        % scvidx = find(pl <= meop_profile(i).dyn_height_pres_BC1 & meop_profile(i).dyn_height_pres_BC1 <= ph);
        % dh_low   = meop_profile(i).dyn_height_anom_BC1(plidx);
        % dh_high  = meop_profile(i).dyn_height_anom_BC1(phidx);
        % dh_peak  = max(meop_profile(i).dyn_height_anom_BC1(scvidx));
        
        %%% Reporting results of test
        % if isempty([dh_low+dh_high+dh_peak])==1 || dh_low >= dh_peak || dh_high >= dh_peak
        %     meop_profile(i).rejected = 1;
        %     meop_profile(i).reason = "Failed DHA Test (Max Not Within SCV Limits)";
        % elseif dh_peak < 0
        %     meop_profile(i).rejected = 1;
        %     meop_profile(i).reason = "Failed DHA Test (Negative Max)";
        % else
        %     meop_profile(i).rejected = 0;
        %     meop_profile(i).reason = strings;
        % end

        meop_profile(i).rejected = 0;
        meop_profile(i).reason = strings;
        
        prof = meop_profile(i);
 
        % % %%% Plot results
        % figure();
        % subplot(121)
        % plot(-meop_profile(i).ref.pmodes(:,1),meop_profile(i).ref.mode_pres,'r','linewidth',2)
        % hold on; grid on; set(gca,'YDir','Reverse')
        % plot(meop_profile(i).ref.pmodes(:,2),meop_profile(i).ref.mode_pres,'b','linewidth',2)
        % plot(meop_profile(i).ref.pmodes(:,3),meop_profile(i).ref.mode_pres,'g','linewidth',2)
        % plot(meop_profile(i).ref.pmodes(:,4),meop_profile(i).ref.mode_pres,'y','linewidth',2)
        % plot(meop_profile(i).ref.pmodes(:,5),meop_profile(i).ref.mode_pres,'color',[0.5 0 0.5],'linewidth',2)
        % title({'\it\bf\fontsize{8}\fontname{Helvetica}Horizontal Velocity','Modes'})
        % set(gca,'XTick',[0])
        % ylabel('Pressure (dbar)')
        % [l,~] = legend('Mode-1','Mode-2','Mode-3','Mode-4','Mode-5','location','southeast');
        % l.Box = 'off';
        % ylim([0 700]);
        % 
        % subplot(122)
        % plot(meop_profile(i).dyn_height_anom,meop_profile(i).dyn_pres,'k','linewidth',2)
        % hold on; grid on; set(gca,'YDir','Reverse')
        % plot(BC1,meop_profile(i).ref.VMD.x_p,':r','linewidth',2)
        % plot(meop_profile(i).dyn_height_anom_BC1,meop_profile(i).dyn_height_pres_BC1,':k','linewidth',2)
        % legend('DH''_{orig}','BC1_{fit}','DH''_{adj}','location','southeast');
        % xlabel('m^2/s^2')
        % title({'\it\bf\fontsize{8}\fontname{Helvetica}Dynamic Height','Anomaly'})
        % set(gca,'YTickLabel',[])
        % ylim([0 700]);
function OW = OWcalculation(dxg, dyg, dxc, dyc, u, v, depths)
    
    dxg = double(dxg);
    dyg = double(dyg);
    dxc = double(dxc);
    dyc = double(dyc);
    u = double(u);
    v = double(v);
    depths = double(depths);
    depths = depths(1:86);

    dudx = NaN(size(dxg));
    dudy = NaN(size(dxg));
    dvdx = NaN(size(dxg));
    dvdy = NaN(size(dxg));

    %%% Getting indices for depth averaging (change as needed)
    ind = (depths < -100) & (depths > -500);
    
    %%% Calculating depth averaged u and v
    u_avg = mean(u(:,:,ind), 3, "omitnan");
    v_avg = mean(v(:,:,ind), 3, "omitnan");

    for i = 2:size(dudx,1)-1
        dudx(i,:) = (u_avg(i+1,:)-u_avg(i,:)) ./ dxg(i,:); % Cell center
        dvdx(i,:) = (v_avg(i,:)-v_avg(i-1,:)) ./ dxc(i,:); % Cell corner
    end

    for j = 2:size(dudx,2)-1
        dudy(:,j) = (u_avg(:,j)-u_avg(:,j-1)) ./ dyc(:,j); % Cell corner
        dvdy(:,j) = (v_avg(:,j+1)-v_avg(:,j)) ./ dyg(:,j); % Cell center
    end

    %%% Calculating normal strain
    s_n_center = dudx - dvdy; % Cell center

    %%% Averaging to get normal strain on the cell corner
    s_n = NaN(size(s_n_center));
    for i = 2:size(s_n_center,1)
        for j = 2:size(s_n_center,2)
            s_n(i,j) = mean([s_n_center(i,j) s_n_center(i-1,j) s_n_center(i-1,j-1) s_n_center(i,j-1)], 'omitnan');
        end
    end
    
    %%% Calculating shear strain
    s_s = dvdx + dudy; % Cell corner
    
    %%% Calculating relative vorticity
    omega = dvdx - dudy; % Cell corner
    
    %%% Calculating Okubo Weiss
    OW = s_n.^2 + s_s.^2 - omega.^2; % Cell corner
    
end
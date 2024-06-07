%%% This script takes the recreates synthetic versions of the MEOP seal
%%% tracks using the LLC4320 snapshots

%%% Load MEOP data
load('qc_ts.mat');

%%% Path's for LLC4320 snapshots
input_path = '/Volumes/Elements/LLCsealdata/Snapshot_';
output_path = '/Volumes/Elements/LLCsealdata/Snapshot_';

%%% LLC4320 snapshot dates
snapshot_dates = {'01-Oct-2011','01-Nov-2011', '01-Dec-2011', '01-Jan-2012', '01-Feb-2012', '01-Mar-2012', '01-Apr-2012', '01-May-2012', ...
    '01-Jun-2012', '01-Jul-2012', '01-Aug-2012', '01-Sep-2012'};

%%% Cycling through snapshots
for i = 1:length(snapshot_dates)
    date = snapshot_dates{i};
    disp(date);
    creatingSyntheticSealTracks(date, qc_ts, input_path, output_path)
end

function creatingSyntheticSealTracks(date, qc_ts, input_path, output_path)

    %%% Loading LLC Snapshot Data
    LLC = cell(4,1);
    sectors = {'LLC_1', 'LLC_2', 'LLC_4', 'LLC_5'};
    for i = 1:4
        load(string(input_path) + string(date) + '/' + string(sectors{i}) + '/lats.mat');
        LLC{i}.lats = lats;
        load(string(input_path) + string(date) + '/' + string(sectors{i}) + '/lons.mat');
        LLC{i}.lons = lons;
        load(string(input_path) + string(date) + '/' + string(sectors{i}) + '/OW.mat');
        LLC{i}.OW = OW;
        load(string(input_path) + string(date) + '/' + string(sectors{i}) + '/vort.mat');
        LLC{i}.vort = vort;
        load(string(input_path) + string(date) + '/' + string(sectors{i}) + '/salt.mat');
        LLC{i}.salt = salt;
        load(string(input_path) + string(date) + '/' + string(sectors{i}) + '/temp.mat');
        LLC{i}.temp = temp;
        load(string(input_path) + string(date) + '/' + string(sectors{i}) + '/depth.mat');
        LLC{i}.depth = depth;
        load(string(input_path) + string(date) + '/' + string(sectors{i}) + '/polygon.mat')
        LLC{i}.polygon = polygon;
        LLC{i}.date = date;
        disp('Snapshot ' + string(i))
        clear lats lons OW vort depth polygon salt temp
    end
    LLC_1 = LLC{1};
    LLC_2 = LLC{2};
    LLC_4 = LLC{3};
    LLC_5 = LLC{4};
    clear LLC i sectors
    
    %%% Seal tags used for interpolation
    test_prof = 1:length(qc_ts);
    
    %%% Pressure grid for interpolation
    depth_grid = 0:800;
    depth_grid = depth_grid';
    
    for tag_no = test_prof

        disp('Tag ' + string(tag_no));
    
        %%% Grabbing data on MEOP seal time series
        LLCsealdata(tag_no).tag = qc_ts(tag_no).tag;
        LLCsealdata(tag_no).cast = qc_ts(tag_no).cast;
        LLCsealdata(tag_no).lat = qc_ts(tag_no).lat;
        LLCsealdata(tag_no).lon = qc_ts(tag_no).lon;
        LLCsealdata(tag_no).time = qc_ts(tag_no).time;
        LLCsealdata(tag_no).date = LLC_1.date;
    
        sector = NaN(1,length(LLCsealdata(tag_no).cast));
        %%% Identifying the LLC sector of each MEOP profile
        for i = 1:length(LLCsealdata(tag_no).cast)
            if isinterior(LLC_1.polygon, geopointshape(LLCsealdata(tag_no).lat(i), LLCsealdata(tag_no).lon(i)))
                sector(i) = 1;
            elseif isinterior(LLC_2.polygon, geopointshape(LLCsealdata(tag_no).lat(i), LLCsealdata(tag_no).lon(i)))
                sector(i) = 2;
            elseif isinterior(LLC_5.polygon, geopointshape(LLCsealdata(tag_no).lat(i), LLCsealdata(tag_no).lon(i)))
                sector(i) = 5;
            else
                sector(i) = 4;
            end
        end
        LLCsealdata(tag_no).sector = sector;
    
        %%% Initializing matrices for interpolated data
        LLCseal_salt = NaN(length(depth_grid), length(LLCsealdata(tag_no).cast));
        LLCseal_temp = NaN(length(depth_grid), length(LLCsealdata(tag_no).cast));
        LLCseal_vort = NaN(length(depth_grid), length(LLCsealdata(tag_no).cast));
        LLCseal_OW = NaN(size(LLC_1.OW, 3), length(LLCsealdata(tag_no).cast));
    
        for i = 1:length(LLCsealdata(tag_no).cast)
    
            %%% Grabbing LLC data for profile sector
            if sector(i) == 1
                LLC = LLC_1;
            elseif sector(i) == 2
                LLC = LLC_2;
            elseif sector(i) == 4
                LLC = LLC_4;
            elseif sector(i) == 5
                LLC = LLC_5;
            end
    
            %%% Finding LLC points close to MEOP profile
            LLC_lats = LLC.lats;
            LLC_lons = LLC.lons;
            nan_ind = LLC_lats > LLCsealdata(tag_no).lat(i) + 0.025 | LLC_lats < LLCsealdata(tag_no).lat(i) - 0.025 | LLC_lons > LLCsealdata(tag_no).lon(i) + 0.04 | LLC_lons < LLCsealdata(tag_no).lon(i) - 0.04;
            LLC_lats(nan_ind) = NaN;
            close_ind = find(~isnan(LLC_lats));
            [rows, cols] = ind2sub(size(LLC_lats), close_ind);
    
            lats_unfmt = NaN(1,length(rows));
            lons_unfmt = NaN(1,length(rows));
            salt = NaN(86, length(rows));
            temp = NaN(86, length(rows));
            vort = NaN(86, length(rows));
            okubo_weiss = NaN(size(LLC.OW, 3), length(rows));
    
            %%% Extracting LLC data close to MEOP profile
            for j = 1:length(rows)
                lats_unfmt(j) = double(LLC.lats(rows(j), cols(j)));
                lons_unfmt(j) = double(LLC.lons(rows(j), cols(j)));
                salt(:,j) = squeeze(LLC.salt(rows(j), cols(j), :));
                temp(:,j) = squeeze(LLC.temp(rows(j), cols(j), :));
                vort(:,j) = squeeze(LLC.vort(rows(j), cols(j), :));
                okubo_weiss(:,j) = squeeze(LLC.OW(rows(j), cols(j),:));
            end
            temp(salt == 0) = NaN;
            vort(salt == 0) = NaN;
            salt(salt == 0) = NaN;
    
            lats = double(lats_unfmt .* ones(size(salt)));
            lons = double(lons_unfmt .* ones(size(salt)));
            depths = LLC.depth(1:size(salt, 1), 1 ) .* ones(size(salt));
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Interpolating LLC data %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
            %%% Salinity
            S = scatteredInterpolant(lats(:), lons(:), depths(:), salt(:),'linear','none');
            S_prof = S(LLCsealdata(tag_no).lat(i).*ones(size(depth_grid)), LLCsealdata(tag_no).lon(i).*ones(size(depth_grid)), -depth_grid);
            if isempty(S_prof)
                LLCseal_salt(:,i) = NaN(size(depth_grid));
            else
                LLCseal_salt(:,i) = S_prof;
            end
    
            %%% Temperature
            T = scatteredInterpolant(lats(:), lons(:), depths(:), temp(:),'linear','none');
            T_prof = T(LLCsealdata(tag_no).lat(i).*ones(size(depth_grid)), LLCsealdata(tag_no).lon(i).*ones(size(depth_grid)), -depth_grid);
            if isempty(T_prof)
                LLCseal_temp(:,i) = NaN(size(depth_grid));
            else
                LLCseal_temp(:,i) = T_prof;
            end
    
            %%% Vorticity
            V = scatteredInterpolant(lats(:), lons(:), depths(:), vort(:),'linear','none');
            V_prof = V(LLCsealdata(tag_no).lat(i).*ones(size(depth_grid)), LLCsealdata(tag_no).lon(i).*ones(size(depth_grid)), -depth_grid);
            if isempty(V_prof)
                LLCseal_vort(:,i) = NaN(size(depth_grid));
            else
                LLCseal_vort(:,i) = V_prof;
            end
    
            %%% Okubo Weiss
            for j = 1:size(LLC.OW, 3)
                OW = scatteredInterpolant(lats_unfmt', lons_unfmt', okubo_weiss(j,:)');
                OW_prof = OW(LLCsealdata(tag_no).lat(i), LLCsealdata(tag_no).lon(i));
                if isempty(OW_prof)
                    LLCseal_OW(j,i) = NaN;
                else
                    LLCseal_OW(j,i) = OW_prof;
                end
            end
    
            ind = isnan(qc_ts(tag_no).temp(:,i));
            LLCseal_salt(ind,i) = NaN;
            LLCseal_temp(ind,i) = NaN;
            LLCseal_vort(ind,i) = NaN;
    
            clear lats_unfmt lons_unfmt salt vort temp okubo_weiss cols rows nan_ind lats lons depths S T V OW j close_ind S_prof T_prof V_prof OW_prof
        end
    
        LLCsealdata(tag_no).salt = LLCseal_salt;
        LLCsealdata(tag_no).temp = LLCseal_temp;
        LLCsealdata(tag_no).vort = LLCseal_vort;
        LLCsealdata(tag_no).OW = LLCseal_OW;
    
        clear sector i LLCseal LLCseal_salt LLCseal_temp LLCseal_vort LLCseal_OW LLC_lats LLC_lons
    end

    save(string(output_path) + string(date) + '/LLCsealdata', 'LLCsealdata',  'depth_grid');

end
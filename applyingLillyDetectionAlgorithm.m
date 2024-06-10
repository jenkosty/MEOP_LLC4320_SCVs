%%% This script applies the jLab detection algorithm to the LLC4320 snapshots 

%%% Path's for LLC4320 snapshots
input_path = '/Volumes/Elements/LLCsealdata/Snapshot_';
output_path = '/Volumes/Elements/LLCsealdata/Snapshot_';

%%% LLC4320 snapshot dates
snapshot_dates = {'01-Oct-2011','01-Nov-2011', '01-Dec-2011', '01-Jan-2012', '01-Feb-2012', '01-Mar-2012', '01-Apr-2012', '01-May-2012', ...
    '01-Jun-2012', '01-Jul-2012', '01-Aug-2012', '01-Sep-2012'};

for u = 1
    
    date = snapshot_dates{u};
    disp(date);

    runningLillyDetectionAlgorithm(date, input_path, output_path)

end

function runningLillyDetectionAlgorithm(date, input_path, output_path)

    %%% Loading LLCsealdata
    load(string(input_path) + string(date) + '/LLCsealdata');

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

    %%% Tags to test
    test_prof = 1:5; %:length(LLCsealdata);
    
    %%% Looping through tags
    clear lilly_data
    for tag_no = test_prof
        disp(tag_no)
    
        %%% Creating structure for lilly detections
        lilly_data(tag_no).date = LLCsealdata(tag_no).date;
        lilly_data(tag_no).scv = zeros(size(LLCsealdata(tag_no).cast));
        lilly_data(tag_no).scv_reason = string(size(LLCsealdata(tag_no).cast));
        lilly_data(tag_no).contourdata = struct('contour_lat', [], 'contour_lon', [], 'ell_lat', [], 'ell_lon', [], 'exp_ell_lat', [], 'exp_ell_lon', [],...
            'ellpt', [], 'ecc', [], 'OW_in_contour', [],'vort_in_contour', [], 'area', [],'area_ratio', [], 'lat', [], 'lon', [], 'OW', [], 'vort', [], 'LLC_bathymetry', [], 'contour_strength', [], 'OW_of_contour', []);
        lilly_data(tag_no).contourdata(length(LLCsealdata(tag_no).cast)) = struct('contour_lat', [], 'contour_lon', [], 'ell_lat', [], 'ell_lon', [], 'exp_ell_lat', [], 'exp_ell_lon', [],...
            'ellpt', [], 'ecc', [], 'OW_in_contour', [], 'vort_in_contour', [], 'area', [],'area_ratio', [], 'lat', [], 'lon', [], 'OW', [], 'vort', [], 'LLC_bathymetry', [], 'contour_strength', [], 'OW_of_contour', []);
    
        %%% Running lilly detection algorithm
        for i = 1:length(LLCsealdata(tag_no).cast)
            [scv, note, contourdata] = lillySCVdetections(LLCsealdata, tag_no, i, LLC_1, LLC_2, LLC_4, LLC_5);
            lilly_data(tag_no).scv(i) = scv;
            lilly_data(tag_no).scv_reason(i) = note;
            lilly_data(tag_no).contourdata(i) = contourdata;
        end
    
    end
    
    %%% Isolating "deep" lilly detections
    for tag_no = 1:length(lilly_data)
        lilly_data(tag_no).scv_deep = lilly_data(tag_no).scv;
        lilly_data(tag_no).scv_reason_deep = lilly_data(tag_no).scv_reason;
        for i = find(lilly_data(tag_no).scv == 1)
            if mean(lilly_data(tag_no).contourdata(i).LLC_bathymetry(:), 'omitnan') < 400
                lilly_data(tag_no).scv_deep(i) = 0;
                lilly_data(tag_no).scv_reason_deep(i) = 'Too Shallow';
            end
        end
    end
    
    %%% Saving lilly data
    %save(string(output_path) + string(date) + '/lilly_data', 'lilly_data');

end

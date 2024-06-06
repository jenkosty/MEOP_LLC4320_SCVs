%%% This script takes the raw LLC4320 snapshot data and formats it into a
%%% usable format

input_path = '/Volumes/Elements/LLC_Snapshots/SO_snapshots_NSF_LLC4320_k1-86_';
grid_path = '/Volumes/Elements/LLC_Snapshots/SO_grids_NSF_LLC4320_k1-86_';
output_path = '/Volumes/Elements/LLCsealdata/Snapshot_';

%%% Date of LLC snapshots
snapshot_dates = {'01-Oct-2011','01-Nov-2011', '01-Dec-2011', '01-Jan-2012', '01-Feb-2012', '01-Mar-2012', '01-Apr-2012', '01-May-2012', ...
    '01-Jun-2012', '01-Jul-2012', '01-Aug-2012', '01-Sep-2012'};

%%% Looping through snapshot dates 
for i = 1:length(snapshot_dates)
    date = snapshot_dates{i};
    disp(date);
    organizingLLCdata(date, input_path, grid_path, output_path);
end

function organizingLLCdata(date, input_path, grid_path, output_path)

    %%%%%%%%%%%%%%%%%%%%%%%%
    %%% Loading LLC data %%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Clearing structures
    clear LLC_1 LLC_2 LLC_4 LLC_5
    
    %%% Loading mat files
    LLC_1.mat = matfile(string(input_path) + 'face1_' + string(date) + '.mat');
    LLC_2.mat = matfile(string(input_path) + 'face2_' + string(date) + '.mat');
    LLC_4.mat = matfile(string(input_path) + 'face4_' + string(date) + '.mat');
    LLC_5.mat = matfile(string(input_path) + 'face5_' + string(date) + '.mat');
    
    %%% Face 1
    LLC_1.edge_lats = [double(LLC_1.mat.yc(1,:)), double(LLC_1.mat.yc(:,end))', flip(double(LLC_1.mat.yc(end,:))), flip(double(LLC_1.mat.yc(:,1)))'];
    LLC_1.edge_lons = [double(LLC_1.mat.xc(1,:)), double(LLC_1.mat.xc(:,end))', flip(double(LLC_1.mat.xc(end,:))), flip(double(LLC_1.mat.xc(:,1)))'];
    LLC_1.polygon = geopolyshape(LLC_1.edge_lats, LLC_1.edge_lons);
    LLC_1.lats = double(LLC_1.mat.yc);
    LLC_1.lons = double(LLC_1.mat.xc);
    LLC_1.salt = double(LLC_1.mat.s);
    LLC_1.temp = double(LLC_1.mat.t);
    LLC_1.vort = double(LLC_1.mat.Ro);
    LLC_1.depth = LLC_1.mat.rc;
    LLC_1.grid = load(string(grid_path) + 'face1.mat');
    LLC_1.date = date;
    
    %%% Face 2
    LLC_2.edge_lats = [double(LLC_2.mat.yc(1,:)), double(LLC_2.mat.yc(:,end))', flip(double(LLC_2.mat.yc(end,:))), flip(double(LLC_2.mat.yc(:,1)))'];
    LLC_2.edge_lons = [double(LLC_2.mat.xc(1,:)), double(LLC_2.mat.xc(:,end))', flip(double(LLC_2.mat.xc(end,:))), flip(double(LLC_2.mat.xc(:,1)))'];
    LLC_2.polygon = geopolyshape(LLC_2.edge_lats, LLC_2.edge_lons);
    LLC_2.lats = double(LLC_2.mat.yc);
    LLC_2.lons = double(LLC_2.mat.xc);
    LLC_2.salt = double(LLC_2.mat.s);
    LLC_2.temp = double(LLC_2.mat.t);
    LLC_2.vort = double(LLC_2.mat.Ro);
    LLC_2.depth = LLC_2.mat.rc;
    LLC_2.grid = load(string(grid_path) + 'face2.mat');
    LLC_2.date = date;
    
    %%% Face 4
    LLC_4.edge_lats = [double(LLC_4.mat.yc(1,:)), double(LLC_4.mat.yc(:,end))', flip(double(LLC_4.mat.yc(end,:))), flip(double(LLC_4.mat.yc(:,1)))'];
    LLC_4.edge_lons = [double(LLC_4.mat.xc(1,:)), double(LLC_4.mat.xc(:,end))', flip(double(LLC_4.mat.xc(end,:))), flip(double(LLC_4.mat.xc(:,1)))'];
    LLC_4.polygon = geopolyshape(LLC_4.edge_lats, LLC_4.edge_lons);
    LLC_4.lats = double(LLC_4.mat.yc);
    LLC_4.lons = double(LLC_4.mat.xc);
    LLC_4.salt = double(LLC_4.mat.s);
    LLC_4.temp = double(LLC_4.mat.t);
    LLC_4.vort = double(LLC_4.mat.Ro);
    LLC_4.depth = LLC_4.mat.rc;
    LLC_4.grid = load(string(grid_path) + 'face4.mat');
    LLC_4.date = date;
    
    %%% Face 5
    LLC_5.edge_lats = [double(LLC_5.mat.yc(1,:)), double(LLC_5.mat.yc(:,end))', flip(double(LLC_5.mat.yc(end,:))), flip(double(LLC_5.mat.yc(:,1)))'];
    LLC_5.edge_lons = [double(LLC_5.mat.xc(1,:)), double(LLC_5.mat.xc(:,end))', flip(double(LLC_5.mat.xc(end,:))), flip(double(LLC_5.mat.xc(:,1)))'];
    LLC_5.polygon = geopolyshape(LLC_5.edge_lats, LLC_5.edge_lons);
    LLC_5.lats = double(LLC_5.mat.yc);
    LLC_5.lons = double(LLC_5.mat.xc);
    LLC_5.salt = double(LLC_5.mat.s);
    LLC_5.temp = double(LLC_5.mat.t);
    LLC_5.vort = double(LLC_5.mat.Ro);
    LLC_5.depth = LLC_5.mat.rc;
    LLC_5.grid = load(string(grid_path) + 'face5.mat');
    LLC_5.date = date;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Calculating the Okubo Weiss parameter %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Depth Averaged Okubo Weiss calculation
    LLC_1.OW = [];
    LLC_2.OW = [];
    LLC_4.OW = [];
    LLC_5.OW = [];
    LLC_1.OW(:,:) = OWcalculation(LLC_1.grid.dxg(:,:), LLC_1.grid.dyg(:,:), LLC_1.grid.dxc(:,:), LLC_1.grid.dyc(:,:), LLC_1.mat.u(:,:,:), LLC_1.mat.v(:,:,:), LLC_1.depth);
    LLC_2.OW(:,:) = OWcalculation(LLC_2.grid.dxg(:,:), LLC_2.grid.dyg(:,:), LLC_2.grid.dxc(:,:), LLC_2.grid.dyc(:,:), LLC_2.mat.u(:,:,:), LLC_2.mat.v(:,:,:), LLC_2.depth);
    LLC_4.OW(:,:) = OWcalculation(LLC_4.grid.dxg(:,:), LLC_4.grid.dyg(:,:), LLC_4.grid.dxc(:,:), LLC_4.grid.dyc(:,:), LLC_4.mat.u(:,:,:), LLC_4.mat.v(:,:,:), LLC_4.depth);
    LLC_5.OW(:,:) = OWcalculation(LLC_5.grid.dxg(:,:), LLC_5.grid.dyg(:,:), LLC_5.grid.dxc(:,:), LLC_5.grid.dyc(:,:), LLC_5.mat.u(:,:,:), LLC_5.mat.v(:,:,:), LLC_5.depth);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Saving LLC Output Data %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    LLC = {LLC_1, LLC_2, LLC_4, LLC_5};
    sectors = {'LLC_1', 'LLC_2', 'LLC_4', 'LLC_5'};
    for i = 1:4
        lats = LLC{i}.lats;
        save(string(output_path) + string(date) + '/' + string(sectors{i}) + '/lats', 'lats');
        lons = LLC{i}.lons;
        save(string(output_path) + string(date) + '/' + string(sectors{i}) + '/lons', 'lons');
        grid = LLC{i}.grid;
        save(string(output_path) + string(date) + '/' + string(sectors{i}) + '/grid', 'grid');
        temp = LLC{i}.temp;
        save(string(output_path) + string(date) + '/' + string(sectors{i}) + '/temp', 'temp');
        salt = LLC{i}.salt;
        save(string(output_path) + string(date) + '/' + string(sectors{i}) + '/salt', 'salt');
        vort = LLC{i}.vort;
        save(string(output_path) + string(date) + '/' + string(sectors{i}) + '/vort', 'vort');
        depth = LLC{i}.depth;
        save(string(output_path) + string(date) + '/' + string(sectors{i}) + '/depth', 'depth');
        OW = LLC{i}.OW;
        save(string(output_path) + string(date) + '/' + string(sectors{i}) + '/OW', 'OW');
        polygon = LLC{i}.polygon;
        save(string(output_path) + string(date) + '/' + string(sectors{i}) + '/polygon', 'polygon');
        save(string(output_path) + string(date) + '/' + string(sectors{i}) + '/date', 'date');
        disp(i);
    
        clear lats lons grid temp salt vort depth OW polygon
    end
    clear LLC

end
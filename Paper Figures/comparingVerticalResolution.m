%%%% Getting vertical grid for LLC4320
input_path = '/Volumes/Elements/LLCsealdata/Snapshot_';
output_path = '/Volumes/Elements/LLCsealdata/Snapshot_';
snapshot_dates = {'01-Oct-2011','01-Nov-2011', '01-Dec-2011', '01-Jan-2012', '01-Feb-2012', '01-Mar-2012', '01-Apr-2012', '01-May-2012', ...
    '01-Jun-2012', '01-Jul-2012', '01-Aug-2012', '01-Sep-2012'};
i = 1;
date = snapshot_dates{i};
LLC = cell(4,1);
sectors = {'LLC_1', 'LLC_2', 'LLC_4', 'LLC_5'};
load(string(input_path) + string(date) + '/' + string(sectors{i}) + '/depth.mat');
LLCdepth = abs(depth);
LLCdepth = LLCdepth(LLCdepth < 800);
idx = abs(LLCdepth) < 100;
shallow_diff = diff(LLCdepth(idx));
disp('LLC Shallow; 25% = ' + string(prctile(shallow_diff, 25)) + ', 75% = ' + string(prctile(shallow_diff, 75)));

deep_diff = diff(LLCdepth(~idx));
disp('LLC Deep; 25% = ' + string(prctile(deep_diff, 25)) + ', 75% = ' + string(prctile(deep_diff, 75)));


%%
%%% Getting vertical resolution for MEOP data
clear
load('qc_ts.mat') 
shallow_diff = [];
deep_diff = [];
tags = [];
u = 1;
for tag_no = 1:length(qc_ts)
    for i = 1:length(qc_ts(tag_no).cast)
        MEOPpres{u} = qc_ts(tag_no).raw_data(i).pres;
        MEOPpres{u} = MEOPpres{u}(abs(MEOPpres{u}) < 800);
        MEOPpreslength(u) = length(MEOPpres{u});
        if MEOPpreslength(u) < 100
            idx = abs(MEOPpres{u}) < 100;
            shallow_diff = vertcat(shallow_diff, diff(MEOPpres{u}(idx)));
            deep_diff = vertcat(deep_diff, diff(MEOPpres{u}(~idx)));
        end
        u = u + 1;
    end
end

%%
figure()
histogram(MEOPpreslength)
mean(MEOPpreslength(MEOPpreslength <= 100))

%%
figure()
tiledlayout(1,2)

nexttile()
histogram(shallow_diff)
title('Shallow; 10% = ' + string(prctile(shallow_diff, 10)) + ', 90% = ' + string(prctile(shallow_diff, 90)))

nexttile()
histogram(deep_diff)
title('Deep; 10% = ' + string(prctile(deep_diff, 10)) + ', 90% = ' + string(prctile(deep_diff, 90)))



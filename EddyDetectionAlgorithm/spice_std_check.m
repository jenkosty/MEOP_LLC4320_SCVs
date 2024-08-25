function spice_std = spice_std_check(ts_data, tag_no)

%%% Creating array to hold profiles that pass MLD check
spice_std = NaN(size(ts_data(tag_no).cast));

for i = 1:length(ts_data(tag_no).cast)

    casts = [ts_data(tag_no).ref_ind{1,i} ts_data(tag_no).ref_ind{2,i}];
    spice = ts_data(tag_no).ps.spice(:,casts);
    pres = ts_data(tag_no).ps.pres(:,casts);
    MLD = ts_data(tag_no).MLD(i);
    for ii = 1:length(casts)
        spice(pres(:,ii) < MLD, ii) = NaN; % can also make 100 dbar, etc.
    end
    spice_std(i) = mean(std(spice, 0, 2), 'omitnan');

end

end
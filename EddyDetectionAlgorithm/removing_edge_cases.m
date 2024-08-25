function [anticyclones, cyclones, rejected, reason] = removing_edge_cases(ts_data, tag_no)

    %%% Extracting arrays that hold profile classification, rejection, and justification data
    anticyclones = ts_data(tag_no).anticyclones;
    cyclones = ts_data(tag_no).cyclones;
    rejected = ts_data(tag_no).rejected;
    reason = ts_data(tag_no).reason;

    %%% Number of edge profiles to exclude from analysis
    edge_limit = 10;

    %%% Removing flagged profiles that are too close to the start/end of
    %%% the time series

    %%% Anticyclones
    ind = find(anticyclones.MLD < edge_limit | anticyclones.MLD > length(ts_data(tag_no).cast) - edge_limit);
    rejected.anticyclones(anticyclones.MLD(ind)) = 1;
    reason.anticyclones(anticyclones.MLD(ind)) = "Too Close to Edge";
    final = anticyclones.MLD;
    final(ind) = [];
    anticyclones.final = final;

    %%% Cyclones
    ind = find(cyclones.MLD < edge_limit | cyclones.MLD > length(ts_data(tag_no).cast) - edge_limit);
    rejected.cyclones(cyclones.MLD(ind)) = 1;
    reason.cyclones(cyclones.MLD(ind)) = "Too Close to Edge";
    final = cyclones.MLD;
    final(ind) = [];
    cyclones.final = final;

end


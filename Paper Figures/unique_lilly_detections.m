snapshot_dates = {'01-Oct-2011','01-Nov-2011', '01-Dec-2011', '01-Jan-2012', '01-Feb-2012', '01-Mar-2012', '01-Apr-2012', '01-May-2012', ...
    '01-Jun-2012', '01-Jul-2012', '01-Aug-2012', '01-Sep-2012'};

for ii = 1:length(snapshot_dates)
    disp(string(snapshot_dates{ii}))
    load('/Volumes/Elements/LLCsealdata/Snapshot_' + string(snapshot_dates{ii}) + '/lilly_data_new');
    date = snapshot_dates{ii};

    u = 1;
    cnt = 1;
    for tag_no = 1:467
        for i = find(lilly_data(tag_no).scv_new == 1 & strcmp([lilly_data(tag_no).scv_reason], "Normal"))
            if lilly_data(tag_no).contourdata(i).vort_in_contour > 0
                cnt = cnt + 1;
                if u == 1
                    contours_unique_poly(u) = polyshape(lilly_data(tag_no).contourdata(i).contour_lon, lilly_data(tag_no).contourdata(i).contour_lat);
                    contours_unique(u).lat = lilly_data(tag_no).contourdata(i).contour_lat;
                    contours_unique(u).lon = lilly_data(tag_no).contourdata(i).contour_lon;
                    contours_unique(u).contourdata = lilly_data(tag_no).contourdata(i);
                    no_overlapping_contours(u) = 1;
                    u = u + 1;
                else
                    contour_being_tested = polyshape(lilly_data(tag_no).contourdata(i).contour_lon, lilly_data(tag_no).contourdata(i).contour_lat);
                    overlapping = [];
                    for uu = 1:length(contours_unique)
                        if mean(lilly_data(tag_no).contourdata(i).contour_lon) <= (mean(contours_unique(uu).lon) + 1)...
                                || mean(lilly_data(tag_no).contourdata(i).contour_lon) >= (mean(contours_unique(uu).lon) - 1)
                            overlapping(uu) = overlaps(contour_being_tested, contours_unique_poly(uu));
                        end
                    end
                    ind = find(overlapping == 1);
                    no_overlapping_contours(ind) = no_overlapping_contours(ind) + 1;

                    if sum(overlapping) == 0
                        contours_unique_poly(u)= contour_being_tested;
                        contours_unique(u).lat = lilly_data(tag_no).contourdata(i).contour_lat;
                        contours_unique(u).lon = lilly_data(tag_no).contourdata(i).contour_lon;
                        contours_unique(u).contourdata = lilly_data(tag_no).contourdata(i);
                        no_overlapping_contours(u) = 1;
                        u = u + 1;
                    end
                end
            end
        end
    end

    unique_scvs_per_snapshot{1,ii} = contours_unique;
    unique_scvs_per_snapshot{2,ii} = cnt;
    clear contour_being_tested contours_unique_poly uu u overlapping contours_unique
end

%%

%%% Loading data
%load('unique_scvs.mat')
load('AntarcticCoastline_rtopo2.mat')
load('rtopo_1080x310.mat')

for ii = 1:length(unique_cyclones_per_snapshot)
    unique_cyclones(ii) = length(unique_cyclones_per_snapshot{1,ii});
    total_cyclones(ii) = unique_cyclones_per_snapshot{2,ii};
    unique_anticyclones(ii) = length(unique_anticyclones_per_snapshot{1,ii});
    total_anticyclones(ii) = unique_anticyclones_per_snapshot{2,ii};
end

cyclone_contour_lat = []; cyclone_contour_lon = [];
for uu = 1:length(unique_cyclones_per_snapshot{1,1})
    cyclone_contour_lat = [cyclone_contour_lat unique_cyclones_per_snapshot{1,1}(uu).contourdata.contour_lat' NaN];
    cyclone_contour_lon = [cyclone_contour_lon unique_cyclones_per_snapshot{1,1}(uu).contourdata.contour_lon' NaN];
end

anticyclone_contour_lat = []; anticyclone_contour_lon = [];
for uu = 1:length(unique_anticyclones_per_snapshot{1,1})
    anticyclone_contour_lat = [anticyclone_contour_lat unique_anticyclones_per_snapshot{1,1}(uu).contourdata.contour_lat' NaN];
    anticyclone_contour_lon = [anticyclone_contour_lon unique_anticyclones_per_snapshot{1,1}(uu).contourdata.contour_lon' NaN];
end

%%
figure('Position', [100 100 1000 800])
tiledlayout(3,2)
fs = 15;

ax = nexttile(1:2);
hold on
plot(1:12, total_cyclones, '-o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'Color', 'r', 'DisplayName', 'Cyclones')
plot(1:12, total_anticyclones, '-o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'Color', 'b', 'DisplayName', 'Anticyclones')
xticks(1:12)
xticklabels(snapshot_dates)
xlim([1 12])
ylabel('Number of Samples', 'FontSize', fs)
ax.XAxis.FontSize = fs;
ax.YAxis.FontSize = fs;
grid on
legend('Location', 'southeast','FontSize', fs)

nexttile([2,1])
axesm('stereo', 'Origin', [-90 0], 'MapLatLimit', [-90 -60]);
axis off; framem on; hold on;
contourm(YC, XC, coastline, [0 0], 'Fill', 'off', 'Color', 'k', 'LineWidth', 2)
plotm(cyclone_contour_lat, cyclone_contour_lon, 'r', 'LineWidth', 3)

nexttile([2,1])
axesm('stereo', 'Origin', [-90 0], 'MapLatLimit', [-90 -60]);
axis off; framem on; hold on;
contourm(YC, XC, coastline, [0 0], 'Fill', 'off', 'Color', 'k', 'LineWidth', 2)
plotm(anticyclone_contour_lat, anticyclone_contour_lon, 'b', 'LineWidth', 3)



function [scv, note, contourdata] = lillySCVdetections(ts_data, tag_no, i, LLC_1, LLC_2, LLC_4, LLC_5)

%%% Creating structure with all fields
contourdata.contour_lat = NaN;  
contourdata.contour_lon = NaN; 
contourdata.ell_lat = NaN; 
contourdata.ell_lon = NaN;
contourdata.exp_ell_lat = NaN; 
contourdata.exp_ell_lon = NaN;
contourdata.ellpt = NaN;
contourdata.ecc = NaN;
contourdata.OW_in_contour = NaN;
contourdata.vort_in_contour = NaN;
contourdata.area = NaN;
contourdata.area_ratio = NaN;
contourdata.lat = NaN;
contourdata.lon = NaN;
contourdata.OW = NaN;
contourdata.vort = NaN;
contourdata.LLC_bathymetry = NaN;
contourdata.contour_strength = NaN;
contourdata.OW_of_contour = NaN;

%%% Minimum area and OW required for closed contour
minarea = 5;
minOW = -0.5;

%%% Extracting profile sector
sector = ts_data(tag_no).sector(i);
if sector == 1
    LLC = LLC_1;
elseif sector == 2
    LLC = LLC_2;
elseif sector == 4
    LLC = LLC_4;
elseif sector == 5
    LLC = LLC_5;
end

%%% Extracting profile location and closest LLC cell
scv_center_lat = ts_data(tag_no).lat(i);
scv_center_lon = ts_data(tag_no).lon(i);
lat_dist = abs(LLC.lats - scv_center_lat);
lon_dist = abs(LLC.lons - scv_center_lon);
[j,k] = find(lat_dist + lon_dist == min(lat_dist(:) + lon_dist(:)));
clear dist scv_center_lat scv_center_lon

%%% Extracting a square grid surrounding the profile
no_cells = 50;
j_ind = j-no_cells:j+no_cells;
if max(j_ind) > 4320 || min(j_ind) < 1
    scv = 0;
    note = 'Grid Issue';
    return
end
k_ind = k-no_cells:k+no_cells;
if max(k_ind) > 3295 || min(k_ind) < 1
    scv = 0;
    note = 'Grid Issue';
    return
end
clear no_cells

%%% Extracting lat/lon/OW/bathymetry data in square grid
lat = LLC.lats(j_ind, k_ind);
lon = LLC.lons(j_ind, k_ind);
OW = LLC.OW(j_ind, k_ind);
OWnans = OW;
OWnans(OW == 0) = NaN; %%% Assigning NaNs where no data exists
vort = LLC.vort(j_ind, k_ind, :);
vort(vort == 0) = NaN;
depth = LLC.depth;
LLC_bathymetry = zeros(size(lat));
for ii = 1:size(lat,1)
    for jj = 1:size(lat,2)
        idx = ~isnan(vort(ii,jj,:));
        if sum(idx) == 0
            LLC_bathymetry(ii,jj) = 0;
        else
            LLC_bathymetry(ii,jj) = max(abs(depth(idx)));
        end
    end
end

%%% Converting lat/lon grid to distance from profile location
[x,y] = latlon2xy(lat, lon, ts_data(tag_no).lat(i), ts_data(tag_no).lon(i));

%%% Calculating closed curves of low Okubo Weiss
contour_level = minOW*std(OWnans(:), 'omitnan');
if isnan(contour_level)
    scv = 0;
    note = 'NaNs in OW field';
    return
end
[xc, yc] = closedcurves(x, y, OW, contour_level);

%%% Returning if no closed curves found
if isempty(xc) || isempty(yc)
    scv = 0;
    note = 'No closed contours';
    return
end

%%% Calculating area of each closed contour
wgs84 = wgs84Ellipsoid("km");
for ii = 1:length(xc)
    [contour_lat, contour_lon] = xy2latlon(xc{ii}, yc{ii}, ts_data(tag_no).lat(i), ts_data(tag_no).lon(i));
    area_all(ii) = areaint(contour_lat,contour_lon,wgs84);
end
ind = (area_all >= minarea);
xc = xc(ind);
yc = yc(ind);

if isempty(xc)
    scv = 0;
    note =  "Not in Contour";
    return
end

%%% Checking if profile lies in one of the closed curves
for ii = 1:length(xc)
    contour_shape = polyshape(xc{ii}, yc{ii});
    [xo,yo,~,~,~,a,b,theta]=curvemoments(xc{ii},yc{ii});
    A = 2 * a;
    B = 2 * b;
    z = 0;
    [contour_lat, contour_lon] = xy2latlon(xc{ii}, yc{ii}, ts_data(tag_no).lat(i), ts_data(tag_no).lon(i));
    area = areaint(contour_lat,contour_lon,wgs84);

    % figure()
    % h = pcolor(x, y, OWnans); %%% Colorplot
    % set(h, 'EdgeColor', 'none');
    % hold on
    % cellplot(xc,yc,'2c') %%% Plotting contours
    % plot(0, 0, 'Marker', 'o', 'MarkerSize', 8, 'Color', 'g', 'MarkerFaceColor', 'g')
    % colormap(cmocean('balance')); colorbar; clim([-1e-9 1e-9]);
    % [k,l]=ab2kl(a,b);
    % ellipseplot(k,l,theta,xo+sqrt(-1)*yo,'2r') %%% Plotting best-fit ellipses
    % [k,l]=ab2kl(A,B);
    % ellipseplot(k,l,theta,xo+sqrt(-1)*yo,'2r--') %%% Plotting best-fit ellipses

    if isinterior(contour_shape, 0, 0) %%% Normal Contour
        ind = ii;
        scv = 1;
        note = 'Normal';
        break
    elseif inellipse(z, sqrt((a^2+b^2)/2),(a^2-b^2)/(a^2+b^2),theta,xo+(sqrt(-1)*yo)) && area > minarea  %%% Fitted Ellipse
        ind = ii;
        scv = 1;
        note = 'Fitted';
        break
    elseif inellipse(z, sqrt((A^2+B^2)/2),(A^2-B^2)/(A^2+B^2),theta,xo+(sqrt(-1)*yo)) && area > minarea  %%% Fitted and Expanded Ellipse
        ind = ii;
        scv = 1;
        note = 'Expanded + Fitted';
        break
    else %%% Not in Contour
        scv = 0;
        note = 'Not in Contour';
    end
end

%%% Skipping if profile is not enclosed in any contour
if scv == 0
    return
end

%%% Calculating ellipticity
[xo,yo,~,~,~,a,b,theta]=curvemoments(xc{ind},yc{ind});
ellpt = b/a;

%%% Converting ellipticity to eccentricity
eccentricity = ecconv(ellpt, 'ell2ecc');

%%% Calculating ratio between area of contour and fitted ellipse
[k,l]=ab2kl(a,b);
[ellxc,ellyc] = ellcurves(k,l,theta,xo + sqrt(-1)*yo,'npoints',64);
[ell_lat, ell_lon] = xy2latlon(ellxc, ellyc, ts_data(tag_no).lat(i), ts_data(tag_no).lon(i));
ellipse_area = areaint(ell_lat,ell_lon,wgs84);
area_ratio = area / ellipse_area;
clear k l ellipse_area

%%% Calculating magnitude of Okubo-Weiss inside contour
points = geopointshape(lat, lon);
[contour_lat, contour_lon] = xy2latlon(xc{ind}, yc{ind}, ts_data(tag_no).lat(i), ts_data(tag_no).lon(i));
incontour = isinterior(geopolyshape(contour_lat, contour_lon), points);
OW_in_contour = OW(incontour);

%%% Calculating magnitude of depth-averaged vorticity inside contour
ind = (LLC_1.depth < -100) & (LLC_1.depth > -500);
vort_2d = mean(vort(:,:,ind), 3, "omitnan");
vort_in_contour = vort_2d(incontour);
clear points incontour ind

%%% Getting lat/lon data for expanded ellipse
[k,l]=ab2kl(2*a,2*b);
[exp_ellxc,exp_ellyc] = ellcurves(k,l,theta,xo + sqrt(-1)*yo,'npoints',64);
[exp_ell_lat, exp_ell_lon] = xy2latlon(exp_ellxc, exp_ellyc, ts_data(tag_no).lat(i), ts_data(tag_no).lon(i));

%%% Saving data
contourdata.contour_lat = contour_lat;  
contourdata.contour_lon = contour_lon; 
contourdata.ell_lat = ell_lat; % ellipse latitude
contourdata.ell_lon = ell_lon; % ellipse longitude
contourdata.exp_ell_lat = exp_ell_lat; % expanded ellipse latitude
contourdata.exp_ell_lon = exp_ell_lon; % expanded ellipse longitude
contourdata.ellpt = ellpt;
contourdata.ecc = eccentricity;
contourdata.OW_in_contour = median(OW_in_contour);
contourdata.vort_in_contour = median(vort_in_contour, 'omitnan');
contourdata.area = area;
contourdata.area_ratio = area_ratio;
contourdata.lat = lat;
contourdata.lon = lon;
contourdata.OW = OW;
contourdata.vort = vort_2d;
contourdata.LLC_bathymetry = LLC_bathymetry;
contourdata.contour_strength = minOW;
contourdata.OW_of_contour = contour_level;

% %%%%%%%%%%%%%%
% %%% Figure %%%
% %%%%%%%%%%%%%%
%
% figure('Position', [100 100 1500 500]);
% tiledlayout(1,3)
% 
% %%% Okubo-Weiss subplot
% nexttile
% h = pcolor(x, y, OW); %%% Colorplot
% set(h, 'EdgeColor', 'none');
% hold on
% plot(0, 0, 'Marker', 'o', 'MarkerSize', 8, 'Color', 'g', 'MarkerFaceColor', 'g')
% colormap(cmocean('balance')); colorbar; clim([-2e-9 2e-9]);
% %title('Depth-Averaged Okubo-Weiss')
% cellplot(xc,yc,'2b') %%% Plotting contours
% %plot(contourdata.centerlat, contourdata.centerlon, 'Marker', 'o', 'MarkerSize', 8, 'Color', 'm', 'MarkerFaceColor', 'm')
% %plot(contourdata.xc_constantlat, contourdata.yc_constantlat, '-o')
% [xo,yo,~,~,~,a,b,theta]=curvemoments(xc,yc);
% [k,l]=ab2kl(a,b);
% ellipseplot(k,l,theta,xo+sqrt(-1)*yo,'2r') %%% Plotting best-fit ellipses
% [k,l]=ab2kl(2*a,2*b);
% ellipseplot(k,l,theta,xo+sqrt(-1)*yo,'2--r') %%% Plotting expanded best-fit ellipses
% daspect('auto')
% title(string(ellpt))
% 
% %%% Bathymetry subplot
% ax2 = nexttile;
% h = pcolor(x, y, LLC.bathymetry); %%% Colorplot
% set(h, 'EdgeColor', 'none');
% hold on
% contour(x, y, LLC.bathymetry, [-1000 -1000], 'k')
% plot(0, 0, 'Marker', 'o', 'MarkerSize', 8, 'Color', 'g', 'MarkerFaceColor', 'g')
% colormap(ax2, flipud(cmocean('deep'))); colorbar;
% 
% %%% Map subplot
% nexttile;
% load coastlines
% axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -57]);
% axis off; framem on;
% geoshow(coastlat, coastlon, 'Color', 'k')
% scatterm(ts_data(tag_no).lat(i), ts_data(tag_no).lon(i), 20, 'ko', 'MarkerFaceColor', 'g')
% 
% %%% Figure title
% sgtitle('Tag #: ' + string(tag_no) + ', Cast: ' + string(i) + ', Mean OW: ' + string(mean(OW_in_contour) * 1e9) + ', Area: ' + string(area) + ', Ratio: ' + string(ratio) + ' STD: ' + string(std(OWnans(:), 'omitnan')));
% 
% drawnow

end
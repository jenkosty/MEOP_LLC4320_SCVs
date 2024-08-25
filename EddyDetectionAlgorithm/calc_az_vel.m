function max_az_vel = calc_az_vel(ts_data, tag_no, i, LLC_1, LLC_2, LLC_4, LLC_5)

%%% Extracting LLC face in which SCV occured
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

%%% Loading LLC u and v data
LLC.u = double(LLC.mat.u);
LLC.v = double(LLC.mat.v);

%%% Extracting region around SCV center
scv_center_lat = ts_data(tag_no).lat(i);
scv_center_lon = ts_data(tag_no).lon(i);
dist = distance(scv_center_lat, scv_center_lon, LLC.lats, LLC.lons);
[j,k] = find(dist == min(dist(:)));
j_ind = j-8:j+8;
k_ind = k-8:k+8;

%%% Extracting data
lat = LLC.lats(j_ind, k_ind);
lon = LLC.lons(j_ind, k_ind);
u = LLC.u(j_ind, k_ind, :);
v = LLC.v(j_ind, k_ind, :);

figure()
pcolor(lon, lat, LLC.OW(j_ind, k_ind))

%%% Calculating theta
for j = 1:size(lat, 1)
    for k = 1:size(lat, 2)
        theta(i,j) = atan2d((lat(j,k) - scv_center_lat), (lon(j,k) - scv_center_lon));
    end
end

%%% Calculating azimuthal velocity
for j = 1:size(lat, 1)
    for k = 1:size(lat, 2)
        vel_az(j,k,:) = -1.*u(j,k,:).*cosd(theta(j,k)) + v(j,k,:).*sind(theta(j,k));
    end
end

%%% Radially averaging azimuthal velocity
[Zr, R] = radialavg(vel_az(:,:,35), 15);
figure()
plot(R, Zr)
    
end


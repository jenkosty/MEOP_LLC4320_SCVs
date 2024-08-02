u = 1;
for tag_no = 1:467
    for i = find(LLCsealdata(tag_no).scv)
        ecc(u) = LLCsealdata(tag_no).contourdata(i).ecc;
        area(u) = LLCsealdata(tag_no).contourdata(i).area;
        u = u + 1;
    end
end

figure()
tiledlayout(1,2)

nexttile
histogram(ecc)
xlabel('Eccentricity', 'FontSize', 15)

nexttile
histogram(area)
xlabel('Area', 'FontSize', 15)
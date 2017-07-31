function study_maculatus_regeneration

  resolution = 7.25;
  nouter = 8;

  data = '/Users/blanchou/Desktop/SB17.zip';
  imgs = '/Users/blanchou/Desktop/SB17.tif';
  ROIs = ReadImageJROI(data);
  [fins, rays] = parse_fin_ROI(ROIs, resolution);
  valids = ~isnan(fins(:,1));

  fins = fins(valids, :);
  rays = rays(valids);

  nfins = length(rays);
  laterals = NaN(nfins, 2*nouter);
  lengths = NaN(nfins, 2*nouter);
  for j=1:nfins
    if (size(rays{j}, 2) > 2*nouter)
      tmp = [rays{j}(:,1:nouter) rays{j}(:,end-nouter+1:end)];
      tmp(1,~tmp(3,:)) = NaN;
      laterals(j,:) = tmp(1,:);
      lengths(j,:) = tmp(2,:);
    end
  end

  ratios = laterals ./ lengths;
  uninj = ratios(1:2:end,:);
  amput = ratios(2:2:end,:);

  figure;boxplot([uninj(:); amput(:)], [ones(numel(uninj),1); 2*ones(numel(amput),1)])

  keyboard

  return;
end

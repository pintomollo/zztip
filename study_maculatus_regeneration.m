function study_maculatus_regeneration

  fname = '/Users/blanchou/Documents/SB17/Brightfield/SB17_Xmaculatus_amputated/SB17';

  data = [fname '.zip'];
  imgs = [fname '.tif'];
  meta = [fname '.txt'];

  [values, headers] = parse_metadata_table(meta);
  types = find(strncmp('type', headers, 4));
  resol = find(strncmp('resol', headers, 5));
  resol = values{resol};

  ROIs = ReadImageJROI(data);
  %[props, coords] = analyze_ROI(ROIs, 'Slice', 'Area', 'Perimeter');
  [fins, rays] = parse_fin_ROI(ROIs, resol);
  %fins(end+1:length(frames{1}), :) = NaN;

  keyboard

  indxs = unique(values{types});
  disp(indxs)
  return;

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

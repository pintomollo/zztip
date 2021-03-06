function study_maculatus_regeneration

  %fname = '/Users/blanchou/Documents/SB17/Brightfield/SB17_Xmaculatus_amputated/SB17';
  fname = '/Users/blanchou/Documents/SB17/SB17';

  resolution = 7.25;
  nouter = 7;

  all_fish = {};
  all_exper = {};
  all_types = NaN(0, 1);
  all_rays = NaN(0, 2*nouter);
  all_lens = NaN(0, 2*nouter);
  all_segs = NaN(0, 2*nouter);
  all_fins = NaN(0, 4);

  [tmp, prefix1, junk] = fileparts(fname);
  [tmp, prefix2, junk] = fileparts(tmp);
  prefix = [prefix2 '_' prefix1 '_'];

  data = [fname '.zip'];
  fimg = [fname '.tif'];
  meta = [fname '.txt'];

  disp(meta)

  fid = fopen(meta, 'r');
  frames = textscan(fid, '%s %s %d %s %s', 'CommentStyle', '#');
  fclose(fid);

  ROIs = ReadImageJROI(data);
  %[props, coords] = analyze_ROI(ROIs, 'Slice', 'Area', 'Perimeter');
  [fins, rays, paths] = parse_fin_ROI(ROIs, resolution);
  fins(end+1:length(frames{1}), :) = NaN;
  rays(end+1:length(frames{1})) = {[]};
  paths(end+1:length(frames{1})) = {[]};
  periods = cell(size(rays));

  size_rays = cellfun(@size, rays, 'UniformOutput', false);
  size_rays = cat(1, size_rays{:});

  indxs = find(~isnan(fins(:,1)) & size_rays(:,2)>nouter);
  ndata = length(indxs);

    imgs = load_data(fimg, indxs);
    paths = paths(indxs);

    periods(indxs) = get_segment_length(imgs, paths, resolution);

    nfins = length(rays);
    laterals = NaN(nfins, 2*nouter);
    lengths = NaN(nfins, 2*nouter);
    segments = NaN(nfins, 2*nouter);
    for j=1:nfins
      if (size(rays{j}, 2) > 2*nouter)
        tmp = [rays{j}(:,1:nouter) rays{j}(:,end-nouter+1:end)];
        tmp(1,~tmp(3,:)) = NaN;
        laterals(j,:) = tmp(1,:);
        lengths(j,:) = tmp(2,:);
        segments(j,:) = [periods{j}(1:nouter).' periods{j}(end-nouter+1:end).'];
      end
    end

    types = [cumsum(isnan(laterals(:,1:nouter)), 2) cumsum(isnan(laterals(:, end:-1:end-nouter+1)), 2)];
    probs = (types(:,nouter) ~= types(:,end));
    laterals(probs, :) = NaN;
    lengths(probs, :) = NaN;
    segments(probs, :) = NaN;

    all_fish = [all_fish; frames{2}];
    all_exper = [all_exper; frames{4}];
    all_types = [all_types; types(:,end)];
    all_rays = [all_rays; laterals];
    all_lens = [all_lens; lengths];
    all_segs = [all_segs; segments];
    all_fins = [all_fins; fins];

  whole = strncmp(all_exper, '-1', 2);
  fulls = (all_types==0 & ~whole);
  amput = (all_types==3);
  tinys = (all_types==5);

  valids = (whole | fulls | amput | tinys);

  whole_bif = all_rays(whole,:);
  regen_bif = all_rays(fulls,:);

  whole_len = all_lens(whole,:);
  regen_len = all_lens(fulls,:);

  whole_seg = all_segs(whole,:);
  regen_seg = all_segs(fulls,:);

  whole_ratio = whole_bif ./ whole_len;
  regen_ratio = regen_bif ./ regen_len;

  whole_nums = whole_bif ./ whole_seg;
  regen_nums = regen_bif ./ regen_seg;

  whole_no_bif = (whole_ratio==1);
  regen_no_bif = (regen_ratio==1);

  whole_bif(whole_no_bif) = NaN;
  regen_bif(regen_no_bif) = NaN;
  whole_ratio(whole_no_bif) = NaN;
  regen_ratio(regen_no_bif) = NaN;
  whole_nums(whole_no_bif) = NaN;
  regen_nums(regen_no_bif) = NaN;

  figure;
  h1=subplot(1,2,1);notBoxPlot(whole_bif);
  h2=subplot(1,2,2);notBoxPlot(regen_bif);
  linkaxes([h1 h2]);
  ylim([0 9000])

  figure;
  h1=subplot(1,2,1);notBoxPlot(whole_ratio);
  h2=subplot(1,2,2);notBoxPlot(regen_ratio);
  linkaxes([h1 h2]);
  ylim([0 1])

  figure;
  h1=subplot(1,2,1);notBoxPlot(whole_seg);
  h2=subplot(1,2,2);notBoxPlot(regen_seg);
  linkaxes([h1 h2]);
  ylim([0 500])

  figure;
  h1=subplot(1,2,1);notBoxPlot(whole_nums);
  h2=subplot(1,2,2);notBoxPlot(regen_nums);
  linkaxes([h1 h2]);

  figure;
  subplot(1,2,1);notBoxPlot([whole_bif(:); regen_bif(:)], [ones(numel(whole_bif), 1); ones(numel(regen_bif(:)),1)*2]);
  ylim([0 9000])
  subplot(1,2,2);notBoxPlot([whole_ratio(:); regen_ratio(:)], [ones(numel(whole_ratio), 1); ones(numel(regen_ratio(:)),1)*2]);
  ylim([0 1])
  keyboard

  return;

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

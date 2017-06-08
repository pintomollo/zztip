function extract_fish_properties

  fpath = '/Users/blanchou/Documents/SB19/Brightfield/modified_data/SB01-04-05';
  if (~exist(fpath, 'dir'))
    fpath = '/Users/blanchou/Documents/SB19/SB01-04-05';
  end

  data = [fpath '.zip'];
  imgs = [fpath '.tif'];
  props = [fpath '.txt'];

  ROIs = ReadImageJROI(data);
  props = analyze_ROI(ROIs, 'Slice', 'Length', 'PCA');

  scales = isnan(props(:,3));
  vals = props(scales, 1:2);
  slices = unique(vals(:,1));
  fish = ismember(props(:,1), slices);

  lengths = props(fish & ~scales,3) ./ vals(:,2);
  [m, s, n] = mymean(lengths)

  keyboard

  return;
end

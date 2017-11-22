function study_maculatus

  fname = '/Users/blanchou/Documents/SB17/SB17';
  resolution = 7.24;
  nouter = 7;

  data = [fname '.zip'];
  imgs = [fname '.tif'];
  meta = [fname '.txt'];

  %fid = fopen(meta, 'r');
  %frames = textscan(fid, '%s %s %d %s', 'CommentStyle', '#');
  %fclose(fid);

  ROIs = ReadImageJROI(data);
  %[props, coords] = analyze_ROI(ROIs, 'Slice', 'Area', 'Perimeter');
  [fins, rays, paths] = parse_fin_ROI(ROIs, resolution);

  nimgs = size_data(imgs);
  imgs = load_data(imgs, [1:nimgs]);

  img = imgs(:,:,1);
  ray = paths{1}{4};

  [v,p,d] = perpendicular_sampling(img, ray, [], [15,1], []);
  avg = mean(v,2);
  davg = differentiator(avg);

  %[ymax, xmax] = find_extrema(max(avg)-avg,11);
  [ymax, xmax] = find_extrema(davg,11);
  dists = abs(bsxfun(@minus, xmax, xmax'));
  [periods] = find_period(dists(:), [1 100], 4);

  f = fft(davg);
  [fmax, indx] = find_extrema(abs(f(1:floor(end/2))),4);
  indx = length(davg)./indx(:);
  [fmax, ii] = sort(fmax, 'descend');
  periods2=[indx(ii(:)), fmax(:)];

  keyboard

  %mins = get_local_minima(avg);
  %dist = all2all
  %count repeats
  %filter
  %extend

  keyboard

  return;
end

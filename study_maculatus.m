function study_maculatus

  %fname = '/Users/blanchou/Documents/SB17/SB17';
  fname = '/Users/blanchou/Documents/SB17/Brightfield/SB17_Xmaculatus_amputated/SB17';
  resolution = 7.24;
  nouter = 7;

  data = [fname '.zip'];
  fimg = [fname '.tif'];
  meta = [fname '.txt'];

  %fid = fopen(meta, 'r');
  %frames = textscan(fid, '%s %s %d %s', 'CommentStyle', '#');
  %fclose(fid);

  ROIs = ReadImageJROI(data);
  %[props, coords] = analyze_ROI(ROIs, 'Slice', 'Area', 'Perimeter');
  [fins, rays, paths] = parse_fin_ROI(ROIs, resolution);

  indxs = find(~isnan(fins(:,1)));
  ndata = length(indxs);

  all_periods = cell(ndata, 1);

  %nimgs = size_data(imgs);
  %imgs = load_data(imgs, [1:nimgs]);
  imgs = load_data(fimg, indxs);
  rays = rays(indxs);
  paths = paths(indxs);

  for i=1:ndata
    figure;h1=subplot(1,2,1);hold on;h2=subplot(1,2,2);hold on;
    img = imgs(:,:,i);
    curr_rays = paths{i};
    nrays = length(curr_rays);
    periods = NaN(nrays,1);

    subplot(h1);
    hold off;
    imagesc(h1, img); colormap(gray)
    hold on;

    subplot(h2);
    hold off;
    plot(1,1);
    hold on;
    for j=1:nrays

      ray = curr_rays{j};
      plot(h1, ray(:,1), ray(:,2), 'r');

      [pimg,full_ray] = perpendicular_sampling(img, ray, [], [15,1], []);
      avg = nanmean(pimg,2);
      %davg = differentiator(avg);
      davg = movmean(avg, 100) - avg;
      plot(h2, davg);

      %[ymax, xmax] = find_extrema(max(avg)-avg,11);
      %[ymax, xmax] = find_extrema(davg,11);
      %dists = abs(bsxfun(@minus, xmax, xmax'));
      %weight = bsxfun(@plus, ymax, ymax');
      %dists = dists(triu(dists,1)~=0);
      %weight = weight(triu(weight,1)~=0);
      %[periods] = find_period(dists(:), weight(:), [10 80], 4);

      f = fft(davg);

      [fmax, indx] = find_extrema(abs(f(1:floor(end/2))),4);
      [junk, ii] = max(fmax);
      period = length(davg)./(indx(ii)-1);
      periods(j) = period;
      %[fmax, ii] = sort(fmax, 'descend');
      %periods2=[fmax(:), indxs(ii(:))];
      %periods2=periods2(1:min(5,end),:);

      freq = indx(ii);

      %f2 = f;
      %f2(1:ii1-1) = 0;
      f(freq+1:end-freq) = 0;
      %f2(end-ii1+2:end) = 0;
      signal = real(ifft(f));
      [ymax, xmax] = find_extrema(signal,period/2);
      goods = (xmax >= period/2 & xmax <= length(signal)-period/2);
      xmax = xmax(goods);
      ymax = ymax(goods);

      plot(h2, signal)
      scatter(h2, xmax, ymax);
      scatter(h1, full_ray(xmax, 1), full_ray(xmax, 2));

    end
    all_periods{i} = periods;
    mper = median(periods);

    wrongs = (abs(periods - mper)>mper/2);

    if (any(wrongs))
      for j=find(wrongs).'
        ray = curr_rays{j};
        [pimg,full_ray] = perpendicular_sampling(img, ray, [], [15,1], []);
        avg = nanmean(pimg,2);
        davg = movmean(avg, 100) - avg;

        %% FFT select around mper !!
        keyboard
      end
    end
  end

  %mins = get_local_minima(avg);
  %dist = all2all
  %count repeats
  %filter
  %extend

  return;
end

function study_maculatus

  %fname = '/Users/blanchou/Documents/SB17/SB17';
  fname = '/Users/blanchou/Documents/SB17/Brightfield/SB17_Xmaculatus_amputated/SB17';
  resolution = 7.24;

  %fname = '/Users/blanchou/Documents/SB01/Brightfield/SB01_Leica_amputated_AB368/modified_data/SB01';
  %resolution = 5.64;

  fname = '/Users/blanchou/Documents/SB05/Brightfield/SB05_Leica_amputated_calcein_SA120814/modified_data/SB05';
  resolution = 5.64;

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

  size_rays = cellfun(@size, rays, 'UniformOutput', false);
  size_rays = cat(1, size_rays{:});

  indxs = find(~isnan(fins(:,1)) & size_rays(:,2)>nouter);
  ndata = length(indxs);

  imgs = load_data(fimg, indxs);
  rays = rays(indxs);
  paths = paths(indxs);

  periods = get_segment_length(imgs, paths);
  keyboard
  %{
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
      f(1:freq-1) = 0;
      f(freq+1:end-freq) = 0;
      f(end-freq+2:end) = 0;

      try
      signal = real(ifft(f));
      [ymax, xmax] = find_extrema(signal,period/2);
      catch
        keyboard
      end
      goods = (xmax >= period/2 & xmax <= length(signal)-period/2);
      xmax = xmax(goods);
      ymax = ymax(goods);

      plot(h2, signal)
      scatter(h2, xmax, ymax);
      scatter(h1, full_ray(xmax, 1), full_ray(xmax, 2));

    end
    all_periods{i} = periods;
  end

  mper = median(cat(1, all_periods{:}));

  for i=1:ndata
    figure;h1=subplot(1,2,1);hold on;h2=subplot(1,2,2);hold on;
    img = imgs(:,:,i);
    curr_rays = paths{i};
    nrays = length(curr_rays);
    periods = all_periods{i};

    subplot(h1);
    hold off;
    imagesc(h1, img); colormap(gray)
    hold on;

    subplot(h2);
    hold off;
    plot(1,1);
    hold on;

    wrongs = (abs(periods - mper)>mper/2);

    if (any(wrongs))
      for j=find(wrongs).'
        ray = curr_rays{j};

        plot(h1, ray(:,1), ray(:,2), 'r');

        [pimg,full_ray] = perpendicular_sampling(img, ray, [], [15,1], []);
        avg = nanmean(pimg,2);
        davg = movmean(avg, 100) - avg;

        %% FFT select around mper !!
        range = 1 + length(davg) ./ ([1.25 0.75]*mper);
%        weight = exp(-([1:floor(length(davg)/2)].'-mean(range)).^2/(8*sqrt(mean(range))));

        f = fft(davg);

%        [fmax, indx] = find_extrema(abs(f(1:floor(end/2))).*weight,4);
        [fmax, indx] = find_extrema(abs(f(1:floor(end/2))),4);
        goods = (indx >= range(1) & indx <= range(2));

        fmax = fmax(goods);
        indx = indx(goods);

        if (~any(goods))
          weight = exp(-([1:floor(length(davg)/2)].'-mean(range)).^2/(8*sqrt(mean(range))));
          [fmax, indx] = find_extrema(abs(f(1:floor(end/2))).*weight,4);
        end


        [junk, ii] = max(fmax);
        period = length(davg)./(indx(ii)-1);
        periods(j) = period;

        freq = indx(ii);

        f(1:freq-1) = 0;
        f(freq+1:end-freq) = 0;
        f(end-freq+2:end) = 0;

        signal = real(ifft(f));
        [ymax, xmax] = find_extrema(signal,period/2);
        goods = (xmax >= period/2 & xmax <= length(signal)-period/2);
        xmax = xmax(goods);
        ymax = ymax(goods);

        plot(h2, signal)
        scatter(h2, xmax, ymax);
        scatter(h1, full_ray(xmax, 1), full_ray(xmax, 2));

      end
    end
  end

  keyboard

  %mins = get_local_minima(avg);
  %dist = all2all
  %count repeats
  %filter
  %extend
  %}

  return;
end

function all_periods = get_segment_length(imgs, rays, params)
% GET_SEGMENT_LENGTH estimates the average length of the ray segments using
%   the Fourier transform of the averaged perpendicular sampling of each ray
%
%   [LENGTHS] = GET_SEGMENT_LENGTH(IMGS, RAYS) computes the LENGTHS of each
%   RAYS on every IMGS. IMGS should be a stack and RAYS a cell list of cells.
%

  width = 15;
  navg = 100;
  if (nargin == 3 && ~isempty(params))
    width = params(1);
    if numel(params) > 1
      navg = params(2);
    end
  end

  ndata = length(rays);
  all_periods = cell(ndata, 1);

  for i=1:ndata
    img = imgs(:,:,i);
    curr_rays = rays{i};
    nrays = length(curr_rays);
    periods = NaN(nrays,1);

    for j=1:nrays

      ray = curr_rays{j};
      periods(j) = get_period(img, ray, [], width, navg);
    end

    all_periods{i} = periods;
  end

  mper = median(cat(1, all_periods{:}));

  for i=1:ndata
    img = imgs(:,:,i);
    curr_rays = rays{i};
    nrays = length(curr_rays);
    periods = all_periods{i};

    wrongs = (abs(periods - mper)>mper/2);

    if (any(wrongs))
      for j=find(wrongs).'
        ray = curr_rays{j};

        periods(j) = get_period(img, ray, mper, width, navg);
      end

      all_periods{i} = periods;
    end
  end

  return;
end

function period = get_period(img, ray, mper, width, navg)

  [pimg,full_ray] = perpendicular_sampling(img, ray, [], [width,1], []);
  avg = nanmean(pimg,2);
  davg = movmean(avg, navg) - avg;

  f = fft(davg);

  [fmax, indx] = find_extrema(abs(f(1:floor(end/2))),4);

  if (~isempty(mper))

    range = 1 + length(davg) ./ ([1.25 0.75]*mper);
    goods = (indx >= range(1) & indx <= range(2));

    fmax = fmax(goods);
    indx = indx(goods);

    if (~any(goods))
      weight = exp(-([1:floor(length(davg)/2)].'-mean(range)).^2/(8*sqrt(mean(range))));
      [fmax, indx] = find_extrema(abs(f(1:floor(end/2))).*weight,4);
    end
  end

  [junk, ii] = max(fmax);

  period = length(davg)./(indx(ii)-1);

  return;
end

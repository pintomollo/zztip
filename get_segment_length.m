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
  imgs = double(imgs);
  imgs = imgs - gaussian_mex(imgs, 10);
  all_sizes = NaN(ndata, 1);

  for i=1:ndata
    img = imgs(:,:,i);
    curr_rays = rays{i};
    nrays = length(curr_rays);
    all_sizes(i) = nrays;
    periods = NaN(nrays,1);

    for j=1:nrays

      ray = curr_rays{j};
      periods(j) = get_period(img, ray, 45, width);
    end

    mper = movmedian(periods, 4);
    mper = movmean((mper + mper(end:-1:1))/2, 4);
    sper = std(mper);

    for j=1:nrays
      ray = curr_rays{j};

      periods(j) = get_period(img, ray, [mper(j) 4*sper], width);
    end

    all_periods{i} = periods;
  end

  max_size = max(all_sizes);
  period_shape = NaN(max_size, ndata);
  for i=1:ndata
    if (all_sizes(i)==max_size)
      period_shape(:,i) = all_periods{i}(:);
    else
      period_shape(:,i) = interp1([0:all_sizes(i)-1].'/(all_sizes(i)-1), all_periods{i}(:), [0:max_size-1]/(max_size-1));
    end
  end
  [all_mper,all_sper] = mymean([period_shape period_shape(end:-1:1, :)], 2);

  for i=1:ndata
    img = imgs(:,:,i);
    curr_rays = rays{i};
    nrays = length(curr_rays);
    periods = all_periods{i};

    if (all_sizes(i)==max_size)
      mper = all_mper;
      sper = all_sper;
    else
      tmp = interp1([0:max_size-1].'/(max_size-1), [all_mper all_sper], [0:all_sizes(i)-1]/(all_sizes(i)-1));
      mper = tmp(:,1);
      sper = tmp(:,2);
    end

    for j=1:nrays
      ray = curr_rays{j};

      periods(j) = get_period(img, ray, [mper(j) sper(j)], width);
    end

    all_periods{i} = periods;
  end

  return;
end

function period = get_period(img, ray, mper, width)

  [pimg] = perpendicular_sampling(img, ray, [], [width,1], []);
  pimg(isnan(pimg)) = 0;

  lavg = size(pimg,1);
  pos = [1:lavg];

  f = fft(pimg);
  if (~isempty(mper))

    if (numel(mper) > 1)
      range = 1 + lavg ./ (mper(1)+[0 +1 -1]*mper(2));
    else
      range = 1 + lavg ./ ([1 1.25 0.75]*mper);
    end

    weight = exp(-(pos(:) - range(1)).^2/(0.5*diff(range(2:3)).^2));
    f = bsxfun(@times, f, weight);
  end
  [fmax, indx] = find_extrema(abs(mean(f(1:floor(lavg/2), :), 2)),2);
  valids = (indx > 3);

  fmax = fmax(valids);
  indx = indx(valids);

  [junk, ii] = max(fmax);

  period = lavg./(indx(ii)-1);

  return;
end

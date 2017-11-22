function [values, perp_path, dpos] = perpendicular_sampling(img, path, perp_path, dpos, opts)

  if (nargin == 2)
    [path, opts] = deal(img, path);
    img = [];
    perp_path = [];
    dpos = [];
  elseif (nargin == 3)
    opts = perp_path;
    perp_path = [];
    dpos = [];
  elseif (nargin == 4)
    if (all(size(path) == size(perp_path)))
      opts = dpos;
      dpos = [];
    else
      opts = dpos;
      dpos = perp_path;
      perp_path = [];
    end
  end

  if (isempty(perp_path))

    [dist, tot_dist] = carth2linear(path, opts);

    if (isempty(dist))
      values = [];
      perp_path = [];
      dpos = [];

      return;
    end

    dist = dist * tot_dist;
    full_dist = [0:tot_dist].';
    path = interp1(dist, path, full_dist);
    dist = full_dist;

    %deriv = differentiator([dist(:), dist(:)], path, 'circular');
    deriv = differentiator([dist(:), dist(:)], path);
    perp_path = [-deriv(:,2), deriv(:,1)];

    nulls = all(perp_path == 0, 2);

    if (any(nulls))
      dpts = diff(path([1:end 1], :));
      perp_path(nulls, :) = dpts(nulls, :);
    end

    perp_path = bsxfun(@rdivide, perp_path, hypot(perp_path(:,1), perp_path(:,2))); 
  end

  if (isempty(dpos))
    bin_dist = 64;
    bin_step = 1;
    dpos = [-bin_dist:bin_step:bin_dist];
  elseif (numel(dpos) == 2 & dpos(1) > dpos(2))
    bin_dist = dpos(1);
    bin_step = dpos(2);
    dpos = [-bin_dist:bin_step:bin_dist];
  end

  if (isempty(img))
    values = perp_path;
    perp_path = dpos;
    dpos = dist;

    return;
  end

  all_pos_x = bsxfun(@plus, perp_path(:,1) * dpos, path(:,1));
  all_pos_y = bsxfun(@plus, perp_path(:,2) * dpos, path(:,2));

  all_pos_x = all_pos_x(:);
  all_pos_y = all_pos_y(:);

  values = imbilinear(double(img), all_pos_x, all_pos_y);
  values = reshape(values, size(path,1), []);

  return;
end

function [linear, indexes, total_distance] = carth2linear(varargin)

  [pts_x, pts_y, ruffles_path, center_posterior, opts] = parse_inputs(varargin{:});

  if (isempty(pts_y))
    if (isempty(pts_x))
      linear = [];
      total_distance = NaN;

      return;
    else
      pts_y = pts_x(:, 2);
      pts_x = pts_x(:, 1);
    end
  end

  %circular = true;

  %[pts_x, pts_y] = insert_ruffles(pts_x, pts_y, ruffles_path);

  %if (pts_x(1) ~= pts_x(end) | pts_y(1) ~= pts_y(end))
  %  pts_x = [pts_x; pts_x(1)];
  %  pts_y = [pts_y; pts_y(1)];

  %  circular = false;
  %end

  distance = sqrt(diff(pts_x).^2 + diff(pts_y).^2);
  cum_dist = [0; cumsum(distance)];

  total_distance = cum_dist(end);

  if (total_distance == 0)
    linear = 0;

    return;
  end

  linear = cum_dist / total_distance;

  %if (~circular)
  %  linear = linear(1:end-1);
  %end

  if (false)
  %if (center_posterior)
    thresh = opts.quantification.pole_threshold;

    goods = (pts_x > 0);
    max_pos = max(pts_x(goods));
    goods = (pts_x > (1-thresh)*max_pos);

    max_r = min(abs(pts_y(goods)));
    post_indx = find(goods & abs(pts_y) == max_r, 1);

    goods = (pts_x < 0);
    max_pos = min(pts_x(goods, 1));
    goods = (pts_x < (1-thresh)*max_pos);

    max_r = min(abs(pts_y(goods)));
    ant_indx = find(goods & abs(pts_y) == max_r, 1);

    linear = linear - linear(post_indx);
    if (ant_indx > post_indx)
      linear = [linear(ant_indx:end, :) - 1; linear(1:ant_indx-1, :)];
    else
      linear = [linear(ant_indx:end, :); linear(1:ant_indx-1, :) + 1];
    end

    %dists = [linear(1), total_dist + linear(1)];
    indexes = [ant_indx:length(linear) 1:ant_indx-1];
  else
    indexes = total_distance;
    total_distance = [];
  end

  return;
end

function [pts_x, pts_y, ruffles_path, center_posterior, opts] = parse_inputs(varargin)

  pts_x = [];
  pts_y = [];
  ruffles_path = {};
  center_posterior = false;
  opts = [];

  for i=1:length(varargin)
    var_type = class(varargin{i});
    switch var_type
      case {'double', 'single', 'int8', 'int16', 'int32', 'int64', 'uint8', 'uint16', 'uint32', 'uint64'}
        if (isempty(pts_x))
          pts_x = varargin{i};
        else
          pts_y = varargin{i};
        end
      case 'logical'
        center_posterior = varargin{i};
      case 'cell'
        ruffles_path = varargin{i};
      case 'struct'
        opts = varargin{i};
    end
  end

  %if (isempty(opts))
  %  opts = get_struct('ASSET');
  %end

  return;
end

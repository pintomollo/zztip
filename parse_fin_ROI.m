function [fins, rays] = parse_fin_ROI(fname, resol)

  ROIs = ReadImageJROI(fname);

  if (~iscell(ROIs))
    ROIs = {ROIs};
  end

  ndata = length(ROIs);
  props = NaN(ndata, 7);
  data = cell(ndata, 1);

  for i=1:ndata
    props(i, 1) = ROIs{i}.nPosition;
    pos = ROIs{i}.mnCoordinates;

    switch ROIs{i}.strType
      case 'PolyLine'
        props(i, 2) = 1;
        props(i, 3) = sum(sqrt(sum(diff(pos).^2, 2)));

        if (pos(1,1) < pos(2,1))
          pos = pos([end:-1:1], :);
        end

        props(i, 4:end) = reshape(pos(1:2, :).', 4, 1).';
      case 'Polygon'
        props(i, 2) = 2;
        props(i, 3) = polyarea(pos(:,1), pos(:,2));
      otherwise
        display(ROIs{i}.strType)
    end

    data{i} = pos;
  end

  ids = unique(props(:,1));
  ids = ids(~isnan(ids));
  nids = max(ids);

  fins = NaN(nids, 4);
  rays = cell(nids, 1);

  for i=1:nids
    curr_id = i;
    curr_rays = (props(:,1) == curr_id & props(:,2) == 1);
    curr_fin = (props(:,1) == curr_id & props(:,2) == 2);

    if (any(curr_fin))

      curr_props = props(curr_rays, 3:end);
      [junk, indx] = unique(curr_props(:,[3 1]), 'rows');
      indx = indx(end:-1:1);
      nrays = length(indx);

      curr_props = curr_props(indx, :);
      curr_data = data(curr_rays);
      curr_data = curr_data(indx);

      ymin = min(cellfun(@(x)(min(x(:,2))), curr_data));
      ymax = max(cellfun(@(x)(max(x(:,2))), curr_data));

      dists = bsxfun(@minus, curr_props(:,2), curr_props(:,2).').^2 + ...
              bsxfun(@minus, curr_props(:,3), curr_props(:,3).').^2;

      thresh = median(min(dists+max(dists(:))*(dists==0)))/5;

      curr_rays = NaN(3, nrays);
      count = 1;
      for j=1:nrays

        hits = dists(j, :) < thresh;
        if (any(hits(j+1:end)))
          l = sort(curr_props(hits, 1));
          curr_rays(1:2, count) = l;
          count = count + 1;
        elseif ~any(hits(1:j-1))
          curr_rays(1:2, count) = curr_props([j j]);
          count = count + 1;
        end
      end

      curr_rays = curr_rays(:, 1:count-1);

      is_full = (curr_rays(2,:) > max(curr_rays(1,:))/2);
      curr_rays(3,:) = is_full;

      fins(curr_id, :) = [props(curr_fin, 3) max(curr_rays(2, :)) min(curr_rays(2, :)) ymax-ymin];
      rays{curr_id} = curr_rays;
    end
  end

  return;
end

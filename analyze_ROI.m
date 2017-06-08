function [props, coords] = analyze_ROI(ROIs, varargin)
% Slice Area Perimeter Centroid BoundingBox PCA

  do_props = false(1,6);
  resol = 1;

  for i=1:nargin-1
    switch varargin{i}
      case 'Slice'
        do_props(1) = true;
      case 'Area'
        do_props(2) = true;
      case {'Perimeter', 'Length'}
        do_props(3) = true;
      case 'BoundingBox'
        do_props(4) = true;
      case 'Centroid'
        do_props(5) = true;
      case 'PCA'
        do_props(6) = true;
      otherwise
        if (isnumeric(varargin{i}))
          resol = varargin{i};
        else
          disp(['Unknown property ''' varargin{i} ''', skipping.']);
        end
    end
  end

  nROIs = length(ROIs);
  props_indx = cumsum(do_props .* [1 1 1 4 2 2]);
  props = NaN(nROIs, props_indx(end));
  coords = cell(nROIs, 1);

  for i=1:nROIs
    type = ROIs{i}.strType;
    if (strncmp(type, 'Line', 4))
      pos = ROIs{i}.vnLinePoints([1 2; 3 4]) * resol;
    else
      pos = ROIs{i}.mnCoordinates * resol;
    end

    if (do_props(1))
      props(i, 1) = ROIs{i}.nPosition;
    end

    if (do_props(2) && strncmp(type, 'Polygon', 7))
      props(i, props_indx(2)) = polyarea(pos(:,1), pos(:,2));
    end

    if (do_props(3))
      props(i, props_indx(3)) = sum(sqrt(sum(diff(pos).^2, 2)));
    end

    if (do_props(4))
      props(i, props_indx(4)-[3:-1:0]) = [min(pos, [], 1) max(pos, [], 1)];
    end

    if (do_props(5))
      props(i, props_indx(5)-[1 0]) = mean(pos, 1);
    end

    if (do_props(6))
      if (size(pos, 1) > size(pos, 2))
        eig = pca(pos);
        rot = pos*eig;

        mins = min(rot, [], 1);
        maxs = max(rot, [], 1);

        props(i, props_indx(6)-[1 0]) = maxs-mins;
      end
    end

    coords{i} = pos;
  end

  return;
end

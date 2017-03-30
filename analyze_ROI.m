function [props, coords] = analyze_ROI(ROIs, varargin)
% Slice Area Perimeter Centroid BoundingBox

  do_props = false(1,5);

  for i=1:nargin-1
    switch varargin{i}
      case 'Slice'
        do_props(1) = true;
      case 'Area'
        do_props(2) = true;
      case 'Perimeter'
        do_props(3) = true;
      case 'BoundingBox'
        do_props(4) = true;
      case 'Centroid'
        do_props(5) = true;
      otherwise
        disp(['Unknown property ''' varargin{i} ''', skipping.']);
    end
  end

  nROIs = length(ROIs);
  props_indx = cumsum(do_props .* [1 1 1 4 2])
  props = NaN(nROIs, props_indx(end));
  coords = cell(nROIs, 1);

  for i=1:nROIs
    type = ROIs{i}.strType;
    pos = ROIs{i}.mnCoordinates;

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

    coords{i} = pos;
  end

  return;
end

% ANALYZE_ROI parses the ImageJ/Fiji ROIs to compute the requested
% properties and return the corresponding path.
%
%   [PROPS, PATHS] = ANALYZE_ROI(ROIs, TYPES) parses the ROIs to calculate
%                     the requested TYPES of properties. Available TYPES are
%                     'Slice', 'Area', 'Perimeter', 'Centroid', 'BoundingBox',
%                     and 'PCA'. TYPES are provided as string arguments. PROPS
%                     is then a NxP matrix with N ROIs and P properties. PATHS
%                     is a Nx1 cell matrix.
%
%   [...] = ANALYZE_ROI(ROIs, RESOL, ...) provides the pixel size to scale the
%                     properties.
%
% Blanchoud Group, UNIFR
% Simon Blanchoud
% 08/06/2020
function [props, coords] = analyze_ROI(ROIs, varargin)

  do_props = false(1,6);
  resol = 1;

  for i=1:nargin-1
    if (isnumeric(varargin{i}))
      resol = varargin{i};
    elseif (ischar(varargin{i}))
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
          disp(['Unknown property ''' varargin{i} ''', skipping.']);
      end
    end
  end

  nROIs = length(ROIs);
  props_indx = cumsum(do_props .* [1 1 1 4 2 2]);
  props = NaN(nROIs, props_indx(end));
  coords = cell(nROIs, 1);
  resol = [resol(:); ones(nROIs - numel(resol), 1)*resol(1)];

  for i=1:nROIs
    type = ROIs{i}.strType;
    if (strncmp(type, 'Line', 4))
      pos = ROIs{i}.vnLinePoints([1 2; 3 4]) * resol(i);
    else
      pos = ROIs{i}.mnCoordinates * resol(i);
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

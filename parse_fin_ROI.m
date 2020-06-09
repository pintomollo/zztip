% PARSE_FIN_ROI parses the ImageJ/Fiji ROIs of segmented fins to sort and
% analyze them based on several properties.
%
% [FINS, RAYS, PATHS] = PARSE_FIN_ROI(ROIs, RESOL) parses all the ROIs of
%                       segmented fins imaged with a pixel size of RESOL.
%                       FINS is a Nx5 matrix, with N segmented fins and measures
%                       of [Area LongestRay ShortestRay FinWidth PaddleWidth].
%                       RAYS is a Nx1 cell matrix with each cell containing a
%                       3xM matrix with M rays as [Bifurcation Length isFull].
%                       PATHS is a Nx1 cell matrix containing the 2D coordinates
%                       of the full ray.
%
% Blanchoud Group, UNIFR
% Simon Blanchoud
% 08/06/2020
function [fins, rays, paths] = parse_fin_ROI(ROIs, resol)

  % Parse the ROIs to extract useful features
  [props, data] = analyze_ROI(ROIs, resol, 'Slice', 'Area', 'Perimeter');

  % Intermediate data
  ndata = length(ROIs);
  coords = NaN(ndata, 2);

  % Make sure the rays are aligned properly
  for i=1:ndata
    pos = data{i};
    if (pos(1,1) < pos(2,1))
      pos = pos([end:-1:1], :);
      data{i} = pos;
    end

    % Keep the position of the base of each ray
    coords(i,:) = pos(1,:);
  end
  
  % Find the maximal number of frames to analyse
  ids = unique(props(:,1));
  ids = ids(~isnan(ids));
  nids = max(ids);

  % Prepare the output data
  fins = NaN(nids, 5);
  rays = cell(nids, 1);
  paths = cell(nids, 1);

  % Loop through each fin
  for i=1:nids

    % Get the indexes of the current ROI data
    curr_id = i;
    curr_rays = (props(:,1) == curr_id & isnan(props(:,2)));
    curr_fin = (props(:,1) == curr_id & ~isnan(props(:,2)));

    % Make sure we have something to analyse
    if (any(curr_fin))

      % Get the actual data
      curr_props = props(curr_rays, :);
      curr_coords = coords(curr_rays, :);

      % Differentiate between rays and bifurcations
      % First we get only single starting points around the peduncle
      [junk, indx] = unique([curr_coords(:, 2), curr_props(:,3)], 'rows');
      indx = indx(end:-1:1);
      nrays = length(indx);

      % Then we measure the distance between the remaining rays to
      % see if two are close enough to be considered the same ray
      curr_props = curr_props(indx, :);
      curr_coords = curr_coords(indx, :);

      dists = bsxfun(@minus, curr_coords(:,1), curr_coords(:,1).').^2 + ...
              bsxfun(@minus, curr_coords(:,2), curr_coords(:,2).').^2;

      thresh = median(diag(dists,1))/5;

      % We then collect the actual segmentations
      curr_data = data(curr_rays);
      curr_data = curr_data(indx);

      % Prepare the intermediate data to collect the proper rays
      curr_rays = NaN(3, nrays);
      curr_paths = cell(1, nrays);
      count = 1;

      % We loop through all rays
      for j=1:nrays

        % We check for very close rays
        hits = dists(j, :) < thresh;

        % If we find one, it might be either the bifurcation or the full ray
        if (any(hits(j+1:end)))

          % We sort and store both lengths
          [l, ii] = sort(curr_props(hits, 3));
          curr_rays(1:2, count) = l(1:2);

          % And we keep the full segmentation only
          tmp_data = curr_data(hits);
          curr_paths{count} = tmp_data{ii(2)};
          count = count + 1;

        % Otherwise, and if it's not been processed perviously, then it's
        % unbifurcated ray
        elseif ~any(hits(1:j-1))

          % In that case, the bifurcation length is identical to the ray length
          curr_rays(1:2, count) = curr_props([j j], 3);
          curr_paths{count} = curr_data{j};
          count = count + 1;
        end
      end

      % We remove the extra empty spaces
      curr_rays = curr_rays(:, 1:count-1);
      curr_paths = curr_paths(:, 1:count-1);

      % We pre-compute whether the ray is bifurcated
      is_full = (curr_rays(2,:) > max(curr_rays(1,:))/2);
      curr_rays(3,:) = is_full;

      % We keep the outer most rays for statistics
      tmp_paths = curr_paths(is_full);

      % We compute the properties of the fin
      fins(curr_id, :) = [props(curr_fin, 2) max(curr_rays(2, :)) ...
                          min(curr_rays(2, is_full)) ...
                          sqrt(sum((tmp_paths{1}(end,:) - tmp_paths{end}(end,:)).^2)) ...
                          sqrt(sum((curr_paths{1}(1,:) - curr_paths{end}(1,:)).^2))];

      % Store the output data into the proper variables
      rays{curr_id} = curr_rays;
      paths{curr_id} = curr_paths;

    end
  end

  return;
end

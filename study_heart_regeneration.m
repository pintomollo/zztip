function [ratios, files] = study_heart_regeneration(do_export)

  if (nargin == 0)
    do_export = false;
  end

  files = find_segmentations('/Users/blanchou/Documents/SB07/Histology/modified_data/');

  title_name = 'SB07';

  ratios = NaN(length(files), 2);
  dpci = NaN(length(files), 1);

  hits = regexp(files, '.*SB\d+_\w(\d+)dpci_.*', 'tokens');

  for i=1:length(files)
    fname = files{i};

    [tmp, prefix1, junk] = fileparts(fname);
    [tmp, prefix2, junk] = fileparts(tmp);
    prefix = [prefix2 '_' prefix1];

    data = [fname '.zip'];
    imgs = [fname '.tif'];

    ROIs = ReadImageJROI(data);
    props = analyze_ROI(ROIs, 'Slice', 'Area');
    [props, indxs] = sortrows(props, [1 -2]);
    ROIs = ROIs(indxs);

    is_injury = ([1; diff(props(:,1))] == 0);
    heart = props(~is_injury, :);
    injury = props(is_injury, :);

    heart = interpolate(heart);
    injury = interpolate(injury);

    largest = sortrows(injury, -2);
    if (size(largest, 1) > 3)
      largest = largest(1:3,:);
    end
    largest = largest(:,1);

    ratios(i, 1) = sum(injury(:,2))/sum(heart(:,2));
    ratios(i, 2) = sum(injury(ismember(injury(:,1), largest),2))/sum(heart(ismember(heart(:,1), largest),2));

    dpci(i) = str2double(hits{i}{1});

    if (do_export)
      export_ROI(prefix, ROIs, imgs);
    end
  end

  ratios(:,end+1) = dpci;

  display_ratios(ratios(:,1), dpci, title_name)
  display_ratios(ratios(:,2), dpci, title_name)

  return;
end

function display_ratios(ratios, dpci, title_name)

  [H, p] = myttest(ratios, dpci);
  ratios = ratios*100;

  ids = unique(dpci(:));
  groups = {};
  pvals = [];
  for i = 1:length(ids)-1
    for j = i+1:length(ids)
      if (p(i,j)<0.05)
        groups{end+1} = ids([i j]);
        pvals(end+1) = p(i,j);
      end
    end
  end

  figure;
  h = axes();
  notBoxPlot(ratios, dpci);
  pos = get(h, 'XTick');
  pos = unique([0, pos]);
  set(h, 'YLim', [0 50], 'XLim', [-1 ids(end)+1], 'XTick', pos, 'XTickLabel', num2str(pos(:)));

  sigstar(groups, pvals);

  ylabel('Fraction of injury (%)')
  xlabel('dpci')

  title(title_name);

  return;
end

function [new] = interpolate(old)

  pts = unique(old(:,1));
  pos = [min(pts):max(pts)].';

  if (numel(pts)==numel(pos))
    new = old;
  else
    new = interp1(old(:,1), old(:,2), pos);
    new = [pos new];
  end

  return;
end

function folders = find_segmentations(fname)

  files = dir(fname);
  folders = {};

  for i=1:length(files)
    if (files(i).isdir && files(i).name(1)~='.')
      fdir = fullfile(fname, files(i).name);

      if (exist([fdir '.tif'], 'file') && exist([fdir '.zip'], 'file'))
        folders{end+1} = fdir;
      else
        subdir = find_segmentations(fdir);
        folders = [folders, subdir];
      end
    end
  end

  return;
end

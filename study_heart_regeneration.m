function [ratios] = study_heart_regeneration

  files = find_segmentations('/Users/blanchou/Documents/SB07/Histology/modified_data/');

  title_name = 'SB07';

  ratios = NaN(length(files), 1);
  dpci = ratios;

  hits = regexp(files, '.*SB\d+_\w(\d+)dpci_.*', 'tokens');

  for i=1:length(files)
    fname = files{i};

    [tmp, prefix1, junk] = fileparts(fname);
    [tmp, prefix2, junk] = fileparts(tmp);
    prefix = [prefix2 '_' prefix1 '_'];

    data = [fname '.zip'];
    imgs = [fname '.tif'];

    ROIs = ReadImageJROI(data);
    props = analyze_ROI(ROIs, 'Slice', 'Area');
    props = sortrows(props);

    is_injury = ([diff(props(:,1)); 1] == 0);
    heart = sum(props(~is_injury, 2));
    injury = sum(props(is_injury, 2));

    ratios(i) = injury/heart;
    dpci(i) = str2double(hits{i}{1});

    export_ROI(prefix, ROIs, imgs);
  end

  [H, p] = myttest(ratios, dpci);
  ratios = ratios*100;

  ids = unique(dpci);
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

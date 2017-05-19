function [ratios, files] = study_heart_regeneration(do_export)

  if (nargin == 0)
    do_export = false;
  end

  %titles = {'SB07', 'SB12', 'SB13', 'SB14'};
  %dirs = {'/Users/blanchou/Documents/SB07/Histology/modified_data/'};
  titles = {'SB07'};

  ratios = {};
  files = {};

  for i=1:length(titles)
    [ratios{end+1}, files{end+1}] = compute_regeneration(titles{i}, do_export);
  end

  return;
end

function [ratios, files] = compute_regeneration(title_name, do_export)

  fpath = ['/Users/blanchou/Documents/' title_name '/Histology/modified_data/'];
  if (~exist(fpath, 'dir'))
    fpath = ['/Users/blanchou/Documents/' title_name '/'];
  end

  files = find_segmentations(fpath);

  meta = fullfile(fpath, [title_name '.txt']);
  if (exist(meta, 'file'))
    fid = fopen(meta, 'r');
    params = textscan(fid, '%s %d %d %f', 'CommentStyle', '#');
    fclose(fid);

    params{2} = double(params{2});
    params{3} = double(params{3});
  else
    params = {''};
  end

  ratios = NaN(length(files), 1);
  volumes = NaN(length(files), 2);
  %ratios = NaN(length(files), 2);
  dpci = NaN(length(files), 1);

  hits = regexp(files, '.*SB\d+_\w(\d+)dpci_.*', 'tokens');

  for i=1:length(files)
    fname = files{i};

    [tmp, prefix1, junk] = fileparts(fname);
    [tmp, prefix2, junk] = fileparts(tmp);
    prefix = [prefix2 '_' prefix1];

    ref = strncmp(params{1}, prefix2, length(prefix));

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
  %  ratios(i, 2) = sum(injury(ismember(injury(:,1), largest),2))/sum(heart(ismember(heart(:,1), largest),2));

    if (any(ref))
      volumes(i,1) = sum(heart(:,2) * params{2}(ref) * params{3}(ref) * params{4}(ref).^2);
      volumes(i,2) = sum(injury(:,2) * params{2}(ref) * params{3}(ref) * params{4}(ref).^2);
    end

    dpci(i) = str2double(hits{i}{1});

    if (do_export)
      export_ROI(prefix, ROIs, imgs);
    end
  end

  ratios(isnan(ratios)) = 0;
  ratios = ratios*100;
  ratios(:,end+1) = dpci;

  o = find_outliers(ratios(:,1), ratios(:,end));

  h1 = display_ratios(ratios(~o,1), dpci(~o), title_name);
  %h2 = display_ratios(ratios(:,2), dpci, [title_name ' - 3 slices']);

  hold on;
  scatter(ratios(o,end), ratios(o,1), [], [0 0 0], 'filled');

  h2 = display_ratios(volumes(~o,1), dpci(~o), [title_name ' - heart']);
  h3 = display_ratios(volumes(~o,2), dpci(~o), [title_name ' - injury']);

  epath = fullfile(pwd, 'export');

  if (~exist(epath, 'dir'))
    mkdir(epath);
  end

  print(h1, '-dpdf', '-noui', '-bestfit', fullfile(epath, [title_name '_ratios.pdf']));
  %print(h2, '-dpdf', '-noui', '-bestfit', fullfile(epath, [title_name '_ratios-3_slices.pdf']));

  delete(h1);
  %delete(h2);

  return;
end

function outliers = find_outliers(vals, groups)

  labels = unique(groups);
  outliers = false(size(vals));

%{ 
C.f. boxplot
'whisker'       Maximum whisker length W.  Default is W=1.5.  Points
                      are drawn as outliers if they are larger than
                      Q3+W*(Q3-Q1) or smaller than Q1-W*(Q3-Q1), where Q1
                      and Q3 are the 25th and 75th percentiles, respectively.
                      The default value 1.5 corresponds to approximately +/-
                      2.7 sigma and 99.3 coverage if the data are normally
                      distributed.
%}
  for i=1:length(labels)
    curr = (groups==labels(i));
    y = vals(curr);
    pc = prctile(y, [25 75]);
    w = 1.5 * diff(pc);

    o = (y < pc(1)-w | y > pc(2)+w);
    outliers(curr) = o;
  end

  return;
end

function hfig = display_ratios(ratios, dpci, title_name)

  [H, p] = myttest(ratios, dpci);

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

  hfig = figure;
  h = axes();
  notBoxPlot(ratios, dpci, 'jitter', 2);
  pos = get(h, 'XTick');
  pos = unique([0, pos]);
  set(h, 'YLim', [0 max(50, ceil(max(ratios)/10)*10)], 'XLim', [-1 ids(end)+1], 'XTick', pos, 'XTickLabel', num2str(pos(:)));

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

      %if (exist([fdir '.tif'], 'file') && exist([fdir '.zip'], 'file'))
      if (exist([fdir '.zip'], 'file'))
        folders{end+1} = fdir;
      else
        subdir = find_segmentations(fdir);
        folders = [folders, subdir];
      end
    end
  end

  return;
end

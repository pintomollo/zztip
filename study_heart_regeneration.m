function [ratios, files] = study_heart_regeneration(do_export)

  if (nargin == 0)
    do_export = false;
  end

  titles = {'SB07_Xhellerii', 'SB12_Xmaculatus', 'SB13_Ptitteya', 'SB14_Dpentazona'};
  %dirs = {'/Users/blanchou/Documents/SB07/Histology/modified_data/'};
  %titles = {'SB07'};

  sizes = get_fish_sizes('SB16');
  groups = unique(sizes(:,end));
  [mvals, svals, ns] = mymean(sizes(:,1), 1, sizes(:,end));

%  figure;hold on;
%  bar(1:length(groups), mvals(:));
%  errorbar(1:length(groups), mvals(:), svals(:));

%  keyboard

  lengths = {};
  volumes = {};

  for i=1:length(titles)
    [lengths{end+1}, volumes{end+1}] = compute_properties(titles{i}, do_export);
  end
  lengths = cat(1, lengths{:});
  [groups, indx, jndx] = unique(lengths(:,2));

  for i=1:length(lengths)
    lengths{i,1}(:,end+1) = jndx(i);
  end
  tmp = lengths(:,1);
  vals = cat(1, tmp{:});
  [mvals, svals, ns] = mymean(vals(:,1:2), 1, vals(:,end));
%  keyboard

%  figure;hold on;
%  bar(1:2*length(groups), mvals(:));
%  errorbar(1:2*length(groups), mvals(:), svals(:));

  for i=1:length(volumes)
    volumes{i}(:,end+1) = i;
  end
  volumes = cat(1, volumes{:});
  [mvals, svals, ns] = mymean(volumes(:,1), 1, volumes(:,end));
%  keyboard

  ratios = {};
  files = {};

  for i=1:length(titles)
    [ratios{end+1}, files{end+1}] = compute_regeneration(titles{i}, do_export);
  end

  return;
end

function [sizes] = get_fish_sizes(title_name)

  fpath = ['/Users/blanchou/Documents/' title_name '/Brightfield/modified_data/'];
  if (~exist(fpath, 'dir'))
    fpath = ['/Users/blanchou/Documents/' title_name '/'];
  end

  files = find_segmentations(fpath);

  meta = fullfile(fpath, [title_name '.txt']);
  if (exist(meta, 'file'))
    fid = fopen(meta, 'r');
    params = textscan(fid, '%s %f %d %f', 'CommentStyle', '#');
    fclose(fid);

    params{2} = double(params{2});
    params{3} = double(params{3});
  else
    params = {''};
  end

  sizes = cell(length(files), 1);

  hits = regexp(files, '.*SB(\d+)_\w_.*', 'tokens');

  for i=1:length(files)
    fname = files{i};

    [tmp, prefix, junk] = fileparts(fname);

    ref = strncmp(params{1}, prefix, length(prefix));

    data = [fname '.zip'];
    imgs = [fname '.tif'];

    if (~isempty(hits{i}))
      ROIs = ReadImageJROI(data);
      props = analyze_ROI(ROIs, 'Length');

      props(:,end+1) = str2double(hits{i}{1});

      if (any(ref))
        props(:,1) = props(:,1) * params{4}(ref);
      end

      sizes{i} = props;
    end
  end

  sizes = sizes(~all(cellfun('isempty', sizes), 2), :);
  sizes = cat(1, sizes{:});

  return;
end

function [lengths, volumes] = compute_properties(title_name, do_export)

  fpath = ['/Users/blanchou/Documents/' title_name '/Histology/AFOG/modified_data/'];
  if (~exist(fpath, 'dir'))
    fpath = ['/Users/blanchou/Documents/' title_name '/'];
  end

  files = find_segmentations(fpath);

  meta = fullfile(fpath, [title_name '.txt']);
  if (exist(meta, 'file'))
    fid = fopen(meta, 'r');
    params = textscan(fid, '%s %f %d %f', 'CommentStyle', '#');
    fclose(fid);

    params{2} = double(params{2});
    params{3} = double(params{3});
  else
    params = {''};
  end

  lengths = cell(length(files), 2);
  volumes = NaN(length(files), 1);

  hits = regexp(files, '.*SB\d+_\w(\d+)dpci_.*', 'tokens');

  for i=1:length(files)
    fname = files{i};

    [tmp, prefix1, junk] = fileparts(fname);
    [tmp, prefix2, junk] = fileparts(tmp);
    prefix = [prefix2 '_' prefix1];

    ref = strncmp(params{1}, prefix2, length(prefix));

    data = [fname '.zip'];
    imgs = [fname '.tif'];

    if (isempty(hits{i}))
      specie = regexp(fname, '.*SB\d+_(\w)uninj_long.*', 'tokens');

      if (~isempty(specie))
        ROIs = ReadImageJROI(data);
        props = analyze_ROI(ROIs, 'Length');
        props = [props(1:2:end, 1), props(2:2:end,1)];

        if (any(ref))
          lengths(i,1) = {props * params{4}(ref)};
          lengths(i,2) = specie{1};
        end

        if (do_export)
          if (any(ref))
            export_ROI(prefix, ROIs, imgs, params{4}(ref));
          else
            export_ROI(prefix, ROIs, imgs);
          end
        end
      end

    elseif (str2double(hits{i}{1}) == 7)
      ROIs = ReadImageJROI(data);
      props = analyze_ROI(ROIs, 'Slice', 'Area');
      [props, indxs] = sortrows(props, [1 -2]);
      ROIs = ROIs(indxs);

      is_injury = ([1; diff(props(:,1))] == 0);
      heart = props(~is_injury, :);
      heart = interpolate(heart);

      if (any(ref))
        volumes(i,1) = sum(heart(:,2) * params{2}(ref) * params{3}(ref) * params{4}(ref).^2);
      end
    end
  end

  lengths = lengths(~all(cellfun('isempty', lengths), 2), :);
  volumes = volumes(~isnan(volumes));

  return;
end


function [ratios, files] = compute_regeneration(title_name, do_export)

  fpath = ['/Users/blanchou/Documents/' title_name '/Histology/AFOG/modified_data/'];
  if (~exist(fpath, 'dir'))
    fpath = ['/Users/blanchou/Documents/' title_name '/'];
  end
  csv_headers = {'Slice index', 'Ventricle area (mm2)', 'Injury area (mm2)'};

  files = find_segmentations(fpath);

  meta = fullfile(fpath, [title_name '.txt']);
  if (exist(meta, 'file'))
    fid = fopen(meta, 'r');
    params = textscan(fid, '%s %f %d %f', 'CommentStyle', '#');
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

    if (~isempty(hits{i}))
      dpci(i) = str2double(hits{i}{1});
    end

    if (do_export)
      data = heart;
      data(:,end+1) = 0;
      if (~isempty(injury))
        data(injury(1,1)-1+[1:size(injury,1)], end) = injury(:,end);
      end

      if (any(ref))
        data(:,2:end) = data(:,2:end) * params{4}(ref).^2;
        export_ROI(prefix, ROIs, imgs, params{4}(ref));
      else
        export_ROI(prefix, ROIs, imgs);
      end

      export_csv(fullfile(tmp, [prefix2 '.csv']), prefix1, data, csv_headers);
    end
  end

  valids = ~isnan(dpci);
  ratios = ratios(valids,:);
  volumes = volumes(valids,:);
  dpci = dpci(valids,:);

  [vals, junk, indxs] = unique(dpci);
  short = (vals(2:end) - vals(1:end-1) < 3);

  if (any(short))
    for i=1:length(short)
      if (short(i))
        vals(i+1) = vals(i);
      end
    end
    dpci = vals(indxs);
  end

  ratios(isnan(ratios)) = 0;
  ratios = ratios*100;
  ratios(:,end+1) = dpci;

  o = find_outliers(ratios(:,1), ratios(:,end));

  tname = strrep(title_name, '_', '\_');

  h1 = display_ratios(ratios(:,1), dpci, o, tname);
  %h2 = display_ratios(ratios(:,2), dpci, [title_name ' - 3 slices']);


  h2 = display_ratios(volumes(:,1), dpci, o, [tname ' - heart'], 'Volume (mm^3)');
  h3 = display_ratios(volumes(:,2), dpci, o, [tname ' - injury'], 'Volume (mm^3)');

  h4 = display_correlations(volumes, dpci, o, tname);

  epath = fullfile(pwd, 'export');

  if (~exist(epath, 'dir'))
    mkdir(epath);
  end

  print(h1, '-dpdf', '-noui', '-bestfit', fullfile(epath, [title_name '_ratios.pdf']));
  %print(h2, '-dpdf', '-noui', '-bestfit', fullfile(epath, [title_name '_ratios-3_slices.pdf']));
  print(h2, '-dpdf', '-noui', '-bestfit', fullfile(epath, [title_name '_heart.pdf']));
  print(h3, '-dpdf', '-noui', '-bestfit', fullfile(epath, [title_name '_injury.pdf']));
  print(h4, '-dpdf', '-noui', '-bestfit', fullfile(epath, [title_name '_correlation.pdf']));

  delete(h1);
  %delete(h2);
  delete(h2);
  delete(h3);
  delete(h4);

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

function hfig = display_ratios(ratios, dpci, outliers, title_name, y_label)

  if (nargin < 5)
    y_label = 'Fraction of injury (%)';
  end

  if (numel(outliers)~=numel(ratios))
    outliers = false(size(ratios));
  end

  or = ratios(outliers);
  od = dpci(outliers);

  ratios = ratios(~outliers);
  dpci = dpci(~outliers);

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

  if (ids(end) < 90)
    ids(end+1) = 90;
  end

  hfig = figure;
  h = axes();
  notBoxPlot(ratios, dpci, 'jitter', 2);
  hold on;
  scatter(od, or, [], [0 0 0], 'filled');
  pos = get(h, 'XTick');
  pos = unique([0, pos]);
  if (all(ratios<10))
    ylims = get(h, 'YLim');
    set(h, 'YLim', [0 ylims(2)],  'XLim', [-1 ids(end)+1], 'XTick', pos, 'XTickLabel', num2str(pos(:)));
  else
    set(h, 'YLim', [0 ceil(max(ratios)/10)*10], 'XLim', [-1 ids(end)+1], 'XTick', pos, 'XTickLabel', num2str(pos(:)));
  end

  sigstar(groups, pvals);

  ylabel(y_label)
  xlabel('dpci')

  title(title_name);

  return;
end

function hfig = display_correlations(values, dpci, outliers, title_name)

  if (numel(outliers)~=numel(values(:,1)))
    outliers = false(size(values,1),1);
  end

  or = values(outliers, :);
  od = dpci(outliers);

  values = values(~outliers, :);
  dpci = dpci(~outliers);

  ids = unique(dpci);
  colors = brewermap(max(92, max(ids)), 'RdYlBu');
  mvals = NaN(length(ids), 2);

  hfig = figure;
  h = axes();
  hold on;
  for i=1:length(ids)
    curr = (dpci==ids(i));
    scatter(values(curr,1), values(curr, 2), [], colors(ids(i),:), 'filled');
    [mvals(i,:), junk] = mymean(values(curr, :));
  end
  scatter(or(:,1), or(:,2), [], [0 0 0], 'filled');
  ylims = get(h, 'YLim');
  xlims = get(h, 'XLim');
  %axis equal
  set(h, 'YLim', [0 ylims(2)],  'XLim', [0 xlims(2)]);

  ylims(1) = 0;
  xlims(1) = 0;

  for i=1:length(ids)
    plot(mvals(i, [1 1]), ylims, 'Color', colors(ids(i),:));
    plot(xlims, mvals(i, [2 2]), 'Color', colors(ids(i),:));
    [mvals(i,:), junk] = mymean(values(curr, :));
  end

  caxis(h, [0 max(92, max(ids))]);
  colormap(h, [0 0 0; colors]);
  colorbar('peer', h);

  ylabel('Injury volume mm^3');
  xlabel('Heart volume mm^3');

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
    if (files(i).name(1)~='.')
      if (files(i).isdir)
        fdir = fullfile(fname, files(i).name);

        %if (exist([fdir '.tif'], 'file') && exist([fdir '.zip'], 'file'))
        if (exist([fdir '.zip'], 'file'))
          folders{end+1} = fdir;
        else
          subdir = find_segmentations(fdir);
          folders = [folders, subdir];
        end
      else
        [tmp, file, ext] = fileparts(files(i).name);
        fdir = fullfile(fname, file);

        if (strncmp(ext, '.zip', 4))
          if (exist([fdir '.tif'], 'file') && ~exist(fdir, 'dir'))
            folders{end+1} = fdir;
          end
        end
      end
    end
  end

  return;
end

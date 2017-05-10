function export_ROI(fname, ROIs, imgs, group, dpi)

  if (nargin < 4)
    group = true;
    dpi = 150;
  elseif (nargin < 5)
    if (islogical(group))
      dpi = 150;
    else
      dpi = group;
      group = true;
    end
  end

  nmax = 1000;

  [fpath, fname, fext] = fileparts(fname);

  if (isempty(fpath))

    fpath = fullfile(pwd, 'export');

    if (~exist(fpath, 'dir'))
      mkdir(fpath);
    end
  end

  if (isempty(fname))
    fname = 'export_roi';
  end

  if (isempty(fext))
    fext = '.pdf';
  end

  if (ischar(imgs))
    nimgs = size_data(imgs);
    imgs = load_data(imgs, [1:nimgs]);
  else
    nimgs = size_data(imgs);
    nimgs = nimgs(end);
  end
  if (ischar(ROIs))
    ROIs = ReadImageJROI(ROIs);
  end

  is_rgb = (ndims(imgs)>3);
  fsize = get( 0, 'Screensize' );
  fsize = [fsize(1) 0 fsize(4)*8.27/11.7 fsize(4)];
  fsize(1) = (fsize(1)-fsize(3))/2;
  fsize = ceil(fsize);

  colors = [0.5 0.5 0.5; 1 0.5 0.5; brewermap(max(length(ROIs), 1), 'Spectral')];

  if (group)

    if (nimgs < 10)
      nj = ceil(sqrt(nimgs));
      ni = ceil(nimgs/nj);
    else
      ni = 5;
      nj = ceil(nimgs/ni);
    end

    hf = figure('position', fsize);
    for i=1:nimgs
      ci = rem(i-1, ni) + 1;
      cj = ceil(i/ni);

      ha = subplot(ni, nj, cj + (ci-1)*nj, 'Parent', hf, 'Visible', 'off');

      if (is_rgb)
        img = imgs(:,:,:,i);
      else
        img = imgs(:,:,i);
      end
      image(img, 'Parent', ha);
      set(ha, 'Visible', 'off', 'NextPlot', 'add', 'LooseInset', get(ha,'TightInset'));
      axis(ha, 'image');

      c = 1;
      for j=1:length(ROIs)
        if (ROIs{j}.nPosition == i)
          switch ROIs{i}.strType
            case 'PolyLine'
              plot(ha, ROIs{j}.mnCoordinates(:,1), ROIs{j}.mnCoordinates(:,2), 'Color', colors(c,:), 'LineWidth', 2);
            case 'Polygon'
              plot(ha, ROIs{j}.mnCoordinates([1:end 1],1), ROIs{j}.mnCoordinates([1:end 1],2), 'Color', colors(c,:), 'LineWidth', 2);
          end
          c = c + 1;
        end
      end
    end

    %print(hf, ['-d' fext(2:end)], ['-r' num2str(dpi)], '-noui', '-bestfit', fullfile(fpath, [fname num2str(count) fext]));
    print(hf, ['-d' fext(2:end)], ['-r' num2str(dpi)], '-noui', '-bestfit', fullfile(fpath, [fname fext]));
  else

    %{
    count = -1;
    for i=1:nmax
      if (~exist(fullfile(fpath, [fname num2str(i) fext]), 'file'))
        count = i;
        break;
      end
    end

    if (count < 0)
      error(['Cannot overwrite ' fname ]);
    end
    %}
    count = 1;

    hf = figure('position', fsize);
    ha = axes('Parent', hf, 'Visible', 'off');
    hi = -1;
    hl = [];
    for i=1:nimgs
      if (is_rgb)
        img = imgs(:,:,:,i);
      else
        img = imgs(:,:,i);
      end
      if (hi < 0)
        hi = image(img, 'Parent', ha);
        set(ha, 'Visible', 'off', 'NextPlot', 'add');
        axis(ha, 'image');
      else
        set(hi, 'CData', img);
      end

      if (any(ishandle(hl)))
        delete(hl)
        hl = [];
      end

      c = 1;
      for j=1:length(ROIs)
        if (ROIs{j}.nPosition == i)
          switch ROIs{i}.strType
            case 'PolyLine'
              hl(end+1) = plot(ha, ROIs{j}.mnCoordinates(:,1), ROIs{j}.mnCoordinates(:,2), 'Color', colors(c,:), 'LineWidth', 2);
            case 'Polygon'
              hl(end+1) = plot(ha, ROIs{j}.mnCoordinates([1:end 1],1), ROIs{j}.mnCoordinates([1:end 1],2), 'Color', colors(c,:), 'LineWidth', 2);
          end
          c = c + 1;
        end
      end

      print(hf, ['-d' fext(2:end)], ['-r' num2str(dpi)], '-noui', '-bestfit', fullfile(fpath, [fname num2str(i+count-1) fext]));
    end
  end

  delete(hf);

  return;
end

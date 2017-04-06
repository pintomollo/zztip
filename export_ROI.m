function export_ROI(fname, ROIs, imgs)

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

  count = -1;
  for i=1:nmax
    if (~exist(fullfile(fpath, [fname num2str(i) fext]), 'file'))
      count = i;
    end
  end

  if (count < 0)
    error(['Cannot overwrite ' fname ]);
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

  hf = figure;
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
    else
      set(hi, 'CData', img);
    end

    if (any(ishandle(hl)))
      delete(hl)
      hl = [];
    end

    for j=1:length(ROIs)
      if (ROIs{j}.nPosition == i)
        switch ROIs{i}.strType
          case 'PolyLine'
            hl(end+1) = plot(ha, ROIs{j}.mnCoordinates(:,1), ROIs{j}.mnCoordinates(:,2));
          case 'Polygon'
            hl(end+1) = plot(ha, ROIs{j}.mnCoordinates([1:end 1],1), ROIs{j}.mnCoordinates([1:end 1],2));
        end
      end
    end

    saveas(hf, fullfile(fpath, [fname num2str(i+count-1) fext]));
  end

  delete(hf);

  return;
end

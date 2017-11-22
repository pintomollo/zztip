function export_shh(folder)

  dpi = 150;
  fext = '.pdf';

  if (nargin==0 || isempty(folder))
    folder = uigetdir();

    if (~ischar(folder))
      return;
    end
  end

  folder = absolutepath(folder);

  if (~exist(folder, 'dir'))
    error([folder ' is not a valid directory.']);
  end

  fi = (1 - 0.1)/6;
  fj = 0.5;

  fsize = get( 0, 'Screensize' );
  fsize = [fsize(1) 0 fsize(4)*8.27/11.7 fsize(4)];
  fsize(1) = (fsize(1)-fsize(3))/2;
  fsize = ceil(fsize);
  hf = figure('position', fsize);
  colormap(gray);
  haxes = NaN(6,2);

  for i=1:6
    hi = subplot('position', [0 (6-i)*fi fj fi], 'Parent', hf, 'Visible', 'off');
    haxes(i,1) = image([], 'Parent', hi, 'CDataMapping', 'scaled');
    set(hi, 'Visible', 'off', 'LooseInset', get(hi,'TightInset'));
    axis(hi, 'image');

    hg = subplot('position', [0.5 (6-i)*fi fj fi], 'Parent', hf, 'Visible', 'off');
    haxes(i,2) = image([], 'Parent', hg, 'CDataMapping', 'scaled');
    set(hg, 'Visible', 'off', 'LooseInset', get(hg,'TightInset'));
    axis(hg, 'image');
  end
  ha = subplot('position', [0 0.9 1 0.01], 'Parent', hf, 'Visible', 'off');
  ht = title(ha, '', 'Visible', 'on', 'FontSize', 24);

  files = dir(folder);
  for i=1:length(files)
    if (files(i).isdir && files(i).name(1)~='.')
      curdir = files(i).name;
      curpath = fullfile(folder, curdir);
      stacks = dir(fullfile(curpath, '*_GFP.tif'));

      for j=1:length(stacks)
        name = stacks(j).name(1:end-8);
        imgs = fullfile(curpath, [name '.tif']);
        gfps = fullfile(curpath, stacks(j).name);

        if (exist(imgs, 'file') == 2)
          nframes = size_data(imgs);
          imgs = load_data(imgs, 1:nframes);
          gfps = load_data(gfps, 1:nframes);

          for k=1:nframes
            ci = rem(k-1, 6) + 1;

            set(haxes(ci,1), 'CData', imgs(:,:,k))
            set(haxes(ci,2), 'CData', gfps(:,:,k))

            if (ci == 6)
              indx = num2str(ceil(k/6));
              set(ht, 'String', [curdir ': ' name ' (' indx ')']);

              print(hf, ['-d' fext(2:end)], ['-r' num2str(dpi)], '-noui', '-bestfit', fullfile(curpath, [name '_' indx fext]));
            end
          end

          if (ci ~= 6)
            for k=ci+1:6
              set(haxes(k,1), 'CData', [])
              set(haxes(k,2), 'CData', [])
            end

            indx = num2str(ceil(nframes/6));
            set(ht, 'String', [curdir ': ' name ' (' indx ')']);

            print(hf, ['-d' fext(2:end)], ['-r' num2str(dpi)], '-noui', '-bestfit', fullfile(curpath, [name '_' indx fext]));
          end
        end
      end
    end
  end

  delete(hf);

  return;
end

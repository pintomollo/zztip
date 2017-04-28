function resize_all(folder)

  folder = absolutepath(folder);

  if (~exist(folder, 'dir'))
    error([folder ' is not a valid directory.']);
  end

  recursive_resize(folder);

  return;
end

function recursive_resize(folder)

  files = dir(folder);
  has_dir = false;

  for i=1:length(files)
    if (files(i).isdir && files(i).name(1)~='.')
      has_dir = true;

      if (~strncmp(files(i).name, '_resized', 8))
        recursive_resize(fullfile(folder, files(i).name));
      end
    end
  end

  if (~has_dir)
    disp(['Resizing ' folder]);
    resize_tif(fullfile(folder, '*.tif'));
  end

  return;
end

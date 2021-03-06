function sections2stack(folder)

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

  folders = recursive_resize(folder);

  for i=1:length(folders)
    fldr = folders{i};

    movefile(fullfile(fldr, ['_resized' filesep '*.tif']), fldr, 'f');

    rmdir(fullfile(fldr, '_resized'));
  end

  create_stacks(folder);

  return;
end

function folders = recursive_resize(folder)

  files = dir(folder);
  has_dir = false;
  folders = {};

  for i=1:length(files)
    if (files(i).isdir && files(i).name(1)~='.')
      has_dir = true;

      if (~exist(fullfile(folder, [files(i).name '.tif']), 'file') && ~strncmp(files(i).name, '_resized', 8))
        new_folders = recursive_resize(fullfile(folder, files(i).name));
        folders = [folders; new_folders];
      end
    end
  end

  if (~has_dir)
    disp(['Resizing ' folder]);
    resize_tif(fullfile(folder, '*.tif'));
    folders{end+1} = folder;
  end

  return;
end

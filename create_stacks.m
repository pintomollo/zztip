function create_stacks(folder)

  if (iscell(folder))

    for i=1:length(folder)
      disp(['Stacking ' folder{i}]);
      stack_folder(folder{i});
    end
  else
    folder = absolutepath(folder);

    if (~exist(folder, 'dir'))
      error([folder ' is not a valid directory.']);
    end

    recursive_stacks(folder);
  end

  return;
end

function recursive_stacks(folder)

  files = dir(folder);
  has_dir = false;

  for i=1:length(files)
    if (files(i).isdir && files(i).name(1)~='.')
      has_dir = true;

      if (~exist(fullfile(folder, [files(i).name '.tif']), 'file'))
        recursive_stacks(fullfile(folder, files(i).name));
      end
    end
  end

  if (~has_dir)
    disp(['Stacking ' folder]);
    stack_folder(folder);
  end

  return;
end

function stack_folder(folder)

  [cur_path, fname, ext] = fileparts(folder);
  files = dir(fullfile(folder, '*.tif'));

  values = struct2cell(files);
  files = values(ismember(fieldnames(files), 'name'), :);
  files = fullfile(folder, files(~strncmp(files, '.', 1)));

  sname = fullfile(cur_path, [fname '.tif']);

  save_data(sname, files);

  return;
end

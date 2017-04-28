function create_stacks(folder)

  folder = absolutepath(folder);

  if (~exist(folder, 'dir'))
    error([folder ' is not a valid directory.']);
  end

  recursive_stacks(folder);

  return;
end

function recursive_stacks(folder)

  files = dir(folder);
  has_dir = false;

  for i=1:length(files)
    if (files(i).isdir && files(i).name(1)~='.')
      has_dir = true;

      recursive_stacks(fullfile(folder, files(i).name));
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
  files = fullfile(folder, values(ismember(fieldnames(files), 'name'), :));

  sname = fullfile(cur_path, [fname '.tif']);

  save_data(sname, files);

  return;
end

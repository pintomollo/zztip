function [ratios] = study_heart_regeneration

  files = {'/Users/blanchou/Documents/SB07/Histology/modified_data/SB07_X7dpci_I_AFOG/heart1', ...
           '/Users/blanchou/Documents/SB07/Histology/modified_data/SB07_X7dpci_I_AFOG/heart2', ...
           '/Users/blanchou/Documents/SB07/Histology/modified_data/SB07_X30dpci_I_AFOG/heart1', ...
           '/Users/blanchou/Documents/SB07/Histology/modified_data/SB07_X30dpci_I_AFOG/heart2', ...
           '/Users/blanchou/Documents/SB07/Histology/modified_data/SB07_X30dpci_I_AFOG/heart3'};
  ratios = NaN(length(files), 1);

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

    export_ROI(prefix, ROIs, imgs);
  end

  return;
end

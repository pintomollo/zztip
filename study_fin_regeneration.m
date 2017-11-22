function study_fin_regeneration

  %files = {'/Users/blanchou/Documents/SB07/Histology/modified_data/SB07_X7dpci_I_AFOG/heart1', ...
  %         '/Users/blanchou/Documents/SB07/Histology/modified_data/SB07_X7dpci_I_AFOG/heart2', ...
  %         '/Users/blanchou/Documents/SB07/Histology/modified_data/SB07_X30dpci_I_AFOG/heart1', ...
  %         '/Users/blanchou/Documents/SB07/Histology/modified_data/SB07_X30dpci_I_AFOG/heart2', ...
  %         '/Users/blanchou/Documents/SB07/Histology/modified_data/SB07_X30dpci_I_AFOG/heart3'};

  %files = {'/Users/blanchou/Documents/Data/Regen/SB04'};
  files = {'/Users/blanchou/Documents/SB01/Brightfield/SB01_Leica_amputated_AB368/modified_data/SB01', ...
           '/Users/blanchou/Documents/SB04/Brightfield/SB04_Leica_amputated_AB371/modified_data/SB04', ...
           '/Users/blanchou/Documents/SB05/Brightfield/SB05_Leica_amputated_calcein_SA120814/modified_data/SB05'};

  %{
  files = {'/Users/blanchou/Documents/Data/Regen/SB01', ...
           '/Users/blanchou/Documents/Data/Regen/SB04', ...
           '/Users/blanchou/Documents/Data/Regen/SB05'};
  %}

  resolution = 5.64;
  nouter = 7;

  all_fish = {};
  all_exper = {};
  all_types = NaN(0, 1);
  all_rays = NaN(0, 2*nouter);
  all_lens = NaN(0, 2*nouter);
  all_fins = NaN(0, 4);

  for i=1:length(files)
    fname = files{i};

    [tmp, prefix1, junk] = fileparts(fname);
    [tmp, prefix2, junk] = fileparts(tmp);
    prefix = [prefix2 '_' prefix1 '_'];

    data = [fname '.zip'];
    imgs = [fname '.tif'];
    meta = [fname '.txt'];

    disp(meta)

    fid = fopen(meta, 'r');
    frames = textscan(fid, '%s %s %d %s', 'CommentStyle', '#');
    fclose(fid);

    ROIs = ReadImageJROI(data);
    %[props, coords] = analyze_ROI(ROIs, 'Slice', 'Area', 'Perimeter');
    [fins, rays] = parse_fin_ROI(ROIs, resolution);
    fins(end+1:length(frames{1}), :) = NaN;
    rays(end+1:length(frames{1})) = {[]};

    nfins = length(rays);
    laterals = NaN(nfins, 2*nouter);
    lengths = NaN(nfins, 2*nouter);
    for j=1:nfins
      if (size(rays{j}, 2) > 2*nouter)
        tmp = [rays{j}(:,1:nouter) rays{j}(:,end-nouter+1:end)];
        tmp(1,~tmp(3,:)) = NaN;
        laterals(j,:) = tmp(1,:);
        lengths(j,:) = tmp(2,:);
      end
    end

    types = [cumsum(isnan(laterals(:,1:nouter)), 2) cumsum(isnan(laterals(:, end:-1:end-nouter+1)), 2)];
    probs = (types(:,nouter) ~= types(:,end));
    laterals(probs, :) = NaN;
    lengths(probs, :) = NaN;

    all_fish = [all_fish; frames{2}];
    all_exper = [all_exper; frames{4}];
    all_types = [all_types; types(:,end)];
    all_rays = [all_rays; laterals];
    all_lens = [all_lens; lengths];
    all_fins = [all_fins; fins];
  end

  whole = strncmp(all_exper, '-1', 2);
  fulls = (all_types==0 & ~whole);
  amput = (all_types==3);
  tinys = (all_types==5);

  valids = (whole | fulls | amput | tinys);
  %figure;boxplot(all_rays(whole, :).', all_types(valids).');

  ratios = [all_lens(:,2:nouter) ./ all_lens(:, 1:nouter-1); all_lens(:,end-1:-1:end-nouter+1) ./ all_lens(:, end:-1:end-nouter+2)];
  figure;boxplot(ratios([valids & whole; valids & whole],:))
  keyboard

  return;
end

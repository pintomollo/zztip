function study_fin_regeneration

  %files = {'/Users/blanchou/Documents/SB07/Histology/modified_data/SB07_X7dpci_I_AFOG/heart1', ...
  %         '/Users/blanchou/Documents/SB07/Histology/modified_data/SB07_X7dpci_I_AFOG/heart2', ...
  %         '/Users/blanchou/Documents/SB07/Histology/modified_data/SB07_X30dpci_I_AFOG/heart1', ...
  %         '/Users/blanchou/Documents/SB07/Histology/modified_data/SB07_X30dpci_I_AFOG/heart2', ...
  %         '/Users/blanchou/Documents/SB07/Histology/modified_data/SB07_X30dpci_I_AFOG/heart3'};

  %files = {'/Users/blanchou/Documents/Data/Regen/SB04'};
  %{
  files = {'/Users/blanchou/Documents/SB01/Brightfield/SB01_Leica_amputated_AB368/modified_data/SB01', ...
           '/Users/blanchou/Documents/SB04/Brightfield/SB04_Leica_amputated_AB371/modified_data/SB04', ...
           '/Users/blanchou/Documents/SB05/Brightfield/SB05_Leica_amputated_calcein_SA120814/modified_data/SB05'};
  %}
  
  files = {'/home/blanchou/Documents/Manuscripts/Fin_regen/Regen/SB01', ...
           '/home/blanchou/Documents/Manuscripts/Fin_regen/Regen/SB04', ...
           '/home/blanchou/Documents/Manuscripts/Fin_regen/Regen/SB05', ...
           '/home/blanchou/Documents/Manuscripts/Fin_regen/Regen/AJ02'};
  

  resolution = 5.64;
  nouter = 7;

  all_fish = {};
  all_exper = {};
  all_types = NaN(0, 1);
  all_rays = NaN(0, 2*nouter);
  all_lens = NaN(0, 2*nouter);
  all_segs = NaN(0, 2*nouter);
  all_tips = NaN(0, 2*nouter, 2);
  all_fins = NaN(0, 4);

  for i=1:length(files)
    fname = files{i};

    [tmp, prefix1, junk] = fileparts(fname);
    [tmp, prefix2, junk] = fileparts(tmp);
    prefix = [prefix2 '_' prefix1 '_'];

    data = [fname '.zip'];
    fimg = [fname '.tif'];
    meta = [fname '.txt'];

    disp(meta)

    fid = fopen(meta, 'r');
    frames = textscan(fid, '%s %s %d %s', 'CommentStyle', '#');
    fclose(fid);

    ROIs = ReadImageJROI(data);
    %[props, coords] = analyze_ROI(ROIs, 'Slice', 'Area', 'Perimeter');
    [fins, rays, paths] = parse_fin_ROI(ROIs, resolution);

    fins(end+1:length(frames{1}), :) = NaN;
    rays(end+1:length(frames{1})) = {[]};
    paths(end+1:length(frames{1})) = {[]};
    periods = cell(size(rays));

    size_rays = cellfun(@size, rays, 'UniformOutput', false);
    size_rays = cat(1, size_rays{:});

    fishes = unique(frames{2});
    for j=1:length(fishes)
      hits = strncmp(frames{2}, fishes{j}, 4) & (size_rays(:,2) > 1);
      nrays = [size_rays(hits, 2) frames{3}(hits)];
      for k=2:size(nrays,1)
        if (nrays(k,1) ~= nrays(1,1))
          disp(['#' num2str(nrays(1,2)) ' (' num2str(nrays(1,1)) ') -> #' num2str(nrays(k,2)) ' (' num2str(nrays(k,1)) ')' ])
        end
      end
    end
    %keyboard

    %indxs = find(~isnan(fins(:,1)) & size_rays(:,2)>nouter);
    %ndata = length(indxs);

    %imgs = load_data(fimg, indxs);
    %paths = paths(indxs);

    %periods(indxs) = get_segment_length(imgs, paths, resolution);

    nfins = length(rays);
    laterals = NaN(nfins, 2*nouter);
    lengths = NaN(nfins, 2*nouter);
    segments = NaN(nfins, 2*nouter);
    tips = NaN(nfins, 2*nouter, 2);
    for j=1:nfins
      if (size(rays{j}, 2) > 2*nouter)
        tmp = [rays{j}(:,1:nouter) rays{j}(:,end-nouter+1:end)];
        tmp(1:2,~tmp(3,:)) = NaN;
        laterals(j,:) = tmp(1,:);
        lengths(j,:) = tmp(2,:);
        
        tmp = [cellfun(@(x)(x(end,1)), paths{j}); cellfun(@(x)(x(end,2)), paths{j})];
        tips(j,:,1) = [tmp(1, 1:nouter) tmp(1, end-nouter+1:end)];
        tips(j,:,2) = [tmp(2, 1:nouter) tmp(2, end-nouter+1:end)];
        %segments(j,:) = [periods{j}(1:nouter).' periods{j}(end-nouter+1:end).'];
      end
    end

    types = [cumsum(isnan(laterals(:,1:nouter+1)), 2) cumsum(isnan(laterals(:, end:-1:end-nouter+1)), 2)];
    probs = (types(:,nouter) ~= types(:,end) & types(:,end)<5);
    laterals(probs, :) = NaN;
    lengths(probs, :) = NaN;
    types(probs,end) = NaN;
    tips(cat(3,isnan(laterals),isnan(laterals))) = NaN;
    %segments(probs, :) = NaN;

    all_fish = [all_fish; frames{2}];
    all_exper = [all_exper; frames{4}];
    all_types = [all_types; types(:,end)];
    all_rays = [all_rays; laterals];
    all_lens = [all_lens; lengths];
    %all_segs = [all_segs; segments];
    all_fins = [all_fins; fins];
    all_tips = [all_tips; tips];
  end

  whole = strncmp(all_exper, '-1', 2) & (all_types==0);
  regro = strncmp(all_exper, '4', 1) & (all_types==0);
  fulls = (all_types==0 & ~whole & ~regro);
  amput = (all_types==3);
  tinys = (all_types==5 | all_types==6);

  length(intersect(unique(all_fish(whole)),unique(all_fish(fulls))))
  length(intersect(unique(all_fish(whole)),unique(all_fish(amput))))
  length(intersect(intersect(unique(all_fish(whole)),unique(all_fish(amput))), unique(all_fish(fulls))))

  valids = (whole | regro | fulls | amput | tinys);

  whole_len = all_lens(whole,:);
  regen_len = all_lens(fulls,:);
  thinn_len = all_lens(amput,:);
  slimm_len = all_lens(tinys,:);

  whole_fin = all_fins(whole,:);
  regen_fin = all_fins(fulls,:);
  thinn_fin = all_fins(amput,:);
  slimm_fin = all_fins(tinys,:);

  whole_bif = all_rays(whole,:);
  regen_bif = all_rays(fulls,:);
  thinn_bif = all_rays(amput,:);
  slimm_bif = all_rays(tinys,:);

  figure;
  subplot(3,3,1);
  boxplot({whole_fin(:, 2), regen_fin(:, 2), thinn_fin(:, 2), slimm_fin(:, 2)});
  [H11,p11] = myttest([whole_fin(:, 2); regen_fin(:, 2); thinn_fin(:, 2); slimm_fin(:, 2)], [ones(size(whole_fin(:, 2))); ones(size(regen_fin(:, 2)))*2; ones(size(thinn_fin(:, 2)))*3; ones(size(slimm_fin(:, 2)))*4]);
  title('Fin length (um)')
  subplot(3,3,2);
  boxplot({whole_fin(:, 2)-whole_fin(:,3), regen_fin(:, 2)-regen_fin(:,3), thinn_fin(:, 2)-thinn_fin(:,3), slimm_fin(:, 2)-slimm_fin(:,3)});
  [H12,p12] = myttest([whole_fin(:, 2)-whole_fin(:,3); regen_fin(:, 2)-regen_fin(:,3); thinn_fin(:, 2)-thinn_fin(:,3); slimm_fin(:, 2)-slimm_fin(:,3)], [ones(size(whole_fin(:, 2))); ones(size(regen_fin(:, 2)))*2; ones(size(thinn_fin(:, 2)))*3; ones(size(slimm_fin(:, 2)))*4]);
  title('Cleft depth (um)')
  subplot(3,3,3);
  boxplot({whole_fin(:, 4), regen_fin(:, 4), thinn_fin(:, 4), slimm_fin(:, 4)});
  [H13,p13] = myttest([whole_fin(:, 4); regen_fin(:, 4); thinn_fin(:, 4); slimm_fin(:, 4)], [ones(size(whole_fin(:, 2))); ones(size(regen_fin(:, 2)))*2; ones(size(thinn_fin(:, 2)))*3; ones(size(slimm_fin(:, 2)))*4]);
  title('Fin width (um)')

  nfins = size(all_fins, 1);
  rel_stats = NaN(nfins, 4);
  rel_bif = NaN(nfins, 2*nouter);
  for i=1:nfins
    if (valids(i))
      ref = find(strncmp(all_fish, all_fish{i}, 4) & whole);
      [len, nindx]= max(all_lens(i,:));
      [outers] = find(isfinite(all_lens(i,:)));
      if (length(outers)==0)
        keyboard
      end

      rel_bif(i,:) = all_rays(i,:) ./ all_rays(ref,:);

      width = sqrt(sum((all_tips(i,outers(1),:) - all_tips(i,outers(end),:)).^2));

      cleft = len - all_fins(i,3);

      rel_stats(i, 1) = whole(i) + 2*fulls(i) + 3*amput(i) + 4*tinys(i);

      if (ref == i)
        rel_stats(i,2) = len/all_fins(i,2);
        rel_stats(i,3) = cleft/(all_fins(i,2)-all_fins(i,3));
        rel_stats(i,4) = width/all_fins(i,4);
      else
        rwidth = sqrt(sum((all_tips(ref,outers(1),:) - all_tips(ref,outers(end),:)).^2));

        rel_stats(i,2) = len/all_lens(ref,nindx);
        rel_stats(i,3) = cleft/(all_lens(ref,nindx)-all_fins(ref,3));
        rel_stats(i,4) = width/rwidth;
      end
    end
  end

  subplot(3,3,4);
  boxplot(rel_stats(:,2), rel_stats(:,1));
  [H21,p21] = myttest([whole_fin(:, 2); regen_fin(:, 2); thinn_fin(:, 2); slimm_fin(:, 2)], [ones(size(whole_fin(:, 2))); ones(size(regen_fin(:, 2)))*2; ones(size(thinn_fin(:, 2)))*3; ones(size(slimm_fin(:, 2)))*4]);
  title('Relative fin length')
  subplot(3,3,5);
  boxplot(rel_stats(:,3), rel_stats(:,1));
  [H22,p22] = myttest([whole_fin(:, 2)-whole_fin(:,3); regen_fin(:, 2)-regen_fin(:,3); thinn_fin(:, 2)-thinn_fin(:,3); slimm_fin(:, 2)-slimm_fin(:,3)], [ones(size(whole_fin(:, 2))); ones(size(regen_fin(:, 2)))*2; ones(size(thinn_fin(:, 2)))*3; ones(size(slimm_fin(:, 2)))*4]);
  title('Relative cleft depth')
  subplot(3,3,6);
  boxplot(rel_stats(:,4), rel_stats(:,1));
  [H23,p23] = myttest([whole_fin(:, 4); regen_fin(:, 4); thinn_fin(:, 4); slimm_fin(:, 4)], [ones(size(whole_fin(:, 2))); ones(size(regen_fin(:, 2)))*2; ones(size(thinn_fin(:, 2)))*3; ones(size(slimm_fin(:, 2)))*4]);
  title('Relative fin width')

  subplot(3,3,7);
  boxplot({whole_bif(:), regen_bif(:), thinn_bif(:), slimm_bif(:)});
  [H31,p31] = myttest([whole_bif(:); regen_bif(:); thinn_bif(:); slimm_bif(:)], [ones(size(whole_bif(:))); ones(size(regen_bif(:)))*2; ones(size(thinn_bif(:)))*3; ones(size(slimm_bif(:)))*4]);
  title('Bifurcation position (um)')
  subplot(3,3,8);
  boxplot({whole_bif(:)./whole_len(:), regen_bif(:)./regen_len(:), thinn_bif(:)./thinn_len(:), slimm_bif(:)./slimm_len(:)});
  [H32,p32] = myttest([whole_bif(:)./whole_len(:); regen_bif(:)./regen_len(:); thinn_bif(:)./thinn_len(:); slimm_bif(:)./slimm_len(:)], [ones(size(whole_bif(:))); ones(size(regen_bif(:)))*2; ones(size(thinn_bif(:)))*3; ones(size(slimm_bif(:)))*4]);
  title('Relative bifurcation position')
  subplot(3,3,9);
  boxplot(rel_bif(:), repmat(rel_stats(:,1), [2*nouter, 1]));
  [H33,p33] = myttest(rel_bif(:), repmat(rel_stats(:,1), [2*nouter, 1]));
  title('Relative bifurcation length')

  indxs = [1:2*nouter]*5;
  full_indxs = [repmat(indxs-4, size(whole_bif,1), 1); ...
                repmat(indxs-3, size(regen_bif,1), 1); ...
                repmat(indxs-2, size(thinn_bif,1), 1); ...
                repmat(indxs-1, size(slimm_bif,1), 1)];
  fulls_bif = [whole_bif; regen_bif; thinn_bif; slimm_bif];
  fulls_len = [whole_len; regen_len; thinn_len; slimm_len];

  figure;
  subplot(3,1,1);
  boxplot(fulls_len(:), full_indxs(:));
  title('Ray length per ray (um)')

  subplot(3,1,2);
  boxplot(fulls_bif(:), full_indxs(:));
  title('Bifurcation position per ray (um)')

  whole_nbi = (whole_bif==whole_len);
  regen_nbi = (regen_bif==regen_len);
  thinn_nbi = (thinn_bif==thinn_len);
  slimm_nbi = (slimm_bif==slimm_len);

  subplot(3,1,3);
  boxplot(fulls_bif(~[whole_nbi; regen_nbi; thinn_nbi; slimm_nbi]), full_indxs(~[whole_nbi; regen_nbi; thinn_nbi; slimm_nbi]));
  title('Bifurcation position per ray - w/o unbifurcated (um)')

  fracs = [mean(whole_nbi); ...
           mean(regen_nbi); ...
           mean(thinn_nbi); ...
           mean(slimm_nbi)];
  figure;bar(fracs.');
  title('Fraction of unbifurcated rays')
  legend('unamputated', 'full', 'thin (-3)', 'slim (-5)')
  keyboard



  %whole_seg = all_segs(whole,:);
  %regen_seg = all_segs(fulls,:);


  %whole_nums = whole_bif ./ whole_seg;
  %regen_nums = regen_bif ./ regen_seg;

  whole_no_bif = (whole_ratio==1);
  regen_no_bif = (regen_ratio==1);

  whole_bif(whole_no_bif) = NaN;
  regen_bif(regen_no_bif) = NaN;
  whole_ratio(whole_no_bif) = NaN;
  regen_ratio(regen_no_bif) = NaN;
  %whole_nums(whole_no_bif) = NaN;
  %regen_nums(regen_no_bif) = NaN;

  keyboard

  figure;
  h1=subplot(1,2,1);notBoxPlot(whole_bif);
  h2=subplot(1,2,2);notBoxPlot(regen_bif);
  linkaxes([h1 h2]);
  ylim([0 9000])

  figure;
  h1=subplot(1,2,1);notBoxPlot(whole_ratio);
  h2=subplot(1,2,2);notBoxPlot(regen_ratio);
  linkaxes([h1 h2]);
  ylim([0 1])

  figure;
  h1=subplot(1,2,1);notBoxPlot(whole_seg);
  h2=subplot(1,2,2);notBoxPlot(regen_seg);
  linkaxes([h1 h2]);
  ylim([0 500])

  %figure;
  %h1=subplot(1,2,1);notBoxPlot(whole_nums);
  %h2=subplot(1,2,2);notBoxPlot(regen_nums);
  %linkaxes([h1 h2]);

  figure;
  subplot(1,2,1);notBoxPlot([whole_bif(:); regen_bif(:)], [ones(numel(whole_bif), 1); ones(numel(regen_bif(:)),1)*2]);
  ylim([0 9000])
  subplot(1,2,2);notBoxPlot([whole_ratio(:); regen_ratio(:)], [ones(numel(whole_ratio), 1); ones(numel(regen_ratio(:)),1)*2]);
  ylim([0 1])
  %figure;boxplot(all_rays(whole, :).', all_types(valids).');

  %ratios = [all_lens(:,2:nouter) ./ all_lens(:, 1:nouter-1); all_lens(:,end-1:-1:end-nouter+1) ./ all_lens(:, end:-1:end-nouter+2)];
  %figure;boxplot(ratios([valids & whole; valids & whole],:))
  %keyboard

  return;
end

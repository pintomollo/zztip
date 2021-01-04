% STUDY_FIN_REGENERATION produces the plots used to analyze fin regeneration
% in the context of the biomechanical study in collaboration with UNIZH
%
% Blanchoud Group, UNIFR
% Simon Blanchoud
% 09/06/2020
function study_fin_regeneration

  tensions = '/home/blanchou/Documents/Manuscripts/Fin_regen/criterion/stress_tension_tip_3.csv';
  T = csvread(tensions);

  c1 = reshape(T(end-8,:), 11, []);
  c2 = reshape(T(end-6,:), 11, []);
  c3 = reshape(T(end-2,:), 11, []);
  c4 = reshape(T(end,:), 11, []);

  c3b = reshape(T(end-3,:), 11, []);
  c4b = reshape(T(end-1,:), 11, []);
  c5 = c3b - c4b;

  e = T([1:2 24:25], :);
  l = T([3:10 16:23], :);
  c = T([11:15], :);

  %figure;
  %subplot(2,2,1);scatter(abs(mean(T(1:25,:))), T(end-9,:));
  %subplot(2,2,2);scatter(abs(mean(l))-abs(mean(c)), T(end-7,:));
  %subplot(2,2,3);scatter(abs(mean([l;c])), T(end-3,:));
  %subplot(2,2,4);scatter(abs(mean([l;c]))-abs(mean(e)), T(end-1,:));

  mc = abs(mean([l; c]));
  dc = abs(mean([l; c])) - abs(mean(e));
  dl = (abs(mean(l)) - abs(mean(e)))/2;
  dc = (abs(mean(l)) - abs(mean(c)))/2;
  me = abs(mean(e));
  c1b = reshape(abs(mean([l;c;e])), 11, []);
  c1b = c1b ./ c1b(:, [1 1 1 1]);
  c2b = reshape(dc, 11, []);
  c2b = c2b ./ c2b(:, [2 2 2 2]);
  c3b = reshape(mc, 11, []);
  c3b = c3b ./ c3b(:, [3 3 3 3]);
  c4b = reshape(dc, 11, []);
  c4b = c4b ./ c4b(:, [2 2 2 2]);
  c4c = reshape(dl, 11, []);
  c4c = c4c ./ c4c(:, [2 2 2 2]);
  c4d = reshape(dl, 11, []);
  c4d = c4d ./ c4d(:, [1 1 1 1]);

  %figure;
  %subplot(2,2,1);scatter(c3b(:), T(end-2,:));
  %subplot(2,2,2);scatter(c4b(:), T(end,:));

  %keyboard

  %figure;
  %subplot(2,2,1);bar(c1(:,1:3))
  %subplot(2,2,2);bar(c2b(:,3:4))
  %subplot(2,2,3);bar(c3(:,2:4))
  %subplot(2,2,4);bar(c4c(:,2:4))

  %figure;
  %subplot(3,2,1);bar(c4b(:,1:4))
  %subplot(3,2,2);bar(c4b(:,3:4))
  %subplot(3,2,3);bar(c4c(:,1:4))
  %subplot(3,2,4);bar(c4c(:,3:4))
  %subplot(3,2,5);bar(c4d(:,1:4))
  %subplot(3,2,6);bar(c4d(:,3:4))

  figure;
  subplot(2,2,1);bar(c1(:,1:3));ylim([0 2])
  subplot(2,2,2);bar(c2(:,3:4));ylim([-5 5])
  subplot(2,2,3);bar(c3(:,2:4));ylim([0 2])
  subplot(2,2,4);bar(c4(:,3:4));ylim([-5 5])

  figure;
  subplot(2,2,1);bar(c1b(:,1:3));ylim([0 2])
  subplot(2,2,3);bar(c2(:,3:4));ylim([-5 5])
  subplot(2,2,2);bar(c3(:,3:4));ylim([0 2])
  subplot(2,2,4);bar(c4d(:,2:4));ylim([-2 2])

  figure;
  subplot(2,2,1);bar(c1b(:,1:3));
  subplot(2,2,3);bar(c2(:,3:4));
  subplot(2,2,2);bar(c3(:,3:4));
  %subplot(2,2,3);bar(c4c(:,2:4));
  subplot(2,2,4);bar(c4d(:,2:4));

  %{
  figure;
  subplot(2,2,1);bar(c1(:,[1 3])-1)
  subplot(2,2,2);bar(c2(:,3:4))
  subplot(2,2,3);bar(c3(:,[2 4])-c3(:,[3 3]))
  subplot(2,2,4);bar(c4(:,3:4))

  figure;
  subplot(2,2,1);bar(c1(:,[1 3])-1);ylim([-1 1])
  subplot(2,2,2);bar(c2(:,3:4));ylim([-5 5]);
  subplot(2,2,3);bar(c3(:,[2 4])-c3(:,[3 3]));ylim([-1 1])
  subplot(2,2,4);bar(c4(:,3:4));ylim([-5 5])
  %}

  %{
  C = [-1.5 -1.0 -1.0 -1.0 -1.0 0; ...
	-0.5 0.0 0.0 0.0 -0.5 1; ...
 	1.0 0.5 2.0 1.5 1.5 1; ...
	0.5 0.5 0 -2 -1.5 0.5; ...
	-0.1 0.0 -0.5 0.0 -0.5 -0.3; ...
	0.01 0.15 0 0.2 0.2 0 ; ...
	1.5 1.0 1.5 3.0 0.5 1.0; ...
	2.5 2.5 4.0 -0.5 -2.0 1.0; ...
	-0.9 0.4 -1.7 1.5 -1.7 -0.4; ...
	-1.05 0 -0.45 0 -0.05 -0.15; ...
	0.1 0.2 0.1 0 0.2 0];

  B = [9.75 9.0 12.0 10 11 10.5 ; ...
	1.0 -1.0 -2.0 -1.75 -0.5 -2; ...
	11.0 11.25 14.0 12.25 13.75 12.5; ...
        0.25 -0.375 0.5 -1.0 -0.375 -0.375; ...
	1.25 1.5 2.75 1.6 2.15 1.75; ...
	0.3 0.325 0.5 0.5 0.6 0.3; ...
	12.75 11.0 10.75 12.75 14.25 12.5; ...
	-0.6 1.875 0.5 0.125 1 0.75; ...
	8.95 10.1 12.65 11.95 12.35 12; ...
	1.175 1.25 1.975 1.4 1.725 1.925; ...
	0.25 0.5 0.65 0.55 0.5 0.35];

  O = [7.5 6.5 9.5 8 8.5 8; ...
	  1 0.5 -1.5 -1.5 -1 0; ...
	  11.5 11 14.05 11.5 12.5 10.5; ...
	  -0.5 0 1 0 -2.5 -0.5;...
	  1 1 1.5 1.3 1.5 1.2; ...
	  0.3 0.5 0.7 0.6 0.6 0.3; ...
	  12.5 10 11 14 13 12; ...
	  -1.5 3.5 1.5 2.5 -3 -1; ...
	  7.5 8.7 9.5 10 10 13; ...
	  1.1 1.1 1.3 0.9 1.5 1.2; ...
	  0.25 0.7 0.9 0.4 0.5 0.4];

  E = [0.1 0.1 1 1 0.2 0.2 1 2 0.2 0.3 0.3].';

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

  A3 = 4* 2*(58 - [77 100 96 97 104 83 64 64 28 28 13 13 18 18 49 49 70 64 61 58 58 44 44 44 44] + 1.4) / 116;
  A2 = 6* 2*(72 - [93 104 104 97 75 75 46 34 34 57 56 56 59 59 105 117 117 112 118 118 85 39 39 23 48] + 1.4) / 145;
  A1 = 6* 2*(72 - [70 97 120 96 89 60 54 53 51 46 56 49 49 45 80 91 91 86 66 54 41 27 16 52 51] + 1.4) / 145;

  indxs = cumsum(ones(size(A1)));
  mirror = [1:13 12:-1:1];

  A3v = cumsum(A3) ./ indxs;
  A3d = cumsum(A3([end:-1:1])) ./ indxs;

  A2v = cumsum(A2) ./ indxs;
  A2d = cumsum(A2([end:-1:1])) ./ indxs;

  A1v = cumsum(A1) ./ indxs;
  A1d = cumsum(A1([end:-1:1])) ./ indxs;

  A3t = A3v(mirror) - A3d(mirror);
  A2t = A2v(mirror) - A2d(mirror);
  A1t = A1v(mirror) - A1d(mirror);

  C = [C; mean([A3t(1:6); A2t(1:6); A1t(1:6)].') - mean([A3t(7:13); A2t(7:13); A1t(7:13)].') NaN(1,3); ...
	  10 0 -5 -15 -20 -10];

  B = [B; mean([A3t; A2t; A1t].') NaN(1,3); ...
  	-10 -2 12 25 25 8];

  O = [O; A3t(1) A2t(1) A1t(1) NaN(1,3);
	  -2 -1 2 2 -1 -1];

  E = [E; 0.1; 2];

  R = B - O;

  %figure;hold on
  %scatter(1:numel(A3), A3t, 'b')
  %scatter(1:numel(A3), A2t, 'k')
  %scatter(1:numel(A3), A1t, 'r')

  S = sign(C(:,1));
  S2 = sign(B(:,1)-B(:,2));
  S3 = sign(R(:,3));

  figure;subplot(2,2,1);
  bar([S2.*(B(:,1)-B(:,2)) S2.*(B(:,2)-B(:,3))])
  title('Slower growth in longer fins')
  subplot(2,2,2);
  bar(S.*(C(:,4)-C(:,3)))
  title('Shallower cleft in narrower fins')
  subplot(2,2,3);
  bar([S3(:,[1 1 1]).*R(:,1:3)]);
  hold on;scatter(1:numel(E), E, 'r');
  title('Bifurcation in A2 or longer');
  subplot(2,2,4);
  bar(S3.*(R(:,4)-R(:,3)));
  title('Distalisation in narrower fins')
  %}

  % Absolute path to the used segmentation files
  files = {'/home/blanchou/Documents/Manuscripts/Fin_regen/revisions/Regen/SB01', ...
           '/home/blanchou/Documents/Manuscripts/Fin_regen/revisions/Regen/SB04', ...
           '/home/blanchou/Documents/Manuscripts/Fin_regen/revisions/Regen/SB05', ...
           '/home/blanchou/Documents/Manuscripts/Fin_regen/revisions/Regen/AJ02', ...
           '/home/blanchou/Documents/Manuscripts/Fin_regen/revisions/Regen/DP03', ...
           '/home/blanchou/Documents/Manuscripts/Fin_regen/revisions/Regen/DP50'};

  % David bifurcation data
  doublets = [0 0 NaN NaN 0 0 NaN NaN 0 0 NaN NaN 0 0; ...
              0 0 NaN NaN 0 0 NaN NaN 0 0 NaN NaN 0 0; ...
              0 0 NaN NaN 0 0 NaN NaN 0 0 NaN NaN 0 0; ...
              0 0 NaN NaN 0 0 NaN NaN 0 0 NaN NaN 0 0; ...
              0 0 NaN NaN 0 0 NaN NaN 0 0 NaN NaN 0 0];

  quadruplets = [0 1 1 0 NaN NaN NaN NaN NaN NaN 0 1 1 0; ...
                 0 1 1 0 NaN NaN NaN NaN NaN NaN 0 1 1 0; ...
                 0 1 1 0 NaN NaN NaN NaN NaN NaN 0 1 1 0; ...
                 0 0 1 1 NaN NaN NaN NaN NaN NaN 0 0 1 0; ...
                 0 0 1 0 NaN NaN NaN NaN NaN NaN 0 1 0 0; ...
                 0 0 0 0 NaN NaN NaN NaN NaN NaN 0 1 0 0; ...
                 0 0 0 0 NaN NaN NaN NaN NaN NaN 0 1 1 0; ...
                 0 0 0 0 NaN NaN NaN NaN NaN NaN 0 0 0 0];

  nquad = size(quadruplets, 1);
  ndoubl = size(doublets, 1);

  % Global variables
  resolution = 5.64;
  nouter = 7;
  bif_thresh = 10;
  plot_segmentation = false;
  nindx = 11;

  % Plot an example of segmentation if requested
  if (plot_segmentation)
    data = [files{1} '.zip'];
    fimg = [files{1} '.tif'];
    ROIs = ReadImageJROI(data);
    [props, data] = analyze_ROI(ROIs, 1, 'Slice', 'Perimeter');
    goods = (props(:,1)==nindx);

    img = load_data(fimg, nindx);
    p = data(goods);

    figure;imshow(img); hold on;
    for i=1:length(p)
      plot(p{i}(:,1), p{i}(:,2))
    end
    keyboard
  end
    
  % Intermediate variables used to store all the data extracted from the segmentations
  all_fish = {};
  all_exper = {};
  all_types = NaN(0, 1);
  all_bifs = NaN(0, 2*nouter);
  all_lens = NaN(0, 2*nouter);
  all_fins = NaN(0, 5);
  all_wids = NaN(0, nouter);
  all_clef = NaN(0, nouter);

  % We loop through all files to gather everything
  for i=1:length(files)

    % The current file name
    fname = files{i};

    % The actual extensions used by the analysis
    data = [fname '.zip'];
    fimg = [fname '.tif'];
    meta = [fname '.txt'];

    % Mention the processing
    disp(meta)

    % Read the custom metadata file
    fid = fopen(meta, 'r');
    frames = textscan(fid, '%s %s %d %s', 'CommentStyle', '#');
    fclose(fid);

    % Read the segmentation file
    ROIs = ReadImageJROI(data);
    [fins, rays, paths] = parse_fin_ROI(ROIs, resolution);

    % The number of segmented fins
    nfins = max(length(rays), length(frames{1}));

    % Allocate the required space for the extracted data
    rays(end+1:nfins) = {[]};
    fins(end+1:nfins, :) = NaN;
    paths(end+1:nfins) = {[]};
    lengths = NaN(nfins, 2*nouter);
    tips = NaN(nfins, 2*nouter, 2);
    bifurcations = NaN(nfins, 2*nouter);
    widths = NaN(nfins, nouter);
    clefts = NaN(nfins, nouter);

    % Extract the number of rays in each fin
    size_rays = cellfun(@size, rays, 'UniformOutput', false);
    size_rays = cat(1, size_rays{:});

    % Make sure that segmentations are consistent by comparing the number of rays per fish
    fishes = unique(frames{2});
    for j=1:length(fishes)
      hits = strncmp(frames{2}, fishes{j}, 4) & (size_rays(:,2) > 1);
      nrays = [size_rays(hits, 2) frames{3}(hits)];

      % If we find the same fish with different number of rays, mention it!
      for k=2:size(nrays,1)
        if (nrays(k,1) ~= nrays(1,1))
          disp(['Variations in ray numbers between timepoints: #' num2str(nrays(1,2)) ' (' num2str(nrays(1,1)) ') -> #' num2str(nrays(k,2)) ' (' num2str(nrays(k,1)) ')' ])
        end
      end
    end

    % Loop through all fins to count the number of amputated rays
    for j=1:nfins

      % Process only the fins with enough rays
      if (size(rays{j}, 2) >= 2*nouter)

        % Keep only the NOUTER outer most rays for analysis
        tmpr = [rays{j}(:,1:nouter) rays{j}(:,end-nouter+1:end)];

        % Extract the position of the tip of each ray
        tmpt = [cellfun(@(x)(x(end,1)), paths{j}); cellfun(@(x)(x(end,2)), paths{j})];

        % Check if the amputation is symmetric
        is_symmetric = sum(tmpr(3,1:nouter) - tmpr(3,end:-1:end-nouter+1));

        % If not, try to fix it
        switch is_symmetric
          case 0
            % Nothing to do regarding the ray lengths

            % Store the position of the NOUTER tips
            tips(j,:,1) = [tmpt(1, 1:nouter) tmpt(1, end-nouter+1:end)];
            tips(j,:,2) = [tmpt(2, 1:nouter) tmpt(2, end-nouter+1:end)];
          case -1
            tmpr = [rays{j}(:,2:nouter+1) rays{j}(:,end-nouter+1:end)];

            tips(j,:,1) = [tmpt(1, 2:nouter+1) tmpt(1, end-nouter+1:end)];
            tips(j,:,2) = [tmpt(2, 2:nouter+1) tmpt(2, end-nouter+1:end)];
          case 1
            tmpr = [rays{j}(:,1:nouter) rays{j}(:,end-nouter:end-1)];

            tips(j,:,1) = [tmpt(1, 1:nouter) tmpt(1, end-nouter:end-1)];
            tips(j,:,2) = [tmpt(2, 1:nouter) tmpt(2, end-nouter:end-1)];
          otherwise
            % Assuming nothing can be done anyway
            tips(j,:,1) = [tmpt(1, 1:nouter) tmpt(1, end-nouter+1:end)];
            tips(j,:,2) = [tmpt(2, 1:nouter) tmpt(2, end-nouter+1:end)];

            disp(['Non-symmetric fin #' num2str(j) ' in file ' files{i}])
        end

        % Discard the amputated rays
        tmpr(1:2,~tmpr(3,:)) = NaN;

        % Extract the bifurcation and the length
        bifurcations(j,:) = tmpr(1,:);
        lengths(j,:) = tmpr(2,:);
        
        % And the corresponding ray-to-ray width
        widths(j,:) = sqrt(sum((tips(j,1:nouter,:) - tips(j,end:-1:end-nouter+1,:)).^2, 3));
        widths(j, isnan(bifurcations(j, 1:nouter))) = NaN;

        % Compute the relative cleft length, the depth will require to subtract the shortest ray
        clefts(j, :) = max(max(triu(lengths(ones(1,nouter)*j,1:nouter)), triu(lengths(ones(1,nouter)*j,end:-1:end-nouter+1))), [],2).';
        clefts(j, isnan(bifurcations(j, 1:nouter))) = NaN;
      end
    end

    % Count the number of amputated rays
    types = [cumsum(isnan(bifurcations(:,1:nouter+1)), 2) cumsum(isnan(bifurcations(:, end:-1:end-nouter+1)), 2)];

    % Discard uneven amputations
    probs = (types(:,nouter) ~= types(:,end) & types(:,end)<5);
    types(probs,end) = NaN;

    % Store all the valid information into our intermediate variables
    all_fish = [all_fish; frames{2}];
    all_exper = [all_exper; frames{4}];
    all_types = [all_types; types(:,end)];
    all_bifs = [all_bifs; bifurcations];
    all_lens = [all_lens; lengths];
    all_fins = [all_fins; fins];
    all_wids = [all_wids; widths];
    all_clef = [all_clef; clefts];
  end

  % Define the 5 types of experiments performed
  whole = strncmp(all_exper, '-1', 2) & (all_types==0);
  bowti = strncmp(all_exper, '4', 1) & (all_types==0);
  dblet = strncmp(all_exper, '5', 1);
  qdlet = strncmp(all_exper, '6', 1);
  treat = strncmp(all_exper, '7', 1);
  regen = (all_types==0 & ~whole & ~bowti & ~treat);
  thinn = (all_types==3 & ~dblet & ~qdlet);
  slimm = ((all_types==5 | all_types==6) & ~dblet & ~qdlet);
  valids = (whole | bowti | regen | thinn | slimm | treat);

  % Get the number of elements in each group
  nwhole = sum(whole);
  nregen = sum(regen);
  nbowti = sum(bowti);
  ntreat = sum(treat);
  nthinn = sum(thinn);
  nslimm = sum(slimm);
  ndblet = sum(dblet);
  nqdlet = sum(qdlet);

  % Prepare the index vector
  indx = [ones(nwhole, 1); ones(nregen, 1)*2; ones(nbowti, 1)*3; ones(ntreat, 1)*4; ones(nthinn, 1)*5; ones(nslimm, 1)*6; ones(ndblet, 1)*7; ones(nqdlet, 1)*8];
  full_indx = [ones(2*nouter*nwhole, 1); ones(2*nouter*nregen, 1)*2; ones(2*nouter*nbowti, 1)*3; ones(2*nouter*ntreat, 1)*4; ones(2*nouter*nthinn, 1)*5; ones(2*nouter*nslimm, 1)*6; ones(2*nouter*ndblet, 1)*7; ones(2*nouter*nqdlet, 1)*8];

  % Split the data per experiment type
  whole_fin = all_fins(whole,:);
  regen_fin = all_fins(regen,:);
  bowti_fin = all_fins(bowti,:);
  treat_fin = all_fins(treat,:);
  thinn_fin = all_fins(thinn,:);
  slimm_fin = all_fins(slimm,:);
  dblet_fin = all_fins(dblet,:);
  qdlet_fin = all_fins(qdlet,:);

  whole_wid = all_wids(whole,:);
  regen_wid = all_wids(regen,:);
  bowti_wid = all_wids(bowti,:);
  treat_wid = all_wids(treat,:);
  thinn_wid = all_wids(thinn,:);
  slimm_wid = all_wids(slimm,:);
  dblet_wid = all_wids(dblet,:);
  qdlet_wid = all_wids(qdlet,:);

  whole_len = all_lens(whole,:);
  regen_len = all_lens(regen,:);
  bowti_len = all_lens(bowti,:);
  treat_len = all_lens(treat,:);
  thinn_len = all_lens(thinn,:);
  slimm_len = all_lens(slimm,:);
  dblet_len = all_lens(dblet,:);
  qdlet_len = all_lens(qdlet,:);

  whole_bif = all_bifs(whole,:);
  regen_bif = all_bifs(regen,:);
  bowti_bif = all_bifs(bowti,:);
  treat_bif = all_bifs(treat,:);
  thinn_bif = all_bifs(thinn,:);
  slimm_bif = all_bifs(slimm,:);
  dblet_bif = all_bifs(dblet,:);
  qdlet_bif = all_bifs(qdlet,:);

  % Figure for the publication
  figure;
  subplot(2,2,1);

  % Figure 1: 
  % We start by analyzing the fin's length

  % Prepare the ray-to-ray averages for relative normalization
  whole_lengths = bsxfun(@rdivide, all_clef(whole,:), whole_fin(:,5));
  avg_lengths = nanmean(whole_lengths, 1);

  % Normalize the data per fish
  whole_norm = whole_fin(:,2) ./ whole_fin(:,5);
  bowti_norm = bowti_fin(:,2) ./ bowti_fin(:,5);
  regen_norm = regen_fin(:,2) ./ regen_fin(:,5);
  treat_norm = treat_fin(:,2) ./ treat_fin(:,5);
  thinn_norm = thinn_fin(:,2) ./ thinn_fin(:,5);
  slimm_norm = slimm_fin(:,2) ./ slimm_fin(:,5);
  qdlet_norm = qdlet_fin(:,2) ./ qdlet_fin(:,5);
  dblet_norm = dblet_fin(:,2) ./ dblet_fin(:,5);

  whole_rel = whole_norm ./ avg_lengths(1);
  bowti_rel = bowti_norm ./ avg_lengths(1);
  regen_rel = regen_norm ./ avg_lengths(1);
  treat_rel = treat_norm ./ avg_lengths(1);
  thinn_rel = thinn_norm ./ avg_lengths(4);
  slimm_rel = slimm_norm ./ avg_lengths(6);
  qdlet_rel = qdlet_norm ./ avg_lengths(1);
  dblet_rel = dblet_norm ./ avg_lengths(1);

  % Plot the two ways to look at the fin's length
  %figure;
  %subplot(1,2,1);
  %boxplot({whole_norm, regen_norm, bowti_norm, thinn_norm, slimm_norm});
  %title('Fin length / peduncle width (a.u.)')
  %subplot(1,2,2);
  %boxplot({whole_rel, regen_rel, bowti_rel, treat_rel, thinn_rel, slimm_rel, qdlet_rel, dblet_rel});
  [mval,sval] = mymean([whole_rel; regen_rel; bowti_rel; treat_rel; thinn_rel; slimm_rel; qdlet_rel; dblet_rel], 1, indx);
  errorbar([1:length(mval)], mval, sval)
  title('Relative Fin length / peduncle width (a.u.)')
  ylim([0.5 1.5])

  % And compute the significance
  %[H11,p11] = myttest([whole_norm; regen_norm; bowti_norm; thinn_norm; slimm_norm; qdlet_norm; dblet_norm], indx);
  [H12,p12] = myttest([whole_rel; regen_rel; bowti_rel; treat_rel; thinn_rel; slimm_rel; qdlet_rel; dblet_rel], indx);

  % Figure 2:
  % Next we look at the depth of the cleft
  whole_clefts = whole_lengths - whole_fin(:,3)./whole_fin(:,5);
  avg_clefts = nanmean(whole_clefts, 1);

  % Same normalization as above
  whole_cle = (whole_fin(:,2) - whole_fin(:,3)) ./ whole_fin(:,5);
  bowti_cle = (bowti_fin(:,2) - bowti_fin(:,3)) ./ bowti_fin(:,5);
  regen_cle = (regen_fin(:,2) - regen_fin(:,3)) ./ regen_fin(:,5);
  treat_cle = (treat_fin(:,2) - treat_fin(:,3)) ./ treat_fin(:,5);
  thinn_cle = (thinn_fin(:,2) - thinn_fin(:,3)) ./ thinn_fin(:,5);
  slimm_cle = (slimm_fin(:,2) - slimm_fin(:,3)) ./ slimm_fin(:,5);
  qdlet_cle = (qdlet_fin(:,2) - qdlet_fin(:,3)) ./ qdlet_fin(:,5);
  dblet_cle = (dblet_fin(:,2) - dblet_fin(:,3)) ./ dblet_fin(:,5);

  % And the plotting
  %figure;
  %subplot(1,2,1);
  %boxplot({whole_cle, regen_cle, bowti_cle, thinn_cle, slimm_cle});
  %title('Cleft depth / peduncle width (a.u.)')
  %subplot(1,2,2);

  subplot(2,2,2);
  %boxplot({whole_cle / avg_clefts(1), regen_cle / avg_clefts(1), bowti_cle / avg_clefts(1), treat_cle / avg_clefts(1), thinn_cle / avg_clefts(4), slimm_cle / avg_clefts(6), qdlet_cle / avg_clefts(1), dblet_cle / avg_clefts(1)});
  [mval,sval] = mymean([whole_cle / avg_clefts(1); regen_cle / avg_clefts(1); bowti_cle / avg_clefts(1); treat_cle / avg_clefts(1); thinn_cle / avg_clefts(4); slimm_cle / avg_clefts(6); qdlet_cle / avg_clefts(1); dblet_cle / avg_clefts(1)], 1, indx);
  errorbar([1:length(mval)], mval, sval)
  title('Relative Cleft depth / peduncle width (a.u.)')
  ylim([0.0 1.5])

  % Stats
  %[H21,p21] = myttest([whole_cle; regen_cle; bowti_cle; thinn_cle; slimm_cle;qdlet_cle; dblet_cle], indx);
  [H22,p22] = myttest([whole_cle / avg_clefts(1); regen_cle / avg_clefts(1); bowti_cle / avg_clefts(1); treat_cle / avg_clefts(1); thinn_cle / avg_clefts(4); slimm_cle / avg_clefts(6); qdlet_cle / avg_clefts(1); dblet_cle / avg_clefts(1)], indx);

  % Figure 3:
  % The width of the fin
  avg_widths = nanmean(bsxfun(@rdivide, whole_wid, whole_fin(:,5)), 1);

  % Same normalization as above
  whole_wid = nanmax(whole_wid, [], 2) ./ whole_fin(:,5);
  bowti_wid = nanmax(bowti_wid, [], 2) ./ bowti_fin(:,5);
  regen_wid = nanmax(regen_wid, [], 2) ./ regen_fin(:,5);
  treat_wid = nanmax(treat_wid, [], 2) ./ treat_fin(:,5);
  thinn_wid = nanmax(thinn_wid, [], 2) ./ thinn_fin(:,5);
  slimm_wid = nanmax(slimm_wid, [], 2) ./ slimm_fin(:,5);
  qdlet_wid = nanmax(qdlet_wid, [], 2) ./ qdlet_fin(:,5);
  dblet_wid = nanmax(dblet_wid, [], 2) ./ dblet_fin(:,5);

  % And the plotting
  %figure;
  %subplot(1,2,1);
  %boxplot({whole_wid, regen_wid, bowti_wid, thinn_wid, slimm_wid});
  %title('Fin width / peduncle width (a.u.)')
  %subplot(1,2,2);

  subplot(2,2,3);
  %boxplot({whole_wid / avg_widths(1), regen_wid / avg_widths(1), bowti_wid / avg_widths(1), treat_wid / avg_widths(1), thinn_wid / avg_widths(4), slimm_wid / avg_widths(6), qdlet_wid / avg_widths(1), dblet_wid / avg_widths(1)});
  [mval,sval] = mymean([whole_wid / avg_widths(1); regen_wid / avg_widths(1); bowti_wid / avg_widths(1); treat_wid / avg_widths(1); thinn_wid / avg_widths(4); slimm_wid / avg_widths(6); qdlet_wid / avg_widths(1); dblet_wid / avg_widths(1)], 1, indx);
  errorbar([1:length(mval)], mval, sval)
  title('Relative fin width / peduncle width (a.u.)')
  ylim([0.0 1.5])

  % Stats
  %[H31,p31] = myttest([whole_wid; regen_wid; bowti_wid; thinn_wid; slimm_wid; qdlet_wid; dblet_wid], indx);
  [H32,p32] = myttest([whole_wid / avg_widths(1); regen_wid / avg_widths(1); bowti_wid / avg_widths(1); treat_wid / avg_widths(1); thinn_wid / avg_widths(4); slimm_wid / avg_widths(6); qdlet_wid / avg_widths(1); dblet_wid / avg_widths(1)], indx);

  % Figure 4:
  % The bifurcation position
  avg_bifs = nanmean(bsxfun(@rdivide, whole_bif, whole_fin(:,5)), 1);

  % Same normalization as above
  whole_norm = bsxfun(@rdivide, whole_bif, whole_fin(:,5));
  bowti_norm = bsxfun(@rdivide, bowti_bif, bowti_fin(:,5));
  regen_norm = bsxfun(@rdivide, regen_bif, regen_fin(:,5));
  treat_norm = bsxfun(@rdivide, treat_bif, treat_fin(:,5));
  thinn_norm = bsxfun(@rdivide, thinn_bif, thinn_fin(:,5));
  slimm_norm = bsxfun(@rdivide, slimm_bif, slimm_fin(:,5));
  qdlet_norm = bsxfun(@rdivide, qdlet_bif, qdlet_fin(:,5));
  dblet_norm = bsxfun(@rdivide, dblet_bif, dblet_fin(:,5));

  whole_rel = bsxfun(@rdivide, whole_norm, avg_bifs);
  bowti_rel = bsxfun(@rdivide, bowti_norm, avg_bifs);
  regen_rel = bsxfun(@rdivide, regen_norm, avg_bifs);
  treat_rel = bsxfun(@rdivide, treat_norm, avg_bifs);
  thinn_rel = bsxfun(@rdivide, thinn_norm, avg_bifs);
  slimm_rel = bsxfun(@rdivide, slimm_norm, avg_bifs);
  qdlet_rel = bsxfun(@rdivide, qdlet_norm, avg_bifs);
  dblet_rel = bsxfun(@rdivide, dblet_norm, avg_bifs);

  % And the plotting
  %figure;
  %subplot(2,2,1);
  %boxplot({whole_norm(:), regen_norm(:), bowti_norm(:), thinn_norm(:), slimm_norm(:)});
  %title('Bifurcations length / peduncle width (a.u.)')
  %subplot(2,2,2);
  %boxplot({nanmean(whole_norm, 2), nanmean(regen_norm, 2), nanmean(bowti_norm, 2), nanmean(thinn_norm, 2), nanmean(slimm_norm, 2)});
  %title('Average bifurcation length / peduncle width (a.u.)')
  %subplot(2,2,3);
  %boxplot({whole_rel(:), regen_rel(:), bowti_rel(:), thinn_rel(:), slimm_rel(:)});
  %title('Relative bifurcation length / peduncle width (a.u.)')
  %subplot(2,2,4);
  %boxplot({nanmean(whole_rel, 2), nanmean(regen_rel, 2), nanmean(bowti_rel, 2), nanmean(thinn_rel, 2), nanmean(slimm_rel, 2)});
  %title('Relative average bifurcation length / peduncle width (a.u.)')

  % Stats
  %[H41,p41] = myttest([whole_norm(:); regen_norm(:); bowti_norm(:); thinn_norm(:); slimm_norm(:); qdlet_norm(:); dblet_norm(:)], full_indx);
  %[H42,p42] = myttest([nanmean(whole_norm, 2); nanmean(regen_norm, 2); nanmean(bowti_norm, 2); nanmean(thinn_norm, 2); nanmean(slimm_norm, 2); nanmean(qdlet_norm, 2); nanmean(dblet_norm, 2)], indx);
  %[H43,p43] = myttest([whole_rel(:); regen_rel(:); bowti_rel(:); thinn_rel(:); slimm_rel(:); qdlet_rel(:); dblet_rel(:)], full_indx);
  %[H44,p44] = myttest([nanmean(whole_rel, 2); nanmean(regen_rel, 2); nanmean(bowti_rel, 2); nanmean(thinn_rel, 2); nanmean(slimm_rel, 2); nanmean(qdlet_rel, 2); nanmean(dblet_rel, 2)], indx);

  % Figure 5:
  % The bifurcation position relative to the ray length
  avg_bifs = nanmean(whole_bif ./ whole_len, 1);

  % Same normalization as above
  whole_norm = whole_bif ./ whole_len;
  regen_norm = regen_bif ./ regen_len;
  bowti_norm = bowti_bif ./ bowti_len;
  treat_norm = treat_bif ./ treat_len;
  thinn_norm = thinn_bif ./ thinn_len;
  slimm_norm = slimm_bif ./ slimm_len;
  qdlet_norm = qdlet_bif ./ qdlet_len;
  dblet_norm = dblet_bif ./ dblet_len;

  whole_rel = bsxfun(@rdivide, whole_norm, avg_bifs);
  regen_rel = bsxfun(@rdivide, regen_norm, avg_bifs);
  bowti_rel = bsxfun(@rdivide, bowti_norm, avg_bifs);
  treat_rel = bsxfun(@rdivide, treat_norm, avg_bifs);
  thinn_rel = bsxfun(@rdivide, thinn_norm, avg_bifs);
  slimm_rel = bsxfun(@rdivide, slimm_norm, avg_bifs);
  qdlet_rel = bsxfun(@rdivide, qdlet_norm, avg_bifs);
  dblet_rel = bsxfun(@rdivide, dblet_norm, avg_bifs);

  % And the plotting
  %figure;
  %subplot(2,2,1);
  %boxplot({whole_norm(:), regen_norm(:), bowti_norm(:), thinn_norm(:), slimm_norm(:)});
  %title('Bifurcations position / peduncle width (a.u.)')
  %subplot(2,2,2);

  subplot(2,2,4);
  %boxplot({nanmean(whole_norm, 2), nanmean(regen_norm, 2), nanmean(bowti_norm, 2), nanmean(treat_norm, 2), nanmean(thinn_norm, 2), nanmean(slimm_norm, 2), nanmean(qdlet_norm, 2), nanmean(dblet_norm, 2)});
  [mval,sval] = mymean([nanmean(whole_norm, 2); nanmean(regen_norm, 2); nanmean(bowti_norm, 2); nanmean(treat_norm, 2); nanmean(thinn_norm, 2); nanmean(slimm_norm, 2); nanmean(qdlet_norm, 2); nanmean(dblet_norm, 2)], 1, indx);
  errorbar([1:length(mval)], mval, sval)
  title('Average bifurcation position / peduncle width (a.u.)')
  ylim([0.5 1.1])
  %subplot(2,2,3);
  %boxplot({whole_rel(:), regen_rel(:), bowti_rel(:), thinn_rel(:), slimm_rel(:)});
  %title('Relative bifurcation position / peduncle width (a.u.)')
  %subplot(2,2,4);
  %boxplot({nanmean(whole_rel, 2), nanmean(regen_rel, 2), nanmean(bowti_rel, 2), nanmean(thinn_rel, 2), nanmean(slimm_rel, 2)});
  %title('Relative average bifurcation position / peduncle width (a.u.)')

  % Stats
  %[H51,p51] = myttest([whole_norm(:); regen_norm(:); bowti_norm(:); thinn_norm(:); slimm_norm(:)], full_indx);
  [H52,p52] = myttest([nanmean(whole_norm, 2); nanmean(regen_norm, 2); nanmean(bowti_norm, 2); nanmean(treat_norm, 2); nanmean(thinn_norm, 2); nanmean(slimm_norm, 2); nanmean(qdlet_norm, 2); nanmean(dblet_norm, 2)], indx);
  %[H53,p53] = myttest([whole_rel(:); regen_rel(:); bowti_rel(:); thinn_rel(:); slimm_rel(:)], full_indx);
  %[H54,p54] = myttest([nanmean(whole_rel, 2); nanmean(regen_rel, 2); nanmean(bowti_rel, 2); nanmean(thinn_rel, 2); nanmean(slimm_rel, 2)], indx);

  display([H12 NaN(8,1) H22 NaN(8,1) H32 NaN(8,1) H52])
  display([p12 NaN(8,1) p22 NaN(8,1) p32 NaN(8,1) p52])
  display([nwhole nregen nbowti ntreat nthinn nslimm nqdlet ndblet])

  keyboard

  % Figure 6: Individual ray length

  % Normalization per peduncle
  whole_norm = bsxfun(@rdivide, whole_len, whole_fin(:,5));
  regen_norm = bsxfun(@rdivide, regen_len, regen_fin(:,5));
  bowti_norm = bsxfun(@rdivide, bowti_len, bowti_fin(:,5));
  treat_norm = bsxfun(@rdivide, treat_len, treat_fin(:,5));
  thinn_norm = bsxfun(@rdivide, thinn_len, thinn_fin(:,5));
  slimm_norm = bsxfun(@rdivide, slimm_len, slimm_fin(:,5));
  qdlet_norm = bsxfun(@rdivide, qdlet_len, qdlet_fin(:,5));
  dblet_norm = bsxfun(@rdivide, dblet_len, dblet_fin(:,5));

  whole_mean = nanmean(whole_norm, 1);
  regen_mean = nanmean(regen_norm, 1);
  bowti_mean = nanmean(bowti_norm, 1);
  treat_mean = nanmean(treat_norm, 1);
  thinn_mean = nanmean(thinn_norm, 1);
  slimm_mean = nanmean(slimm_norm, 1);
  qdlet_mean = nanmean(qdlet_norm, 1);
  dblet_mean = nanmean(dblet_norm, 1);
  whole_std = nanstd(whole_norm, 1);
  regen_std = nanstd(regen_norm, 1);
  bowti_std = nanstd(bowti_norm, 1);
  treat_std = nanstd(treat_norm, 1);
  thinn_std = nanstd(thinn_norm, 1);
  slimm_std = nanstd(slimm_norm, 1);
  qdlet_std = nanstd(qdlet_norm, 1);
  dblet_std = nanstd(dblet_norm, 1);

  peduncle_indx = [1:2*nouter]/(2*nouter);
  figure;
  subplot(1,2,1);
  %plot(peduncle_indx, [whole_mean; regen_mean; bowti_mean; thinn_mean; slimm_mean]);
  plot(repmat(peduncle_indx, 1, 3), [whole_mean whole_mean+whole_std whole_mean-whole_std; regen_mean regen_mean+regen_std regen_mean-regen_std; bowti_mean bowti_mean+bowti_std bowti_mean-bowti_std; treat_mean treat_mean+treat_std treat_mean-treat_std; thinn_mean thinn_mean+thinn_std thinn_mean-thinn_std; slimm_mean slimm_mean+slimm_std slimm_mean-slimm_std; qdlet_mean qdlet_mean+qdlet_std qdlet_mean-qdlet_std; dblet_mean dblet_mean+dblet_std dblet_mean-dblet_std]);
  %errorbar(repmat(peduncle_indx, 5, 1).', [whole_mean; regen_mean; bowti_mean; thinn_mean; slimm_mean].', [whole_std; regen_std; bowti_std; thinn_std; slimm_std].');
  %errorbar(peduncle_indx, whole_mean, whole_std);
  %axis('equal')
  legend({'whole', 'flat', 'step', 'RA', 'narrow', 'trimm', '4-let', '2-let'})

  H61 = cell(1, 2*nouter);
  p61 = cell(1, 2*nouter);
  H62 = NaN(8, 2*nouter);
  p62 = NaN(8, 2*nouter);
  %part_indx = [ones(nwhole, 1); ones(nregen, 1)*2; ones(nbowti, 1)*3; ones(nthinn, 1)*4; ones(nslimm, 1)*5; ones(n)];
  for i=1:2*nouter
    [H61{i},p61{i}] = myttest([whole_norm(:, i); regen_norm(:, i); bowti_norm(:, i); treat_norm(:, i); thinn_norm(:, i); slimm_norm(:, i); qdlet_norm(:, i); dblet_norm(:, i)], indx);
    H62(:,i) = H61{i}(1,:).';
    p62(:,i) = p61{i}(1,:).';
  end

  whole_norm = [whole_norm(:,1:nouter); whole_norm(:, end:-1:end-nouter+1)];
  regen_norm = [regen_norm(:,1:nouter); regen_norm(:, end:-1:end-nouter+1)];
  bowti_norm = [bowti_norm(:,1:nouter); bowti_norm(:, end:-1:end-nouter+1)];
  treat_norm = [treat_norm(:,1:nouter); treat_norm(:, end:-1:end-nouter+1)];
  thinn_norm = [thinn_norm(:,1:nouter); thinn_norm(:, end:-1:end-nouter+1)];
  slimm_norm = [slimm_norm(:,1:nouter); slimm_norm(:, end:-1:end-nouter+1)];

  whole_mean = nanmean(whole_norm, 1);
  regen_mean = nanmean(regen_norm, 1);
  bowti_mean = nanmean(bowti_norm, 1);
  treat_mean = nanmean(treat_norm, 1);
  thinn_mean = nanmean(thinn_norm, 1);
  slimm_mean = nanmean(slimm_norm, 1);
  whole_std = nanstd(whole_norm, 1);
  regen_std = nanstd(regen_norm, 1);
  bowti_std = nanstd(bowti_norm, 1);
  treat_std = nanstd(treat_norm, 1);
  thinn_std = nanstd(thinn_norm, 1);
  slimm_std = nanstd(slimm_norm, 1);

  %subplot(1,2,2);
  %plot(repmat(peduncle_indx(1:nouter), 1, 3), [whole_mean whole_mean+whole_std whole_mean-whole_std]);
  %hold on;
  %plot(repmat(peduncle_indx(end:-1:end-nouter+1), 1, 3), [regen_mean regen_mean+regen_std regen_mean-regen_std; bowti_mean bowti_mean+bowti_std bowti_mean-bowti_std; thinn_mean thinn_mean+thinn_std thinn_mean-thinn_std; slimm_mean slimm_mean+slimm_std slimm_mean-slimm_std]);

  %Figure 7: fraction of unbifurcated rays
  whole_norm = (whole_len - whole_bif < bif_thresh)+0;
  regen_norm = (regen_len - regen_bif < bif_thresh)+0;
  bowti_norm = (bowti_len - bowti_bif < bif_thresh)+0;
  treat_norm = (treat_len - treat_bif < bif_thresh)+0;
  thinn_norm = (thinn_len - thinn_bif < bif_thresh)+0;
  slimm_norm = (slimm_len - slimm_bif < bif_thresh)+0;
  qdlet_norm = (qdlet_len - qdlet_bif < bif_thresh)+0;
  dblet_norm = (dblet_len - dblet_bif < bif_thresh)+0;
  %qdlet_norm = 1-quadruplets;
  %dblet_norm = 1-doublets;

  whole_norm(isnan(whole_len) | isnan(whole_bif)) = NaN;
  regen_norm(isnan(regen_len) | isnan(regen_bif)) = NaN;
  bowti_norm(isnan(bowti_len) | isnan(bowti_bif)) = NaN;
  treat_norm(isnan(treat_len) | isnan(treat_bif)) = NaN;
  thinn_norm(isnan(thinn_len) | isnan(thinn_bif)) = NaN;
  slimm_norm(isnan(slimm_len) | isnan(slimm_bif)) = NaN;
  qdlet_norm(isnan(qdlet_len) | isnan(qdlet_bif)) = NaN;
  dblet_norm(isnan(dblet_len) | isnan(dblet_bif)) = NaN;

  whole_mean = nanmean(whole_norm, 1);
  regen_mean = nanmean(regen_norm, 1);
  bowti_mean = nanmean(bowti_norm, 1);
  treat_mean = nanmean(treat_norm, 1);
  thinn_mean = nanmean(thinn_norm, 1);
  slimm_mean = nanmean(slimm_norm, 1);
  qdlet_mean = nanmean(qdlet_norm, 1);
  dblet_mean = nanmean(dblet_norm, 1);
  whole_std = nanstd(whole_norm, 1);
  regen_std = nanstd(regen_norm, 1);
  bowti_std = nanstd(bowti_norm, 1);
  treat_std = nanstd(treat_norm, 1);
  thinn_std = nanstd(thinn_norm, 1);
  slimm_std = nanstd(slimm_norm, 1);
  qdlet_std = nanstd(qdlet_norm, 1);
  dblet_std = nanstd(dblet_norm, 1);

  peduncle_indx = [1:2*nouter]/(2*nouter);
  %figure;
  %subplot(1,2,1);
  subplot(1,2,2);
  %plot(peduncle_indx, [whole_mean; regen_mean; bowti_mean; thinn_mean; slimm_mean]);
  plot(repmat(peduncle_indx, 1, 3), [whole_mean whole_mean+whole_std whole_mean-whole_std; regen_mean regen_mean+regen_std regen_mean-regen_std; bowti_mean bowti_mean+bowti_std bowti_mean-bowti_std; treat_mean treat_mean+treat_std treat_mean-treat_std; thinn_mean thinn_mean+thinn_std thinn_mean-thinn_std; slimm_mean slimm_mean+slimm_std slimm_mean-slimm_std; qdlet_mean qdlet_mean+qdlet_std qdlet_mean-qdlet_std; dblet_mean dblet_mean+dblet_std dblet_mean-dblet_std]);
  %errorbar(repmat(peduncle_indx, 5, 1).', [whole_mean; regen_mean; bowti_mean; thinn_mean; slimm_mean].', [whole_std; regen_std; bowti_std; thinn_std; slimm_std].');
  %errorbar(peduncle_indx, whole_mean, whole_std);
  %axis('equal')
  legend({'whole', 'flat', 'step', 'RA', 'narrow', 'trimm', 'qdlet', 'dblet'})

  H71 = cell(1, 2*nouter);
  p71 = cell(1, 2*nouter);
  H72 = NaN(8, 2*nouter);
  p72 = NaN(8, 2*nouter);
  p73 = NaN(8, 2*nouter);
  %part_indx = [ones(nwhole, 1); ones(nregen, 1)*2; ones(nbowti, 1)*3; ones(nthinn, 1)*4; ones(nslimm, 1)*5; ones(nquad, 1)*6; ones(ndoubl, 1)*7];
  for i=1:2*nouter
    [H71{i},p71{i}] = myttest([whole_norm(:, i); regen_norm(:, i); bowti_norm(:, i); treat_norm(:, i); thinn_norm(:, i); slimm_norm(:, i); qdlet_norm(:,i); dblet_norm(:,i)], indx);
    %H72(1:size(H71{i},1),i) = H71{i}(1,:).';
    %p72(1:size(p71{i},1),i) = p71{i}(1,:).';
    H72(:,i) = H71{i}(1,:).';
    p72(:,i) = p71{i}(1,:).';
  end

  p73 = NaN(8, 2*nouter);
  whole_counts = nansum(whole_norm, 1);
  regen_counts = nansum(regen_norm, 1);
  bowti_counts = nansum(bowti_norm, 1);
  treat_counts = nansum(treat_norm, 1);
  thinn_counts = nansum(thinn_norm, 1);
  slimm_counts = nansum(slimm_norm, 1);
  dblet_counts = nansum(dblet_norm, 1);
  qdlet_counts = nansum(qdlet_norm, 1);

  whole_sizes = nansum(~isnan(whole_norm), 1);
  regen_sizes = nansum(~isnan(regen_norm), 1);
  bowti_sizes = nansum(~isnan(bowti_norm), 1);
  treat_sizes = nansum(~isnan(treat_norm), 1);
  thinn_sizes = nansum(~isnan(thinn_norm), 1);
  slimm_sizes = nansum(~isnan(slimm_norm), 1);
  dblet_sizes = nansum(~isnan(dblet_norm), 1);
  qdlet_sizes = nansum(~isnan(qdlet_norm), 1);

  pkg load nan
  all_counts = cat(1, whole_counts, regen_counts, bowti_counts, treat_counts, thinn_counts, slimm_counts, dblet_counts, qdlet_counts);
  all_sizes =  cat(1, whole_sizes, regen_sizes, bowti_sizes, treat_sizes, thinn_sizes, slimm_sizes, dblet_sizes, qdlet_sizes);

  for i=1:2*nouter
    for j=1:size(all_sizes, 1)
      p73(j,i) = fishers_exact_test(all_counts(j,i), all_sizes(j,i) - all_counts(j,i), all_counts(2,i), all_sizes(2,i) - all_counts(2,i));
      if all_sizes(j,i)==0
        p73(j,i) = NaN;
      end
    end
  end
  H73 = (p73 < 0.05) + (p73 < 0.01) + (p73 < 0.001);
  H73(isnan(p73)) = NaN;

  pkg unload nan
  keyboard

  whole_norm = [whole_norm(:,1:nouter); whole_norm(:, end:-1:end-nouter+1)];
  regen_norm = [regen_norm(:,1:nouter); regen_norm(:, end:-1:end-nouter+1)];
  bowti_norm = [bowti_norm(:,1:nouter); bowti_norm(:, end:-1:end-nouter+1)];
  treat_norm = [treat_norm(:,1:nouter); treat_norm(:, end:-1:end-nouter+1)];
  thinn_norm = [thinn_norm(:,1:nouter); thinn_norm(:, end:-1:end-nouter+1)];
  slimm_norm = [slimm_norm(:,1:nouter); slimm_norm(:, end:-1:end-nouter+1)];
  qdlet_norm = [qdlet_norm(:,1:nouter); qdlet_norm(:, end:-1:end-nouter+1)];
  dblet_norm = [dblet_norm(:,1:nouter); dblet_norm(:, end:-1:end-nouter+1)];

  whole_mean = nanmean(whole_norm, 1);
  regen_mean = nanmean(regen_norm, 1);
  bowti_mean = nanmean(bowti_norm, 1);
  treat_mean = nanmean(treat_norm, 1);
  thinn_mean = nanmean(thinn_norm, 1);
  slimm_mean = nanmean(slimm_norm, 1);
  qdlet_mean = nanmean(qdlet_norm, 1);
  dblet_mean = nanmean(dblet_norm, 1);
  whole_std = nanstd(whole_norm, 1);
  regen_std = nanstd(regen_norm, 1);
  bowti_std = nanstd(bowti_norm, 1);
  treat_std = nanstd(treat_norm, 1);
  thinn_std = nanstd(thinn_norm, 1);
  slimm_std = nanstd(slimm_norm, 1);
  qdlet_std = nanstd(qdlet_norm, 1);
  dblet_std = nanstd(dblet_norm, 1);

  %subplot(1,2,2);
  %plot(repmat(peduncle_indx(1:nouter), 1, 3), [whole_mean whole_mean+whole_std whole_mean-whole_std]);
  %hold on;
  %plot(repmat(peduncle_indx(end:-1:end-nouter+1), 1, 3), [regen_mean regen_mean+regen_std regen_mean-regen_std; bowti_mean bowti_mean+bowti_std bowti_mean-bowti_std; thinn_mean thinn_mean+thinn_std thinn_mean-thinn_std; slimm_mean slimm_mean+slimm_std slimm_mean-slimm_std; qdlet_mean qdlet_mean+qdlet_std qdlet_mean-qdlet_std]);

  keyboard

  % General pictures
  figure;
  subplot(2,3,1);
  boxplot({whole_fin(:, 2), regen_fin(:, 2), bowti_fin(:,2), thinn_fin(:, 2), slimm_fin(:, 2)});
  [H11,p11] = myttest([whole_fin(:, 2); regen_fin(:, 2); bowti_fin(:,2); thinn_fin(:, 2); slimm_fin(:, 2)], [ones(size(whole_fin(:, 2))); ones(size(regen_fin(:, 2)))*2; ones(size(bowti_fin(:,2)))*3; ones(size(thinn_fin(:, 2)))*4; ones(size(slimm_fin(:, 2)))*5]);
  title('Fin length (um)')

  subplot(2,3,2);
  boxplot({whole_fin(:, 2)-whole_fin(:,3), regen_fin(:, 2)-regen_fin(:,3), bowti_fin(:, 2)-bowti_fin(:,3), thinn_fin(:, 2)-thinn_fin(:,3), slimm_fin(:, 2)-slimm_fin(:,3)});
  [H12,p12] = myttest([whole_fin(:, 2)-whole_fin(:,3); regen_fin(:, 2)-regen_fin(:,3); bowti_fin(:, 2)-bowti_fin(:,3); thinn_fin(:, 2)-thinn_fin(:,3); slimm_fin(:, 2)-slimm_fin(:,3)], [ones(size(whole_fin(:, 2))); ones(size(regen_fin(:, 2)))*2; ones(size(bowti_fin(:, 2)))*3; ones(size(thinn_fin(:, 2)))*4; ones(size(slimm_fin(:, 2)))*5]);
  title('Cleft depth (um)')

  subplot(2,3,3);
  boxplot({whole_fin(:, 4), regen_fin(:, 4), bowti_fin(:,4), thinn_fin(:, 4), slimm_fin(:, 4)});
  [H13,p13] = myttest([whole_fin(:, 4); regen_fin(:, 4); bowti_fin(:, 4); thinn_fin(:, 4); slimm_fin(:, 4)], [ones(size(whole_fin(:, 2))); ones(size(regen_fin(:, 2)))*2; ones(size(bowti_fin(:, 2)))*3; ones(size(thinn_fin(:, 2)))*4; ones(size(slimm_fin(:, 2)))*5]);
  title('Fin width (um)')

  subplot(2,3,4);
  boxplot({whole_fin(:, 2)./whole_fin(:, 5), regen_fin(:, 2)./regen_fin(:, 5), bowti_fin(:,2)./bowti_fin(:, 5), thinn_fin(:, 2)./thinn_fin(:,5), slimm_fin(:, 2)./slimm_fin(:,5)});
  [H21,p21] = myttest([whole_fin(:, 2)./whole_fin(:,5); regen_fin(:, 2)./regen_fin(:,5); bowti_fin(:,2)./bowti_fin(:,5); thinn_fin(:, 2)./thinn_fin(:,5); slimm_fin(:, 2)./slimm_fin(:,5)], [ones(size(whole_fin(:, 2))); ones(size(regen_fin(:, 2)))*2; ones(size(bowti_fin(:,2)))*3; ones(size(thinn_fin(:, 2)))*4; ones(size(slimm_fin(:, 2)))*5]);
  title('Fin length (um)')

  subplot(2,3,5);
  boxplot({(whole_fin(:, 2)-whole_fin(:,3))./whole_fin(:,5), (regen_fin(:, 2)-regen_fin(:,3))./regen_fin(:,5), (bowti_fin(:, 2)-bowti_fin(:,3))./bowti_fin(:,5), (thinn_fin(:, 2)-thinn_fin(:,3))./thinn_fin(:,5), (slimm_fin(:, 2)-slimm_fin(:,3))./slimm_fin(:,5)});
  [H22,p22] = myttest([(whole_fin(:, 2)-whole_fin(:,3))./whole_fin(:,5); (regen_fin(:, 2)-regen_fin(:,3))./regen_fin(:,5); (bowti_fin(:, 2)-bowti_fin(:,3))./bowti_fin(:,5); (thinn_fin(:, 2)-thinn_fin(:,3))./thinn_fin(:,5); (slimm_fin(:, 2)-slimm_fin(:,3))./slimm_fin(:,5)], [ones(size(whole_fin(:, 2))); ones(size(regen_fin(:, 2)))*2; ones(size(bowti_fin(:, 2)))*3; ones(size(thinn_fin(:, 2)))*4; ones(size(slimm_fin(:, 2)))*5]);
  title('Cleft depth (um)')

  subplot(2,3,6);
  boxplot({whole_fin(:, 4)./whole_fin(:,5), regen_fin(:, 4)./regen_fin(:,5), bowti_fin(:,4)./bowti_fin(:,5), thinn_fin(:, 4)./thinn_fin(:,5), slimm_fin(:, 4)./slimm_fin(:,5)});
  [H23,p23] = myttest([whole_fin(:, 4)./whole_fin(:,5); regen_fin(:, 4)./regen_fin(:,5); bowti_fin(:, 4)./bowti_fin(:,5); thinn_fin(:, 4)./thinn_fin(:,5); slimm_fin(:, 4)./slimm_fin(:,5)], [ones(size(whole_fin(:, 2))); ones(size(regen_fin(:, 2)))*2; ones(size(bowti_fin(:, 2)))*3; ones(size(thinn_fin(:, 2)))*4; ones(size(slimm_fin(:, 2)))*5]);
  title('Fin width (um)')

  keyboard

  nfins = size(all_fins, 1);
  rel_stats = NaN(nfins, 5);
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

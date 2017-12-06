function study_stitching

  trial1 = [0 5;
    0 8;
    0 17;
    0 39;
    1 79];

  [x,f] = get_survival(trial1);
  figure;
  stairs(x,f)
  ylim([0 1.1]);


  latex = [1 15;
      0 8;
      1 22;
      1 22;
      1 22;
      0 10;
      0 3;
      0 2;
      0 4;
      0 5];

  [x,f] = get_survival(latex);
  %figure;
  hold on;
  stairs(x,f)
  %ylim([0 1.1]);

  pdms = [1 1;
      1 20;
      0 2;
      0 4;
      0 4;
      1 22;
      0 1;
      0 1;
      0 4;
      0 4];

  [x,f] = get_survival(pdms);
  %figure;
  stairs(x,f)
  %ylim([0 1.1]);

  keyboard

  return;
end

function [x,f] = get_survival(vals)

  [f,x] = ecdf(vals(:,2), 'censoring', vals(:,1), 'function','survivor');
  f = [1; f; f(end)];
  x = [0; x; max(vals(:,2))];

  return;
end

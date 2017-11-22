function [periods] = find_period(vals, range, tol)

  vals = sort(vals);
  mval = max(vals);

  tests = [range(1):range(2)].';

  nperiods = length(tests);
  periods = NaN(nperiods, 1);

  if (tol < 1)
    thresh = 1 + [-tol +tol];
  else
    thresh = [-tol tol];
  end

  for i=1:nperiods
    p = tests(i);
    nmult = floor(mval/p);
    nhits = zeros(nmult, 1);
    for n=1:nmult
      curr_per = n*p;

      if (thresh(1)<0)
        currn = sum(vals >= curr_per+thresh(1) & vals <= curr_per+thresh(2));
      else
        currn = sum(vals >= curr_per*thresh(1) & vals <= curr_per*thresh(2));
      end

      nhits(n) = currn;
      if nhits(1)==0
        break;
      end
    end

    if (nhits(1)>0)
      nexpec = nhits(1)./[1:nmult].';
      cc = cov(nhits, nexpec);

      periods(i) = cc(1,2);
    end
  end

  if (tol<1)
    [vals, indxs] = local_extrema(periods, tol*mean(tests));
  else
    [vals, indxs] = local_extrema(periods, tol);
  end

  [vals, rev] = sort(vals, 'descend');
  indxs = indxs(rev);

  periods = [vals(:), tests(indxs(:))];

  return;
end

function [periods, hits] = find_period(vals, weight, range, tol)

  mval = max(vals);

  tests = [range(1):range(2)].';

  nperiods = length(tests);
  periods = NaN(nperiods, 2);

  if (tol < 1)
    thresh = 1 + [-tol +tol];
  else
    thresh = [-tol tol];
  end

  for i=1:nperiods
    p = tests(i);
    nmult = floor(mval/p);
    nhits = zeros(nmult, 2);
    for n=1:nmult
      curr_per = n*p;

      if (thresh(1)<0)
        goods = (vals >= curr_per+thresh(1) & vals <= curr_per+thresh(2));
        currn = sum(goods);
        w = sum(weight(goods));
      else
        goods = (vals >= curr_per*thresh(1) & vals <= curr_per*thresh(2));
        currn = sum(goods);
        w = sum(weight(goods));
      end

      nhits(n,1) = currn;
      nhits(n,2) = w;
      if nhits(1,1)==0
        break;
      end
    end

    if (nhits(1,1)>0)
      nexpec = nhits(1,1)./[1:nmult].';
      cc = cov(nhits(:,1), nexpec)/(std(nhits(:,1))*std(nexpec));
      height = (nhits(:,2)./nhits(:,1));
      height(isnan(height)) = 0;

      periods(i,1) = cc(1,2);
      periods(i,2) = mean(height);
    end
  end

  if (tol<1)
    [vals, indxs] = local_extrema(periods(:,1).*periods(:,2), tol*mean(tests));
  else
    [vals, indxs] = local_extrema(periods(:,1).*periods(:,2), tol);
  end

  [vals, rev] = sort(vals, 'descend');
  indxs = indxs(rev);

  periods = [vals(:), tests(indxs(:))];

  return;
end

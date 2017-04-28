function [img, params] = scaled_cast(img, params, class_type)
% Simon Blanchoud
% 15.05.2014

  orig_type = class(img);

  if (class_type(1) == 'd')
    img = double(img);

    return;
  elseif (orig_type(1) == 'd')
    img = cast(img, class_type);

    return;
  end

  if (isempty(params))
    omin = double(intmin(orig_type));
    omax = double(intmax(orig_type));
    oran = omax - omin;

    cmin = double(intmin(class_type));
    cmax = double(intmax(class_type));
    cran = cmax - cmin;

    params = [cran/oran cmin - omin*cran/oran];
  end

  img = double(img)*params(1) + params(2);
  img = cast(img, class_type);

  return;
end

function img = imwhitebalance(img, prct)

  if (size(img, 3) ~= 3 || isfloat(img))
    error('Need RGB UINT images');
  end

  type = class(img);
  maxval = double(intmax(type));

  img = double(img) / maxval;

  if ndims(img) > 3
    ei = prctile(reshape(permute(img, [1 2 4 3]), [], 3), 100-prct);
  else
    ei = prctile(reshape(img, [], 3), 100-prct);
  end
  img = maxval * bsxfun(@times, img, 1./ permute(ei, [1 3 2]));

  img = cast(img, type);

  return;
end

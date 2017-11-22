function pixels = extract_ROI_pixels(imgs, ROIs)

  if (ischar(imgs))
    [nimgs,ssize] = size_data(imgs);
    imgs = load_data(imgs, [1:nimgs]);
    imgs = imwhitebalance(imgs, 45);
  else
    [nimgs, ssize] = size_data(imgs);
    nimgs = nimgs(end);
  end

  nROIs = length(ROIs);
  pixels = cell(nROIs, 1);

  for i=1:nROIs
    type = ROIs{i}.strType;

    if (strncmp(type, 'Polygon', 7))
      slice = ROIs{i}.nPosition;
      pos = ROIs{i}.mnCoordinates;

      if (any(pos(1,:) ~= pos(end,:)))
        pos(end+1,:) = pos(1,:);
      end

      mask = poly2mask(pos(:,1),pos(:,2),ssize(1),ssize(2));
      img = squeeze(imgs(:,:,:,slice));
      pix = img(repmat(mask, [1 1 ssize(3)]));
      pixels{i} = reshape(pix, [], ssize(3));
    end
  end

  return;
end

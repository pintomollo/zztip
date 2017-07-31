function export_csv(fname, title, data, headers)

  fid = fopen(fname, 'a');

  if (fid==-1)
    error(['Cannot create the file ' fname]);
  end
  nfields = size(data,2);

  fprintf(fid, ['"%s"' repmat(',', 1, nfields) '\n'], title);
  if (length(headers)==nfields)
    for i=1:nfields
      fprintf(fid, '"%s",', headers{i});
    end
    fprintf(fid, '\n');
  end
  for i=1:size(data,1)
    fprintf(fid, '%.4e,', data(i,:));
    fprintf(fid, '\n');
  end
  fprintf(fid, '\n');
  fclose(fid);

  return;
end

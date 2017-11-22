function [values, fields] = parse_metadata_table(fname)

  fid = fopen(fname, 'r');
  if (ischar(fid))
    fid = fname;
  elseif (fid == -1)
    error([fname ' cannot be opened']);
  end

  header = fgetl(fid);
  if (header(1) ~= '#')
    error(['Header line in ' fname ' is missing']);
  end
  fields = strsplit(header(2:end));

  ncolumns = length(fields);
  parsing = '';
  for (i=1:ncolumns)
    tmp = strsplit(fields{i}, '%');

    if (length(tmp)<2)
      type = 's';
    else
      type = tmp{2};
    end
    fields{i} = tmp{1};
    parsing = [parsing '%' type ' '];
  end

  values = textscan(fid, parsing, 'CommentStyle', '#');
  fclose(fid);

  return;
end

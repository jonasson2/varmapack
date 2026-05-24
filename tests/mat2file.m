function mat2file(fid, A, name)
  if nargin < 3
    name = inputname(2);
    if isempty(name)
      error('When A is a compound expression, filename must be specified');
    end
  end
  if isscalar(A)
    dims = [];
  elseif isvector(A)
    dims = numel(A);
  else
    dims = size(A);
  end
  fprintf(fid, "var %s %d", name, numel(dims));
  for i=1:numel(dims)
    fprintf(fid, " %d", dims(i));
  end
  fprintf(fid, "\n");
  for i = 1:numel(A)
    if i > 1, fprintf(fid, " "); end
    fprintf(fid, "%.17g", A(i));
  end
  fprintf(fid, "\n");
end

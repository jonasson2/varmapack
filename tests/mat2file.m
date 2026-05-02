function mat2file(fid, A, name)
  if nargin < 3
    name = inputname(2);
    if isempty(name)
      error('When A is a compound expression, filename must be specified');
    end
  end
  [m, n] = size(A);
  fprintf(fid, "%s,%d,%d\n", name, m, n);
  for i=1:numel(A)
    fprintf(fid, "%.17g,", A(i));
  end
  fprintf(fid, "\n");
end
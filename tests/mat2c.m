function mat2c(A, name)
  if nargin < 2
    name = inputname(1);
    if isempty(name)
      error('When A is a compound expression, filename must be specified');
    end
    if isempty(A)
      error('A must not be empty');
    end
  end
  [m, n] = size(A);
  filename = name + ".txt";
  f = fopen(filename, 'w');
  fprintf(f, "%s,%d,%d\n%.17g", name, m, n, A(1));
  for i=2:numel(A)
    fprintf(f, ",%.17g", A(i));
  end
  fprintf(f, ",\n");
  fclose(f);
end
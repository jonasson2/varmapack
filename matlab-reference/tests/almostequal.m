function close = almostequal(x, y, tol)
  if nargin < 3, tol = 5e-14; end
  if iscell(x)
    if ~iscell(y), close = false; return, end
    x = cell2mat(x);
    y = cell2mat(y);
  end
  M = max([1; abs(x(:)); abs(y(:))]);
  close = isequal(size(x), size(y)) && all(abs(x(:) - y(:)) < tol*M);
end

function assert_equal(x, y, tol)
  if nargin < 3, tol = 0; end
  if iscell(x), x = cell2mat(x); end
  if iscell(y), y = cell2mat(y); end
  ascertain(isequal(size(x), size(y)), 'Sizes differ');
  M = max([1; abs(x(:)); abs(y(:))]);
  d = max(abs(x(:) - y(:)))/M;
  ascertain(d <= tol, sprintf('Maximum difference %.3e exceeds %.3e', d, tol));
end

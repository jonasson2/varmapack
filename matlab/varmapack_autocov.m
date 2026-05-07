function C = varmapack_autocov(X, maxlag, norm)
%VARMAPACK_AUTOCOV  Sample autocovariance of a multivariate time series.
  [~, n] = size(X);
  if nargin < 2 || isempty(maxlag), maxlag = n - 1; end
  if nargin < 3 || isempty(norm), norm = "ML"; end
  norm = char(norm);
  C = varmapack_autocov_gateway(X, maxlag, norm);
end

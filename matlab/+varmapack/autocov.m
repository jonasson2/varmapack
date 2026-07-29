function C = autocov(X, maxlag, norm)
%VARMAPACK.AUTOCOV  Sample autocovariance of a multivariate time series.
%
%   C = VARMAPACK.AUTOCOV(X) computes autocovariance matrices for all lags
%   of the series X, where x(t) is stored in column t.
%
%   C = VARMAPACK.AUTOCOV(X,MAXLAG) computes lags 0,...,MAXLAG, where the
%   lag-k matrix is Cov(x(t), x(t-k)).
%
%   C = VARMAPACK.AUTOCOV(X,MAXLAG,NORM) uses NORM = "ML" for 1/n
%   normalization, or "C", "corr", or "corrected" for 1/(n-k).
%
%   See also VARMAPACK.ACVF.
  [~, n] = size(X);
  if nargin < 2 || isempty(maxlag), maxlag = n - 1; end
  if nargin < 3 || isempty(norm), norm = "ML"; end
  norm = char(norm);
  C = varmapack_autocov_gateway(X, maxlag, norm);
end

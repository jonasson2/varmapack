function C = ref_varma_autocov(X, maxlag, norm)
%REF_VARMA_AUTOCOV  Sample autocovariance of a multivariate time series.
%
%   C = REF_VARMA_AUTOCOV(X) computes autocovariance matrices for all lags of
%   the series X, where x(t) is stored in column t.
%
%   C = REF_VARMA_AUTOCOV(X, MAXLAG) computes lags 0,...,MAXLAG.
%
%   C = REF_VARMA_AUTOCOV(X, MAXLAG, NORM) uses NORM = "ML" for 1/n
%   normalization, or "C", "corr", or "corrected" for 1/(n-k).

  [r, n] = size(X);
  if nargin < 2 || isempty(maxlag), maxlag = n - 1; end
  if nargin < 3 || isempty(norm), norm = "ML"; end
  if maxlag < 0 || maxlag >= n
    error('maxlag must satisfy 0 <= maxlag < n')
  end
  norm = lower(string(norm));
  ML = norm == "ml";
  CORRECTED = ismember(norm, ["c", "corr", "corrected"]);
  if ~ML && ~CORRECTED
    error('norm must be "ML", "C", "corr", or "corrected"')
  end
  X = X - mean(X, 2);
  C = zeros(r, r, maxlag + 1);
  for k = 0:maxlag
    if ML, f = n; else, f = n - k; end
    Y = X(:, 1:end-k);
    Z = X(:, 1+k:end);
    C(:, :, k+1) = Y*Z'/f;
  end
end

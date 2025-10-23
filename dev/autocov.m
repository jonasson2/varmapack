function S = autocov(X,k)
  %AUTOCOV  Sample autocovariance of multivariate series X at lag k
  % X is n×r (rows=time, cols=vars). Returns S = Cov(X_t, X_{t-k}) (r×r).

  [n,r] = size(X);
  if k < 0 || k >= n, error('k must satisfy 0 <= k < n'); end
  N = n - k;
  X1 = X(1+k:end,:);
  X0 = X(1:end-k,:);
  X1 = X1 - mean(X1,1);
  X0 = X0 - mean(X0,1);
  if N > 1
    S = (X1.' * X0) / (N - 1);   % unbiased
  else
    S = NaN(r,r);
  end
end
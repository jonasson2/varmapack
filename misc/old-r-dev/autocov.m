function S = autocov(X,k)
  %AUTOCOV  Sample autocovariance of multivariate series X at lag k
  % X is r×n (rows=vars, cols=time). Returns S = Cov(X_t, X_{t-k}) (r×r).

  [n,r] = size(X);
  if k < 0 || k >= n, error('k must satisfy 0 <= k < n'); end
  
  X1 = X(:, 1+k:end);
  X0 = X(:, 1:end-k);
  X0 = X0 - mean(X0,2);
  X1 = X1 - mean(X1,2);
  if n > 1
    S = (X1*X0')/(n - 1);   % unbiased
  else
    S = NaN(r,r);
  end
end

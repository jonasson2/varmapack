function [X, E, ok] = varmapack_sim(A, B, Sig, n, mu, M, x0)
%VARMAPACK_SIM  Simulate a VARMA time series with the C varmapack engine.
%
%   X = VARMAPACK_SIM(A,B,Sig,n) simulates one zero-mean VARMA series.
%   [X,E] = VARMAPACK_SIM(...) also returns the shocks. A and B are r-by-r*p
%   and r-by-r*q matrices, with [] allowed for p=0 or q=0.

  r = size(Sig, 1);
  if isempty(A), A = zeros(r, 0); end
  if isempty(B), B = zeros(r, 0); end
  if nargin < 5 || isempty(mu)
    mu = [];
  else
    if isscalar(mu), mu = repmat(mu, r, 1); end
    mu = mu(:);
  end
  if nargin < 6 || isempty(M), M = 1; end
  if nargin < 7, x0 = []; end
  global varmapack_rng
  rng = varmapack_rng;
  if isempty(rng), error('varmapack_sim: no Randompack RNG specified'); end
  [X, E, ok] = varmapack_sim_gateway(A, B, Sig, n, mu, M, x0, rng);
  if r == 1 && M == 1
    X = reshape(X, 1, n);
    E = reshape(E, 1, n);
  elseif r == 1
    X = reshape(X, n, M);
    E = reshape(E, n, M);
  else
    X = reshape(X, r, n, M);
    E = reshape(E, r, n, M);
  end
end

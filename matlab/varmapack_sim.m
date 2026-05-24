function [X, E] = varmapack_sim(A, B, Sig, n, mu, M, X0, rng)
%VARMAPACK_SIM  Simulate a VARMA time series with the C varmapack engine.
%
%   X = VARMAPACK_SIM(A,B,Sig,n) simulates one zero-mean VARMA series.
%   [X,E] = VARMAPACK_SIM(...) also returns the shocks. A and B are r-by-r*p
%   and r-by-r*q matrices, with [] allowed for p=0 or q=0.
%   VARMAPACK_SIM(...,X0,rng) uses the Randompack RNG handle rng.

  r = size(Sig, 1);
  if isempty(A), A = zeros(r, 0); end
  if isempty(B), B = zeros(r, 0); end
  if nargin < 5 || isempty(mu)
    mu = [];
  else
    if isscalar(mu), mu = repmat(mu, r, 1); end
    if size(mu, 1) ~= r, error('varmapack_sim: mu must have r rows'); end
  end
  if nargin < 6 || isempty(M), M = 1; end
  if nargin < 7, X0 = []; end
  if nargin < 8 || isempty(rng)
    global varmapack_rng
    rng = varmapack_rng;
  end
  if isempty(rng), error('varmapack_sim: no Randompack RNG specified'); end
  if nargout < 2
    X = varmapack_sim_gateway(A, B, Sig, n, mu, M, X0, rng);
    E = [];
  else
    [X, E] = varmapack_sim_gateway(A, B, Sig, n, mu, M, X0, rng);
  end
  if r == 1 && M == 1
    X = reshape(X, 1, n);
    if nargout >= 2, E = reshape(E, 1, n); end
  elseif r == 1
    X = reshape(X, n, M);
    if nargout >= 2, E = reshape(E, n, M); end
  else
    X = reshape(X, r, n, M);
    if nargout >= 2, E = reshape(E, r, n, M); end
  end
end

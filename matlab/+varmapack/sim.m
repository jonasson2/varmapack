function [X, E] = sim(A, B, Sig, mu, n, M, X0, rng)
%VARMAPACK.SIM  Simulate a VARMA time series with the C varmapack engine.
%
%   [X,E] = VARMAPACK.SIM(A,B,Sig,mu,n,M,X0,rng) simulates M independent
%   VARMA(p,q) series of length n:
%
%     x(t) - mu(t) = sum Ai*(x(t-i) - mu(t-i)) + eps(t) + sum Bi*eps(t-i),
%
%   where eps(t) is N(0,Sig). A and B have r rows and contain the coefficient
%   blocks [A1,...,Ap] and [B1,...,Bq]; Sig is r-by-r. Use [] for A or B to
%   select a pure VMA or VAR model, respectively.
%
%   mu is [] or an r-by-nmu sequence of series means. A scalar is expanded to
%   all r components, and the last supplied column repeats through time n. M
%   defaults to 1.
%
%   With X0=[], a stationary model is initialized by drawing the startup states
%   and shocks from their exact joint distribution, so no burn-in is needed. X0
%   may instead be r-by-nX0 for a common fixed startup path or r-by-nX0-by-M
%   for replicate-specific paths, where max(p,q) <= nX0 <= n. For stationary
%   models, the resulting paths have the exact conditional distribution given
%   X0. A nonstationary model requires X0. With MA terms, its startup shocks are
%   drawn conditionally on the residual equations implied by X0; Sig must then
%   be positive definite.
%
%   rng is a VARMAPACK.RNG object. If it is omitted or empty, a temporary
%   randomized generator is created for this call.
%
%   For r > 1, X and E have shape r-by-n-by-M. Scalar models return 1-by-n for
%   M=1 and n-by-M otherwise. E is returned only when requested.

  r = size(Sig, 1);
  if isempty(A), A = zeros(r, 0); end
  if isempty(B), B = zeros(r, 0); end
  if nargin < 4 || isempty(mu)
    mu = [];
  else
    if isscalar(mu), mu = repmat(mu, r, 1); end
    if size(mu, 1) ~= r, error('varmapack_sim: mu must have r rows'); end
  end
  if nargin < 5 || isempty(n), error('varmapack_sim: n must be specified'); end
  if nargin < 6 || isempty(M), M = 1; end
  if nargin < 7, X0 = []; end
  temporary_rng = [];
  if nargin < 8 || isempty(rng)
    temporary_rng = varmapack.Rng();
    cleanup = onCleanup(@() delete(temporary_rng));
    rng = temporary_rng;
  end
  if ~isa(rng, 'varmapack.Rng') || ~isscalar(rng)
    error('varmapack:sim:rng', 'rng must be a varmapack.Rng object');
  end
  rng = rng.mex_handle();
  if nargout < 2
    X = varmapack_sim_gateway(A, B, Sig, mu, n, M, X0, rng);
    E = [];
  else
    [X, E] = varmapack_sim_gateway(A, B, Sig, mu, n, M, X0, rng);
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

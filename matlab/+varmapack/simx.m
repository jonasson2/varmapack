function [X, E] = simx(A, B, C, z, Sig, n, M, X0, h, rng)
%VARMAPACK.SIMX  Simulate a VARMAX time series with the C varmapack engine.
%
%   [X,E] = VARMAPACK.SIMX(A,B,C,z,Sig,n,M,X0,h,rng) simulates
%
%     x(t) = eps(t) + sum Ai*x(t-i) + sum Bi*eps(t-i) + sum Ci*z(t-i+1),
%
%   using X0 as a known startup segment. C is r-by-s, with C(:,i) multiplying
%   the scalar sequence z(t-i+1); C=[] specifies no exogenous terms. z has at
%   least n rows and is either n-by-M for replicate-specific data or n-by-1 to
%   broadcast the same sequence to all replicates. X0 is r-by-h or r-by-h-by-M.
%
%   The formula time axis is zero-based: X0(:,1) is x(0), and z(1,:) is z(0).
%   MATLAB column h+1 is therefore formula time h.
%
%   M defaults to 1 and h defaults to size(X0,2). X0 is required, n must be at
%   least h, and h must be greater than max(p,s-1).
%
%   rng is a VARMAPACK.RNG object. If it is omitted or empty, a temporary
%   randomized generator is created for this call.
%
%   For r > 1, X and E have shape r-by-n-by-M. Scalar models return 1-by-n for
%   M=1 and n-by-M otherwise. E is returned only when requested.

  r = size(Sig, 1);
  if isempty(A), A = zeros(r, 0); end
  if isempty(B), B = zeros(r, 0); end
  if isempty(C), C = zeros(r, 0); end
  if nargin < 7 || isempty(M), M = 1; end
  if nargin < 8 || isempty(X0)
    error('varmapack_simx: X0 must be specified');
  end
  if nargin < 9 || isempty(h), h = size(X0, 2); end
  temporary_rng = [];
  if nargin < 10 || isempty(rng)
    temporary_rng = varmapack.Rng();
    cleanup = onCleanup(@() delete(temporary_rng));
    rng = temporary_rng;
  end
  if ~isa(rng, 'varmapack.Rng') || ~isscalar(rng)
    error('varmapack:simx:rng', 'rng must be a varmapack.Rng object');
  end
  rng = rng.mex_handle();
  if nargout < 2
    X = varmapack_simx_gateway(A, B, C, z, Sig, n, M, X0, h, rng);
    E = [];
  else
    [X, E] = varmapack_simx_gateway(A, B, C, z, Sig, n, M, X0, h, rng);
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

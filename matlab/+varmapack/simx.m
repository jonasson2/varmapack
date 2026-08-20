function [X, E] = simx(A, B, C, z, Sig, n, M, X0, h, rng)
%VARMAPACK.SIMX  Simulate a VARMAX time series with the C varmapack engine.
%
%   [X,E] = VARMAPACK.SIMX(A,B,C,z,Sig,n,M,X0,h,rng) simulates
%
%     x(t) = eps(t) + sum Ai*x(t-i) + sum Bi*eps(t-i) + sum Ci*z(t-i+1),
%
%   using X0 as a known startup segment. C consists of r-by-d blocks C_i and
%   has shape r-by-(d*s). Each z(t) is a d-vector. z has shape d-by-n to
%   broadcast a common path or d-by-n-by-M for replicate-specific paths. The
%   historical scalar forms n-by-1 and n-by-M are also accepted. C=[] specifies
%   no exogenous terms. X0 is r-by-h or r-by-h-by-M. It may be empty when
%   max(p,q,s-1) is zero.
%
%   The formula time axis is zero-based: X0(:,1) is x(0), and z(1,:) is z(0).
%   MATLAB column h+1 is therefore formula time h.
%
%   M defaults to 1 and h defaults to size(X0,2). The minimum startup length is
%   max(p,q,s-1), and n must be at least h.
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
  if nargin < 8 || isempty(X0), X0 = zeros(r, 0); end
  if nargin < 9 || isempty(h), h = size(X0, 2); end
  [C, z] = normalize_exogenous(C, z, r, n, M);
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

function [C, z] = normalize_exogenous(C, z, r, n, M)
  if isempty(C)
    C = zeros(r, 0);
    z = zeros(0, n);
    return
  end
  if size(C, 1) ~= r || ndims(C) ~= 2
    error('varmapack:simx:C', 'C must be an r-by-(d*s) matrix');
  end
  if ndims(z) == 2 && (size(z, 2) == 1 || size(z, 2) == M) && size(z, 1) >= n
    z = reshape(z, 1, size(z, 1), size(z, 2));
  elseif ndims(z) == 2
    z = reshape(z, size(z, 1), size(z, 2), 1);
  elseif ndims(z) ~= 3
    error('varmapack:simx:z', 'z must have shape d-by-n or d-by-n-by-M');
  end
  d = size(z, 1);
  if d == 0 || size(z, 2) < n || (size(z, 3) ~= 1 && size(z, 3) ~= M)
    error('varmapack:simx:z', 'z must have shape d-by-n or d-by-n-by-M');
  end
  if mod(size(C, 2), d) ~= 0
    error('varmapack:simx:C', 'C must have d columns per term');
  end
end

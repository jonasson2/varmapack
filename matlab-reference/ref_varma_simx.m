%REF_VARMA_SIMX  Reference simulation of a VARMAX model with fixed exogenous data
%
%  [X,E] = REF_VARMA_SIMX(A,B,C,z,Sig,n,M,X0,h,e0,rng) simulates
%
%    x(t) = eps(t) + sum Ai*x(t-i) + sum Bi*eps(t-i) + sum Ci*z(t-i+1),
%
%  using X0 as a known startup segment. The exogenous sequence z is M by n for
%  replicate-specific data, or 1 by n to broadcast the same sequence to all
%  replicates. The startup segment X0 is r by h, or r by h by M.
%
%  The formula time axis is zero-based: X0(:,1) is x(0), and z(:,1) is z(0).
%  MATLAB column h+1 is therefore formula time h.
%
%  The number of exogenous terms is s = size(C,2), with s = 0 when C is empty.
%
%  The e0 parameter exists to test ref_varma_simx against ref_varma_sim. It is
%  not present in the corresponding varmapack_simx C function.

function [X, E] = ref_varma_simx(A, B, C, z, Sig, n, M, X0, h, e0, rng)
  r = size(Sig, 1);
  if isempty(A), A = zeros(r,0); end
  if isempty(B), B = zeros(r,0); end
  if nargin < 7 || isempty(M), M = 1; end
  if nargin < 8 || isempty(X0), error('X0 must be specified'); end
  if nargin < 9 || isempty(h), h = size(X0, 2); end
  if nargin < 10, e0 = []; end
  if nargin < 11, rng = []; end
  [p, q, r] = get_dimensions(A, B, Sig);
  if size(X0, 1) ~= r || size(X0, 2) ~= h || ...
      (ndims(X0) == 3 && size(X0, 3) ~= M) || ndims(X0) > 3
    error('X0 must be r by h or r by h by M');
  end
  if n < h, error('n must be at least h'); end
  if isempty(C)
    s = 0;
    Cc = {};
    z = zeros(1, n);
  else
    if size(C, 1) ~= r, error('C must have r rows'); end
    s = size(C, 2);
    Cc = cell(1, s);
    for i = 1:s
      Cc{i} = C(:, i);
    end
  end
  if size(z, 2) < n, error('z must have at least n columns'); end
  if size(z, 1) ~= 1 && size(z, 1) ~= M
    error('z must be 1 by n or M by n');
  end
  zlag = max(s - 1, 0);
  t0 = max(p, zlag);
  if h <= t0, error('h must be greater than max(p,s-1)'); end
  first_active_shock = t0 - q;
  first_shock_time = min(first_active_shock, 0);
  he = h - first_shock_time;
  if ~isempty(e0) && (size(e0, 1) ~= r || size(e0, 2) ~= he)
    error('e0 must be r by h - min(max(p,s-1)-q,0)');
  end
  m = h - t0;
  A = makecell(A);
  B = makecell(B);
  E = zeros(r*n, M);
  Eall = zeros(r*(n - first_shock_time), M);
  X = zeros(r*n, M);
  if isempty(e0)
    Eall(1:r*he,:) = startup_shocks(A, B, Cc, z, Sig, X0, p, q, s, r, h, ...
      t0, first_shock_time, first_active_shock, M, rng);
  else
    Eall(1:r*he,:) = repmat(e0(:), 1, M);
  end
  for j = 1:M
    X0j = X0rep(X0, j);
    X(1:r*h,j) = X0j(:);
  end
  for j = 1:M
    zj = zrep(z, j);
    if n > h
      I = r*(h - first_shock_time) + 1:r*(n - first_shock_time);
      Eall(I,j) = reshape(randnm(n - h, Sig, "T", rng), r*(n - h), 1);
    end
    for t = h+1:n
      Xt = Eall(r*(t - first_shock_time - 1)+(1:r), j);
      for i = 1:p
        Xt = Xt + A{i}*X(r*(t-i-1)+(1:r), j);
      end
      for i = 1:q
        Xt = Xt + B{i}*Eall(r*(t - i - first_shock_time - 1)+(1:r), j);
      end
      for i = 1:s
        Xt = Xt + Cc{i}*zj(t-i+1);
      end
      X(r*(t-1)+(1:r), j) = Xt;
    end
  end
  E(:,:) = Eall(r*(-first_shock_time)+1:r*(n - first_shock_time),:);
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

function E0 = startup_shocks(A, B, C, z, Sig, X0, p, q, s, r, h, t0, ...
                             first_shock_time, first_active_shock, M, rng)
  he = h - first_shock_time;
  E0 = zeros(r*he, M);
  nPrefix = max(first_active_shock, 0);
  if nPrefix > 0
    for j = 1:M
      E0(1:r*nPrefix,j) = reshape(randnm(nPrefix, Sig, "T", rng), r*nPrefix, 1);
    end
  end
  activeOffset = r*(first_active_shock - first_shock_time);
  if q == 0
    for j = 1:M
      E0(activeOffset+1:end,j) = startup_residual(A, C, zrep(z, j), X0rep(X0, j), ...
        p, s, r, h, t0);
    end
    return
  end
  m = h - t0;
  na = h - first_active_shock;
  H = build_H(B, r, q, h, t0, first_active_shock);
  HD = postmultiply_sig(H, Sig, na);
  Wlag = find_W(cell2mat(B), Sig);
  W = W_build(Wlag, m);
  K = W\HD;
  R = blockdiag_sig(Sig, na) - HD'*K;
  for j = 1:M
    rvec = startup_residual(A, C, zrep(z, j), X0rep(X0, j), p, s, r, h, t0);
    Ehat = HD'*(W\rvec);
    E0(activeOffset+1:end,j) = Ehat + randnm(1, R, "T", rng);
  end
end

function H = build_H(B, r, q, h, t0, first_shock_time)
  m = h - t0;
  ne = h - first_shock_time;
  H = zeros(r*m, r*ne);
  for t = t0:h-1
    I = r*(t - t0) + (1:r);
    for ell = first_shock_time:h-1
      lag = t - ell;
      J = r*(ell - first_shock_time) + (1:r);
      if lag == 0
        H(I,J) = eye(r);
      elseif 1 <= lag && lag <= q
        H(I,J) = B{lag};
      end
    end
  end
end

function HD = postmultiply_sig(H, Sig, nblocks)
  r = size(Sig, 1);
  HD = H;
  for j = 1:nblocks
    J = r*(j-1) + (1:r);
    HD(:,J) = H(:,J)*Sig;
  end
end

function D = blockdiag_sig(Sig, nblocks)
  r = size(Sig, 1);
  D = zeros(r*nblocks);
  for j = 1:nblocks
    J = r*(j-1) + (1:r);
    D(J,J) = Sig;
  end
end

function rvec = startup_residual(A, C, z, X0, p, s, r, h, t0)
  m = h - t0;
  rvec = zeros(r*m, 1);
  for t = t0:h-1
    rt = X0(:, t+1);
    for i = 1:p
      rt = rt - A{i}*X0(:, t+1-i);
    end
    for i = 1:s
      rt = rt - C{i}*z(t+2-i);
    end
    rvec(r*(t-t0)+(1:r)) = rt;
  end
end

function zj = zrep(z, j)
  if size(z, 1) == 1
    zj = z(1,:);
  else
    zj = z(j,:);
  end
end

function X0j = X0rep(X0, j)
  if ndims(X0) < 3
    X0j = X0;
  else
    X0j = X0(:,:,j);
  end
end

function x = randnm(n, Sig, transpose, rng)
  TRANSP = nargin > 2 && startsWith(transpose, {'t', 'T'});
  r = size(Sig, 1);
  if n == 0, x = zeros(0, r); return, end
  if nargin < 4 || isempty(rng)
    global varmapack_rng
    rng = varmapack_rng;
  end
  if isempty(rng), error('Randompack RNG is not initialized'); end
  if TRANSP
    x = randompack_mvn(rng, Sig, n, zeros(r, 1));
  else
    x = randompack_mvn(rng, Sig, n, zeros(r, 1))';
  end
end

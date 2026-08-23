%REF_VARMA_SIM  Reference simulation of a VARMA time series model.
%
%  [X,E] = REF_VARMA_SIM(A,B,Sig,mu,n,M,X0,rng,e0) simulates M independent
%  VARMA(p,q) series of length n:
%
%    x(t) - mu(t) = sum Ai*(x(t-i) - mu(t-i)) + eps(t) + sum Bi*eps(t-i),
%
%  where eps(t) is N(0,Sig). A and B have r rows and contain the coefficient
%  blocks [A1,...,Ap] and [B1,...,Bq]; Sig is r by r. Use [] for A or B to
%  select a pure VMA or VAR model, respectively.
%
%  mu is [] or an r by nmu mean path. A scalar is expanded to all r components,
%  and the last supplied column repeats through time n. M defaults to 1.
%
%  With X0=[], a stationary model is initialized by drawing the startup states
%  and shocks from their exact joint distribution, so no burn-in is needed. X0
%  may instead be an r by nX0 common fixed startup path or an r by nX0 by M
%  array of replicate-specific paths, where max(p,q) <= nX0 <= n. For
%  stationary models, the resulting paths have the exact conditional
%  distribution given X0. A nonstationary model requires X0. With MA terms,
%  startup shocks are drawn conditionally on the residual equations implied by
%  X0; Sig must then be positive definite.
%
%  rng is a Randompack RNG handle. e0 is for the AgainstMatlab comparison tests
%  only: an r by nE0 fixed startup shock path used for exact comparisons with
%  REF_VARMA_SIMX. It must have max(p,q) <= nE0 <= n columns and, with X0,
%  nE0 must equal nX0.
%
%  For r > 1, X and E have shape r by n by M. Scalar models return 1 by n for
%  M=1 and n by M otherwise. E is returned only when requested.

function [X, E] = ref_varma_sim(A, B, Sig, mu, n, M, X0, rng, e0)
  r = size(Sig, 1);
  if isempty(A), A = zeros(r,0); end
  if isempty(B), B = zeros(r,0); end
  [~, p, q, h] = getdims(A, B, Sig);
  Aflp = flipmat(A);
  Bflp = flipmat(B);
  if nargin < 4, mu = []; end
  if nargin < 5 || isempty(n), error('n must be specified'); end
  if n<h, error('Too short series'); end
  if isempty(mu)
    mu = zeros(r,0);
  elseif isscalar(mu)
    mu = repmat(mu, r, 1);
  elseif isvector(mu) && r > 1 && numel(mu) == r
    mu = mu(:);
  elseif size(mu, 1) ~= r
    error('mu must have r rows')
  end
  if size(mu, 2) > n, error('mu cannot have more columns than n'); end
  if nargin < 6 || isempty(M), M=1; end
  if nargin < 7, X0 = []; end
  if nargin < 8, rng = []; end
  if nargin < 9, e0 = []; end
  if isempty(X0) && ref_varma_specrad(A) >= 1
    error("Cannot run varma_sim with unspecified X0 and rho(A) ≥ 1");
  end
  returnE = nargout >= 2;
  rollingE = ~returnE;

  % Check size of provided start vectors, set h to their size if ok
  if ~isempty(X0)
    nX0 = size(X0, 2);
    assert(size(X0, 1) == r)
    assert(ndims(X0) <= 3 && (size(X0, 3) == 1 || size(X0, 3) == M))
    assert(h <= nX0 && nX0 <= n)
    h = nX0;
  end
  if ~isempty(e0)
    ne0 = size(e0, 2);
    assert(size(e0, 1) == r)
    assert(max(p, q) <= ne0 && ne0 <= n)
    if ~isempty(X0), assert(nX0 == ne0); end
    h = ne0;
  end
  if ~isempty(X0) && q > 0 && ref_varma_specrad(A) >= 1
    if ~isempty(e0)
      error('e0 is not supported for nonstationary models with MA terms');
    end
    Mu0 = meanpath(mu, r, h);
    if ndims(X0) < 3
      X0bar = X0 - Mu0;
    else
      X0bar = X0 - repmat(Mu0, [1, 1, M]);
    end
    if returnE
      [X, E] = ref_varma_simx(A, B, zeros(r, 0), zeros(0, n), Sig, n, M, ...
                               X0bar, h, rng);
    else
      X = ref_varma_simx(A, B, zeros(r, 0), zeros(0, n), Sig, n, M, ...
                          X0bar, h, rng);
      E = [];
    end
    Mu = meanpath(mu, r, n);
    if r==1 && M==1
      X = reshape(X, 1, n) + Mu;
    elseif r==1
      X = reshape(X, n, M) + repmat(Mu', 1, M);
    else
      X = reshape(X, r, n, M) + repmat(Mu, [1, 1, M]);
    end
    return
  end
  e0 = e0(:);

  [~, G] = find_CG(A, B, Sig);
  PLU = vyw_factorize(A);
  assert(isempty(PLU) || isempty(PLU{1}) || PLU{1}(1) ~= 0)  % vyw_factorize ok
  S = vyw_solve(A, PLU, G);

  SS = S_build(S, A, G, h);
  if rollingE
    E = zeros(r*(h + 1), M);
  else
    E = zeros(r*n, M);
  end

  % Build theoretical covariance of xt
  if isempty(X0)  % Generate x{1:h} or draw x{1:h}|eps{1:h}
    if isempty(e0)
      for j = 1:M
        E(1:r*h,j) = reshape(randnm(h, Sig, "T", rng), r*h, 1);
      end
    else
      E(1:r*h,:) = repmat(e0, 1, M);
    end
    Psi = find_Psi(A, B, h);
    Psi_hat = find_Psi_hat(Psi, Sig);
    R = SS - Psi_hat*Psi_hat';
    e = Psi*E(1:r*h, :);
    Wrk = randnm(M, R, "T", rng);
    X1 = e + Wrk;
    h = h;
  elseif ~isempty(e0)  % Fixed x{1:h} and eps{1:h}
    X0bar = centered_startup(X0, mu, r, h, M);
    E(1:r*h,:) = repmat(e0, 1, M);
    X1 = X0bar;
  elseif q == 0 && ref_varma_specrad(A) >= 1  % Fixed-history pure AR path.
    X0bar = centered_startup(X0, mu, r, h, M);
    X1 = X0bar;
    R = Sig;
  else  % X0 given
    SS = S_build(S, A, G, h);
    C = find_C(A, B, Sig, h);
    CC = CC_build(A, C, h);
    LS = chol(SS, 'lower'); % TODO: Check this
    Chat = LS\CC;
    X0bar = centered_startup(X0, mu, r, h, M);
    e = Chat'*(LS\X0bar);
    R = -Chat'*Chat;
    J = 1:r;
    for j = 1:h
      R(J,J) = R(J,J) + Sig;
      J = J + r;
    end

    % Draw eps{1:h} and fill x{1:h}
    E(1:r*h, :) = e + randnm(M, R, "T", rng);
    X1 = X0bar;
  end
  X2 = zeros(n*r - h*r, M);
  X = [X1; X2];

  % Generate x{h+1:n}
  I = r*h + (1:r);
  J = r*(h-p)+1 : r*h;
  K = r*(h-q)+1 : r*h;
  if rollingE
    Ehist = zeros(r*q, M);
    for t = h+1:n
      slot = mod(t - 1, h + 1) + 1;
      iE = r*(slot - 1) + (1:r);
      E(iE,:) = randnm(M, Sig, "T", rng);
      X(I,:) = E(iE,:);
      if p > 0
        X(I,:) = X(I,:) + Aflp*X(J,:);
      end
      if q > 0
        for j = 1:q
          slotj = mod(t - q + j - 2, h + 1) + 1;
          Ehist(r*(j - 1) + (1:r), :) = E(r*(slotj - 1) + (1:r), :);
        end
        X(I,:) = X(I,:) + Bflp*Ehist;
      end
      I = I+r;
      J = J+r;
    end
    E = [];
  else
    for t = h+1:n
      E(I,:) = randnm(M, Sig, "T", rng);
      X(I,:) = E(I,:);
      if p > 0
        X(I,:) = X(I,:) + Aflp*X(J,:);
      end
      if q > 0
        X(I,:) = X(I,:) + Bflp*E(K,:);
      end
      I = I+r;
      J = J+r;
      K = K+r;
    end
  end

  % Reshape as appropriate for ARMA or VARMA
  Mu = meanpath(mu, r, n);
  if r==1 && M==1  %  one ARMA sequence:
    X = reshape(X,1,n) + Mu;
    if returnE, E = reshape(E,1,n); end
  elseif r==1      %  several ARMA sequences:
    X = reshape(X,n,M) + repmat(Mu', 1, M);
    if returnE, E = reshape(E,n,M); end
  else             %   one or more VARMA sequences in r×n×M array:
    X = reshape(X,r,n,M) + repmat(Mu, [1,1,M]);
    if returnE, E = reshape(E,r,n,M); end
  end
end

function X0bar = centered_startup(X0, mu, r, h, M)
  Mu = meanpath(mu, r, h);
  if size(X0, 3) == 1
    X0bar = repmat(reshape(X0(:, :, 1) - Mu, r*h, 1), 1, M);
  else
    X0bar = reshape(X0 - repmat(Mu, [1, 1, M]), r*h, M);
  end
end

function Mu = meanpath(mu, r, n)
  if isempty(mu)
    Mu = zeros(r, n);
    return
  end
  nmu = size(mu, 2);
  Mu = zeros(r, n);
  if nmu >= n
    Mu = mu(:, 1:n);
    return
  end
  Mu(:, 1:nmu) = mu;
  Mu(:, nmu+1:n) = repmat(mu(:, nmu), 1, n - nmu);
end

function x = randnm(n, Sig, transpose, rng)
  % Multivariate normal random vectors. In this new (2025) version the
  % eigendecomposition of Sig = U·Lam·U' is used to find the linear
  % transformation applied to the independent random variates when its Cholesky
  % decompositions fails due to it being only postivesemidefinite. The former
  % version added small numbers to its diagonal instead. See e.g. the R function
  % MASS::mvrnorm. To have this happen use ref_varma_testcase{ 5, 6 or 11. Another new
  % feature is the transpose parameter, to allow compatibility with C.
  TRANSP = nargin > 2 && startsWith(transpose, {'t', 'T'});
  r = size(Sig, 1);
  if n==0, x = zeros(0, r); return, end
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

function Psi_hat = find_Psi_hat(Psi, Sig)
  LSig = chol(Sig, 'lower');
  r = size(Sig, 1);
  h = size(Psi, 1)/r;
  I = 1:r;
  Psi_hat = zeros(r*h, r*h);
  for k = 1:h
    Psi_hat(:,I) = Psi(:,I)*LSig;
    I = I+r;
  end
end

% NEW_VARMA_SIM  Simulate an ARMA or a VARMA time series model
% 
%   ARMA MODELS:
%     X = NEW_VARMA_SIM(A,B,sigma,n) with scalar sigma generates a zero-mean
%     time series of length n using the ARMA(p,q) model:
% 
%              x(t) = A1·x(t-1) + ... + Ap·x(t-p) + y(t)                    (1a)
%     where:
%              y(t) = eps(t) + B1·eps(t-1) + ... + Bq·eps(t-q)              (1b)
% 
%     and eps(t) is N(0,sigma). A should be the row vector [A1,...,Ap] and B
%     should be [B1,...,Bq]. X returns the series as a row vector. To start the
%     simulation, (x(1),...,x(h),eps(1),...,eps(h)), with h=max(p,q), are drawn
%     from N(0,SCE) where SCE = [SS CC; CC' EE] is obtained by solving the
%     modified Yule-Walker equations (see vyw_factorize, vyw_solve and S_build),
%     so there is no need to throw away the first values to avoid spin-up
%     effects. The model must be stationary. 
% 
%     X = NEW_VARMA_SIM(A,B,sigma,n,mu) returns a time series with mean mu using
%     y(t) as above, and x(t) given by:
% 
%             x(t) - mu = A1·(x(t-1) - mu) + ... + Ap·(x(t-p) - mu) + y(t)   (2)
% 
%     In this case x and eps are drawn from N([mu; 0], SCE).
% 
%     X = NEW_VARMA_SIM(A,B,sigma,n,mu,M) generates M sequences simultaneously
%     and returns the i-th one in the i-th column of X. Use mu=[] to obtain
%     zero-mean.
% 
%     X = NEW_VARMA_SIM(A,B,sigma,n,mu,M,X0) sets x(1)...x(h) to X0(end-h+1:end)
%     - mu and draws eps(1)...eps(h) from the conditional distribution of
%     eps(1), ..., eps(h) given x(1),...,x(h). For this option the model need
%     not be stationary.
% 
%     [X,EPS] = NEW_VARMA_SIM(...) returns also the shock series, eps(t), in the
%     t-th row of EPS.
% 
%   VARMA MODELS:
%     X = NEW_VARMA_SIM(A,B,Cov,n), where Cov is an r×r matrix uses a VARMA(p,q)
%     model given by (1), but now x(t) is r-dimensional, eps(t) is r-variate
%     normal with mean 0 and covariance Cov, and Cov and the Ai's and Bi's are
%     r×r matrices. A and B should contain A = [A1 A2...Ap] (an r × r·p matrix)
%     and B = [B1 B2 ... Bq] (an r × r·q matrix). The r×n matrix X returns x(t)
%     in its t-th column. As in the scalar case the simulated series is spin-up
%     free, the starting values x(1),...,x(h) and eps(1),...,eps(h) being drawn
%     from the correct distribution. The model must be stationary.
%  
%     X = NEW_VARMA_SIM(...,mu) uses (2) for x(t) instead of (1a).
% 
%     X = NEW_VARMA_SIM(...,mu,M) returns M sequences in an r×n×M
%     multidimensional X (when r > 1) with the i-th sequence in X(1:r, 1:n, i).
%     Use mu = [] to obtain zero-mean.
% 
%     X = NEW_VARMA_SIM(...,mu,M,X0) initializes x(1),...,x(h) with the last h
%     columns of X0 instead of drawing from N(mu,SS). The model need not be
%     stationary.
% 
%     [X,EPS] = NEW_VARMA_SIM(...) returns also the shock series; the i-th
%     eps(t) is returned in EPS(:, t, i).
% 
%   For both ARMA and VARMA, use NEW_VARMA_SIM(A,[],...) for a pure
%   autoregressive model, and NEW_VARMA_SIM([],B,...) for a pure moving average
%   model.
% 
%   The method used is described in [3]. It is an improved version of varma_sim
%   in the original Vauto package aka ACM TOMS Algorithm 878 as described in [1]
%   and [2]. The primary difference is that in this new function the first h
%   shocks (eps_t) are drawn first and afterwards the first h states (x_t) are
%   drawn using the conditional distribution of x_t|eps_t, whereas the original
%   version did the opposite (eps_t after x_t).
%
%   [1] K Jonasson and SE Ferrando 2008. Evaluating exact VARMA likelihood and
%       its gradient when data are incomplete. ACM Trans. Math. Softw. (TOMS),
%       35(1)
%
%   [2] K Jonasson 2008. Algorithm 878: Exact VARMA likelihood and its gradient
%       for complete and incomplete data with Matlab. ACM Trans. Math. Softw.
%       (TOMS), 35(1)
%
%   [3] K Jonasson 2025. Burn-in free simulation of VARMA time series.
%       Manuscript in preparation.
%
%   (C) Kristján Jónasson, Dept. of Computer Science, University of Iceland,
%   2025. jonasson@hi.is.

function [X, E, condR] = ref_varma_sim(A, B, Sig, n, mu, M, x0, rng)
  r = size(Sig, 1);
  if isempty(A), A = zeros(r,0); end
  if isempty(B), B = zeros(r,0); end
  [~, p, q, h] = getdims(A, B, Sig);
  Aflp = flipmat(A);
  Bflp = flipmat(B);
  if n<h, error('Too short series'); end
  if nargin < 5 || isempty(mu)
    mu = zeros(r,1); 
  else
    if isscalar(mu), mu = repmat(mu, r, 1); end
    mu = mu(:); 
  end
  if nargin < 6 || isempty(M), M=1; end
  if nargin < 7, x0 = []; end
  if nargin < 8, rng = []; end
  if isempty(x0) && ref_varma_specrad(A) >= 1
    error("Cannot run varma_sim with unspecified x0 and rho(A) ≥ 1");
  end
  returnE = nargout >= 2;
  rollingE = ~returnE;
  [~, G] = find_CG(A, B, Sig);
  PLU = vyw_factorize(A);
  assert(isempty(PLU) || isempty(PLU{1}) || PLU{1}(1) ~= 0)  % vyw_factorize ok
  S = vyw_solve(A, PLU, G);

  % Check size of provided start vector, set h to its size if ok
  if ~isempty(x0)
    nx0 = size(x0, 2);
    x0 = x0(:);
    assert(h <= nx0 && nx0 <= n)
    h = nx0;
  end

  SS = S_build(S, A, G, h);
  if rollingE
    E = zeros(r*(h + 1), M);
  else
    E = zeros(r*n, M);
  end

  % Build theoretical covariance of xt
  if isempty(x0)  % Generate x{1:h}
    for j = 1:M
      E(1:r*h,j) = reshape(randnm(h, Sig, "T", rng), r*h, 1);
    end
    Psi = find_Psi(A, B);
    Psi_hat = find_Psi_hat(Psi, Sig);
    R = SS - Psi_hat*Psi_hat';
    e = Psi*E(1:r*h, :);
    Wrk = randnm(M, R, "T", rng);
    X1 = e + Wrk;
    h = h;
  elseif q == 0 && ref_varma_specrad(A) >= 1  % Fixed-history pure AR path.
    x0bar = x0(:) - repmat(mu, h, 1);
    X1 = repmat(x0bar,1,M);
    R = Sig;
  else  % x0 given
    SS = S_build(S, A, G, h);
    C = find_C(A, B, Sig, h);
    CC = CC_build(A, C, h);
    LS = chol(SS, 'lower'); % TODO: Check this
    Chat = LS\CC;
    x0bar = x0(:) - repmat(mu, h, 1);
    e = Chat'*(LS\x0bar);
    R = -Chat'*Chat;
    J = 1:r;
    for j = 1:h
      R(J,J) = R(J,J) + Sig;
      J = J + r;
    end

    % Draw eps{1:h} and fill x{1:h}
    E(1:r*h, :) = e + randnm(M, R, "T", rng);
    X1 = repmat(x0bar,1,M);
  end
  condR = cond(R);
  condSig = cond(Sig);
  rho = ref_varma_specrad(A);
  % fprintf("cond(R) = %.2e, cond(Sig) = %.2e, rho = %.4f\n", condR, condSig, rho);
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
    for j = 1:M
      E(r*h+1:r*n,j) = reshape(randnm(n - h, Sig, "T", rng), r*(n - h), 1);
    end
    for t = h+1:n
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
  if r==1 && M==1  %  one ARMA sequence:
    X = reshape(X,1,n) + mu;
    if returnE, E = reshape(E,1,n); end
  elseif r==1      %  several ARMA sequences:
    X = reshape(X,n,M) + mu;
    if returnE, E = reshape(E,n,M); end
  else             %   one or more VARMA sequences in r×n×M array:
    X = reshape(X,r,n,M) + repmat(mu,[1,n,M]);
    if returnE, E = reshape(E,r,n,M); end
  end
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

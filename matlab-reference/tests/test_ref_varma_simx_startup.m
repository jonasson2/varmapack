function zmax = test_ref_varma_simx_startup(cases, M, n, seed)
  fprintf('TESTING REF_VARMA_SIMX STARTUP DRAWS... ');
  if nargin < 1 || isempty(cases), cases = 1:ref_varma_testcasex('count'); end
  if nargin < 2 || isempty(M), M = 5000; end
  if nargin < 3 || isempty(n), n = 10; end
  if nargin < 4 || isempty(seed), seed = 2000; end
  rng = randompack_create();
  cleanup = onCleanup(@() randompack_free(rng));
  zmax = 0;
  for icase = cases
    [A, B, C, Sig, z, p, q, s, r, name] = ref_varma_testcasex(icase, n);
    h = max([p, q, s]) + 2;
    if n < h, error('n must be at least h'); end
    StartX0 = startup_values(r, h, icase);
    [muE, covE] = startup_theory(A, B, C, Sig, z, p, q, s, r, h, StartX0);
    randompack_seed(rng, seed + icase);
    [X, E] = ref_varma_simx(A, B, C, z, Sig, n, M, StartX0, h, [], rng);
    E0 = reshape_startup_shocks(E, r, h, M);
    X0drawn = reshape_startup_values(X, r, h, M);
    ascertain(almostequal(X0drawn, repmat(StartX0(:), 1, M), 1e-12), ...
              ['simx startup X mismatch in ' name]);
    z = max(startup_zscores(E0, muE, covE, M));
    zmax = max(zmax, z);
  end
  ascertain(zmax < 3.5, sprintf('simx startup draw z-score %.2f is too high', zmax));
  fprintf('OK\n');
end

function [muE, covE] = startup_theory(A, B, C, Sig, z, p, q, s, r, h, X0)
  t0 = max(p, s - 1);
  a = min(t0 - q, 0);
  mE = h - a;
  mR = h - t0;
  H = zeros(r*mR, r*mE);
  rho = zeros(r*mR, 1);
  for t = t0:h-1
    I = r*(t - t0) + (1:r);
    rt = X0(:, t+1);
    for i = 1:p
      rt = rt - A(:, r*(i-1)+(1:r))*X0(:, t-i+1);
    end
    for i = 1:s
      rt = rt - C(:, i)*z(t-i+2);
    end
    rho(I) = rt;
    for i = 0:q
      J = r*(t - i - a) + (1:r);
      if i == 0
        H(I,J) = eye(r);
      else
        H(I,J) = B(:, r*(i-1)+(1:r));
      end
    end
  end
  D = kron(eye(mE), Sig);
  W = H*D*H';
  K = D*H'/W;
  muStar = K*rho;
  covStar = D - K*H*D;
  I0 = r*(-a) + (1:r*h);
  muE = muStar(I0);
  covE = covStar(I0, I0);
end

function X0 = startup_values(r, h, icase)
  X0 = zeros(r, h);
  for t = 0:h-1
    X0(:, t+1) = 0.07*icase + 0.20*(1:r)' - 0.11*t;
  end
end

function E0 = reshape_startup_shocks(E, r, h, M)
  if r == 1
    E = reshape(E, 1, [], M);
  end
  E0 = reshape(E(:, 1:h, :), r*h, M);
end

function X0 = reshape_startup_values(X, r, h, M)
  if r == 1
    X = reshape(X, 1, [], M);
  end
  X0 = reshape(X(:, 1:h, :), r*h, M);
end

function z = startup_zscores(E0, muE, covE, M)
  mud = mean(E0, 2);
  covd = cov(E0');
  zmean = vector_mean_zscore(mud, muE, covE, M);
  zcov = covariance_zscore(covd, covE, M);
  z = [zmean zcov];
end

function z = vector_mean_zscore(mud, mu, Sigma, M)
  v = diag(Sigma);
  keep = v > zero_cov_tolerance(Sigma);
  se = sqrt(v(keep)/M);
  z = max(abs(mud(keep) - mu(keep))./se(keep));
end

function z = covariance_zscore(Shat, S, M)
  v = diag(S);
  SE = sqrt((v*v' + S.^2)/(M - 1));
  tol = zero_cov_tolerance(S);
  keep = v*v' > tol^2;
  Z = abs(Shat(keep) - S(keep))./SE(keep);
  z = max(Z(:));
end

function tol = zero_cov_tolerance(S)
  tol = eps(max(1, max(abs(S(:)))));
end

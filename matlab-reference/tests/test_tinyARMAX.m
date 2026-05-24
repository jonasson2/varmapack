function maxAbsRel = test_tinyARMAX(varargin)
  % Moment formulas used here are derived in ../docs/tiny-armax-moments.txt.
  showTable = false;
  if nargin >= 1 && (ischar(varargin{1}) || isstring(varargin{1})) && ...
      strcmpi(varargin{1}, 'table')
    showTable = true;
    varargin = varargin(2:end);
  end
  if length(varargin) >= 1, M = varargin{1}; else, M = []; end
  if isempty(M), M = 200; end
  fprintf('TESTING TINY ARMAX REFERENCE SIMULATION... ');
  rng = randompack_create();
  cleanup = onCleanup(@() randompack_free(rng));
  n = 10;
  h = 2;
  [A, B, C, Sig, z, p, q, s, r, name] = ref_varma_testcasex(1, n);
  ascertain(strcmp(name, 'tinyARMAX'));
  randompack_seed(rng, 41);
  maxAbsRel = check_simx_sim_agreement(A, B, C, z, Sig, n, h, rng);
  randompack_seed(rng, 42);
  check_stationary_start(A, B, C, z, Sig, n, h, M, rng, showTable);
  randompack_seed(rng, 43);
  check_fixed_zero_start(A, B, C, z, Sig, n, h, max(50*M, M), rng, showTable);
  fprintf('OK\n');
end

function maxAbsRel = check_simx_sim_agreement(A, B, C, z, Sig, n, h, rng)
  M = 7;
  mu = zeros(1, n);
  mu(1) = 0.3;
  for t = 2:n
    mu(t) = A*mu(t-1) + C*z(t);
  end
  y0 = [0.2 -0.4];
  x0 = y0 + mu(1:h);
  e0 = [-0.7 0.4];
  randompack_seed(rng, 41);
  [Xx, Ex] = ref_varma_simx(A, B, C, z, Sig, n, M, x0, h, e0, rng);
  randompack_seed(rng, 41);
  [Xs, Es] = ref_varma_sim(A, B, Sig, mu, n, M, x0, e0, rng);
  maxX = absrel_difference(Xx, Xs);
  maxE = absrel_difference(Ex, Es);
  maxAbsRel = max(maxX, maxE);
  ascertain(maxX < 1e-12, 'tinyARMAX simx/sim X mismatch');
  ascertain(maxE < 1e-12, 'tinyARMAX simx/sim E mismatch');
end

function d = absrel_difference(x, y)
  absd = abs(x(:) - y(:));
  scale = max(abs(x(:)), abs(y(:)));
  reld = absd./scale;
  reld(scale == 0) = absd(scale == 0);
  d = max(min(absd, reld));
end

function check_stationary_start(A, B, C, z, Sig, n, h, M, rng, showTable)
  mu0 = 0.5*(-1).^(0:h-1);
  V0 = [2 1.4; 1.4 2];
  X = zeros(M, n);
  for j = 1:M
    x0 = randompack_mvn(rng, V0, 1, mu0')';
    X(j, :) = ref_varma_simx(A, B, C, z, Sig, n, 1, x0, h, [], rng);
  end
  mu = 0.5*(-1).^(0:n-1);
  V = 2*ones(1, n);
  C1 = 1.4*ones(1, n-1);
  check_moments(X, mu, V, C1, showTable, 'A');
end

function check_fixed_zero_start(A, B, C, z, Sig, n, h, M, rng, showTable)
  x0 = zeros(1, h);
  X = ref_varma_simx(A, B, C, z, Sig, n, M, x0, h, [], rng)';
  mu = [0 0 62/65 zeros(1, n-3)];
  for t = 4:n
    mu(t) = 0.6*mu(t-1) + 0.8*(-1)^(t-1);
  end
  V = zeros(1, n);
  V(3:n) = 2 - 649/650*0.36.^(0:n-3);
  C1 = zeros(1, n-1);
  C1(3:n-1) = 1.4 - 1947/3250*0.36.^(0:n-4);
  check_moments(X, mu, V, C1, showTable, 'B');
end

function check_moments(X, mu, V, C1, showTable, label)
  M = size(X, 1);
  mud = mean(X, 1);
  Vd = var(X, 0, 1);
  C1d = zeros(1, length(C1));
  for t = 1:length(C1)
    x = X(:, t) - mud(t);
    y = X(:, t+1) - mud(t+1);
    C1d(t) = sum(x.*y)/(M - 1);
  end
  Z = moment_zscores(mu, mud, V, Vd, C1, C1d, M);
  if showTable
    print_table(label, M, mu, mud, V, Vd, C1, C1d, Z);
  end
  ascertain(max(Z(:)) < 3, 'tinyARMAX z-score check failed');
end

function Z = moment_zscores(mu, mud, V, Vd, C1, C1d, M)
  Z = zeros(length(mu), 3);
  nz = V > 0;
  Z(nz, 1) = abs(mud(nz) - mu(nz))./sqrt(V(nz)/M);
  Z(nz, 2) = abs(Vd(nz) - V(nz))./sqrt(2*V(nz).^2/(M - 1));
  for t = 1:length(C1)
    if V(t) > 0 && V(t+1) > 0
      se = sqrt((V(t)*V(t+1) + C1(t)^2)/(M - 1));
      Z(t, 3) = abs(C1d(t) - C1(t))/se;
    end
  end
end

function print_table(label, M, mu, mud, V, Vd, C1, C1d, Z)
  fprintf('\nCase (%s), M = %d, max-z = %.2f\n', label, M, max(Z(:)));
  fprintf('  t     mean obs-mean       var  obs-var      cov1 obs-cov1    max-z\n');
  for t = 1:length(mu)
    if t < length(mu)
      fprintf('%3d %8.3f %8.3f  %8.3f %8.3f  %8.3f %8.3f  %7.2f\n', ...
              t-1, mu(t), mud(t), V(t), Vd(t), C1(t), C1d(t), max(Z(t, :)));
    else
      fprintf('%3d %8.3f %8.3f  %8.3f %8.3f                     %7.2f\n', ...
              t-1, mu(t), mud(t), V(t), Vd(t), max(Z(t, :)));
    end
  end
end

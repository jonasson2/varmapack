function maxAbsRel = test_ref_varma_simx_forward(cases, n, M, tol)
  fprintf('TESTING REF_VARMA_SIMX FORWARD RECURSION... ');
  if nargin < 1 || isempty(cases), cases = 1:ref_varma_testcasex('count'); end
  if nargin < 2 || isempty(n), n = 10; end
  if nargin < 3 || isempty(M), M = 7; end
  if nargin < 4 || isempty(tol), tol = 5*eps; end
  rng = randompack_create();
  cleanup = onCleanup(@() randompack_free(rng));
  maxAbsRel = 0;
  for icase = cases
    [A, B, C, Sig, z, p, q, s, r, name] = ref_varma_testcasex(icase, n);
    h = max([p, q, s]) + 2;
    if n < h, error('n must be at least h'); end
    [mu, x0, e0x, e0] = fixed_start(A, B, C, z, p, q, s, r, h, n, icase);
    randompack_seed(rng, 2000 + icase);
    [Xx, Ex] = ref_varma_simx(A, B, C, z, Sig, n, M, x0, h, e0x, rng);
    randompack_seed(rng, 2000 + icase);
    [Xs, Es] = ref_varma_sim(A, B, Sig, mu, n, M, x0, e0, rng);
    dx = absrel_difference(Xx, Xs);
    de = absrel_difference(Ex, Es);
    d = max(dx, de);
    maxAbsRel = max(maxAbsRel, d);
    ascertain(d < tol, ['simx/sim forward mismatch in ' name]);
  end
  fprintf('OK\n');
end

function [mu, x0, e0x, e0] = fixed_start(A, B, C, z, p, q, s, r, h, n, icase)
  first_shock_time = min(max(p, s - 1) - q, 0);
  he = h - first_shock_time;
  mu = zeros(r, n);
  for t = 0:h-1
    mu(:, t+1) = 0.04*icase + 0.03*(1:r)' - 0.02*t;
  end
  for t = h:n-1
    mt = zeros(r, 1);
    for i = 1:p
      mt = mt + A(:, r*(i-1)+(1:r))*mu(:, t-i+1);
    end
    for i = 1:s
      mt = mt + C(:, i)*z(t-i+2);
    end
    mu(:, t+1) = mt;
  end
  y0 = zeros(r, h);
  for t = 0:h-1
    y0(:, t+1) = 0.07*icase - 0.05*(1:r)' + 0.01*t;
  end
  x0 = y0 + mu(:, 1:h);
  e0x = zeros(r, he);
  for t = first_shock_time:h-1
    e0x(:, t - first_shock_time + 1) = 0.11*icase + 0.02*(1:r)' - 0.03*t;
  end
  e0 = e0x(:, 1 - first_shock_time:end);
end

function d = absrel_difference(x, y)
  absd = abs(x(:) - y(:));
  scale = max(abs(x(:)), abs(y(:)));
  reld = absd./scale;
  reld(scale == 0) = absd(scale == 0);
  d = max(min(absd, reld));
end

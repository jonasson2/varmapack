function test_ref_varma_sim(cases, M, n, showTable)
  fprintf('TESTING REF_VARMA_SIM... ');
  rng = randompack_create();
  cleanup = onCleanup(@() randompack_free(rng));
  D = zeros(length(cases), 4);
  Z = zeros(length(cases), 8);
  if showTable
    fprintf('\nPercentage difference:\n');
    fprintf(' case      mean    startup      first       last     max z\n');
  end
  for i = cases
    k = find(cases == i, 1);
    [A, B, Sig, p, q, r] = ref_varma_testcase(i);
    h = max(p, q);
    mu = (1:r)'/10;
    % Smoke test: generated states and shocks should be finite.
    randompack_seed(rng, 42);
    [X, E] = ref_varma_sim(A, B, Sig, [], h + 4, 3, [], [], rng);
    ascertain(all(isfinite(X(:))));
    ascertain(all(isfinite(E(:))));
    % Constant mean: using the same shocks should only shift X by mu.
    randompack_seed(rng, 42);
    [Xmu, Emu] = ref_varma_sim(A, B, Sig, mu, h + 4, 3, [], [], rng);
    ascertain(almostequal(Emu, E));
    % Mean path: supplied means should shift X while leaving shocks unchanged.
    mupath = [mu, 2*mu, 3*mu];
    Mupath = [mupath, repmat(mupath(:, end), 1, h + 1)];
    randompack_seed(rng, 42);
    [Xpath, Epath] = ref_varma_sim(A, B, Sig, mupath, h + 4, 3, [], [], rng);
    ascertain(almostequal(Epath, E));
    if r == 1
      X = reshape(X, h + 4, 3);
      Xmu = reshape(Xmu, h + 4, 3);
      ascertain(almostequal(Xmu - X, repmat(mu, h + 4, 3), 1e-12));
      Xpath = reshape(Xpath, h + 4, 3);
      ascertain(almostequal(Xpath - X, repmat(Mupath', 1, 3), 1e-12));
    else
      ascertain(almostequal(Xmu - X, repmat(mu, [1, h + 4, 3]), 1e-12));
      ascertain(almostequal(Xpath - X, repmat(Mupath, [1, 1, 3]), 1e-12));
    end
    % Fixed startup: supplied X0 should be copied into every replicate.
    X0 = 2 + repmat((1:r)', 1, h);
    randompack_seed(rng, 42);
    Xstart = ref_varma_sim(A, B, Sig, mu, h + 3, 2, X0, [], rng);
    if r == 1
      Xstart = reshape(Xstart, 1, h + 3, 2);
    end
    ascertain(almostequal(Xstart(:, 1:h, 1), X0, 1e-12));
    ascertain(almostequal(Xstart(:, 1:h, 2), X0, 1e-12));
    [D(k, :), Z(k, :)] = check_statistics(A, B, Sig, p, q, r, n, M, rng);
    if showTable
      fprintf('%5d %9.2f %10.2f %10.2f %10.2f %9.2f\n', ...
              i, 100*D(k, :), max(Z(k, :)));
    end
  end
  if showTable
    fprintf('  max %9.2f %10.2f %10.2f %10.2f %9.2f\n', ...
            100*max(D, [], 1), max(Z(:)));
  end
  ascertain(max(Z(:)) < 4);
  fprintf('OK\n');
end

function [d, z] = check_statistics(A, B, Sig, p, q, r, n, M, rng)
  h = max(p, q);
  mu = (1:r)';
  nCov = h + 3;
  d = zeros(1, 4);
  z = zeros(1, 4);
  [C, G] = find_CG(A, B, Sig);
  PLU = vyw_factorize(A);
  S = vyw_solve(A, PLU, G);
  % z(1): long-run sample mean agrees with the supplied mean.
  randompack_seed(rng, 42);
  X = ref_varma_sim(A, B, Sig, mu, n, M, [], [], rng);
  X = asarray3(X, r, n, M);
  mud = mean(mean(X, 3), 2);
  d(1) = relabsdiff(mu, mud);
  SSmean = S_build(S, A, G, n);
  z(1) = mean_zscore(mud, mu, SSmean, r, n, M);
  % z(2): covariance of the startup/early block agrees with S_build.
  randompack_seed(rng, 42);
  X = ref_varma_sim(A, B, Sig, [], nCov, M, [], [], rng);
  X = reshape(asarray3(X, r, nCov, M), r*nCov, M);
  SSd = cov(X');
  SS = S_build(S, A, G, nCov);
  d(2) = relabsdiff(SS, SSd);
  z(2) = cov_zscore(SSd, SS, M);
  % z(5): covariance of the first two values agrees with S_build.
  nLag = max(h, 2);
  randompack_seed(rng, 42);
  X = ref_varma_sim(A, B, Sig, [], nLag, M, [], [], rng);
  X = asarray3(X, r, nLag, M);
  X = reshape(X(:, 1:2, :), 2*r, M);
  SS2d = cov(X');
  SS2 = S_build(S, A, G, 2);
  z(5) = cov_zscore(SS2d, SS2, M);
  if h == 0, return, end
  % z(3) and z(4): state-shock cross-covariances agree at start and end.
  CC = CC_build(A, makecell(C), h);
  SSh = S_build(S, A, G, h);
  EE = kron(eye(h), Sig);
  randompack_seed(rng, 42);
  [X, E] = ref_varma_sim(A, B, Sig, [], 2*h, M, [], [], rng);
  X = asarray3(X, r, 2*h, M);
  E = asarray3(E, r, 2*h, M);
  X1 = reshape(X(:, 1:h, :), r*h, []);
  E1 = reshape(E(:, 1:h, :), r*h, []);
  X2 = reshape(X(:, end-h+1:end, :), r*h, []);
  E2 = reshape(E(:, end-h+1:end, :), r*h, []);
  CC1 = cov([X1' E1']);
  CC2 = cov([X2' E2']);
  CC1 = CC1(1:r*h, r*h+1:end);
  CC2 = CC2(1:r*h, r*h+1:end);
  d(3) = relabsdiff(CC, CC1);
  d(4) = relabsdiff(CC, CC2);
  z(3) = crosscov_zscore(CC1, CC, SSh, EE, M);
  z(4) = crosscov_zscore(CC2, CC, SSh, EE, M);
  % z(6:8): conditional one- and two-step forecasts from fixed X0 agree.
  [z(6), z(7), z(8)] = forecast_zscores(A, B, Sig, C, SSh, h, r, mu, M, rng);
end

function X = asarray3(X, r, n, M)
  if r == 1
    X = reshape(X, 1, n, M);
  end
end

function d = relabsdiff(a, b)
  d = norm(a(:) - b(:), inf)/max([1; abs(a(:)); abs(b(:))]);
end

function z = mean_zscore(mud, mu, SS, r, n, M)
  z = 0;
  for i = 1:r
    idx = i:r:r*n;
    se = sqrt(sum(sum(SS(idx, idx)))/(n*n*M));
    z = max(z, abs(mud(i) - mu(i))/se);
  end
end

function z = cov_zscore(Shat, S, M)
  v = diag(S);
  SE = sqrt((v*v' + S.^2)/(M - 1));
  Z = abs(Shat - S)./SE;
  z = max(Z(:));
end

function z = crosscov_zscore(CChat, CC, SS, EE, M)
  vx = diag(SS);
  ve = diag(EE);
  SE = sqrt((vx*ve' + CC.^2)/(M - 1));
  Z = abs(CChat - CC)./SE;
  z = max(Z(:));
end

function [z1, z2, z3] = forecast_zscores(A, B, Sig, C, SS, h, r, mu, M, rng)
  X0 = 2 + repmat((1:r)', 1, h);
  randompack_seed(rng, 43);
  X = ref_varma_sim(A, B, Sig, mu, h + 2, M, X0, [], rng);
  X = asarray3(X, r, h + 2, M);
  X1 = reshape(X(:, h+1, :), r, M);
  X2 = reshape(X(:, h+2, :), r, M);
  [EX1, EX2, VX1] = forecast_moments(A, B, Sig, C, SS, X0, mu);
  z1 = vector_mean_zscore(mean(X1, 2), EX1, VX1, M);
  VX2 = cov(X2');
  z2 = vector_mean_zscore(mean(X2, 2), EX2, VX2, M);
  z3 = cov_zscore(cov(X1'), VX1, M);
end

function [EX1, EX2, VX1] = forecast_moments(A, B, Sig, C, SS, X0, mu)
  r = size(Sig, 1);
  t = size(X0, 2) + 1;
  C = makecell(C);
  CC = CC_build(A, C, t - 1);
  A = makecell(A);
  B = makecell(B);
  X0 = X0 - repmat(mu, 1, t - 1);
  p = length(A);
  q = length(B);
  EX1 = zeros(r, 1);
  EX2 = zeros(r, 1);
  xp = reshape(X0, [], 1);
  for i = 1:p
    EX1 = EX1 + A{i}*X0(:, t - i);
  end
  Eeps = reshape(CC'*(SS\xp), r, t - 1);
  BB = zeros(r, 0);
  EE = [];
  for i = 1:q
    EX1 = EX1 + B{i}*Eeps(:, t - i);
    BB = [B{i} BB];
    EE = blkdiag(EE, Sig);
  end
  VX1 = CC'*(SS\CC);
  VX1 = Sig + BB*(EE - VX1(end-r*q+1:end, end-r*q+1:end))*BB';
  for i = 1:p
    if i == 1, EX2 = EX2 + A{i}*EX1; end
    if i > 1, EX2 = EX2 + A{i}*X0(:, t - i + 1); end
  end
  for i = 2:q
    EX2 = EX2 + B{i}*Eeps(:, t - i + 1);
  end
  EX1 = EX1 + mu;
  EX2 = EX2 + mu;
end

function z = vector_mean_zscore(mud, mu, S, M)
  z = max(abs(mud - mu)./sqrt(diag(S)/M));
end

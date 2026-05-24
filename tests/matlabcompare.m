function matlabcompare(cases, Ms, n)
  % MATLABCOMPARE writes tests/matlabcompare.txt, the reference data consumed by
  % TestAgainstMatlab.c. The text format is:
  %
  %   # comment
  %   var <name> <ndim> <dim1> ... <dimN>
  %   <values in Matlab column-major order>
  %
  % Scalars use rank 0, for example "var n 0"; vectors use rank 1, for example
  % "var cases 1 10"; matrices use rank 2; arrays use their natural rank.
  % Values are whitespace separated. See mat2file.m for the writer.
  here = fileparts(mfilename('fullpath'));
  addpath(fullfile(here, "matlab"));
  addpath(fullfile(here, "..", "matlab-reference"));
  if nargin < 3, n = 20; end
  if nargin < 2, Ms = 1; end
  if nargin < 1 || isempty(cases), cases = 1:ref_varma_testcase('count'); end
  fid = fopen(fullfile(here, "matlabcompare.txt"), "w");
  fidRolling = fopen(fullfile(here, "matlabcompare-rolling.txt"), "w");
  fprintf(fid, "# varmapack Matlab reference data v1\n");
  fprintf(fid, "# Format documented in tests/matlabcompare.m\n");
  fprintf(fidRolling, "# varmapack Matlab rolling-E reference data v1\n");
  fprintf(fidRolling, "# Format documented in tests/matlabcompare.m\n");
  mat2file(fid, length(Ms), "#Ms");
  mat2file(fid, length(cases), "#cases");
  mat2file(fid, n);
  mat2file(fid, Ms);
  mat2file(fid, cases);
  mat2file(fidRolling, length(Ms), "#Ms");
  mat2file(fidRolling, length(cases), "#cases");
  mat2file(fidRolling, n);
  mat2file(fidRolling, Ms);
  mat2file(fidRolling, cases);
  reducedCases = [3, 7, 12, 15];
  reducedM = 3;
  mat2file(fid, length(reducedCases), "#reducedCases");
  mat2file(fid, reducedM, "reducedM");
  mat2file(fid, reducedCases, "reducedCases");
  for k = cases
    icase = k + 1;
    [A, B, Sig, p, q, r] = ref_varma_testcase(k);    
    StartX0 = startmat(r, max(p,q));
    mu = (1:r)'/10;
    muPath = [mu, 2*mu, 3*mu];
    mat2file(fid, p, "p" + icase);
    mat2file(fid, q, "q" + icase);
    mat2file(fid, r, "r" + icase);
    mat2file(fid, A, "A" + icase);
    mat2file(fid, B, "B" + icase);
    mat2file(fid, Sig, "Sig" + icase);
    mat2file(fid, StartX0, "X0" + icase);
    mat2file(fid, mu, "mu" + icase);
    mat2file(fid, muPath, "muPath" + icase);
    mat2file(fidRolling, StartX0, "X0" + icase);
    mat2file(fidRolling, mu, "mu" + icase);
    if p == 0
      rho = 0;
    else
      rho = ref_varma_specrad(A);
    end
    maRho = ref_varma_ma_specrad(B);
    mat2file(fid, rho, "specrad-" + k);
    mat2file(fid, maRho, "maSpecrad-" + k);
    acvfMaxlag = max(n - 1, p);
    Gamma = ref_varma_acvf(A, B, Sig, acvfMaxlag);
    Psi = ref_varma_psi(A, B, acvfMaxlag);
    Theta = ref_varma_irf(A, B, Sig, acvfMaxlag);
    mat2file(fid, acvfMaxlag, "acvfMaxlag-" + k);
    mat2file(fid, Gamma, "Gamma-" + k);
    mat2file(fid, Psi, "Psi-" + k);
    mat2file(fid, Theta, "Theta-" + k);
    for M = Ms
      rng = randompack_create();
      cleanup = onCleanup(@() randompack_free(rng));
      randompack_seed(rng, 42);
      [X, E]        = ref_varma_sim(A, B, Sig, 0, n, M, [], [], rng);
      randompack_seed(rng, 42);
      [X0, E0]      = ref_varma_sim(A, B, Sig, 0, n, M, StartX0, [], rng);
      randompack_seed(rng, 42);
      [Xmu, Emu]    = ref_varma_sim(A, B, Sig, mu, n, M, [], [], rng);
      randompack_seed(rng, 42);
      [X0mu, E0mu]  = ref_varma_sim(A, B, Sig, mu, n, M, StartX0, [], rng);
      randompack_seed(rng, 42);
      [Xpath, Epath] = ref_varma_sim(A, B, Sig, muPath, n, M, [], [], rng);
      Cml = ref_varma_autocov(X(:, :, 1), n - 1, "ML");
      Cc = ref_varma_autocov(X(:, :, 1), n - 1, "C");
      mat2file(fid, Cml, "AutoML" + M + "-" + k);
      mat2file(fid, Cc, "AutoC" + M + "-" + k);
      mat2file(fid, X, "X" + M + "-" + k);
      mat2file(fid, E, "E" + M + "-" + k);
      mat2file(fid, X0, "X0" + M + "-" + k);
      mat2file(fid, E0, "E0" + M + "-" + k);
      mat2file(fid, Xmu, "Xmu" + M + "-" + k);
      mat2file(fid, Emu, "Emu" + M + "-" + k);
      mat2file(fid, X0mu, "X0mu" + M + "-" + k);
      mat2file(fid, E0mu, "E0mu" + M + "-" + k);
      mat2file(fid, Xpath, "Xpath" + M + "-" + k);
      mat2file(fid, Epath, "Epath" + M + "-" + k);
      randompack_seed(rng, 42);
      Xroll = ref_varma_sim(A, B, Sig, 0, n, M, [], [], rng);
      randompack_seed(rng, 42);
      X0roll = ref_varma_sim(A, B, Sig, 0, n, M, StartX0, [], rng);
      randompack_seed(rng, 42);
      Xmuroll = ref_varma_sim(A, B, Sig, mu, n, M, [], [], rng);
      randompack_seed(rng, 42);
      X0muroll = ref_varma_sim(A, B, Sig, mu, n, M, StartX0, [], rng);
      mat2file(fidRolling, Xroll, "X" + M + "-" + k);
      mat2file(fidRolling, X0roll, "X0" + M + "-" + k);
      mat2file(fidRolling, Xmuroll, "Xmu" + M + "-" + k);
      mat2file(fidRolling, X0muroll, "X0mu" + M + "-" + k);
    end
  end
  for k = reducedCases
    [A, B, Sig, p, q, r] = ref_varma_testcase(k);
    StartX0 = startmat(r, max(p,q));
    mu = (1:r)'/10;
    M = reducedM;
    rng = randompack_create();
    cleanup = onCleanup(@() randompack_free(rng));
    randompack_seed(rng, 42);
    [X, E] = ref_varma_sim(A, B, Sig, 0, n, M, [], [], rng);
    randompack_seed(rng, 42);
    [X0, E0] = ref_varma_sim(A, B, Sig, 0, n, M, StartX0, [], rng);
    randompack_seed(rng, 42);
    [Xmu, Emu] = ref_varma_sim(A, B, Sig, mu, n, M, [], [], rng);
    randompack_seed(rng, 42);
    [X0mu, E0mu] = ref_varma_sim(A, B, Sig, mu, n, M, StartX0, [], rng);
    if r == 1
      Xrep2 = X(:, 2)';
    else
      Xrep2 = X(:, :, 2);
    end
    CmlT = ref_varma_autocov(Xrep2, n - 1, "ML");
    CcT = ref_varma_autocov(Xrep2, n - 1, "C");
    mat2file(fid, CmlT, "RedAutoMLT" + M + "-" + k);
    mat2file(fid, CcT, "RedAutoCT" + M + "-" + k);
    mat2file(fid, X, "RedX" + M + "-" + k);
    mat2file(fid, E, "RedE" + M + "-" + k);
    mat2file(fid, X0, "RedX0" + M + "-" + k);
    mat2file(fid, E0, "RedE0" + M + "-" + k);
    mat2file(fid, Xmu, "RedXmu" + M + "-" + k);
    mat2file(fid, Emu, "RedEmu" + M + "-" + k);
    mat2file(fid, X0mu, "RedX0mu" + M + "-" + k);
    mat2file(fid, E0mu, "RedE0mu" + M + "-" + k);
  end
  xcases = 1:ref_varma_testcasex('count');
  xMultiM = 3;
  mat2file(fid, length(xcases), "#xcases");
  mat2file(fid, xMultiM, "xMultiM");
  mat2file(fid, xcases, "xcases");
  for k = xcases
    [A, B, C, Sig, z, p, q, s, r] = ref_varma_testcasex(k, n);
    h = max([p, q, s]) + 1;
    StartX0 = startmat(r, h);
    mat2file(fid, p, "xp" + k);
    mat2file(fid, q, "xq" + k);
    mat2file(fid, s, "xs" + k);
    mat2file(fid, r, "xr" + k);
    mat2file(fid, h, "xh" + k);
    mat2file(fid, A, "xA" + k);
    mat2file(fid, B, "xB" + k);
    mat2file(fid, C, "xC" + k);
    mat2file(fid, Sig, "xSig" + k);
    mat2file(fid, z, "xz" + k);
    mat2file(fid, StartX0, "xStartX0" + k);
    zM = zpaths(z, xMultiM);
    StartX0M = startpaths(r, h, xMultiM);
    mat2file(fid, zM', "xzM" + xMultiM + "-" + k);
    mat2file(fid, StartX0M, "xStartX0M" + xMultiM + "-" + k);
    for M = Ms
      rng = randompack_create();
      cleanup = onCleanup(@() randompack_free(rng));
      randompack_seed(rng, 42);
      [X, E] = ref_varma_simx(A, B, C, z, Sig, n, M, StartX0, h, [], rng);
      randompack_seed(rng, 42);
      XnoE = ref_varma_simx(A, B, C, z, Sig, n, M, StartX0, h, [], rng);
      mat2file(fid, X, "Xx" + M + "-" + k);
      mat2file(fid, E, "Ex" + M + "-" + k);
      mat2file(fid, XnoE, "XxNoE" + M + "-" + k);
    end
    rng = randompack_create();
    cleanup = onCleanup(@() randompack_free(rng));
    randompack_seed(rng, 42);
    [X, E] = ref_varma_simx(A, B, C, zM, Sig, n, xMultiM, StartX0M, h, [], rng);
    mat2file(fid, X, "XxM" + xMultiM + "-" + k);
    mat2file(fid, E, "ExM" + xMultiM + "-" + k);
  end
  fclose(fid);
  fclose(fidRolling);
end

function X0 = startmat(r, h)
  X0 = repmat([-2,-1,0,1,2], 1, ceil(r*h/5));
  X0 = reshape(X0(1:r*h), r, h);
end

function X0 = startpaths(r, h, M)
  X0 = zeros(r, h, M);
  base = startmat(r, h);
  for j = 1:M
    X0(:,:,j) = base + (j - 1)/10;
  end
end

function zM = zpaths(z, M)
  n = length(z);
  zM = zeros(M, n);
  for j = 1:M
    zM(j,:) = z + (j - 1)/20;
  end
end

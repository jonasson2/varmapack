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
  if nargin < 3, n = 8; end
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
  reducedCases = [3, 7, 12];
  reducedM = 3;
  mat2file(fid, length(reducedCases), "#reducedCases");
  mat2file(fid, reducedM, "reducedM");
  mat2file(fid, reducedCases, "reducedCases");
  for k = cases
    icase = k + 1;
    [A, B, Sig, p, q, r] = ref_varma_testcase(k);    
    x0 = startmat(r, max(p,q));
    mu = (1:r)'/10;
    mat2file(fid, p, "p" + icase);
    mat2file(fid, q, "q" + icase);
    mat2file(fid, r, "r" + icase);
    mat2file(fid, A, "A" + icase);
    mat2file(fid, B, "B" + icase);
    mat2file(fid, Sig, "Sig" + icase);
    mat2file(fid, x0, "x0" + icase);
    mat2file(fid, mu, "mu" + icase);
    mat2file(fidRolling, x0, "x0" + icase);
    mat2file(fidRolling, mu, "mu" + icase);
    if p == 0
      rho = 0;
    else
      rho = ref_varma_specrad(A);
    end
    mat2file(fid, rho, "specrad-" + k);
    acvfMaxlag = max(n - 1, p);
    Gamma = ref_varma_acvf(A, B, Sig, acvfMaxlag);
    mat2file(fid, acvfMaxlag, "acvfMaxlag-" + k);
    mat2file(fid, Gamma, "Gamma-" + k);
    for M = Ms
      rng = randompack_create();
      cleanup = onCleanup(@() randompack_free(rng));
      randompack_seed(rng, 42);
      [X, E, condR] = ref_varma_sim(A, B, Sig, n, 0, M, [], rng);
      randompack_seed(rng, 42);
      [X0, E0]      = ref_varma_sim(A, B, Sig, n, 0, M, x0, rng);
      randompack_seed(rng, 42);
      [Xmu, Emu]    = ref_varma_sim(A, B, Sig, n, mu, M, [], rng);
      randompack_seed(rng, 42);
      [X0mu, E0mu]  = ref_varma_sim(A, B, Sig, n, mu, M, x0, rng);
      if M == Ms(1), mat2file(fid, condR, "condR-" + k); end
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
      randompack_seed(rng, 42);
      Xroll = ref_varma_sim(A, B, Sig, n, 0, M, [], rng);
      randompack_seed(rng, 42);
      X0roll = ref_varma_sim(A, B, Sig, n, 0, M, x0, rng);
      randompack_seed(rng, 42);
      Xmuroll = ref_varma_sim(A, B, Sig, n, mu, M, [], rng);
      randompack_seed(rng, 42);
      X0muroll = ref_varma_sim(A, B, Sig, n, mu, M, x0, rng);
      mat2file(fidRolling, Xroll, "X" + M + "-" + k);
      mat2file(fidRolling, X0roll, "X0" + M + "-" + k);
      mat2file(fidRolling, Xmuroll, "Xmu" + M + "-" + k);
      mat2file(fidRolling, X0muroll, "X0mu" + M + "-" + k);
    end
  end
  for k = reducedCases
    [A, B, Sig, p, q, r] = ref_varma_testcase(k);
    x0 = startmat(r, max(p,q));
    mu = (1:r)'/10;
    M = reducedM;
    rng = randompack_create();
    cleanup = onCleanup(@() randompack_free(rng));
    randompack_seed(rng, 42);
    [X, E] = ref_varma_sim(A, B, Sig, n, 0, M, [], rng);
    randompack_seed(rng, 42);
    [X0, E0] = ref_varma_sim(A, B, Sig, n, 0, M, x0, rng);
    randompack_seed(rng, 42);
    [Xmu, Emu] = ref_varma_sim(A, B, Sig, n, mu, M, [], rng);
    randompack_seed(rng, 42);
    [X0mu, E0mu] = ref_varma_sim(A, B, Sig, n, mu, M, x0, rng);
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
  fclose(fid);
  fclose(fidRolling);
end

function x0 = startmat(r, h)
  x0 = repmat([-2,-1,0,1,2], 1, ceil(r*h/5));
  x0 = reshape(x0(1:r*h), r, h);
end

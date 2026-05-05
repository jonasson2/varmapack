function matlabcompare(cases, Ms, n)
  % MATLABCOMPARE writes matlabcompare.txt, the reference data consumed by
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
  if nargin < 3, n = 5; end
  if nargin < 2, Ms = 1; end
  if nargin < 1 || isempty(cases), cases = 1:ref_varma_testcase('number'); end
  fid = fopen("../matlabcompare.txt", "w");
  fprintf(fid, "# varmapack Matlab reference data v1\n");
  fprintf(fid, "# Format documented in tests/matlabcompare.m\n");
  mat2file(fid, length(Ms), "#Ms");
  mat2file(fid, length(cases), "#cases");
  mat2file(fid, n);
  mat2file(fid, Ms);
  mat2file(fid, cases);
  % folder = """" + pwd() + filesep + """";
  % cmd = "  strcpy(matmatfolder, " + folder + ");\n";
  % fprintf(fid, cmd);
  % char *matmatfolder[] = """ + pwd() + filesep + """");
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
      global varmapack_rng
      varmapack_rng = randompack_create();
      cleanup = onCleanup(@() randompack_free(varmapack_rng));
      randompack_seed(varmapack_rng, 42);
      [X, E, condR] = ref_varma_sim(A, B, Sig, n, 0, M);
      randompack_seed(varmapack_rng, 42);
      [X0, E0]      = ref_varma_sim(A, B, Sig, n, 0, M, x0);
      randompack_seed(varmapack_rng, 42);
      [Xmu, Emu]    = ref_varma_sim(A, B, Sig, n, mu, M);
      randompack_seed(varmapack_rng, 42);
      [X0mu, E0mu]  = ref_varma_sim(A, B, Sig, n, mu, M, x0);
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
    end
  end
  fclose(fid);
end

function x0 = startmat(r, h)
  x0 = repmat([-2,-1,0,1,2], 1, ceil(r*h/5));
  x0 = reshape(x0(1:r*h), r, h);
end

function TimeSimulateVar(varargin)
  parser = inputParser;
  addParameter(parser, "t", 0.1);
  addParameter(parser, "w", 0.1);
  addParameter(parser, "n", 100);
  addParameter(parser, "M", 1000);
  addParameter(parser, "burnin", 200);
  addParameter(parser, "selected", false);
  addParameter(parser, "platform", false);
  parse(parser, varargin{:});
  opts = parser.Results;
  if opts.t <= 0 || opts.w < 0 || opts.n < 1 || opts.M < 1 || opts.burnin < 0
    error("TimeSimulateVar:options", "timings and dimensions must be valid");
  end
  if opts.selected && opts.platform
    error("TimeSimulateVar:options", "selected and platform cannot both be true");
  end
  if exist("varm", "file") ~= 2
    error("TimeSimulateVar:varm", ...
          "TimeSimulateVar requires Econometrics Toolbox (varm)");
  end
  time_problem_table(opts);
end

function time_problem_table(opts)
  rngC = varmapack.Rng();
  cleanupC = onCleanup(@() delete(rngC));
  warm_cpu(opts.w);
  warm_simulations(rngC, opts);
  fprintf("Varmapack MATLAB benchmark and comparison with Econometrics Toolbox VAR\n");
  fprintf("Benchmark unit:          ns/retained value\n");
  fprintf("Length per series:       %d\n", opts.n);
  fprintf("Replicates per call:     %d\n", opts.M);
  fprintf("Discarded burn-in:        %d per path\n", opts.burnin);
  fprintf("Benchmark time per case: %.1f s\n\n", opts.t);
  if opts.selected
    fprintf("Models:                  selected compatible paper set\n\n");
  end
  if opts.platform
    fprintf("Models:                  selected platform set\n\n");
  end
  fprintf("%-12s %2s %2s %2s %5s %10s %10s %6s\n", ...
          "Testcase", "p", "q", "r", "rho", "Varmapack", "varm", "Ratio");
  names = {"tinyAR", "smallAR1", "smallAR2", "mediumAR", "largeAR"};
  if opts.selected
    names = {"tinyAR", "tinyARMA", "smallAR1", "smallARMA2", ...
             "mediumAR", "mediumARMA2", "largeAR", "largeARMA"};
  end
  if opts.platform
    names = {"tinyARMA", "smallAR2", "mediumAR", "mediumARMA1", "largeAR"};
  end
  for i = 1:numel(names)
    [A, B, Sig, p, q, r] = varmapack.testcase(names{i});
    time_case(names{i}, A, B, Sig, p, q, r, opts, rngC);
  end
  if ~opts.selected && ~opts.platform
    [A, B, Sig, p, q, r] = varmapack.testcase(3, 0, 10, .98);
    time_case("largeVAR", A, B, Sig, p, q, r, opts, rngC);
  end
  clear cleanupC
end

function time_case(name, A, B, Sig, p, q, r, opts, rngC)
  rho = varmapack.specrad(A);
  values = r*opts.n*opts.M;
  total = opts.n + opts.burnin;
  rngC.seed(12345);
  varmapack_simulate(A, B, Sig, total, opts.M, opts.burnin, rngC);
  rngC.seed(12345);
  varmapack_time = time_one(@() ...
      varmapack_simulate(A, B, Sig, total, opts.M, opts.burnin, rngC), values, opts.t);
  varmapack_ns = 1e9*varmapack_time/values;
  if q == 0
    mdl = matlab_var(A, Sig, p, r);
    rng(12345);
    varm_simulate(mdl, total, opts.M, opts.burnin);
    rng(12345);
    varm_time = time_one(@() varm_simulate(mdl, total, opts.M, opts.burnin), ...
                         values, opts.t);
    varm_ns = 1e9*varm_time/values;
    fprintf("%-12s %2d %2d %2d %5.3f %10.1f %10.1f %6.1f\n", ...
            name, p, q, r, rho, varmapack_ns, varm_ns, varm_ns/varmapack_ns);
  else
    fprintf("%-12s %2d %2d %2d %5.3f %10.1f %10s %6s\n", ...
            name, p, q, r, rho, varmapack_ns, "-", "-");
  end
end

function mdl = matlab_var(A, Sig, p, r)
  AR = cell(1, p);
  for i = 1:p
    AR{i} = A(:, (i - 1)*r + 1:i*r);
  end
  mdl = varm(r, p);
  mdl.Constant = zeros(r, 1);
  mdl.AR = AR;
  mdl.Covariance = Sig;
end

function varm_simulate(mdl, total, M, burnin)
  X = simulate(mdl, total, NumPaths = M);
  if size(X, 1) - burnin < 1
    error("TimeSimulateVar:burnin", "burn-in leaves no retained values");
  end
end

function varmapack_simulate(A, B, Sig, total, M, burnin, rngC)
  X = varmapack.sim(A, B, Sig, [], total, M, [], rngC);
  if size(X, 2) - burnin < 1
    error("TimeSimulateVar:burnin", "burn-in leaves no retained values");
  end
end

function seconds = time_one(fun, values, timingTarget)
  repsPerCheck = max(1, floor(1e6/values));
  reps = 0;
  elapsed = 0;
  tic;
  while elapsed < timingTarget
    for i = 1:repsPerCheck
      fun();
    end
    reps = reps + repsPerCheck;
    elapsed = toc;
  end
  seconds = elapsed/reps;
end

function warm_simulations(rngC, opts)
  [A, B, Sig] = varmapack.testcase("mediumAR");
  mdl = matlab_var(A, Sig, 1, 3);
  rng(12345);
  varm_simulate(mdl, opts.n + opts.burnin, 10, opts.burnin);
  rngC.seed(12345);
  varmapack_simulate(A, B, Sig, opts.n + opts.burnin, 10, opts.burnin, rngC);
end

function warm_cpu(seconds)
  A = rand(96);
  elapsed = 0;
  while elapsed < seconds
    tic;
    A = A*A;
    A = A/norm(A, "fro");
    elapsed = elapsed + toc;
  end
end

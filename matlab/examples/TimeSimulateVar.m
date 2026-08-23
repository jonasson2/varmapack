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
  rngVarmapack = varmapack.Rng();
  cleanupVarmapack = onCleanup(@() delete(rngVarmapack));
  benchmark_warm_cpu(opts.w);
  fprintf("Varmapack MATLAB benchmark and comparison with Econometrics Toolbox VAR\n");
  fprintf("Benchmark unit:          ns/retained value\n");
  fprintf("Length per series:       %d\n", opts.n);
  fprintf("Replicates per call:     %d\n", opts.M);
  fprintf("varm discarded burn-in:  %d per path\n", opts.burnin);
  fprintf("Target time per method:  %.1f s per testcase\n", opts.t);
  if opts.selected
    fprintf("Models:                  selected paper set\n");
  end
  if opts.platform
    fprintf("Models:                  selected platform set\n");
  end
  fprintf("\n");
  fprintf("%-12s %2s %2s %2s %5s %10s %10s %6s\n", ...
          "Testcase", "p", "q", "r", "rho", "Varmapack", "varm", "Ratio");
  names = benchmark_testcases(opts.selected, opts.platform, "var");
  for i = 1:numel(names)
    [A, B, Sig, p, q, r] = varmapack.testcase(names{i});
    time_case(names{i}, A, B, Sig, p, q, r, opts, rngVarmapack);
  end
  if ~opts.selected && ~opts.platform
    [A, B, Sig, p, q, r] = varmapack.testcase(3, 0, 10, .98);
    time_case("largeVAR", A, B, Sig, p, q, r, opts, rngVarmapack);
  end
end

function time_case(name, A, B, Sig, p, q, r, opts, rngVarmapack)
  rho = varmapack.specrad(A);
  values = r*opts.n*opts.M;
  rngVarmapack.seed(12345);
  varmapackNs = benchmark_time( ...
      @() varmapack_simulate(A, B, Sig, opts.n, opts.M, rngVarmapack), ...
      values, opts.t);
  if q == 0
    mdl = matlab_var(A, Sig, p, r);
    rng(12345);
    varmNs = benchmark_time( ...
        @() varm_simulate(mdl, opts.n + opts.burnin, opts.M), values, opts.t);
    fprintf("%-12s %2d %2d %2d %5.3f %10.1f %10.1f %6.1f\n", ...
            name, p, q, r, rho, varmapackNs, varmNs, varmNs/varmapackNs);
  else
    fprintf("%-12s %2d %2d %2d %5.3f %10.1f %10s %6s\n", ...
            name, p, q, r, rho, varmapackNs, "-", "-");
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

function varm_simulate(mdl, total, M)
  simulate(mdl, total, NumPaths = M);
end

function varmapack_simulate(A, B, Sig, n, M, rngVarmapack)
  varmapack.sim(A, B, Sig, [], n, M, [], rngVarmapack);
end

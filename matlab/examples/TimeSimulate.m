function TimeSimulate(varargin)
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
    error("TimeSimulate:options", "timings and dimensions must be valid");
  end
  if opts.selected && opts.platform
    error("TimeSimulate:options", "selected and platform cannot both be true");
  end
  time_problem_table(opts, exist("varm", "file") == 2);
end

function time_problem_table(opts, hasVarm)
  rngRef = randompack_create();
  rngVarmapack = varmapack.Rng();
  cleanupRef = onCleanup(@() randompack_free(rngRef));
  cleanupVarmapack = onCleanup(@() delete(rngVarmapack));
  benchmark_warm_cpu(opts.w);
  fprintf("Varmapack MATLAB benchmark and comparison with reference simulation");
  if hasVarm
    fprintf(" and varm\n");
  else
    fprintf("\n");
  end
  fprintf("Benchmark unit:          ns/retained value\n");
  fprintf("Length per series:       %d\n", opts.n);
  fprintf("Replicates per call:     %d\n", opts.M);
  fprintf("Target time per method:  %.1f s per testcase\n", opts.t);
  if hasVarm
    fprintf("varm discarded burn-in:  %d per path\n", opts.burnin);
  else
    fprintf("varm:                    unavailable\n");
  end
  if opts.selected
    fprintf("Models:                  selected paper set\n");
  end
  if opts.platform
    fprintf("Models:                  selected platform set\n");
  end
  fprintf("\n");
  fprintf("%-12s %2s %2s %2s %5s %10s %10s %10s %6s %7s\n", ...
          "Testcase", "p", "q", "r", "rho", "Varmapack", "Reference", ...
          "varm", "Ref/VP", "varm/VP");
  names = benchmark_testcases(opts.selected, opts.platform, "all");
  for i = 1:numel(names)
    [A, B, Sig, p, q, r] = varmapack.testcase(names{i});
    time_case(names{i}, A, B, Sig, p, q, r, opts, hasVarm, rngRef, ...
              rngVarmapack);
  end
end

function time_case(name, A, B, Sig, p, q, r, opts, hasVarm, rngRef, ...
                   rngVarmapack)
  if p == 0
    rho = 0;
  else
    rho = varmapack.specrad(A);
  end
  values = r*opts.n*opts.M;
  rngVarmapack.seed(12345);
  varmapackNs = benchmark_time( ...
      @() call_varmapack_sim(A, B, Sig, opts.n, opts.M, rngVarmapack), ...
      values, opts.t);
  randompack_seed(rngRef, 12345);
  referenceNs = benchmark_time( ...
      @() call_ref_sim(A, B, Sig, opts.n, opts.M, rngRef), values, opts.t);
  if hasVarm && q == 0
    mdl = matlab_var(A, Sig, p, r);
    rng(12345);
    varmNs = benchmark_time( ...
        @() varm_simulate(mdl, opts.n + opts.burnin, opts.M), values, opts.t);
    fprintf("%-12s %2d %2d %2d %5.3f %10.1f %10.0f %10.1f %6.1f %7.1f\n", ...
            name, p, q, r, rho, varmapackNs, referenceNs, varmNs, ...
            referenceNs/varmapackNs, varmNs/varmapackNs);
  else
    fprintf("%-12s %2d %2d %2d %5.3f %10.1f %10.0f %10s %6.1f %7s\n", ...
            name, p, q, r, rho, varmapackNs, referenceNs, "-", ...
            referenceNs/varmapackNs, "-");
  end
end

function call_ref_sim(A, B, Sig, n, M, rng)
  ref_varma_sim(A, B, Sig, [], n, M, [], rng);
end

function call_varmapack_sim(A, B, Sig, n, M, rng)
  varmapack.sim(A, B, Sig, [], n, M, [], rng);
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

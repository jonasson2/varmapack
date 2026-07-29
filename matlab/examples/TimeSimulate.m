function TimeSimulate(varargin)
  parser = inputParser;
  addParameter(parser, "t", 0.1);
  addParameter(parser, "w", 0.1);
  addParameter(parser, "n", 100);
  addParameter(parser, "M", 1000);
  parse(parser, varargin{:});
  opts = parser.Results;
  if opts.t <= 0 || opts.w < 0 || opts.n < 1 || opts.M < 1
    error("TimeSimulate:options", "t and dimensions must be positive");
  end
  time_problem_table(opts);
end

function time_problem_table(opts)
  rngRef = randompack_create();
  rngC = varmapack.Rng();
  cleanupRef = onCleanup(@() randompack_free(rngRef));
  cleanupC = onCleanup(@() delete(rngC));
  warm_cpu(opts.w);
  warm_simulations(rngRef, rngC);
  fprintf("Varmapack benchmark and comparison with MATLAB reference " + ...
          "simulation\n");
  fprintf("Benchmark unit:          ns/value\n");
  fprintf("length per series:       %d\n", opts.n);
  fprintf("replicates per call:     %d\n", opts.M);
  fprintf("benchmark time per case: %.1f s\n", opts.t);
  fprintf("\n");
  fprintf("%-12s %2s %2s %2s %5s %10s %10s %6s\n", ...
          "Testcase", "p", "q", "r", "rho", "Varmapack", "Reference", "Ratio");
  names = { ...
    "tinyAR", "tinyMA", "tinyARMA", "smallAR1", "smallAR2", "smallMA1", ...
    "smallMA2", "smallARMA1", "smallARMA2", "mediumAR", "mediumMA1", ...
    "mediumARMA1", "mediumARMA2", "mediumMA2", "largeAR"};
  for i = 1:numel(names)
    [A, B, Sig, p, q, r] = varmapack.testcase(names{i});
    time_case(names{i}, A, B, Sig, p, q, r, opts, rngRef, rngC);
  end
  [A, B, Sig, p, q, r] = varmapack.testcase(3, 3, 10, .98);
  time_case("Unnamed", A, B, Sig, p, q, r, opts, rngRef, rngC);
  clear cleanupRef cleanupC
end

function time_case(name, A, B, Sig, p, q, r, opts, rngRef, rngC)
  if p == 0
    rho = 0;
  else
    rho = varmapack.specrad(A);
  end
  values = r*opts.n*opts.M;
  randompack_seed(rngRef, 12345);
  refTime = time_one(@() call_ref_sim(A, B, Sig, opts.n, opts.M, rngRef), ...
                     values, opts.t);
  rngC.seed(12345);
  cTime = time_one(@() call_c_sim(A, B, Sig, opts.n, opts.M, rngC), ...
                   values, opts.t);
  refNs = 1e9*refTime/values;
  cNs = 1e9*cTime/values;
  fprintf("%-12s %2d %2d %2d %5.2f %10.1f %10.0f %6.1f\n", ...
          name, p, q, r, rho, cNs, refNs, refNs/cNs);
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

function warm_simulations(rngRef, rngC)
  [A, B, Sig] = varmapack.testcase("mediumARMA1");
  randompack_seed(rngRef, 12345);
  call_ref_sim(A, B, Sig, 100, 100, rngRef);
  rngC.seed(12345);
  call_c_sim(A, B, Sig, 100, 100, rngC);
end

function call_ref_sim(A, B, Sig, n, M, rng)
  X = ref_varma_sim(A, B, Sig, [], n, M, [], rng);
end

function call_c_sim(A, B, Sig, n, M, rng)
  X = varmapack.sim(A, B, Sig, [], n, M, [], rng);
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

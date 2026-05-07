function TimeVarmaSim(varargin)
  parser = inputParser;
  addParameter(parser, "t", 0.2);
  addParameter(parser, "w", 0.1);
  addParameter(parser, "d", 2);
  addParameter(parser, "p", [1 3 5]);
  addParameter(parser, "q", [1 3 5]);
  addParameter(parser, "r", [2 5 16]);
  addParameter(parser, "n", [0 100]);
  addParameter(parser, "M", [1 100]);
  addParameter(parser, "rho", [.5 .95 .995]);
  parse(parser, varargin{:});
  opts = parser.Results;
  fprintf("TimeVarmaSim options\n");
  fprintf("  t   = %.16g\n", opts.t);
  fprintf("  w   = %.16g\n", opts.w);
  fprintf("  d   = %.16g\n", opts.d);
  print_list("p", opts.p);
  print_list("q", opts.q);
  print_list("r", opts.r);
  print_list("n", opts.n);
  print_list("M", opts.M);
  print_list("rho", opts.rho);
  fprintf("\n");
  time_problem_table(opts);
end

function print_list(name, x)
  fprintf("  %-3s = ", name);
  fprintf("%.16g ", x);
  fprintf("\n");
end

function time_problem_table(opts)
  cases = benchmark_cases(opts);
  rngRef = randompack_create();
  rngC = randompack_create();
  cleanupRef = onCleanup(@() randompack_free(rngRef));
  cleanupC = onCleanup(@() randompack_free(rngC));
  warm_cpu(max(opts.w, 0.1));
  fmt = "%3s %3s %3s %5s %5s %8s %10s %12s %12s %8s\n";
  fprintf(fmt, "p", "q", "r", "n", "M", "rho", "base", "ref ns/val", ...
          "C ns/val", "speedup");
  for i = 1:size(cases, 1)
    p = cases(i, 1);
    q = cases(i, 2);
    r = cases(i, 3);
    n = cases(i, 4);
    M = cases(i, 5);
    rho = cases(i, 6);
    if n == 0
      n = max(p, q);
    end
    [A, B, Sig, baseRho] = benchmark_problem(p, q, r, rho);
    achieved = varmapack_specrad(A);
    if achieved >= 1
      error("TimeVarmaSim:rho", "rho must be below 1 for stationary timing");
    end
    randompack_seed(rngRef, 12345);
    refTime = time_one(@() call_ref_sim(A, B, Sig, n, M, rngRef), opts.t);
    randompack_seed(rngC, 12345);
    cTime = time_one(@() call_c_sim(A, B, Sig, n, M, rngC), opts.t);
    scale = 1e9/(r*n*M);
    fprintf("%3d %3d %3d %5d %5d %8.4f %10.4f %12.*f %12.*f %8.2f\n", ...
            p, q, r, n, M, achieved, baseRho, opts.d, refTime*scale, ...
            opts.d, cTime*scale, refTime/cTime);
  end
  clear cleanupRef cleanupC
end

function cases = benchmark_cases(opts)
  base = [3 3 5 100 100 .95];
  cases = base;
  cases = add_axis(cases, base, 1, opts.p);
  cases = add_axis(cases, base, 2, opts.q);
  cases = add_axis(cases, base, 3, opts.r);
  cases = add_axis(cases, base, 4, opts.n);
  cases = add_axis(cases, base, 5, opts.M);
  cases = add_axis(cases, base, 6, opts.rho);
  cases = unique(cases, "rows", "stable");
end

function cases = add_axis(cases, base, k, values)
  for i = 1:numel(values)
    row = base;
    row(k) = values(i);
    cases(end + 1, :) = row;
  end
end

function [A, B, Sig, baseRho] = benchmark_problem(p, q, r, targetRho)
  Ahi = benchmark_Ahi(p, r);
  baseRho = varmapack_specrad(Ahi);
  [A, B, Sig] = varmapack_testcase(p, q, r, targetRho);
end

function seconds = time_one(fun, timingTarget)
  reps = 0;
  elapsed = 0;
  while elapsed < timingTarget
    tic;
    fun();
    elapsed = elapsed + toc;
    reps = reps + 1;
  end
  seconds = elapsed/reps;
end

function call_ref_sim(A, B, Sig, n, M, rng)
  [~, ~] = ref_varma_sim(A, B, Sig, n, [], M, [], rng);
end

function call_c_sim(A, B, Sig, n, M, rng)
  [~, ~] = varmapack_sim(A, B, Sig, n, [], M, [], rng);
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

function A = benchmark_Ahi(p, r)
  A = zeros(r, r*p);
  if p == 0
    return
  end
  [I, J] = ndgrid(1:r, 1:(p*r));
  A = 2*min(I, J)./max(I, J)/(r*p^(1/3));
end

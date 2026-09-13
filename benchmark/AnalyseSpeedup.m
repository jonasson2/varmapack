function AnalyseSpeedup(varargin)
%ANALYSESPEEDUP  Controlled Varmapack/MATLAB varm comparison.

  parser = inputParser;
  addParameter(parser, "R", 5);
  addParameter(parser, "t", 0.1);
  addParameter(parser, "w", 0.1);
  addParameter(parser, "n", 100);
  addParameter(parser, "M", 1000);
  addParameter(parser, "burnin", 200);
  parse(parser, varargin{:});
  opts = parser.Results;
  if opts.R < 1 || opts.t <= 0 || opts.w < 0 || opts.n < 1 || opts.M < 2 || ...
      opts.burnin < 0
    error("AnalyseSpeedup:options", "timings and dimensions must be valid");
  end
  benchmark_warm_cpu(opts.w);
  names = {"tinyAR", "smallAR1", "mediumAR", "largeAR"};
  runs = struct([]);
  index = 1;
  for trial = 1:opts.R
    fprintf(2, "MATLAB interface run %d of %d\n", trial, opts.R);
    for i = 1:numel(names)
      row = time_case(names{i}, opts);
      row.trial = trial;
      if index == 1
        runs = row;
      else
        runs(index) = row;
      end
      index = index + 1;
    end
  end
  print_table(summarize(median_timings(runs, names), opts), opts);
end

function row = time_case(name, opts)
  [A, B, Sig, p, q, r] = varmapack.testcase(name);
  if q ~= 0
    error("AnalyseSpeedup:model", "varm only supports VAR models");
  end
  startup = p;
  if opts.n < startup
    error("AnalyseSpeedup:length", "series length is shorter than startup");
  end
  mdl = matlab_var(A, Sig, p, r);
  rngVarmapack = varmapack.Rng();
  cleanup = onCleanup(@() delete(rngVarmapack));
  rngVarmapack.seed(12345);
  row = struct( ...
      "model", name, "p", p, "q", q, "r", r, "startup", startup, ...
      "vpGeneration", time_vp_generation(A, B, Sig, opts.n, opts.t, ...
          rngVarmapack), ...
      "vpTotal1", time_vp(A, B, Sig, opts.n, 1, opts.t, rngVarmapack), ...
      "vpSetupBatch", time_vp(A, B, Sig, startup, opts.M, opts.t, ...
          rngVarmapack), ...
      "vpBatch", time_vp(A, B, Sig, opts.n, opts.M, opts.t, ...
          rngVarmapack), ...
      "varmSingle", time_varm(mdl, opts.n, 1, opts.t), ...
      "varmBatch", time_varm(mdl, opts.n, opts.M, opts.t), ...
      "varmBurnin", time_varm(mdl, opts.n + opts.burnin, opts.M, opts.t));
end

function ns = time_vp(A, B, Sig, n, M, target, rngVarmapack)
  ns = time_call(@() varmapack.sim(A, B, Sig, [], n, M, [], rngVarmapack), ...
                 target);
end

function ns = time_varm(mdl, n, M, target)
  rng(12345);
  ns = time_call(@() simulate(mdl, n, NumPaths=M), target);
end

function ns = time_call(fun, target)
  result = fun();
  if ~all(isfinite(result), "all")
    error("AnalyseSpeedup:simulation", "simulation returned non-finite values");
  end
  clear result
  calls = 0;
  elapsed = 0;
  timer = tic;
  while elapsed < target
    fun();
    calls = calls + 1;
    elapsed = toc(timer);
  end
  ns = 1e9*elapsed/calls;
end

function ns = time_vp_generation(A, B, Sig, n, target, rngVarmapack)
  lengths = n*[10, 20, 40, 80];
  times = zeros(size(lengths));
  for i = 1:numel(lengths)
    times(i) = time_vp(A, B, Sig, lengths(i), 1, target, rngVarmapack);
  end
  coefficient = polyfit(lengths, times, 1);
  ns = coefficient(1);
  if ~isfinite(ns) || ns <= 0
    error("AnalyseSpeedup:fit", "invalid fitted Varmapack generation time");
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

function rows = median_timings(runs, names)
  dimensions = {"p", "q", "r", "startup"};
  timings = {"vpGeneration", "vpTotal1", "vpSetupBatch", "vpBatch", ...
             "varmSingle", "varmBatch", "varmBurnin"};
  rows = struct([]);
  for i = 1:numel(names)
    selected = runs(arrayfun(@(item) strcmp(item.model, names{i}), runs));
    row = struct("model", names{i});
    for j = 1:numel(dimensions)
      field = dimensions{j};
      row.(field) = selected(1).(field);
    end
    for j = 1:numel(timings)
      field = timings{j};
      row.(field) = median([selected.(field)]);
    end
    if i == 1
      rows = row;
    else
      rows(i) = row;
    end
  end
end

function rows = summarize(rows, opts)
  summary = struct([]);
  for i = 1:numel(rows)
    row = rows(i);
    intermediate = row.vpSetupBatch/opts.M + ...
                   (opts.n - row.startup)*row.vpGeneration;
    row.overall = row.varmBurnin/row.vpBatch;
    row.burnin = row.varmBurnin/row.varmBatch;
    row.setupAmortization = row.vpTotal1/intermediate;
    row.bulkBatch = (intermediate/(row.vpBatch/opts.M))/ ...
                    (row.varmSingle/(row.varmBatch/opts.M));
    row.implementation = row.varmSingle/row.vpTotal1;
    if i == 1
      summary = row;
    else
      summary(i) = row;
    end
  end
  rows = summary;
end

function print_table(rows, opts)
  fprintf("Runs: %d; target: %.1f s; length: %d; batch: %d.\n\n", ...
          opts.R, opts.t, opts.n, opts.M);
  fprintf("Speed-up factor due to:\n\n");
  fprintf("| Model     | Overall | Burn-in | Setup/amort. | Bulk/batch | Implementation |\n");
  fprintf("| --------- | ------: | ------: | -----------: | ----------: | -------------: |\n");
  for i = 1:numel(rows)
    row = rows(i);
    fprintf("| %-9s | %7.0f | %7.1f | %12.1f | %11.1f | %14.1f |\n", ...
            row.model, row.overall, row.burnin, row.setupAmortization, ...
            row.bulkBatch, row.implementation);
  end
end

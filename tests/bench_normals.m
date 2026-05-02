function bench_normals
  tic
  warmup_randn(0.1);
  toc
  n = 2e7;
  trials = 7;

  rng(123,'twister');
  x = zeros(n,1,'double');

  % Warm-up (~0.1 s total)
  for k = 1:5
    x = randn(2e6,1);
  end

  % Timed trials (fill existing array shape)
  t = zeros(trials,1);
  for k = 1:trials
    tic;
    x = randn(size(x));
    t(k) = toc;
  end

  best = min(t);
  ns_per = best * 1e9 / numel(x);

  fprintf('MATLAB randn: best %.3f s  =>  %.3f ns/value (n=%g)\n', ...
          best, ns_per, n);
end

function warmup_randn(seconds)
  if nargin < 1
    seconds = 0.1;
  end

  n = 2e6;                 % chunk size (large enough to stay compute-bound)
  x = zeros(n,1,'double');

  t0 = tic;
  while toc(t0) < seconds
    x = randn(size(x));    % fill existing shape
  end
end
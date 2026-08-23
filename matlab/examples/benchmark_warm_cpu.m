function benchmark_warm_cpu(seconds)
%BENCHMARK_WARM_CPU  Warm up the CPU for the requested time.
  if seconds <= 0
    return
  end
  A = rand(96);
  elapsed = 0;
  while elapsed < seconds
    timer = tic;
    A = A*A;
    A = A/norm(A, "fro");
    elapsed = elapsed + toc(timer);
  end
end

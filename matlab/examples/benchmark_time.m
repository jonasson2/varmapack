function ns = benchmark_time(fun, values, target)
%BENCHMARK_TIME  Time repeated calls and return nanoseconds per value.
  fun();
  calls = 0;
  elapsed = 0;
  timer = tic;
  while elapsed < target
    fun();
    calls = calls + 1;
    elapsed = toc(timer);
  end
  ns = 1e9*elapsed/(calls*values);
end

function names = benchmark_testcases(selected, platform, defaultSet)
%BENCHMARK_TESTCASES  Return the requested ordered set of testcase names.
  if selected
    setName = "paper";
  elseif platform
    setName = "platform";
  else
    setName = defaultSet;
  end
  switch setName
    case "all"
      names = { ...
        "tinyAR", "tinyMA", "tinyARMA", "smallAR1", "smallAR2", ...
        "smallMA1", "smallMA2", "smallARMA1", "smallARMA2", "mediumAR", ...
        "mediumMA1", "mediumARMA1", "mediumARMA2", "mediumMA2", ...
        "largeAR", "largeARMA"};
    case "var"
      names = {"tinyAR", "smallAR1", "smallAR2", "mediumAR", "largeAR"};
    case "paper"
      names = {"tinyAR", "tinyARMA", "smallAR1", "smallARMA2", ...
               "mediumAR", "mediumARMA2", "largeAR", "largeARMA"};
    case "platform"
      names = {"tinyARMA", "smallAR2", "mediumAR", "mediumARMA1", "largeAR"};
    otherwise
      error("benchmark_testcases:set", "Unknown default testcase set");
  end
end

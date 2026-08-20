function test_varmapack_cov2corr
  cov = cat(3, [4, 3; 3, 9], [2, -3; 6, 1.5]);
  original = cov;
  expected = cat(3, [1, 0.5; 0.5, 1], [0.5, -0.5; 1, 1/6]);
  corr = varmapack.cov2corr(cov);
  assert(isequal(size(corr), size(cov)))
  assert(norm(corr(:) - expected(:), inf) < 1e-15)
  assert(isequal(cov, original))
  corr0 = varmapack.cov2corr(cov(:, :, 1));
  assert(isequal(size(corr0), [2, 2]))
  assert(norm(corr0 - expected(:, :, 1), inf) < 1e-15)
  fprintf('test_varmapack_cov2corr passed\n');
end

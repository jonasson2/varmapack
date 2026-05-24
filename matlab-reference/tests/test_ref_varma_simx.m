function [shortDiff, longDiff, startupZ] = test_ref_varma_simx()
  fprintf('TESTING REF_VARMA_SIMX... ');
  fprintf('\n');
  % Tiny scalar ARMAX: checks closed-form moments and a scalar simx/sim match.
  test_tinyARMAX;
  % Short deterministic comparison: catches indexing errors near startup.
  shortDiff = test_ref_varma_simx_forward(1:15, 10, 7, 5*eps);
  % Long deterministic comparison: catches accumulated recursion/indexing drift.
  longDiff = test_ref_varma_simx_forward(1:15, 1000, 7, 50*eps);
  % Startup shocks: checks the conditional Gaussian distribution of eps0...eps{h-1}.
  startupZ = test_ref_varma_simx_startup;
  % Multipath startup: checks r by h by M startup paths and M by n z paths.
  test_ref_varma_simx_multipath;
  fprintf('REF_VARMA_SIMX TESTS OK\n');
end

function test_ref_varma_simx_multipath()
  n = 8; M = 3; h = 2;
  [A, B, C, Sig, z] = ref_varma_testcasex(1, n);
  z = repmat(z, M, 1) + (0:M-1)'/20;
  X0 = zeros(1, h, M);
  for j = 1:M
    X0(:,:,j) = [j/10, -j/10];
  end
  rng = randompack_create();
  cleanup = onCleanup(@() randompack_free(rng));
  randompack_seed(rng, 42);
  X = ref_varma_simx(A, B, C, z, Sig, n, M, X0, h, [], rng);
  for j = 1:M
    ascertain(almostequal(X(1:h,j)', X0(:,:,j), 1e-12));
  end
end

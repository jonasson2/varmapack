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
  fprintf('REF_VARMA_SIMX TESTS OK\n');
end

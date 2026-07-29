function [shortDiff, longDiff, startupZ] = test_ref_varma_simx()
  fprintf('TESTING REF_VARMA_SIMX... ');
  fprintf('\n');
  % Tiny scalar ARMAX: checks closed-form moments and a scalar simx/sim match.
  test_tinyARMAX;
  % Short deterministic comparison: catches indexing errors near startup.
  shortDiff = test_ref_varma_simx_forward(1:15, 10, 7, 5*eps);
  % Long deterministic comparison: catches accumulated recursion/indexing drift.
  longDiff = test_ref_varma_simx_forward(1:15, 1000, 7, 100*eps);
  % Startup shocks: checks the conditional Gaussian distribution of eps0...eps{h-1}.
  startupZ = test_ref_varma_simx_startup;
  % Multipath startup: checks r by h by M startup paths and scalar z paths.
  test_ref_varma_simx_multipath;
  % Vector exogenous input: checks d by n and d by n by M z paths.
  test_ref_varma_simx_vector;
  fprintf('REF_VARMA_SIMX TESTS OK\n');
end

function test_ref_varma_simx_multipath()
  n = 8; M = 3; h = 2;
  [A, B, Sig, p, q, r] = ref_varma_testcase(1);
  s = 2;
  [C, z] = ref_varma_testcasex(s, 1, r, n);
  z = repmat(z, 1, 1, M);
  for j = 1:M
    z(:,:,j) = z(:,:,j) + (j - 1)/20;
  end
  X0 = zeros(1, h, M);
  for j = 1:M
    X0(:,:,j) = [j/10, -j/10];
  end
  rng = randompack_create();
  cleanup = onCleanup(@() randompack_free(rng));
  randompack_seed(rng, 42);
  X = ref_varma_simx(A, B, C, z, Sig, n, M, X0, h, rng);
  for j = 1:M
    ascertain(almostequal(X(1:h,j)', X0(:,:,j), 1e-12));
  end
end

function test_ref_varma_simx_vector()
  n = 6; M = 2; h = 2;
  A = [0.3, 0.1; -0.1, 0.2];
  B = [0.2, 0; 0.1, -0.2];
  C = [0.5, -0.2, 0.1, 0.3; -0.4, 0.2, 0.2, 0.1];
  Sig = [1, 0.2; 0.2, 0.8];
  z = [0.6, -0.2, 0.8, 0, -1, 0.4; -0.2, 0.8, 0, -1, 0.4, 0.2];
  zmulti = cat(3, z, z + 0.1);
  X0 = zeros(2, h, M);
  rng = randompack_create();
  cleanup = onCleanup(@() randompack_free(rng));
  randompack_seed(rng, 42);
  [Xcommon, Ecommon] = ref_varma_simx(A, B, C, z, Sig, n, M, X0, h, rng);
  for j = 1:M
    for t = h+1:n
      expected = Ecommon(:,t,j) + A*Xcommon(:,t-1,j) + B*Ecommon(:,t-1,j) + ...
        C(:,1:2)*z(:,t) + C(:,3:4)*z(:,t-1);
      ascertain(almostequal(Xcommon(:,t,j), expected, 1e-12));
    end
  end
  randompack_seed(rng, 42);
  [X, E] = ref_varma_simx(A, B, C, zmulti, Sig, n, M, X0, h, rng);
  for j = 1:M
    for t = h+1:n
      expected = E(:,t,j) + A*X(:,t-1,j) + B*E(:,t-1,j) + ...
        C(:,1:2)*zmulti(:,t,j) + C(:,3:4)*zmulti(:,t-1,j);
      ascertain(almostequal(X(:,t,j), expected, 1e-12));
    end
  end
end

function test_varmapack_analysis
  test_scalar_models
  test_multivariate_responses
  test_sample_autocovariances
  fprintf('test_varmapack_analysis passed\n');
end

function test_scalar_models
  A = 0.5;
  Sig = 1;
  Gamma = varmapack.acvf(A, [], Sig, 4);
  expectedGamma = [4/3, 2/3, 1/3, 1/6, 1/12];
  assert(norm(squeeze(Gamma)' - expectedGamma, inf) < 1e-12)
  assert(abs(varmapack.specrad(A) - 0.5) < 1e-12)
  assert(varmapack.specrad([]) == 0)
  expectedPsi = [1, 0.5, 0.25, 0.125, 0.0625];
  assert(norm(squeeze(varmapack.psi(A, [], 4))' - expectedPsi, inf) < 1e-12)
  assert(norm(squeeze(varmapack.irf(A, [], Sig, 4))' - expectedPsi, inf) < 1e-12)
  B = 0.25;
  Gamma = varmapack.acvf([], B, 4, 3);
  assert(norm(squeeze(Gamma)' - [4.25, 1, 0, 0], inf) < 1e-12)
  assert(abs(varmapack.ma_specrad(B) - 0.25) < 1e-12)
  assert(varmapack.ma_specrad([]) == 0)
end

function test_multivariate_responses
  A = [0.1, 0.2; 0.3, 0.4];
  B = [0.5, 0.6; 0.7, 0.8];
  Sig = [1, 0.25; 0.25, 2];
  expected = cat(3, eye(2), [0.6, 0.8; 1, 1.2], ...
                 [0.26, 0.32; 0.58, 0.72]);
  Psi = varmapack.psi(A, B, 2);
  assert(isequal(size(Psi), [2, 2, 3]))
  assert(norm(Psi(:) - expected(:), inf) < 1e-12)
  Theta = varmapack.irf(A, B, Sig, 2);
  L = chol(Sig, 'lower');
  for k = 1:3
    assert(norm(Theta(:, :, k) - expected(:, :, k)*L, inf) < 1e-12)
  end
  singularSig = [1, 1; 1, 1];
  Theta = varmapack.irf(A, B, singularSig, 0);
  assert(norm(Theta*Theta' - singularSig, inf) < 1e-12)
end

function test_sample_autocovariances
  X = [1, 3, 5, 7];
  ml = varmapack.autocov(X, 2, "ML");
  corrected = varmapack.autocov(X, 2, "C");
  assert(norm(squeeze(ml)' - [5, 1.25, -1.5], inf) < 1e-12)
  assert(norm(squeeze(corrected)' - [5, 5/3, -3], inf) < 1e-12)
  X = [1:8; 2*(1:8)];
  Gamma = varmapack.autocov(X, 1);
  Xc = X - mean(X, 2);
  expected = Xc(:, 2:8)*Xc(:, 1:7)'/8;
  assert(isequal(size(Gamma), [2, 2, 2]))
  assert(norm(Gamma(:, :, 2) - expected, inf) < 1e-12)
end

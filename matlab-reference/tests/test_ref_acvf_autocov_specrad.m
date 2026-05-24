function test_ref_acvf_autocov_specrad
  fprintf('TESTING ACVF, AUTOCOV AND SPECRAD... ');
  Gamma = ref_varma_acvf(0.5, [], 0.8, 3);
  assert_equal(Gamma(:), [0.8/(1 - 0.25); 0.5*0.8/(1 - 0.25); ...
                          0.25*0.8/(1 - 0.25); 0.125*0.8/(1 - 0.25)], ...
               1e-14);
  Gamma = ref_varma_acvf([], 0.5, 0.8, 3);
  assert_equal(Gamma(:), [1; 0.4; 0; 0], 1e-14);
  [A, B, Sig] = ref_varma_testcase(8);
  Gamma = ref_varma_acvf(A, B, Sig, 6);
  ascertain(isequal(size(Gamma), [2, 2, 7]));
  ascertain(almostequal(Gamma(:, :, 1), Gamma(:, :, 1)'));
  X = [1 2 4 7 11; 3 5 9 15 23];
  C = ref_varma_autocov(X, 2, "ML");
  assert_equal(C, direct_autocov(X, 2, true), 1e-14);
  C = ref_varma_autocov(X, 2, "corrected");
  assert_equal(C, direct_autocov(X, 2, false), 1e-14);
  ascertain(ref_varma_specrad(zeros(2, 0)) == 0);
  ascertain(almostequal(ref_varma_specrad(0.5), 0.5));
  A = [0.6 0.1; 0.2 0.3];
  ascertain(almostequal(ref_varma_specrad(A), max(abs(eig(A)))));
  ascertain(ref_varma_ma_specrad(zeros(2, 0)) == 0);
  ascertain(almostequal(ref_varma_ma_specrad(0.7), 0.7));
  ascertain(almostequal(ref_varma_ma_specrad(-1.2), 1.2));
  B = [0.4 -0.12];
  ascertain(almostequal(ref_varma_ma_specrad(B), max(abs(roots([1 B])))));
  B = diag([0.3 -0.8]);
  ascertain(almostequal(ref_varma_ma_specrad(B), 0.8));
  for icase = 1:ref_varma_testcase('count')
    [A, B] = ref_varma_testcase(icase);
    rho = ref_varma_ma_specrad(B);
    ascertain(isfinite(rho) && rho >= 0);
  end
  fprintf('OK\n');
end

function C = direct_autocov(X, maxlag, ML)
  [r, n] = size(X);
  X = X - mean(X, 2);
  C = zeros(r, r, maxlag + 1);
  for k = 0:maxlag
    if ML, f = n; else, f = n - k; end
    C(:, :, k+1) = X(:, 1:n-k)*X(:, k+1:n)'/f;
  end
end

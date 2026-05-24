function test_ref_psi_irf
  fprintf('TESTING PSI AND IRF... ');
  Psi = ref_varma_psi(0.5, 0.2, 4);
  assert_equal(reshape(Psi, 1, []), [1 0.7 0.35 0.175 0.0875], 1e-15);
  B = reshape(cat(3, [1 2; 3 4], [5 6; 7 8]), 2, 4);
  Psi = ref_varma_psi([], B, 3);
  assert_equal(Psi(:, :, 1), eye(2));
  assert_equal(Psi(:, :, 2), B(:, 1:2));
  assert_equal(Psi(:, :, 3), B(:, 3:4));
  assert_equal(Psi(:, :, 4), zeros(2));
  A = [0.4 0.1; 0.2 0.3];
  B = [0.1 0.2; 0.3 0.4];
  Psi = ref_varma_psi(A, B, 3);
  assert_equal(Psi(:, :, 1), eye(2));
  assert_equal(Psi(:, :, 2), A + B, 1e-15);
  assert_equal(Psi(:, :, 3), A*(A + B), 1e-15);
  assert_equal(Psi(:, :, 4), A*A*(A + B), 1e-15);
  Sig = [2 0.5; 0.5 1];
  Theta = ref_varma_irf(A, B, Sig, 2);
  L = chol(Sig, 'lower');
  for j = 1:3
    assert_equal(Theta(:, :, j), Psi(:, :, j)*L, 1e-15);
  end
  Sigs = cat(3, [1 1; 1 1], [1 2; 2 4]);
  for i = 1:2
    Theta = ref_varma_irf([], [], Sigs(:, :, i), 0);
    assert_equal(Theta(:, :, 1)*Theta(:, :, 1)', Sigs(:, :, i), 1e-12);
  end
  fprintf('OK\n');
end

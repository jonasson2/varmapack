function test_ref_varma_sim_nonstationary
  rng = randompack_create();
  cleanup = onCleanup(@() randompack_free(rng));
  A = 1.25;
  B = 0.5;
  Sig = 1;
  X0 = [2, 3, 4];
  [X, E] = ref_varma_sim(A, B, Sig, [], 6, 1, X0, rng);
  assert(norm(X(:,1:3) - X0, inf) < 1e-12)
  for t = 2:6
    expected = A*X(:,t-1) + E(:,t) + B*E(:,t-1);
    assert(norm(X(:,t) - expected, inf) < 1e-12)
  end
  fprintf('test_ref_varma_sim_nonstationary passed\n');
end

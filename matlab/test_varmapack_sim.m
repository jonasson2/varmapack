function test_varmapack_sim
  A = 1.25;
  B = 0.5;
  Sig = 1;
  X0 = [2, 3, 4];
  rng = varmapack.Rng(46);
  cleanup = onCleanup(@() delete(rng));
  [X, E] = varmapack.sim(A, B, Sig, [], 6, 1, X0, rng);
  assert(norm(X(:,1:3) - X0, inf) < 1e-12)
  for t = 2:6
    expected = A*X(:,t-1) + E(:,t) + B*E(:,t-1);
    assert(norm(X(:,t) - expected, inf) < 1e-12)
  end
  fprintf('test_varmapack_sim passed\n');
end

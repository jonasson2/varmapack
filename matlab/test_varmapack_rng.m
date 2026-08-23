function test_varmapack_rng
  %TEST_VARMAPACK_RNG  Check public MATLAB RNG ownership and simulation use.

  A = [0.4, 0.1; 0, 0.3];
  B = [0.1, 0; 0, 0.2];
  C = [0.4; 0.1];
  Sig = [1, 0.2; 0.2, 0.8];
  rng1 = varmapack.Rng(123);
  rng2 = varmapack.Rng(123);
  cleanup1 = onCleanup(@() delete(rng1));
  cleanup2 = onCleanup(@() delete(rng2));
  X1 = varmapack.sim(A, B, Sig, [], 12, 2, [], rng1);
  X2 = varmapack.sim(A, B, Sig, [], 12, 2, [], rng2);
  assert(isequal(X1, X2))
  rng1.seed(456);
  X1 = varmapack.sim(A, B, Sig, [], 12, 2, [], rng1);
  rng1.seed(456);
  X2 = varmapack.sim(A, B, Sig, [], 12, 2, [], rng1);
  assert(isequal(X1, X2))
  X3 = varmapack.sim(A, B, Sig, [], 12, 2, [], rng1);
  X4 = varmapack.sim(A, B, Sig, [], 12, 2, [], rng1);
  assert(~isequal(X3, X4))
  z = cos((0:11)'/5);
  X = varmapack.simx(A, B, C, z, Sig, 12, 2, zeros(2), 2, rng1);
  assert(isequal(size(X), [2, 12, 2]))
  X = varmapack.sim(A, [], Sig, [], 12);
  assert(isequal(size(X), [2, 12]))
  X = varmapack.simx(A, B, C, z, Sig, 12, 2, zeros(2), 2);
  assert(isequal(size(X), [2, 12, 2]))
  fprintf('test_varmapack_rng passed\n');
end

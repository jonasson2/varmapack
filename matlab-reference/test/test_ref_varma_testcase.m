function test_ref_varma_testcase
  n = ref_varma_testcase('count');
  assert(n == 15);
  [A, B, Sig, p, q, r, name] = ref_varma_testcase(8);
  assert(strcmp(name, 'smallARMA1'));
  assert(isequal([p q r], [1 1 2]));
  [A2, B2, Sig2, p2, q2, r2, icase] = ref_varma_testcase('smallARMA1');
  assert(icase == 8);
  assert(isequal([p2 q2 r2], [1 1 2]));
  assert(max(abs(A(:) - A2(:))) == 0);
  assert(max(abs(B(:) - B2(:))) == 0);
  assert(max(abs(Sig(:) - Sig2(:))) == 0);
  [A, B, Sig, p, q, r] = ref_varma_testcase(2, 1, 3);
  assert(isequal([p q r], [2 1 3]));
  assert(all(A(:) == 0.5/(2*3)));
  assert(all(B(:) == 1/(1*3)));
  assert(max(abs(Sig(:) - reshape(hilb(3) + 0.2*eye(3), [], 1))) == 0);
  rng = randompack_create();
  cleanup = onCleanup(@() randompack_free(rng));
  randompack_seed(rng, 123);
  [A1, B1, Sig1] = ref_varma_testcase(0, 2, 2, rng);
  randompack_seed(rng, 123);
  [A2, B2, Sig2] = ref_varma_testcase(0, 2, 2, rng);
  assert(isempty(A1) && isempty(A2));
  assert(max(abs(B1(:) - B2(:))) == 0);
  assert(max(abs(Sig1(:) - Sig2(:))) == 0);
end

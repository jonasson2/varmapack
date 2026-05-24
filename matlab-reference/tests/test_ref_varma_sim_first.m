function test_ref_varma_sim_first
  fprintf('TESTING REF_VARMA_SIM (PRELIMINARY)... ');
  rng = randompack_create();
  cleanup = onCleanup(@() randompack_free(rng));
  for i = 1:ref_varma_testcase('count')
    [A, B, Sig, p, q, r] = ref_varma_testcase(i);
    randompack_seed(rng, 42);
    X1 = ref_varma_sim(A, B, Sig, [], 10, 1, [], [], rng);
    randompack_seed(rng, 42);
    X2 = ref_varma_sim(A, B, Sig, [], 10, 1, [], [], rng);
    randompack_seed(rng, 42);
    X3 = ref_varma_sim(A, B, Sig, [], 10, 5, [], [], rng);
    randompack_seed(rng, 42);
    [X4, E4] = ref_varma_sim(A, B, Sig, [], 10, 1, [], [], rng);
    randompack_seed(rng, 42);
    [X5, E5] = ref_varma_sim(A, B, Sig, [], 10, 5, [], [], rng);
    assertsize(X1, r, 10);
    assertsize(X2, r, 10);
    assertsize(X4, r, 10);
    assertsize(E4, r, 10);
    if r == 1
      assertsize(X3, 10, 5);
      assertsize(X5, 10, 5);
      assertsize(E5, 10, 5);
    else
      assertsize(X3, r, 10, 5);
      assertsize(X5, r, 10, 5);
      assertsize(E5, r, 10, 5);
    end
  end
  fprintf('OK\n');
end

function assertsize(A, n, m, k)
  if nargin == 3 || k == 1
    ascertain(isequal(size(A), [n, m]));
  else
    ascertain(isequal(size(A), [n, m, k]));
  end
end

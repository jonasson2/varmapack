function test_ref_varma_sim
  fprintf('TESTING REF_VARMA_SIM... ');
  rng = randompack_create();
  cleanup = onCleanup(@() randompack_free(rng));
  cases = [1 3 5 7 12];
  for i = cases
    [A, B, Sig, p, q, r] = ref_varma_testcase(i);
    h = max(p, q);
    mu = (1:r)'/10;
    randompack_seed(rng, 42);
    [X, E] = ref_varma_sim(A, B, Sig, h + 4, [], 3, [], rng);
    ascertain(all(isfinite(X(:))));
    ascertain(all(isfinite(E(:))));
    randompack_seed(rng, 42);
    [Xmu, Emu] = ref_varma_sim(A, B, Sig, h + 4, mu, 3, [], rng);
    ascertain(almostequal(Emu, E));
    if r == 1
      X = reshape(X, h + 4, 3);
      Xmu = reshape(Xmu, h + 4, 3);
      ascertain(almostequal(Xmu - X, repmat(mu, h + 4, 3), 1e-12));
    else
      ascertain(almostequal(Xmu - X, repmat(mu, [1, h + 4, 3]), 1e-12));
    end
    X0 = 2 + repmat((1:r)', 1, h);
    randompack_seed(rng, 42);
    Xstart = ref_varma_sim(A, B, Sig, h + 3, mu, 2, X0, rng);
    if r == 1
      Xstart = reshape(Xstart, 1, h + 3, 2);
    end
    ascertain(almostequal(Xstart(:, 1:h, 1), X0, 1e-12));
    ascertain(almostequal(Xstart(:, 1:h, 2), X0, 1e-12));
  end
  fprintf('OK\n');
end

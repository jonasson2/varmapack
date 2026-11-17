function test_varmapack_cov2sig
  assert(abs(varmapack.cov2sig(0.5, [], 4) - 3) < 1e-12)
  assert(abs(varmapack.cov2sig([], 0.25, 2.125) - 2) < 1e-12)
  assert(abs(varmapack.cov2sig(0.5, 0.25, 3.5) - 2) < 1e-12)
  A = [0.4, 0.2; -0.1, 0.3];
  cov = [2, 0.3; 0.3, 1];
  expected = cov - A*cov*A';
  assert(norm(varmapack.cov2sig([], [], cov) - cov, inf) < 1e-12)
  cov(1, 2) = NaN;
  original = cov;
  Sig = varmapack.cov2sig(A, [], cov);
  assert(isequal(size(Sig), [2, 2]))
  assert(norm(Sig - expected, inf) < 1e-12)
  assert(isequaln(cov, original))
  for value = [0, 1]
    cov = value*ones(2);
    assert(norm(varmapack.cov2sig(0.5*eye(2), [], cov) - 0.75*cov, inf) < 1e-12)
  end
  for index = 1:varmapack.testcase('count')
    [A, B, Sig] = varmapack.testcase(index);
    cov = varmapack.acvf(A, B, Sig, 0);
    original = cov;
    recovered = varmapack.cov2sig(A, B, cov);
    assert(norm(recovered - Sig, inf) < 1e-9)
    roundtrip = varmapack.acvf(A, B, recovered, 0);
    assert(norm(roundtrip - cov, inf) < 1e-9*max(1, norm(cov, inf)))
    assert(isequal(cov, original))
  end
  varmapack.cov2sig([], [], eye(2)); % No output is also supported.
  fprintf('test_varmapack_cov2sig passed\n');
end

function test_varmapack_testcase
  assert(varmapack.testcase('count') == 16)
  [A, B, Sig, p, q, r, index] = varmapack.testcase('smallARMA1');
  [indexedA, indexedB, indexedSig, indexedP, indexedQ, indexedR, name] = ...
    varmapack.testcase(index);
  assert(index == 8 && strcmp(name, 'smallARMA1'))
  assert(isequal([p, q, r], [1, 1, 2]))
  assert(isequal([indexedP, indexedQ, indexedR], [p, q, r]))
  assert(isequal(A, indexedA) && isequal(B, indexedB) && isequal(Sig, indexedSig))
  check_catalog_entry('tinyAR', 1, 1, 0, 1)
  check_catalog_entry('largeAR', 15, 5, 0, 7)
  check_catalog_entry('largeARMA', 16, 3, 3, 7)
  [A, B, Sig] = varmapack.testcase(2, 1, 2);
  assert(isequal(size(A), [2, 4]))
  assert(isequal(size(B), [2, 2]))
  assert(isequal(size(Sig), [2, 2]))
  [A, ~, ~] = varmapack.testcase(3, 1, 2, 0.5);
  assert(abs(varmapack.specrad(A) - 0.5) < 2e-6)
  rng1 = varmapack.Rng(321);
  rng2 = varmapack.Rng(321);
  cleanup1 = onCleanup(@() delete(rng1));
  cleanup2 = onCleanup(@() delete(rng2));
  [A1, B1, Sig1] = varmapack.testcase(2, 1, 2, rng1);
  [A2, B2, Sig2] = varmapack.testcase(2, 1, 2, rng2);
  assert(isequal(A1, A2) && isequal(B1, B2) && isequal(Sig1, Sig2))
  summary = evalc("varmapack.testcase('summary')");
  assert(contains(summary, 'tinyAR') && contains(summary, 'largeARMA'))
  fprintf('test_varmapack_testcase passed\n');
end

function check_catalog_entry(name, expectedIndex, expectedP, expectedQ, expectedR)
  [~, ~, ~, p, q, r, index] = varmapack.testcase(name);
  assert(isequal([index, p, q, r], ...
                 [expectedIndex, expectedP, expectedQ, expectedR]))
end

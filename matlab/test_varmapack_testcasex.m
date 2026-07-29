function test_varmapack_testcasex
  [C, z] = varmapack.testcasex(3, 1, 2, 20);
  expectedC = [0.3 -0.3 0.2; -0.5 0 0.5];
  expectedZ = [0.6 -0.2 0.8 0 -1 0.4 0.2 -0.6 1 -0.4 0 0.8 -0.8 0.2 ...
               -0.2 0.4 -0.6 0.6 -0.2 0.8];
  assert(isequal(C, expectedC))
  assert(isequal(z, expectedZ))
  [C, z] = varmapack.testcasex(0, 0, 3, 0);
  assert(isequal(size(C), [3 0]))
  assert(isempty(z))
  [C, z] = varmapack.testcasex(2, 2, 2, 3);
  expectedC = [0.3 -0.1 -0.3 0.4; -0.5 0.2 0 -0.4];
  expectedZ = [0.6 -0.2 0.8; -0.2 0.8 0];
  assert(isequal(size(C), [2 4]))
  assert(isequal(C, expectedC))
  assert(isequal(z, expectedZ))
  fprintf('test_varmapack_testcasex passed\n');
end

function test_varmapack_testcasex
  [C, z] = varmapack.testcasex(3, 2, 20);
  expectedC = [0.3 -0.3 0.2; -0.5 0 0.5];
  expectedZ = [0.6 -0.2 0.8 0 -1 0.4 0.2 -0.6 1 -0.4 0 0.8 -0.8 0.2 ...
               -0.2 0.4 -0.6 0.6 -0.2 0.8]';
  assert(isequal(C, expectedC))
  assert(isequal(z, expectedZ))
  [C, z] = varmapack.testcasex(0, 3, 0);
  assert(isequal(size(C), [3 0]))
  assert(isempty(z))
  fprintf('test_varmapack_testcasex passed\n');
end

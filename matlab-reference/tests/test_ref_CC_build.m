function test_ref_CC_build
  fprintf('TESTING CC_BUILD... ');
  A = [0.5 0.3];
  C = {4, 5, 6};
  CC = CC_build(A, C, 2);
  assert_equal(CC, [4 0; 5 4]);
  CC = CC_build(A, C, 3);
  assert_equal(CC, [4 0 0; 5 4 0; 6 5 4]);
  A = 0.5;
  CC = CC_build(A, C, 3);
  assert_equal(CC, [4 0 0; 5 4 0; 6 5 4]);
  A = [0.5 0.3];
  C = {2};
  CC = CC_build(A, C, 3);
  assert_equal(CC, [2 0 0; 1 2 0; 1.1 1 2], 1e-15);
  C0 = ones(2);
  C1 = [5 6; 7 8];
  C = {C0, C1};
  A = [0.4 0.3 0.2 0.1; 0.1 0.1 0.2 0.3];
  CC = CC_build(A, C, 2);
  assert_equal(CC, [1 1 0 0; 1 1 0 0; 5 6 1 1; 7 8 1 1]);
  CC = CC_build(A, C, 3);
  CCok = [1 1 0 0 0 0; 1 1 0 0 0 0; 5 6 1 1 0 0; 7 8 1 1 0 0;
          4.4 5.1 5 6 1 1; 1.7 1.9 7 8 1 1];
  assert_equal(CC, CCok, 1e-15);
  fprintf('OK\n');
end

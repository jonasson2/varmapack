function test_ref_S_build
  fprintf('TESTING S_BUILD... ');
  S = {2, 3, 1};
  G = {4, 5, 0};
  A = [0.5 0.3];
  SS = S_build(S, A, G, 3);
  SSok = [2, 3, 1; 3, 2, 3; 1, 3, 2];
  assert_equal(SS, SSok);
  S0 = eye(2);
  S1 = [1 2; 3 4];
  G0 = ones(2);
  G1 = [5 6; 7 8];
  S = {S0, S1};
  G = {G0, G1};
  A = [0.4 0.3; 0.1 0.1];
  SS = S_build(S, A, G, 3);
  S2 = [1.3 2; 0.4 0.6];
  SSok = [S0 S1' S2'; S1 S0 S1'; S2 S1 S0];
  assert_equal(SS, SSok, 1e-15);
  fprintf('OK\n');
end

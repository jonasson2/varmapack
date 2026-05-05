function SmokeTest
  rng = randompack_create();
  cleanup = onCleanup(@() randompack_free(rng));
  randompack_seed(rng, 42);
  u = randompack_u01(rng, 3);
  z = randompack_normal(rng, 3, 2, 3);
  Sig = [2, 1; 1, 2];
  X = randompack_mvn(rng, Sig, 3, [0; 0]);
  disp('u01 =');
  disp(u');
  disp('normal N(2,3) =');
  disp(z');
  disp('mvn Sigma=[2 1; 1 2] =');
  disp(X);
end

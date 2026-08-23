function test_varmapack_simx
  A = [0.3, 0.1; -0.1, 0.2];
  B = [0.2, 0; 0.1, -0.2];
  C = [0.5, -0.2, 0.1, 0.3; -0.4, 0.2, 0.2, 0.1];
  Sig = [1, 0.2; 0.2, 0.8];
  z = [0.6, -0.2, 0.8, 0, -1, 0.4; -0.2, 0.8, 0, -1, 0.4, 0.2];
  zmulti = cat(3, z, z + 0.1);
  n = 6; M = 2; h = 2;
  rng = varmapack.Rng(42);
  cleanup = onCleanup(@() delete(rng));
  X0 = zeros(2, h, M);
  [Xcommon, Ecommon] = varmapack.simx(A, B, C, z, Sig, n, M, X0, h, rng);
  difference = Xcommon(:, 1:h, :) - X0;
  assert(norm(difference(:), inf) < 1e-12)
  for j = 1:M
    for t = h+1:n
      expected = Ecommon(:,t,j) + A*Xcommon(:,t-1,j) + B*Ecommon(:,t-1,j) + ...
        C(:,1:2)*z(:,t) + C(:,3:4)*z(:,t-1);
      assert(norm(Xcommon(:,t,j) - expected, inf) < 1e-12)
    end
  end
  rng.seed(42)
  [X, E] = varmapack.simx(A, B, C, zmulti, Sig, n, M, X0, h, rng);
  assert(isequal(size(X), [2, n, M]))
  assert(isequal(size(E), [2, n, M]))
  difference = X(:, 1:h, :) - X0;
  assert(norm(difference(:), inf) < 1e-12)
  for j = 1:M
    for t = h+1:n
      expected = E(:,t,j) + A*X(:,t-1,j) + B*E(:,t-1,j) + ...
        C(:,1:2)*zmulti(:,t,j) + C(:,3:4)*zmulti(:,t-1,j);
      assert(norm(X(:,t,j) - expected, inf) < 1e-12)
    end
  end
  rng.seed(43)
  hminimum = 1;
  X0minimum = zeros(2, hminimum, M);
  [Xminimum, Eminimum] = varmapack.simx(A, B, C, z, Sig, n, M, X0minimum, ...
                                         hminimum, rng);
  for j = 1:M
    for t = 2:n
      expected = Eminimum(:,t,j) + A*Xminimum(:,t-1,j) + ...
        B*Eminimum(:,t-1,j) + C(:,1:2)*z(:,t) + C(:,3:4)*z(:,t-1);
      assert(norm(Xminimum(:,t,j) - expected, inf) < 1e-12)
    end
  end
  rng.seed(44)
  Czero = 0.5*eye(2);
  [Xzero, Ezero] = varmapack.simx([], [], Czero, z, Sig, n, M, [], 0, rng);
  for j = 1:M
    assert(norm(Xzero(:,:,j) - Ezero(:,:,j) - Czero*z, inf) < 1e-12)
  end
  rng.seed(45)
  Cscalar = [0.4, -0.2];
  zscalar = (0:n-1)';
  [Xscalar, Escalar] = varmapack.simx([], [], Cscalar, zscalar, 1, n, 1, ...
                                      0, 1, rng);
  assert(isequal(size(Xscalar), [1, n]))
  for t = 2:n
    expected = Escalar(t) + Cscalar(1)*zscalar(t) + Cscalar(2)*zscalar(t-1);
    assert(abs(Xscalar(t) - expected) < 1e-12)
  end
  zscalar = [zscalar, zscalar + 1];
  [Xscalar, Escalar] = varmapack.simx([], [], Cscalar, zscalar, 1, n, M, ...
                                      0, 1, rng);
  assert(isequal(size(Xscalar), [n, M]))
  for j = 1:M
    for t = 2:n
      expected = Escalar(t,j) + Cscalar(1)*zscalar(t,j) + ...
        Cscalar(2)*zscalar(t-1,j);
      assert(abs(Xscalar(t,j) - expected) < 1e-12)
    end
  end
  fprintf('test_varmapack_simx passed\n');
end

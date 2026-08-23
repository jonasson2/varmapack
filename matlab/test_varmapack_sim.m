function test_varmapack_sim
  test_white_noise_and_mean
  test_starting_values
  test_nonstationary_varma
  fprintf('test_varmapack_sim passed\n');
end

function test_white_noise_and_mean
  Sig = [1, 0.25; 0.25, 2];
  rng = varmapack.Rng(45);
  cleanup = onCleanup(@() delete(rng));
  [X, E] = varmapack.sim([], [], Sig, [], 5, 3, [], rng);
  assert(isequal(size(X), [2, 5, 3]))
  assert(norm(X(:) - E(:), inf) < 1e-12)
  [X, E] = varmapack.sim([], [], 1, [], 5, 1, [], rng);
  assert(isequal(size(X), [1, 5]) && isequal(size(E), [1, 5]))
  [X, E] = varmapack.sim([], [], 1, [], 5, 3, [], rng);
  assert(isequal(size(X), [5, 3]) && isequal(size(E), [5, 3]))
  mu = [10, 11, 12; 20, 21, 22];
  [X, E] = varmapack.sim([], [], Sig, mu, 5, 1, [], rng);
  expected = [mu, mu(:, end), mu(:, end)];
  assert(isequal(size(X), [2, 5]))
  assert(norm(X - E - expected, inf) < 1e-12)
  singularSig = [1, 1; 1, 1];
  [X, E] = varmapack.sim([], [], singularSig, [], 5, 2, [], rng);
  assert(norm(X(:) - E(:), inf) < 1e-12)
  difference = X(1, :, :) - X(2, :, :);
  assert(norm(difference(:), inf) < 1e-12)
end

function test_starting_values
  A = [0.1, 0.2; 0.3, 0.4];
  B = [0.5, 0.6; 0.7, 0.8];
  Sig = [1, 0.25; 0.25, 2];
  rng = varmapack.Rng(46);
  cleanup = onCleanup(@() delete(rng));
  X0 = [1, 2; 3, 4];
  X = varmapack.sim(A, B, Sig, [], 5, 2, X0, rng);
  assert(norm(X(:, 1:2, 1) - X0, inf) < 1e-12)
  assert(norm(X(:, 1:2, 2) - X0, inf) < 1e-12)
  X0 = cat(3, [1; 2], [3; 4], [5; 6]);
  X = varmapack.sim(A, B, Sig, [], 5, 3, X0, rng);
  assert(norm(reshape(X(:, 1, :), 2, 3) - reshape(X0, 2, 3), inf) < 1e-12)
end

function test_nonstationary_varma
  A = 1.25;
  B = 0.5;
  Sig = 1;
  X0 = [2, 3, 4];
  rng = varmapack.Rng(46);
  cleanup = onCleanup(@() delete(rng));
  [X, E] = varmapack.sim(A, B, Sig, [], 6, 1, X0, rng);
  assert(norm(X(:,1:3) - X0, inf) < 1e-12)
  for t = 2:6
    expected = A*X(:,t-1) + E(:,t) + B*E(:,t-1);
    assert(norm(X(:,t) - expected, inf) < 1e-12)
  end
end

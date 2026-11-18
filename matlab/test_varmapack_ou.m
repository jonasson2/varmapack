function test_varmapack_ou
  rng = varmapack.Rng(423);
  model = varmapack.OUModel([1; 2], [1 .3; .3 2], [2; -1]);
  t = [-3, -2.9, -1, 3];
  M = 40000;
  X = model.sim(t, M, rng);
  assert(isequal(size(X), [2 4 M]));
  for k = 1:4
    x = reshape(X(:,k,:), 2, M) - model.mu;
    assert(max(abs(mean(x, 2))) < .03);
    for h = 1:k
      y = reshape(X(:,h,:), 2, M) - model.mu;
      expected = exp(-(t(k)-t(h))./model.tau).*model.cov;
      assert(max(abs(x*y'/M-expected), [], 'all') < .05);
    end
  end
  rng.seed(423);
  assert(isequal(X, model.sim(t, M, rng)));
  assert(isequal(size(model.sim(0)), [2 1]));
  model = varmapack.OUModel([1, 1], ones(2));
  X = model.sim([0 1e-15 10000], 10, rng);
  assert(max(abs(X(1,:,:)-X(2,:,:)), [], 'all') < 1e-12);
  model = varmapack.OUModel([1 2], zeros(2), [3 4]);
  X = model.sim(t);
  assert(isequal(X, repmat([3;4], 1, 4)));
  model = varmapack.OUModel([1 2], [1 NaN; .2 1]);
  assert(isequal(model.cov, [1 .2; .2 1]));
  must_fail(@() varmapack.OUModel([1 2], ones(2)));
  must_fail(@() varmapack.OUModel(0, 1));
  must_fail(@() varmapack.OUModel(Inf, 1));
  must_fail(@() varmapack.OUModel(1, -1));
  must_fail(@() varmapack.OUModel([1 2], 1));
  K = [1 -1; 0 1];
  model = varmapack.OUModel(K, eye(2), [], 'drift');
  assert(isempty(model.tau) && isequal(model.drift, K));
  X = model.sim([0 .5], 20000, rng);
  empirical = squeeze(X(:,2,:))*squeeze(X(:,1,:))'/20000;
  assert(max(abs(empirical-exp(-.5)*[1 .5; 0 1]), [], 'all') < .04);
  must_fail(@() varmapack.OUModel([1 3; 0 1], eye(2), [], 'drift'));
  must_fail(@() model.sim([0 0]));
  must_fail(@() model.sim([0 NaN]));
  must_fail(@() model.sim([0 1], 0));
  for shape = [.5 1 4 16]
    rng.seed(11);
    t = varmapack.gamma_times(50001, shape, 2, rng);
    d = diff(t);
    assert(t(1) == 0 && all(d > 0));
    assert(abs(mean(d)-2) < 7*sqrt(4/shape/numel(d)));
    assert(abs(var(d)-4/shape) < .08*4/shape);
    rng.seed(11);
    assert(isequal(t, varmapack.gamma_times(numel(t), shape, 2, rng)));
  end
  assert(isequal(varmapack.gamma_times(5, Inf, 2), (0:2:8)'));
  assert(isequal(varmapack.gamma_times(1), 0));
  must_fail(@() varmapack.gamma_times(0));
  must_fail(@() varmapack.gamma_times(3, 0));
  must_fail(@() varmapack.gamma_times(3, NaN));
  must_fail(@() varmapack.gamma_times(3, 4, 0));
  fprintf('OU and Gamma-time tests passed\n');
end

function must_fail(f)
  failed = false;
  try
    f();
  catch
    failed = true;
  end
  assert(failed, 'Expected an error');
end

function test_varmapack_slepian
  theta = [0 .5; .25 0];
  cov = diag([1.25 1.0625]);
  mu = [2; -1];
  model = varmapack.SlepianModel(theta, cov, mu);
  rng = varmapack.Rng(423);
  t = [0 .5 1 2];
  M = 20000;
  X = model.sim(t, M, rng);
  assert(isequal(size(X), [2 4 M]));
  assert(isequal(model.theta, theta) && isequal(model.cov, cov));
  expected = {cov, (cov+theta)/2, theta, zeros(2)};
  first = squeeze(X(:,1,:)) - mu;
  for k = 1:numel(t)
    empirical = (squeeze(X(:,k,:))-mu)*first'/M;
    delta = empirical - expected{k};
    assert(max(abs(delta(:))) < .035);
  end
  rng.seed(423);
  assert(isequal(X, model.sim(t, M, rng)));
  must_fail(@() model.sim([0 0]));
  must_fail(@() varmapack.SlepianModel([0 1; 0 0], diag([1 2])));
  model = varmapack.SlepianModel(0, 0, 3);
  assert(isequal(model.sim(0), 3));
  fprintf('Slepian model tests passed\n');
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

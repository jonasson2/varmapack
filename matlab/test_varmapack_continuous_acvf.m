function test_varmapack_continuous_acvf
  ou = varmapack.OUModel([1 -1; 0 1], eye(2), [], 'drift');
  theta = [0 .5; .25 0];
  cov = eye(2) + theta*theta';
  ma = varmapack.SlepianModel(theta, cov);
  check_lags(ou, exp(-.5)*[1 .5; 0 1], eye(2));
  assert(abs(ou.specrad(.5) - exp(-.5)) < 1e-12);
  check_lags(ma, (cov+theta)/2, cov);
  assert(isequal(ma.acvf(2), zeros(2)));
  fprintf('Continuous-time autocovariance tests passed\n');
end

function check_lags(model, half, cov)
  G = model.acvf([.5; -.5; 0]);
  assert(isequal(size(G), [2 2 3]));
  assert(max(abs(G - cat(3, half, half', cov)), [], 'all') < 1e-12);
  assert(isequal(size(model.acvf([])), [2 2 0]));
  failed = false;
  try
    model.acvf(Inf);
  catch
    failed = true;
  end
  assert(failed, 'Expected nonfinite lags to be rejected');
end

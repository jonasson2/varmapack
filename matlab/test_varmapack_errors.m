function test_varmapack_errors
  test_analysis_errors
  test_simulation_errors
  test_testcase_errors
  test_rng_errors
  fprintf('test_varmapack_errors passed\n');
end

function test_analysis_errors
  assert_error(@() varmapack.acvf(zeros(2, 3), [], eye(2), 1), ...
               'varmapack:acvf:A')
  assert_error(@() varmapack.acvf([], [], eye(2), -1), ...
               'varmapack:acvf:argument')
  assert_error(@() varmapack.psi(zeros(2, 3), [], 1), ...
               'varmapack:psi:argument')
  assert_error(@() varmapack.irf([], [], [1, 2; 2, 1], 0), ...
               'varmapack:irf:error')
  assert_error(@() varmapack.specrad(zeros(2, 3)), ...
               'varmapack:specrad:A')
  assert_error(@() varmapack.ma_specrad(zeros(2, 3)), ...
               'varmapack:ma_specrad:B')
  assert_error(@() varmapack.autocov(zeros(2, 2, 2), 1), ...
               'varmapack:autocov:argument')
  assert_error(@() varmapack.autocov(zeros(2, 3), 3), ...
               'varmapack:autocov:error')
  assert_error(@() varmapack.autocov(zeros(2, 3), 1, "bad"), ...
               'varmapack:autocov:error')
  assert_error(@() varmapack.cov2corr([0, 0; 0, 1]), ...
               'varmapack:cov2corr:error')
end

function test_simulation_errors
  rng = varmapack.Rng(1);
  cleanup = onCleanup(@() delete(rng));
  assert_error(@() varmapack.sim(zeros(2, 3), [], eye(2), [], 4, 1, [], rng), ...
               'varmapack:sim:A')
  assert_error(@() varmapack.sim(0.4, [], 1, [], 0, 1, [], rng), ...
               'varmapack:sim:error')
  assert_error(@() varmapack.sim(0.4, [], 1, [], 4, 0, [], rng), ...
               'varmapack:sim:error')
  assert_error(@() varmapack.sim(0.4, [], 1, [], 4, 1, [], struct), ...
               'varmapack:sim:rng')
  assert_error(@() varmapack.sim(1.1, [], 1, [], 4, 1, [], rng), ...
               'varmapack:sim:error')
  assert_error(@() varmapack.simx(0.2, [], 0.3, [], 1, 4, 1, 0, 1, rng), ...
               'varmapack:simx:z')
  assert_error(@() varmapack.simx(0.2, [], 0.3, zeros(2, 4), 1, 4, 1, ...
                                  0, 1, rng), 'varmapack:simx:C')
end

function test_testcase_errors
  assert_error(@() three_outputs('does-not-exist'), 'varmapack:testcase:error')
  assert_error(@() three_outputs(-1), 'varmapack:testcase:argument')
  assert_error(@() rho_outputs(-1), 'varmapack:testcase:argument')
  assert_error(@() rho_outputs(NaN), 'varmapack:testcase:error')
end

function test_rng_errors
  assert_error(@() varmapack.Rng(1.5), 'varmapack:rng:seed')
  rng = varmapack.Rng(1);
  delete(rng)
  assert_error(@() rng.seed(1), 'MATLAB:class:InvalidHandle')
end

function three_outputs(which)
  [~, ~, ~] = varmapack.testcase(which);
end

function rho_outputs(rho)
  [~, ~, ~] = varmapack.testcase(1, 1, 1, rho);
end

function assert_error(f, identifier)
  try
    f();
  catch exception
    assert(strcmp(exception.identifier, identifier), ...
           'Expected %s, got %s', identifier, exception.identifier)
    return
  end
  error('Expected error %s', identifier)
end

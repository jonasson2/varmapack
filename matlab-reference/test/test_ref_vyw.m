function test_ref_vyw
  fprintf('TESTING VYW_FACTORIZE AND VYW_SOLVE... ');
  cases = {1, 3, 5, 7, 12, [1 1 2], [2 1 2], [2 3 2]};
  for i = 1:length(cases)
    if isscalar(cases{i})
      [A, B, Sig] = ref_varma_testcase(cases{i});
    else
      dims = num2cell(cases{i});
      rng = randompack_create();
      cleanup = onCleanup(@() randompack_free(rng));
      randompack_seed(rng, 42);
      [A, B, Sig] = ref_varma_testcase(dims{:}, rng);
    end
    [~, G] = find_CG(A, B, Sig);
    PLU = vyw_factorize(A);
    S = vyw_solve(A, PLU, G);
    check_solution(A, G, S);
  end
  fprintf('OK\n');
end

function check_solution(A, G, S)
  q = size(G, 2)/size(G, 1) - 1;
  r = size(G, 1);
  p = size(A, 2)/r;
  A = makecell(A);
  if ~iscell(S), S = makecell(S); end
  G = makecell(G);
  S0 = S{1};
  S = S(2:end);
  G0 = G{1};
  G = G(2:end);
  maxres = 0;
  for i = 0:p
    if i == 0
      sum = S0 - G0;
    elseif i <= q
      sum = S{i} - G{i};
    else
      sum = S{i};
    end
    for j = 1:p
      if j < i, sum = sum - A{j}*S{i-j}; end
      if j == i, sum = sum - A{j}*S0; end
      if j > i, sum = sum - A{j}*S{j-i}'; end
    end
    maxres = max(maxres, norm(sum));
  end
  ascertain(maxres < 1e-12);
end

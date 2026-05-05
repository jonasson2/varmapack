function rho = ref_varma_specrad(A)
  r = size(A, 1);
  p = size(A, 2)/r;
  if p == 0
    rho = 0;
    return
  end
  n = r*p;
  if n <= 50
    Ac = zeros(n, n);
    Ac(1:r, :) = A;
    Ac(r+1:n, 1:n-r) = eye(n-r);
    rho = max(abs(eig(Ac)));
  else
    rho = abs(eigs(@(x) companion_multiply(A, r, n, x), n, 1, "largestabs"));
  end
end

function y = companion_multiply(A, r, n, x)
  y = zeros(n, 1);
  y(1:r) = A*x;
  y(r+1:n) = x(1:n-r);
end

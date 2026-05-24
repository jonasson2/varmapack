function rho = ref_companion_specrad(P, sign)
  if isempty(P)
    rho = 0;
    return
  end
  r = size(P, 1);
  k = size(P, 2)/r;
  if k == 0
    rho = 0;
    return
  end
  n = r*k;
  if n <= 50
    Ac = zeros(n, n);
    Ac(1:r, :) = sign*P;
    Ac(r+1:n, 1:n-r) = eye(n-r);
    rho = max(abs(eig(Ac)));
  else
    rho = abs(eigs(@(x) companion_multiply(P, sign, r, n, x), ...
                   n, 1, "largestabs"));
  end
end

function y = companion_multiply(P, sign, r, n, x)
  y = zeros(n, 1);
  y(1:r) = sign*P*x;
  y(r+1:n) = x(1:n-r);
end

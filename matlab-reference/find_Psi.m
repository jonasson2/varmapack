function Psi = find_Psi(A, B)
  % Get dimensions and prepare
  [r, p, q] = getdims(A, B);
  h = max(p, q);
  if h == 0, Psi = []; return, end
  Aflp = flipmat(A);
  Psi = zeros(h*r, h*r);

  % Compute first block colu
  I = 1:r;
  J = 1:r;
  Psi(I, J) = eye(r);
  for i = 1:h-1
    l = i*r;
    k = p*r - i*r + 1;
    I = I + r;
    if i <= p, Psi(I, J) = Aflp(:, k:end)*Psi(1:l, J); end
    if i <= q, Psi(I, J) = Psi(I, J) + B(:, I - r); end
  end

  % Copy to remaining lower Psi blocks
  J = 1:r;
  k = h*r;
  l = 1;
  for j = 2:h
    J = J + r;
    k = k - r;
    l = l + r;
    Psi(l:end, J) = Psi(1:k, 1:r);
  end
end
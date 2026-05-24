function Psi = ref_varma_psi(A, B, maxlag)
%REF_VARMA_PSI  VARMA impulse-response coefficient matrices.
%
%   Psi = REF_VARMA_PSI(A,B,MAXLAG) returns Psi(:,:,j+1) for j=0,...,MAXLAG,
%   where x_t = sum_j Psi(:,:,j+1)*eps_{t-j}.

  if size(A, 1) > 0
    r = size(A, 1);
  elseif size(B, 1) > 0
    r = size(B, 1);
  else
    r = 1;
  end
  p = size(A, 2)/r;
  q = size(B, 2)/r;
  if maxlag < 0, error('maxlag must be nonnegative'); end
  if isempty(A), A = zeros(r, 0); end
  if isempty(B), B = zeros(r, 0); end
  Psi = zeros(r, r, maxlag + 1);
  Psi(:, :, 1) = eye(r);
  for j = 1:maxlag
    if j <= q
      Psi(:, :, j+1) = B(:, (j-1)*r+1:j*r);
    end
    for i = 1:min(p, j)
      Ai = A(:, (i-1)*r+1:i*r);
      Psi(:, :, j+1) = Psi(:, :, j+1) + Ai*Psi(:, :, j-i+1);
    end
  end
end

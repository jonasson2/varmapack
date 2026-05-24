function Theta = ref_varma_irf(A, B, Sig, maxlag)
%REF_VARMA_IRF  Orthogonalized VARMA impulse responses.
%
%   Theta = REF_VARMA_IRF(A,B,Sig,MAXLAG) returns Theta_j = Psi_j*L, where
%   L*L' = Sig.

  r = size(Sig, 1);
  if isempty(A), A = zeros(r, 0); end
  if isempty(B), B = zeros(r, 0); end
  Psi = ref_varma_psi(A, B, maxlag);
  L = psdfactor(Sig);
  Theta = zeros(size(Psi));
  for j = 1:maxlag+1
    Theta(:, :, j) = Psi(:, :, j)*L;
  end
end

function L = psdfactor(Sig)
  [L, flag] = chol(Sig, 'lower');
  if flag == 0, return, end
  Sig = (Sig + Sig')/2;
  [U, D] = eig(Sig);
  d = diag(D);
  if min(d) < -1e-12
    error('Sig must be positive semidefinite')
  end
  d = max(d, 0);
  L = U*diag(sqrt(d));
end

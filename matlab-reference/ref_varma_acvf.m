function Gamma = ref_varma_acvf(A, B, Sig, maxlag)
%REF_VARMA_ACVF  Theoretical autocovariance function of a VARMA model.
%
%   Gamma = REF_VARMA_ACVF(A,B,Sig,MAXLAG) returns Gamma(:,:,k+1) =
%   Cov(x_t, x_{t-k}) for k = 0,...,MAXLAG.

  r = size(Sig, 1);
  if isempty(A), A = zeros(r, 0); end
  if isempty(B), B = zeros(r, 0); end
  [p, q] = get_dimensions(A, B, Sig);
  [~, G] = find_CG(A, B, Sig);
  PLU = vyw_factorize(A);
  S = vyw_solve(A, PLU, G);
  if p == 0
    Gcell = makecell(G);
    Gamma = zeros(r, r, maxlag + 1);
    for k = 0:min(q, maxlag)
      Gamma(:, :, k+1) = Gcell{k+1};
    end
  else
    Scol = S_extend(A, makecell(G), S, maxlag + 1);
    Gamma = zeros(r, r, maxlag + 1);
    for k = 0:maxlag
      Gamma(:, :, k+1) = Scol(k*r+1:(k+1)*r, :);
    end
  end
end

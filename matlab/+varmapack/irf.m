function Theta = irf(A, B, Sig, maxlag)
%VARMAPACK.IRF  Orthogonalized VARMA impulse responses.
%
%   Theta = VARMAPACK.IRF(A,B,Sig,MAXLAG) returns Theta(:,:,j+1) =
%   Psi(:,:,j+1)*L for j=0,...,MAXLAG, where L*L' = Sig.
%
%   A = [A1 ... Ap] and B = [B1 ... Bq] contain the r-by-r autoregressive
%   and moving-average coefficient matrices. Use [] for A or B when p=0 or
%   q=0, respectively. Sig is the r-by-r shock covariance matrix.
%
%   See also VARMAPACK.PSI.
  r = size(Sig, 1);
  if isempty(A), A = zeros(r, 0); end
  if isempty(B), B = zeros(r, 0); end
  Theta = varmapack_irf_gateway(A, B, Sig, maxlag);
end

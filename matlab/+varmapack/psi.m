function Psi = psi(A, B, maxlag)
%VARMAPACK.PSI  VARMA impulse-response coefficient matrices.
%
%   Psi = VARMAPACK.PSI(A,B,MAXLAG) returns Psi(:,:,j+1) for j=0,...,MAXLAG,
%   where x_t = sum_j Psi(:,:,j+1)*eps_{t-j}.
%
%   A = [A1 ... Ap] and B = [B1 ... Bq] contain the r-by-r autoregressive
%   and moving-average coefficient matrices. Use [] for A or B when p=0 or
%   q=0, respectively.
%
%   See also VARMAPACK.IRF.
  if ~isempty(A)
    r = size(A, 1);
  elseif ~isempty(B)
    r = size(B, 1);
  else
    r = 1;
  end
  if isempty(A), A = zeros(r, 0); end
  if isempty(B), B = zeros(r, 0); end
  Psi = varmapack_psi_gateway(A, B, maxlag);
end

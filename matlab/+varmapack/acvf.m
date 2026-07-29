function Gamma = acvf(A, B, Sig, maxlag)
%VARMAPACK.ACVF  Theoretical autocovariance function of a VARMA model.
%
%   Gamma = VARMAPACK.ACVF(A,B,Sig,MAXLAG) returns Gamma(:,:,k+1) =
%   Cov(x_t, x_{t-k}) for k = 0,...,MAXLAG.
%
%   A = [A1 ... Ap] and B = [B1 ... Bq] contain the r-by-r autoregressive
%   and moving-average coefficient matrices. Use [] for A or B when p=0 or
%   q=0, respectively. Sig is the r-by-r shock covariance matrix.
%
%   See also VARMAPACK.AUTOCOV, VARMAPACK.SPECRAD.
  r = size(Sig, 1);
  if isempty(A), A = zeros(r, 0); end
  if isempty(B), B = zeros(r, 0); end
  Gamma = varmapack_acvf_gateway(A, B, Sig, maxlag);
end

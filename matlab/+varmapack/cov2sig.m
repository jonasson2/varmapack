function Sig = cov2sig(A, B, cov)
%VARMAPACK.COV2SIG  Recover innovation covariance from model covariance.
%
%   Sig = VARMAPACK.COV2SIG(A,B,cov) recovers the r-by-r innovation
%   covariance for a specified r-by-r lag-zero model covariance.
%   A = [A1 ... Ap] and B = [B1 ... Bq]; use [] for absent terms.
%   Only the lower triangle of cov is used. Sig is returned symmetric.
%
%   The AR model must be stationary. Both covariances must be positive
%   semidefinite; singular covariances are allowed. Incompatible targets
%   and numerically unreliable inverse problems raise errors.
%
%   See also VARMAPACK.ACVF, VARMAPACK.COV2CORR.
  r = size(cov, 1);
  if isempty(A), A = zeros(r, 0); end
  if isempty(B), B = zeros(r, 0); end
  Sig = varmapack_cov2sig_gateway(A, B, cov);
end

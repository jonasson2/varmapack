function Gamma = varmapack_acvf(A, B, Sig, maxlag)
%VARMAPACK_ACVF  Theoretical autocovariance function of a VARMA model.
  r = size(Sig, 1);
  if isempty(A), A = zeros(r, 0); end
  if isempty(B), B = zeros(r, 0); end
  p = size(A, 2)/r;
  maxlagC = max(maxlag, p);
  Gamma = varmapack_acvf_gateway(A, B, Sig, maxlagC);
  Gamma = Gamma(:, :, 1:maxlag+1);
end

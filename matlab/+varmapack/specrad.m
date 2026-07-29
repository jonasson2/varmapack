function rho = specrad(A)
%VARMAPACK.SPECRAD  Spectral radius of a VAR model companion matrix.
%
%   rho = VARMAPACK.SPECRAD(A) returns the spectral radius of the companion
%   matrix for the autoregressive coefficient matrices in A.
%
%   A = [A1 A2...Ap] is r-by-(r*p), formed by horizontally concatenating the
%   r-by-r coefficient matrices A1,...,Ap. An empty A represents p = 0.
%
%   rho is the largest absolute eigenvalue of the model companion matrix. A VAR
%   model is stationary exactly when rho < 1. VARMAPACK.SPECRAD([]) returns 0.
%
%   See also VARMAPACK.MA_SPECRAD, VARMAPACK.ACVF, VARMAPACK.TESTCASE.
  if isempty(A)
    rho = 0;
  else
    rho = varmapack_specrad_gateway(A);
  end
end

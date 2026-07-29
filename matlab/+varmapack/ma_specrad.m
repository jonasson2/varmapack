function rho = ma_specrad(B)
%VARMAPACK.MA_SPECRAD  Spectral radius of an MA companion matrix.
%
%   rho = VARMAPACK.MA_SPECRAD(B) returns the spectral radius of the companion
%   matrix for the moving-average coefficient matrices in B.
%
%   B = [B1 B2...Bq] is r-by-(r*q), formed by horizontally concatenating the
%   r-by-r coefficient matrices B1,...,Bq. An empty B represents q = 0.
%
%   rho is the largest absolute eigenvalue of the MA companion matrix, whose
%   top block row is [-B1 -B2 ... -Bq]. A VARMA model is invertible exactly
%   when rho < 1. VARMAPACK.MA_SPECRAD([]) returns 0.
%
%   See also VARMAPACK.SPECRAD, VARMAPACK.ACVF, VARMAPACK.TESTCASE.
  if isempty(B)
    rho = 0;
    return
  end
  rho = varmapack_ma_specrad_gateway(B);
end

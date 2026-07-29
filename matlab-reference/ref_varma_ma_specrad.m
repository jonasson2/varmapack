function rho = ref_varma_ma_specrad(B)
%REF_VARMA_MA_SPECRAD  Spectral radius of an MA companion matrix.
%
%   rho = REF_VARMA_MA_SPECRAD(B) returns the spectral radius of the companion
%   matrix for the moving-average coefficient matrices in B.
%
%   B = [B1 B2...Bq] is r-by-(r*q), formed by horizontally concatenating the
%   r-by-r coefficient matrices B1,...,Bq. An empty B represents q = 0.
%
%   rho is the largest absolute eigenvalue of the MA companion matrix, whose
%   top block row is [-B1 -B2 ... -Bq]. A VARMA model is invertible exactly
%   when rho < 1. REF_VARMA_MA_SPECRAD([]) returns 0.
%
%   See also REF_VARMA_SPECRAD, REF_VARMA_ACVF, REF_VARMA_TESTCASE.
  rho = ref_companion_specrad(B, -1);
end

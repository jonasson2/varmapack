function rho = ref_varma_specrad(A)
%REF_VARMA_SPECRAD  Spectral radius of a VAR model companion matrix.
%
%   rho = REF_VARMA_SPECRAD(A) returns the spectral radius of the companion
%   matrix for the autoregressive coefficient matrices in A.
%
%   A = [A1 A2...Ap] is r-by-(r*p), formed by horizontally concatenating the
%   r-by-r coefficient matrices A1,...,Ap. An empty A represents p = 0.
%
%   rho is the largest absolute eigenvalue of the model companion matrix. A VAR
%   model is stationary exactly when rho < 1. REF_VARMA_SPECRAD([]) returns 0.
%
%   See also REF_VARMA_MA_SPECRAD, REF_VARMA_ACVF, REF_VARMA_TESTCASE.
  rho = ref_companion_specrad(A, 1);
end

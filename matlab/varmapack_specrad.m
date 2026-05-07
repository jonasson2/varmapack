function rho = varmapack_specrad(A)
%VARMAPACK_SPECRAD  Spectral radius of a VAR companion matrix.
  if isempty(A)
    rho = 0;
    return
  end
  rho = varmapack_specrad_gateway(A);
end

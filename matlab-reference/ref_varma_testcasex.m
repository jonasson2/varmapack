function [C, z] = ref_varma_testcasex(s, r, n)
%REF_VARMA_TESTCASEX  Create deterministic exogenous data for a VARMAX testcase.
%
%   [C,z] = REF_VARMA_TESTCASEX(s,r,n) returns r-by-s coefficient vectors C
%   and an n-by-1 exogenous sequence z. s is the number of exogenous terms.
%
%   The returned values are deterministic and suitable for testing
%   REF_VARMA_SIMX. C is empty when s is zero.

  if s < 0 || r <= 0 || n < 0
    error('s and n must be nonnegative, and r must be positive')
  end
  C = zeros(r, s);
  for j = 1:s
    for i = 1:r
      C(i, j) = (mod(3*i + 5*j, 11) - 5)/10;
    end
  end
  pattern = [3 -1 4 0 -5 2 1 -3 5 -2 0 4 -4 1 -1 2 -3]/5;
  z = pattern(mod(0:n-1, length(pattern)) + 1)';
end

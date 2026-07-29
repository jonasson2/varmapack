function [C, z] = ref_varma_testcasex(s, d, r, n)
%REF_VARMA_TESTCASEX  Create deterministic exogenous data for a VARMAX testcase.
%
%   [C,z] = REF_VARMA_TESTCASEX(s,d,r,n) returns r-by-(d*s) coefficient matrix C
%   and a d-by-n exogenous sequence z. s is the number of exogenous terms and
%   d is the dimension of each exogenous value.
%
%   The returned values are deterministic and suitable for testing
%   REF_VARMA_SIMX. C is empty when s is zero.

  if s < 0 || d < 0 || r <= 0 || n < 0 || (s == 0 && d ~= 0) || ...
      (s > 0 && d == 0)
    error('s and d must be nonnegative, r must be positive, and s and d must agree')
  end
  C = zeros(r, d*s);
  for k = 1:s
    for j = 1:d
      for i = 1:r
        C(i, j + (k - 1)*d) = (mod(3*i + 5*k + 7*(j - 1), 11) - 5)/10;
      end
    end
  end
  pattern = [3 -1 4 0 -5 2 1 -3 5 -2 0 4 -4 1 -1 2 -3]/5;
  z = zeros(d, n);
  for j = 1:n
    for i = 1:d
      z(i, j) = pattern(mod(i + j - 2, length(pattern)) + 1);
    end
  end
end

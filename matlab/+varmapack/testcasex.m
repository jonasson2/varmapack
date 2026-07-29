function [C, z] = testcasex(s, r, n)
%VARMAPACK.TESTCASEX  Create deterministic exogenous data for a VARMAX testcase.
%
%   [C,z] = VARMAPACK.TESTCASEX(s,r,n) returns r-by-s coefficient vectors C
%   and an n-by-1 exogenous sequence z. s is the number of exogenous terms.
%
%   The returned values are deterministic and suitable for testing
%   VARMAPACK.SIMX. C is empty when s is zero.
%
%   See also VARMAPACK.TESTCASE, VARMAPACK.SIMX.
  [C, z] = varmapack_testcasex_gateway(s, r, n);
end

function [C, z] = testcasex(s, d, r, n)
%VARMAPACK.TESTCASEX  Create deterministic exogenous data for a VARMAX testcase.
%
%   [C,z] = VARMAPACK.TESTCASEX(s,d,r,n) returns r-by-(d*s) coefficient matrix C
%   and a d-by-n exogenous sequence z. s is the number of exogenous terms and
%   d is the dimension of each exogenous value.
%
%   The returned values are deterministic and suitable for testing
%   VARMAPACK.SIMX. C is empty when s is zero.
%
%   See also VARMAPACK.TESTCASE, VARMAPACK.SIMX.
  [C, z] = varmapack_testcasex_gateway(s, d, r, n);
end

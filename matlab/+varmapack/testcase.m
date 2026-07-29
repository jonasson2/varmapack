function varargout = testcase(varargin)
%VARMAPACK.TESTCASE  Create testcases for VARMA calculations.
%
%   [A,B,Sig] = VARMAPACK.TESTCASE(NAME) returns a named testcase in matrices
%   A = [A1...Ap], B = [B1...Bq], and Sig, suitable for VARMA calculations.
%   VARMAPACK.TESTCASE('summary') prints a list of the named testcases with
%   possible values for NAME.
%
%   [A,B,Sig,p,q,r] = VARMAPACK.TESTCASE(NAME) also returns the dimensions of
%   the case.
%
%   [A,B,Sig,name] = VARMAPACK.TESTCASE(i) returns the one-based i-th named
%   testcase. To also get dimensions, use [A,B,Sig,p,q,r,name] =
%   VARMAPACK.TESTCASE(i).
%
%   n = VARMAPACK.TESTCASE('COUNT') returns the number of named testcases.
%
%   [A,B,Sig] = VARMAPACK.TESTCASE(p,q,r) returns a deterministic unnamed
%   testcase with dimensions p, q, r.
%
%   [A,B,Sig] = VARMAPACK.TESTCASE(p,q,r,rng) returns a random unnamed
%   testcase using the supplied VARMAPACK.RNG object.
%
%   See also VARMAPACK.TESTCASEX.
  if nargin == 4 && isa(varargin{4}, 'varmapack.Rng') && isscalar(varargin{4})
    varargin{4} = varargin{4}.mex_handle();
  end
  if nargin == 1 && (isstring(varargin{1}) || ischar(varargin{1}))
    key = char(varargin{1});
    if strcmpi(key, 'summary')
      print_summary;
      return
    end
    if strcmpi(key, 'count')
      key = 'count';
    end
    varargin{1} = key;
  end
  if nargout == 0
    varmapack_testcase_gateway(varargin{:});
  else
    [varargout{1:nargout}] = varmapack_testcase_gateway(varargin{:});
  end
end

function print_summary
  n = varmapack.testcase('count');
  fprintf('No. Name          p  q  r  spec.rad\n');
  for i = 1:n
    [A, ~, ~, p, q, r, name] = varmapack.testcase(i);
    if p == 0
      rho = 0;
    else
      rho = varmapack.specrad(A);
    end
    fprintf('%-2i  %-12s %2d %2d %2d   %6.4f\n', i, name, p, q, r, rho);
  end
end

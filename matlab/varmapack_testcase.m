function varargout = varmapack_testcase(varargin)
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
  n = varmapack_testcase('count');
  fprintf('No. Name          p  q  r  spec.rad\n');
  for i = 1:n
    [A, ~, ~, p, q, r, name] = varmapack_testcase(i);
    if p == 0
      rho = 0;
    else
      rho = varmapack_specrad(A);
    end
    fprintf('%-2i  %-12s %2d %2d %2d   %6.4f\n', i, name, p, q, r, rho);
  end
end

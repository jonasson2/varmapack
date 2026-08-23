function test_varmapack
  root = fileparts(mfilename('fullpath'));
  files = dir(fullfile(root, 'test_varmapack_*.m'));
  names = sort(erase({files.name}, '.m'));
  for i = 1:numel(names)
    feval(names{i})
  end
  fprintf('All MATLAB tests passed\n');
end

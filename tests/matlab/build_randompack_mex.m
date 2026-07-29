function build_randompack_mex(randompack_root)
  if nargin < 1 || isempty(randompack_root)
    randompack_root = fullfile(getenv("HOME"), "randompack");
  end
  if ~isfile(fullfile(randompack_root, "src", "randompack.h"))
    error("build_randompack_mex: Randompack source not found; supply randompack_root")
  end
  inc = fullfile(randompack_root, "src");
  lib = fullfile(randompack_root, "build", "src");
  files = [
    "randompack_create.c"
    "randompack_free.c"
    "randompack_seed.c"
    "randompack_u01.c"
    "randompack_mvn.c"
  ];
  linkargs = ["-L" + lib, "-lrandompack"];
  if ~ispc
    linkargs = [linkargs, "LDFLAGS=$LDFLAGS -Wl,-rpath," + lib];
  end
  for k = 1:numel(files)
    mex("-I" + inc, linkargs{:}, files(k));
  end
end

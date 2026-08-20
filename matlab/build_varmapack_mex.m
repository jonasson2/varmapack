function build_varmapack_mex(randompack_root)
  root = fileparts(fileparts(mfilename('fullpath')));
  if nargin < 1 || isempty(randompack_root)
    randompack_root = fullfile(getenv("HOME"), "randompack");
  end
  if ~isfile(fullfile(randompack_root, "src", "randompack.h"))
    error("build_varmapack_mex: Randompack source not found; supply randompack_root")
  end
  inc = fullfile(root, "src");
  lib = build_library_dir(root, "Varmapack");
  rpinc = fullfile(randompack_root, "src");
  rplib = build_library_dir(randompack_root, "Randompack");
  gateways = [
    "varmapack_sim_gateway"
    "varmapack_simx_gateway"
    "varmapack_specrad_gateway"
    "varmapack_ma_specrad_gateway"
    "varmapack_psi_gateway"
    "varmapack_irf_gateway"
    "varmapack_acvf_gateway"
    "varmapack_autocov_gateway"
    "varmapack_cov2corr_gateway"
    "varmapack_testcase_gateway"
    "varmapack_testcasex_gateway"
  ];
  linkargs = ["-L" + lib, "-lvarmapack", "-L" + rplib, "-lrandompack"];
  if ~ispc
    linkargs = [linkargs, "LDFLAGS=$LDFLAGS -Wl,-rpath," + lib + ...
                " -Wl,-rpath," + rplib];
  end
  for k = 1:numel(gateways)
    mex("-I" + inc, "-I" + rpinc, linkargs{:}, ...
        fullfile(root, "matlab", gateways(k) + ".c"), ...
        "-outdir", fullfile(root, "matlab"));
  end
  rngargs = ["-L" + rplib, "-lrandompack"];
  if ~ispc
    rngargs = [rngargs, "LDFLAGS=$LDFLAGS -Wl,-rpath," + rplib];
  end
  mex("-I" + rpinc, rngargs{:}, fullfile(root, "matlab", "varmapack_rng_gateway.c"), ...
      "-outdir", fullfile(root, "matlab"));
end

function lib = build_library_dir(root, name)
  build = fullfile(root, "build", "src");
  if isfolder(build)
    lib = build;
  else
    error("build_varmapack_mex: " + name + " has not been built")
  end
end

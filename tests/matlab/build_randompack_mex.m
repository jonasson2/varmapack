function build_randompack_mex
  rp = "/Users/jonasson/randompack";
  inc = fullfile(rp, "src");
  lib = fullfile(rp, "release", "src");
  files = [
    "randompack_create.c"
    "randompack_free.c"
    "randompack_seed.c"
    "randompack_u01.c"
    "randompack_normal.c"
    "randompack_mvn.c"
    "randompack_set_state.c"
    "randompack_squares_set_key.c"
    "randompack_philox_set_key.c"
  ];
  for k = 1:numel(files)
    mex("-I" + inc, "-L" + lib, "-lrandompack", ...
        "LDFLAGS=$LDFLAGS -Wl,-rpath," + lib, files(k));
  end
end

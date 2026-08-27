#!/bin/sh
set -eu

usage() {
  cat <<'EOF'
Usage: scripts/compiler-check.sh -c CC -f FC -b BLAS [-s]

Configure, build, and test Varmapack with an explicit compiler pair and BLAS
provider. The named build directory is recreated on each run.

  -c CC    C compiler command, for example clang or gcc
  -f FC    Fortran compiler command, for example gfortran or ifx
  -b BLAS  BLAS provider: accelerate, openblas, or mkl
  -s       Run a debug build with C AddressSanitizer and UndefinedBehaviorSanitizer
  -h       Show this help message

Without -s, the script performs a release build. All builds treat C warnings
as errors; warnings from the bundled SLICOT Fortran source are suppressed.
EOF
}

show_compile_command() {
  awk -v file="$1" '
    /"command":/ { command = $0 }
    index($0, "\"file\": \"" file "\"") {
      sub(/^[[:space:]]*"command": "/, "", command)
      sub(/",$/, "", command)
      print command
      exit
    }
  ' "$builddir/compile_commands.json"
}

cc=
fc=
blas=
sanitizers=0
while getopts 'c:f:b:sh' option; do
  case "$option" in
    c) cc=$OPTARG ;;
    f) fc=$OPTARG ;;
    b) blas=$OPTARG ;;
    s) sanitizers=1 ;;
    h)
      usage
      exit 0
      ;;
    *)
      usage >&2
      exit 2
      ;;
  esac
done
shift $((OPTIND - 1))
if [ "$#" -ne 0 ] || [ -z "$cc" ] || [ -z "$fc" ] || [ -z "$blas" ]; then
  usage >&2
  exit 2
fi

case "$blas" in
  accelerate) blas_name=Accelerate ;;
  openblas) blas_name=OpenBLAS ;;
  mkl) blas_name=MKL ;;
  *)
    echo "Unsupported BLAS provider: $blas" >&2
    usage >&2
    exit 2
    ;;
esac

if ! command -v "$cc" >/dev/null 2>&1; then
  echo "C compiler not found: $cc" >&2
  exit 1
fi
if ! command -v "$fc" >/dev/null 2>&1; then
  echo "Fortran compiler not found: $fc" >&2
  exit 1
fi

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(dirname "$script_dir")
cc_name=$(basename "$cc")
fc_name=$(basename "$fc")
builddir="$repo_root/build-$cc_name-$fc_name-$blas"
buildtype=release
c_std=c11
if [ "$cc_name" = nvc ]; then
  # Meson 1.3 does not advertise NVC's supported C standard selections.
  c_std=none
fi
set -- "-Dbuildtype=$buildtype" "-Dc_std=$c_std" -Dwarning_level=3 -Dwerror=true \
  "-Dblas_order=$blas"
if [ "$cc_name" = clang ] && [ "$(uname -s)" = Linux ]; then
  set -- "$@" -Dc_link_args=-no-pie
fi
if [ "$sanitizers" -eq 1 ]; then
  builddir="${builddir}-san"
  buildtype=debug
  set -- "-Dbuildtype=$buildtype" "-Dc_std=$c_std" -Dwarning_level=3 -Dwerror=true \
    -Denable_c_sanitizers=true "-Dblas_order=$blas" -Db_sanitize=none \
    -Db_lundef=false
fi

export CC=$cc
export FC=$fc

echo "C compiler: $CC"
"$CC" --version 2>&1 | sed -n '1p' || :
echo "Fortran compiler: $FC"
"$FC" --version 2>&1 | sed -n '1p' || :
echo "Build directory: $builddir"
echo "Build type: $buildtype"

cd "$repo_root"
# This is the compiler/Fortran/BLAS-specific build directory derived above.
rm -rf "$builddir"
meson setup "$builddir" "$@"
if ! grep -F "Using BLAS library: $blas_name" \
  "$builddir/meson-logs/meson-log.txt" >/dev/null; then
  echo "Meson did not select requested BLAS provider: $blas_name" >&2
  exit 1
fi
echo "C compile command:"
show_compile_command ../src/FindW.c
echo "Fortran compile command:"
show_compile_command ../src/sb03md-complete.F
meson compile -C "$builddir"
meson test -C "$builddir" --no-rebuild --print-errorlogs

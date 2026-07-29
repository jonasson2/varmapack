#!/bin/sh
set -e

root_dir=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
t=0.1
w=0.1
n=100
M=1000

usage()
{
  cat <<EOF
Usage: TimeSimulate.sh [options]

Options:
  -h          show this help
  -t seconds  timing target per case (default: $t)
  -w seconds  CPU warmup time (default: $w)
  -n length   length of each series (default: $n)
  -M count    replicates per call (default: $M)
EOF
}

while getopts "ht:w:n:M:" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    t) t=$OPTARG ;;
    w) w=$OPTARG ;;
    n) n=$OPTARG ;;
    M) M=$OPTARG ;;
    *)
      usage >&2
      exit 2
      ;;
  esac
done

shift $((OPTIND - 1))
if [ "$#" -ne 0 ]; then
  usage >&2
  exit 2
fi

cmd="addpath('$root_dir/matlab/examples'); addpath('$root_dir/matlab'); "
cmd="$cmd addpath('$root_dir/matlab-reference'); addpath('$root_dir/tests/matlab'); "
cmd="$cmd TimeSimulate('t', $t, 'w', $w, 'n', $n, 'M', $M)"

matlab -batch "$cmd"

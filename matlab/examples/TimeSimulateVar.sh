#!/bin/sh
set -e

root_dir=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
t=0.1
w=0.1
n=100
M=1000
burnin=200
selected=false
platform=false

usage()
{
  cat <<EOF
Usage: TimeSimulateVar.sh [options]

Options:
  -h          show this help
  -t seconds  timing target per case (default: $t)
  -w seconds  CPU warmup time (default: $w)
  -n length   retained length of each series (default: $n)
  -M count    replicates per call (default: $M)
  -b count    discarded burn-in values per path (default: $burnin)
  -5          time the five selected platform models only
  -8          time the selected compatible paper models only
EOF
}

while getopts "58ht:w:n:M:b:" opt; do
  case "$opt" in
    5) platform=true ;;
    8) selected=true ;;
    h)
      usage
      exit 0
      ;;
    t) t=$OPTARG ;;
    w) w=$OPTARG ;;
    n) n=$OPTARG ;;
    M) M=$OPTARG ;;
    b) burnin=$OPTARG ;;
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
cmd="$cmd TimeSimulateVar('t', $t, 'w', $w, 'n', $n, 'M', $M, 'burnin', $burnin, "
cmd="$cmd 'selected', $selected, 'platform', $platform)"

matlab -batch "$cmd"

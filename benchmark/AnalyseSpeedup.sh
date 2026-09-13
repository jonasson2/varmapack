#!/bin/sh
set -e

root_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
runs=5
target=0.1
warmup=0.1
length=100
paths=1000
burnin=200

usage()
{
  cat <<EOF
Usage: AnalyseSpeedup.sh [options]

Options:
  -h          show this help
  -R count    independent timing runs (default: $runs)
  -t seconds  target per individual timing (default: $target)
  -w seconds  CPU warmup time (default: $warmup)
  -n length   length of each series (default: $length)
  -M count    paths in a batch (default: $paths)
  -b count    discarded varm burn-in values (default: $burnin)
EOF
}

while getopts "hR:t:w:n:M:b:" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    R) runs=$OPTARG ;;
    t) target=$OPTARG ;;
    w) warmup=$OPTARG ;;
    n) length=$OPTARG ;;
    M) paths=$OPTARG ;;
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

command="addpath('$root_dir/benchmark'); addpath('$root_dir/matlab/examples'); "
command="$command addpath('$root_dir/matlab'); AnalyseSpeedup('R', $runs, "
command="$command 't', $target, 'w', $warmup, 'n', $length, 'M', $paths, "
command="$command 'burnin', $burnin)"
matlab -batch "$command"

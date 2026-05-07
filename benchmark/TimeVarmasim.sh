#!/bin/sh
set -e

root_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
t=0.2
w=0.1
d=2
p_list=1,3,5
q_list=1,3,5
r_list=2,5,16
n_list=0,100
M_list=1,100
rho_list=.5,.95,.995

usage()
{
  cat <<EOF
Usage: TimeVarmasim.sh [options]

Options:
  -h          show this help
  -t seconds  timing target per case (default: $t)
  -w seconds  warmup target per case (default: $w)
  -d digits   printed digits (default: $d)
  -p list     p values (default: $p_list)
  -q list     q values (default: $q_list)
  -r list     r values (default: $r_list)
  -n list     n values, with 0 meaning max(p,q) (default: $n_list)
  -M list     M values (default: $M_list)
  -R list     rho values (default: $rho_list)
EOF
}

while getopts "ht:w:d:p:q:r:n:M:R:" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    t) t=$OPTARG ;;
    w) w=$OPTARG ;;
    d) d=$OPTARG ;;
    p) p_list=$OPTARG ;;
    q) q_list=$OPTARG ;;
    r) r_list=$OPTARG ;;
    n) n_list=$OPTARG ;;
    M) M_list=$OPTARG ;;
    R) rho_list=$OPTARG ;;
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

cmd="addpath('$root_dir/benchmark'); addpath('$root_dir/matlab'); "
cmd="$cmd addpath('$root_dir/matlab-reference'); addpath('$root_dir/tests/matlab'); "
cmd="$cmd TimeVarmaSim("
cmd="$cmd't', $t, 'w', $w, 'd', $d, "
cmd="$cmd'p', [$p_list], 'q', [$q_list], 'r', [$r_list], "
cmd="$cmd'n', [$n_list], 'M', [$M_list], 'rho', [$rho_list])"

matlab -batch "$cmd"

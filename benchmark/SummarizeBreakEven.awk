# Summarize TimeBreakEven output.
#
# Usage:
#   awk -f benchmark/SummarizeBreakEven.awk time_breakeven.out
#   release/benchmark/TimeBreakEven -t 0.05 | awk -f benchmark/SummarizeBreakEven.awk

NR > 4 && $1 ~ /^[0-9]+$/ {
  r = $1
  p = $2
  q = $3
  ratio[p "," q "," r] = $8
  if (!(p in pseen)) {
    pseen[p] = 1
    ps[np++] = p
  }
  if (!(q in qseen)) {
    qseen[q] = 1
    qs[nq++] = q
  }
  if (!(p "," q in minr) || r < minr[p "," q]) minr[p "," q] = r
  if (!(p "," q in maxr) || r > maxr[p "," q]) maxr[p "," q] = r
}

END {
  sort_ints(ps, np)
  sort_ints(qs, nq)
  print_table("first winning r", "first")
  print ""
  print_table("stable break-even r", "stable")
  print ""
  print_table("combined rc", "combined")
  print ""
  print "reported gaps r1-r0 >= 5"
  any = 0
  for (ip=0; ip<np; ip++) {
    p = ps[ip]
    for (iq=0; iq<nq; iq++) {
      q = qs[iq]
      r0 = first_win(p, q)
      r1 = stable_win(p, q)
      if (r0 != "-" && r1 != "-" && r1 - r0 >= 5) {
        printf "p=%d q=%d r0=%d r1=%d gap=%d\n", p, q, r0, r1, r1 - r0
        any = 1
      }
    }
  }
  if (!any) print "none"
}

function print_table(title, kind, ip, iq, p, q, v) {
  print title
  printf "     "
  for (iq=0; iq<nq; iq++) printf " q=%-3d", qs[iq]
  print ""
  for (ip=0; ip<np; ip++) {
    p = ps[ip]
    printf "p=%-2d ", p
    for (iq=0; iq<nq; iq++) {
      q = qs[iq]
      if (kind == "first") v = first_win(p, q)
      else if (kind == "stable") v = stable_win(p, q)
      else v = combined_rc(p, q)
      printf "%5s ", v
    }
    print ""
  }
}

function first_win(p, q, r, key) {
  key = p "," q
  if (!(key in minr)) return "-"
  for (r=minr[key]; r<=maxr[key]; r++) {
    if (ratio[p "," q "," r] < 1) return r
  }
  return "-"
}

function stable_win(p, q, r, s, ok, key) {
  key = p "," q
  if (!(key in minr)) return "-"
  for (r=minr[key]; r<=maxr[key]; r++) {
    ok = 1
    for (s=r; s<=maxr[key]; s++) {
      if (ratio[p "," q "," s] >= 1) ok = 0
    }
    if (ok) return r
  }
  return "-"
}

function combined_rc(p, q, r0, r1, d) {
  r0 = first_win(p, q)
  r1 = stable_win(p, q)
  if (r0 == "-" || r1 == "-") return "-"
  d = r1 - r0
  if (d == 0) return r0
  if (d%2 == 0) return int((r0 + r1)/2)
  if (ratio[p "," q "," (r0 + 1)] < 1) return r0 + 1
  return r0 + 2
}

function sort_ints(a, n, i, j, x) {
  for (i=1; i<n; i++) {
    x = a[i]
    j = i
    while (j > 0 && x < a[j-1]) {
      a[j] = a[j-1]
      j--
    }
    a[j] = x
  }
}

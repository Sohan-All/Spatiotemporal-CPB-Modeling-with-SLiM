"""Settle the n_real / n_big disagreement with replicates (CLAUDE.md 7.9.12F).

For each arm and gap, the F_c drift term is averaged over independent replicates at
POPMULT 2000 and 4000 (K=40 fixed), and the elasticity in N is reported with a CI:

    elasticity = ln(mean_4000 / mean_2000) / ln 2

  0.0  = the statistic does not read deme size at all
 -1.0  = it reads deme size exactly (a persistent deme measures -0.996, CLAUDE.md 7.9.4)

Every replicate re-rolls SLiM, recapitation, the mutation overlay and the subsample, so
these sds are a REAL replicate floor, not the fixed-tree subsample spread the single
points carried.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

recs = [json.loads(l) for l in open(sys.argv[1], encoding="utf-8") if l.strip()]

# Guard against pooling different statistics -- the lesson of CLAUDE.md 7.9.8D.
specs = {r["fc_spec"] for r in recs}
assert len(specs) == 1, f"records disagree on fc_spec: {specs} -- refusing to pool"
for k in ("recomb", "anc_ne", "mu", "total_migration", "num_clusters", "deme_choice"):
    vals = {r[k] for r in recs}
    assert len(vals) == 1, f"records disagree on {k}: {vals} -- refusing to pool"
print(f"spec {specs.pop()}, constants agree across {len(recs)} records")


def exc(res, t):
    return res[t]["fc_mean"] - res[t]["pedestal_mean"]


def grab(arm, pop, fn):
    return np.array([fn(r["res"]) for r in recs
                     if r["arm"] == arm and r["popmult"] == pop])


def line(label, a, b):
    ma, mb = a.mean(), b.mean()
    sa, sb = a.std(ddof=1) / math.sqrt(len(a)), b.std(ddof=1) / math.sqrt(len(b))
    el = math.log(mb / ma) / math.log(2.0)
    # delta method on ln(mb/ma)
    se = math.sqrt((sb / mb) ** 2 + (sa / ma) ** 2) / math.log(2.0)
    print(f"  {label:<22} {ma:.5f} +/- {sa:.5f}   {mb:.5f} +/- {sb:.5f}   "
          f"{el:+.3f}  [{el - 1.96 * se:+.3f}, {el + 1.96 * se:+.3f}]")
    return el, se


print("\n" + "=" * 100)
print("REPLICATED N-SENSITIVITY UNDER RE-FOUNDING (K=40 fixed, kernel 5e-5)")
print("=" * 100)
for arm in ("n_big", "n_real"):
    a8, b8 = grab(arm, 2000, lambda r: exc(r, "t8")), grab(arm, 4000, lambda r: exc(r, "t8"))
    if len(a8) == 0 or len(b8) == 0:
        continue
    print(f"\n--- {arm} arm   (n={len(a8)} at POPMULT 2000, n={len(b8)} at 4000) ---")
    print(f"  {'quantity':<22} {'POPMULT 2000':>18}   {'POPMULT 4000':>18}   "
          f"elasticity in N [95% CI]")
    line("F_c drift, t=8", a8, b8)
    line("F_c drift, t=16",
         grab(arm, 2000, lambda r: exc(r, "t16")), grab(arm, 4000, lambda r: exc(r, "t16")))
    line("F_st", grab(arm, 2000, lambda r: r["fst_hudson"]),
         grab(arm, 4000, lambda r: r["fst_hudson"]))

    ga = grab(arm, 2000, lambda r: exc(r, "t16") - exc(r, "t8"))
    gb = grab(arm, 4000, lambda r: exc(r, "t16") - exc(r, "t8"))
    print(f"  {'growth (t16 - t8)':<22} {ga.mean():+.5f} +/- "
          f"{ga.std(ddof=1)/math.sqrt(len(ga)):.5f}   {gb.mean():+.5f} +/- "
          f"{gb.std(ddof=1)/math.sqrt(len(gb)):.5f}     (observed -0.00151)")

print("\n" + "=" * 100)
print("REFERENCE  persistent deme, same statistic: elasticity in N = -0.996 (CLAUDE.md 7.9.4)")
print("           an elasticity whose CI excludes -1 but contains 0 means the statistic has")
print("           stopped reading deme size.")

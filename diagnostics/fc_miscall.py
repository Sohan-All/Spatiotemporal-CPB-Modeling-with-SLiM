"""Does the genotype heterozygote deficit move simulated F_c enough to matter? And does purging
relatives? -- CLAUDE.md 7.2.2F. Seconds per rep on a saved tree.

Every empirical sample carries a per-individual het-miscall rate e (median ~0.24): about a quarter
of true heterozygotes read as homozygous for a random allele. The simulated side has clean calls.
Per individual that is marginally an inbreeding coefficient F = e, which inflates the sampling
variance of a site frequency by ~(1+F) -- naively +0.034 on F_c's pedestal at n=7. The MAF filter is
an ascertainment at n=7 (fc_common.MIN_MAF), so the real size has to be MEASURED, not derived.

Separately, fc_common.EXCLUDED_SAMPLES purges relatives from the empirical side. But a random sample
of a small simulated deme contains relatives too, and those are part of the drift signal
(Waples & Anderson 2017). The symmetric alternative is to draw simulated samples with no shared
parent -- SLiM records pedigree_p1/p2, so relatives are known exactly. This measures what that does.

FOUR CONDITIONS, all on the same draws per rep so the differences are paired:
    clean              : a random draw, genotypes as simulated          (what production does now)
    miscall            : the same draw, with each individual given the e of one of that field's
                         real retained samples (random order); each heterozygous site is, with
                         probability e, replaced by a homozygote for one of its two alleles
    unrelated          : the same draw with individuals swapped out until no two share a parent
    unrelated_miscall  : both -- the candidate production rule
F_c is pooled by generation gap exactly as fc_loss does. Compare each shift with the POPMULT signal
between two trees (--summarize).

Usage (from diagnostics/; the 7.9.4 trees are already recapitated+simplified+mutated):
    python fc_miscall.py --load-ts ../out/ld_probe_p2000.trees --popmult 2000 --reps 16
    python fc_miscall.py --load-ts ../out/ld_probe_p5000.trees --popmult 5000 --reps 16
    python fc_miscall.py --summarize
"""
import argparse
import json
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import tskit

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import temporal_fc as tfc        # noqa: E402  matched_pairs (nearest-deme attribution)
import kinship_correct as kc     # noqa: E402  per-sample miscall rates
fcc = tfc.fcc

OUT = HERE.parent / "out" / "fc_miscall.jsonl"
VERSION = 2
CONDITIONS = ("clean", "miscall", "unrelated", "unrelated_miscall")


def _individuals(ts, deme, tval):
    """{individual id: [node, node]} for one deme at one sampling time."""
    by = {}
    for nd in ts.samples(population=deme, time=tval):
        by.setdefault(ts.node(nd).individual, []).append(int(nd))
    return by


def _freq(G):
    """G: (sites, n, 2) allele indices -> (frequency of allele 1, any allele index > 1)."""
    return (G == 1).sum(axis=(1, 2)) / (2.0 * G.shape[1]), (G > 1).any(axis=(1, 2))


def run(args):
    pairs = tfc.matched_pairs("nearest", verbose=False)
    rates = kc.miscall_rates()
    members = {y: fcc.popfile_members(tfc._pop_path(y)) for y in fcc.YEAR_TIME}

    ts = tskit.load(args.load_ts)
    print(f"fc_common spec {fcc.spec_hash()}; {len(pairs)} field pairs; tree {ts.num_samples} "
          f"samples, {ts.num_sites} sites", flush=True)

    fields = {}                                   # (year, site) -> (deme, n, real e values)
    for p in pairs:
        for side in ("a", "b"):
            key = (p[f"year_{side}"], p[f"site_{side}"])
            ev = np.clip([rates[s][2] for s in members[key[0]][key[1]]], 0.0, 1.0)
            if len(ev) != p[f"n_{side}"]:
                raise RuntimeError(f"{key}: {len(ev)} retained samples but n={p[f'n_{side}']}")
            fields[key] = (p["cluster"], p[f"n_{side}"], ev)
    pools = {(d, fcc.YEAR_TIME[y]): _individuals(ts, d, fcc.YEAR_TIME[y])
             for (y, _), (d, _, _) in fields.items()}

    parents = {}

    def _par(i):
        if i not in parents:
            md = ts.individual(i).metadata
            parents[i] = {q for q in (md["pedigree_p1"], md["pedigree_p2"]) if q >= 0}
        return parents[i]

    def _related_pair(chosen):
        for a, b in combinations(range(len(chosen)), 2):
            pa, pb = _par(chosen[a]), _par(chosen[b])
            if pa & pb:
                return a, b, pa == pb
        return None

    def _unrelated(chosen, pool, rng):
        """Swap later members of shared-parent pairs for fresh individuals until none remain."""
        chosen = list(chosen)
        spare = [int(i) for i in rng.permutation(sorted(pool)) if int(i) not in set(chosen)]
        k = 0
        while (hit := _related_pair(chosen)) is not None:
            chosen[hit[1]] = spare[k]
            k += 1
        return chosen, k

    reps = []
    sib = {"pairs": 0, "shared_parent": 0, "full_sib": 0, "groups": 0, "groups_with": 0, "swaps": 0}
    for rep in range(args.reps):
        t0 = time.time()
        rng = np.random.default_rng(args.seed * 1000 + rep)
        draw, draw_u = {}, {}
        for key, (deme, n, _) in fields.items():
            pool = pools[(deme, fcc.YEAR_TIME[key[0]])]
            draw[key] = [int(i) for i in rng.choice(sorted(pool), size=n, replace=False)]
            sib["groups"] += 1
            hit = False
            for i, j in combinations(draw[key], 2):
                pi, pj = _par(i), _par(j)
                sib["pairs"] += 1
                if pi & pj:
                    sib["shared_parent"] += 1
                    sib["full_sib"] += int(pi == pj)
                    hit = True
            sib["groups_with"] += int(hit)
            draw_u[key], swaps = _unrelated(draw[key], pool, rng)
            sib["swaps"] += swaps

        def _nodes(d):
            return {nd for key in d for i in d[key]
                    for nd in pools[(fields[key][0], fcc.YEAR_TIME[key[0]])][i]}

        allnodes = np.array(sorted(_nodes(draw) | _nodes(draw_u)), dtype=np.int32)
        work, nmap = ts.simplify(samples=allnodes, filter_populations=False, map_nodes=True)
        G = work.genotype_matrix()

        freq = {}
        for key, (deme, n, ev) in fields.items():
            pool = pools[(deme, fcc.YEAR_TIME[key[0]])]
            e_perm = rng.permutation(ev)
            Gr = G[:, np.array([[nmap[nd] for nd in pool[i]] for i in draw[key]])]
            Gu = G[:, np.array([[nmap[nd] for nd in pool[i]] for i in draw_u[key]])]
            freq[key] = {"clean": _freq(Gr), "miscall": _freq(fcc.apply_miscall(Gr, e_perm, rng)),
                         "unrelated": _freq(Gu),
                         "unrelated_miscall": _freq(fcc.apply_miscall(Gu, e_perm, rng))}
        del G

        res = {}
        for cond in CONDITIONS:
            rows = []
            for p in pairs:
                x, mx = freq[(p["year_a"], p["site_a"])][cond]
                y, my = freq[(p["year_b"], p["site_b"])][cond]
                fc, _, _, _ = fcc.waples_fc(x, y, p["n_a"], p["n_b"], p["t"], valid=~mx & ~my)
                rows.append({"t": p["t"], "fc": fc})
            res[cond] = {k: v["fc_mean"] for k, v in fcc.pool_by_gap(rows).items()}
        reps.append(res)
        print(f"  rep {rep:>2}: " + "  ".join(
            f"{g}: clean {res['clean'][g]:.5f} miscall {res['miscall'][g] - res['clean'][g]:+.5f} "
            f"unrel {res['unrelated'][g] - res['clean'][g]:+.5f}" for g in ("t8", "t16"))
            + f"  [{time.time() - t0:.0f} s]", flush=True)

    out = {}
    for g in ("t8", "t16", "all"):
        clean = np.array([r["clean"][g] for r in reps])
        out[g] = {"clean": float(clean.mean()), "clean_sd": float(clean.std(ddof=1))}
        for cond in CONDITIONS[1:]:
            v = np.array([r[cond][g] for r in reps])
            out[g][cond] = float(v.mean())
            out[g][f"{cond}_shift"] = float((v - clean).mean())
            out[g][f"{cond}_shift_se"] = float((v - clean).std(ddof=1) / np.sqrt(len(reps)))
    rec = {"version": VERSION, "popmult": args.popmult, "ts": args.load_ts, "reps": args.reps,
           "seed": args.seed, "fc_spec": fcc.spec_hash(),
           "when": time.strftime("%Y-%m-%dT%H:%M:%S"), "res": out, "sib": sib}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "a", encoding="utf-8") as f:
        f.write(json.dumps(rec) + "\n")
    return rec


def summarize():
    recs = [json.loads(l) for l in open(OUT, encoding="utf-8") if l.strip()]
    n_all = len(recs)
    # Only records under the CURRENT spec are pooled -- never mix sample sets (CLAUDE.md 7.9.8D).
    recs = [r for r in recs if r.get("version") == VERSION and r["fc_spec"] == fcc.spec_hash()]
    if not recs:
        raise SystemExit(f"no version-{VERSION} records under spec {fcc.spec_hash()} in {OUT}")
    recs.sort(key=lambda r: r["popmult"])
    print(f"fc_common spec {fcc.spec_hash()}: {len(recs)} of {n_all} records   "
          f"(shift = condition minus clean, paired; +- is its SE)")
    print(f"{'POPMULT':>8} {'gap':>4} {'clean':>9} " + " ".join(f"{c + ' shift':>24}" for c in CONDITIONS[1:]))
    for r in recs:
        for g in ("t8", "t16", "all"):
            v = r["res"][g]
            print(f"{r['popmult']:>8.0f} {g:>4} {v['clean']:>9.5f} " + " ".join(
                f"{v[c + '_shift']:>+14.5f} +-{v[c + '_shift_se']:.5f}" for c in CONDITIONS[1:]))
        s = r["sib"]
        print(f"{'':>8} random simulated draws: {s['shared_parent']}/{s['pairs']} individual pairs share "
              f"a parent ({s['full_sib']} full sibs); {s['groups_with']}/{s['groups']} field samples "
              f"contain one; {s['swaps']} swaps to purge")
    if len(recs) >= 2:
        lo, hi = recs[0], recs[-1]
        print(f"\nPOPMULT {lo['popmult']:.0f} -> {hi['popmult']:.0f}: F_c signal under each rule, "
              f"and each shift as a fraction of the CLEAN signal")
        for g in ("t8", "t16", "all"):
            sig = {c: (lo["res"][g]["clean"] if c == "clean" else lo["res"][g][c])
                   - (hi["res"][g]["clean"] if c == "clean" else hi["res"][g][c]) for c in CONDITIONS}
            print(f"  {g:>4}: " + "  ".join(f"{c} {sig[c]:+.5f}" for c in CONDITIONS))
            for c in CONDITIONS[1:]:
                shift = 0.5 * (lo["res"][g][c + "_shift"] + hi["res"][g][c + "_shift"])
                print(f"        mean {c} shift {shift:+.5f} = {shift / sig['clean']:+.2f} x clean signal")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--load-ts")
    ap.add_argument("--popmult", type=float)
    ap.add_argument("--reps", type=int, default=16)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--summarize", action="store_true")
    args = ap.parse_args()
    if args.summarize:
        return summarize()
    if not (args.load_ts and args.popmult):
        raise SystemExit("need --load-ts and --popmult")
    run(args)


if __name__ == "__main__":
    main()

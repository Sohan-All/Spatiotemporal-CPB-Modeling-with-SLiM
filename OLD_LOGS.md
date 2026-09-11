# OLD_LOGS.md — archived detail from CLAUDE.md

**Do NOT read this file unless a session explicitly asks for it, or `CLAUDE.md` sends you here by
name.** It is not project context. It is the long-form record behind conclusions that `CLAUDE.md`
now states in a paragraph — kept because the arguments, measurements and caveats are real and
re-deriving them would cost days.

Archived 2026-09-10. Everything below is verbatim as it stood in `CLAUDE.md`; section numbers,
cross-references and dates are as-written and are **not** maintained. Where a claim here disagrees
with `CLAUDE.md`, `CLAUDE.md` wins.

Contents:

- **§A — the LD thread** (old §7.5, §7.6, §7.7, §7.8, plus the two LD reference blocks from §11).
  Retired as a fitted statistic; `CLAUDE.md` §7.5 states the reason.
- **§B — resolved defects** (old §6.3, §6.4, §6.5, §6.6, §6.7, §6.7.1). All fixed in code;
  `CLAUDE.md` §6.3–6.7 keeps the standing rules.
- **§C — the 99-deme cap** (old §7.9.5). Lifted; `CLAUDE.md` §7.9.5 keeps the stub.
- **§D — the per-file code table** (old §9).

---

# §A — the LD thread [ARCHIVED 2026-09-10]

`ld_loss` is retired as a fitted statistic. The single-sentence reason is that the empirical target
pools 17 chromosomes whose effective recombination rates span >=10x, and a single-locus single-`r`
simulation cannot be a mixture. See `CLAUDE.md` §7.5. The full account follows.

### 7.5 LD decay — built and measured, then RETIRED [see §7.8 for the verdict]

> **STATUS 2026-09-08 — LD IS RETIRED AS A FITTED STATISTIC. Read §7.8 first.** Everything in
> §7.5 is built, verified and measured, and `ld_loss` is the single most `pop`-informative
> statistic in the project (§7.6.1, unique R² 0.258). It is still **not fitted**, and that is now
> a settled conclusion, not a pending question: §7.6.2 showed it takes over `D` and F_st stops
> contributing; §7.7 showed freeing `r` cannot help; §6.8.1 pinned `r` from a linkage map; and
> **§7.8 measured why none of that rescues it.** The material below remains correct and is worth
> keeping — the empirical run, the cost model, the subsampling requirement — but do not restart
> the fitting effort from it.

**Why.** §7.4.2: F_st ≈ `1/(1+4Nm)` identifies only the *product*, and π cannot help because §6.2
shows ~7/8 of it is set by the fixed ancestral phase. LD gives a **second equation**: it depends on
`N·r`, so pinning N from LD would free F_st to speak about migration.

> **Two corrections to that premise, both measured after the fact.**
> **(1) "Migration does not enter it" is FALSE** — `ld_loss` loads 0.186 on `total_migration` and
> 0.115 on `m` (§7.6.1), together more than its 0.258 on `pop`. LD is a better handle on N than
> F_st, not a clean independent equation.
> **(2) The second equation only exists if `r` is fixed.** LD constrains `4·N·r` exactly as π
> constrains `4·Ne·μ`; with `r` free it is a ridge, not an equation (§7.7). Precedent: **Boitard et al. 2016 (PopSizeABC)** does exactly this — binned LD as
an ABC summary statistic (§11). Time-depth rationale: **Hayes et al. 2003** — long-range LD
reflects *recent* `N_e`, which is precisely the window the 324-generation forward phase controls.

**Statistics already checked and rejected as alternatives:** d_xy (96% noise, R²=0.043 across
batch 1), relatedness (+0.777 correlated with `fst_loss` — same signal, another vote on Nm not a
second equation), IBD slope (null in all three years, §7.1), and the temporal method
(§11 Waples 1989 — 13 resampled fields exist but n=5–8 makes the sampling correction 11–21× the
drift signal).

#### 7.5.1 The empirical run landed — the decay is at ~60 bp [VERIFIED 2026-09-07]

`LD_WORKERS=9 python CalculateLD.py` on the Beagle machine: 17/17 chromosomes, **7.7 h wall**,
spec `4d1d1d92b25b` on both sides. Log in `ldCalcOut.txt` (repo root); targets in
`data/empiricalStats/averaged_ldDecay_{year}.csv`; per-chromosome files kept in `data/ld_per_chr/`
(640 KB — the only route to a between-chromosome jackknife, and to re-pooling without another
7.7 h run, so **do not delete them**).

**The curve.** Pooled over demes, mean r² runs 0.435 → 0.107 (2015), 0.399 → 0.086 (2019),
0.407 → 0.087 (2023) from 1 bp to 1 Mb, essentially flat past ~10 kb.

**1. The plateau is the `1/n` floor — Hill 1981 confirmed directly, and this is the strongest
possible argument for §7.5a.** Fitting `1/(1+Cd) + F` per population, the fitted floor `F` tracks
`1/n_hap` at **r = 0.995 / 0.984 / 0.769** by year. Within 2015 it runs **0.41 at the two n=2
sites down to 0.080 at `H53-2015` (n=19)** — a 5× spread inside one year, from sample size alone.
So: the simulated side **must** subsample to real `n_i` (already implemented), and `ld_loss`
**must** apply `get_keep_mask` — the n≤3 sites are pure floor.

**2. Thinning is unbiased — verified for free.** At the 10 kb `STAGES` boundary the pair count
drops 625× (thin 1 → 25) and mean r² steps only **−2.3 / −2.5 / −2.7%**, which is the decay
continuing rather than a discontinuity. This is the empirical confirmation of the claim
`ld_common.STAGES` rests on.

**3. Half-decay of the excess-over-plateau is ~60 bp** (between the 32–56 and 56–100 bp bins in
all three years). §6.8 is settled by measurement and **neither reading of `ρ_HAN` was right** —
see the note there. `ρ_obs ≈ 0.017 per bp`.

**4. This is the awkward branch, but not the fatal one.** §7.5(c) predicted that a decay in the
first few bins would put LD "off the table in this form", because sub-100 bp resolution needs
pre-thin=1 at **11.6 h/trial**. 60 bp is indeed below thin=25's ~98 bp spacing. But the curve does
not end there: **34% of the dynamic range remains above 100 bp and 13–14% above 562 bp**, backed
by 20–186 billion pairs per bin, i.e. very precisely measured.

**5. The cut: fit bins with `bin_lo ≥ 562`.** Two independent arguments land on the same place,
which is what makes it more than a convenience:

- **Cost.** Measured at **20.6 s/trial** for all three years at `thin=25` (§7.5.2). Even full
  resolution is only 5.4 min, so **cost is no longer the binding argument for this cut** — the
  time-depth argument below is, and it stands on its own.
- **Time depth.** At `r = 2.75e-6` the 324-generation forward window controls distances
  `d ≥ 1/(2·324·r) ≈ 561 bp`. Everything shorter reflects coalescence in the **fixed** ancestral
  phase (Ne=6700) and therefore carries **no POPMULT signal at all** (Hayes et al. 2003, §11).

So the region we cannot afford is also the region that could not have informed POPMULT anyway.
**Do not "improve" the statistic by adding short bins** — that is 10× the cost for a region the
forward simulation does not control.

**6. What LD does and does not buy — state this carefully.** LD half-decay measures **`4·N_e·r`**,
structurally the *same* confound as π's `4·N_e·μ`. It breaks §7.4.2 only because **`m` is absent
from it**, not because it is assumption-free: any N it yields is conditional on `r` exactly as π's
is on μ. And `r` is disputed (§6.8). The favourable part is that the current `r = 2.75e-6` puts
the observed 60 bp at **POPMULT ≈ 15,000** — inside the prior and within ~20% of F_st's 12,469
(§6.7) and batch 1's median. Two independent statistics agreeing to 20% is real corroboration;
it is **not** independent confirmation of `r`, which is what §6.8's note now spells out.

**7. [OPEN] Still unquantified: Beagle imputation inflates empirical LD.** It pushes `1/ρ_obs`
*longer*, so the true decay is at ≤60 bp, i.e. the affordable region is if anything thinner than
measured here. Subsampling does not fix it (unlike the `1/n` bias). Unchanged from the original
plan — treat the LD *level* as suspect and lean on the decay *shape*.

**Next measurement: `diagnostics/ld_probe.py`** — does the ≥562 bp region actually move with
POPMULT? If it does not, LD does not break the N/m confound and steps 4–6 below should not be
built.

#### 7.5.2 The per-trial cost — the bottleneck was extraction, not pair enumeration [VERIFIED 2026-09-07]

**§7.5c had the cost model wrong, and it was wrong in a way that would have killed the statistic.**
It timed `ld_common.decay_one_pop` in isolation and concluded `thin` was the cost knob. In the
real call it is not even close. Measured per deme at POPMULT=5000 (`diagnostics/ld_probe.py`,
POPMULT=5000 tree, 2015):

| step | thin=1 | thin=25 |
|---|---|---|
| `ts.genotype_matrix(samples=nodes)` | **51.3 s** | 51.3 s (unchanged) |
| `decay_one_pop` (pair enumeration) | 1.2 s | 0.02 s |
| **per-deme total** | **52.5 s** | **51.3 s — a 2% saving** |

**Extraction cost tracks SITES AND TREES, not sample count.** The deme benchmarked has *four*
haplotypes and still took 51 s. So thinning inside the deme loop saves nothing, and the old code
paid the full extraction **once per deme — 61 times (24+17+20) over the same tree sequence.**

**The fix: one `genotype_matrix` call per year, then slice columns per deme.** Measured, 3 years:

| approach | 3 years | overhead on a ~1.2 h trial |
|---|---|---|
| per-deme extraction, thin=1 *(the old code)* | **53 min** | 74% |
| per-deme extraction, thin=25 | ~52 min | 72% |
| **batched extraction, thin=1** | **5.4 min** | 7.5% |
| **batched extraction, thin=25** | **20.6 s** | **0.5%** |

**This does not contradict §7.5b.** "Never call `genotype_matrix` on the full tree" is about
255,690 sites × 70,078 samples = 1.8e10. The batched call uses only the ~334 **subsampled** nodes
of one year (`Σ 2·n_i`): 342 MB at thin=1, 14 MB at thin=25. The site thin also moved out of the
loop into a single `ts.delete_sites()`.

**Verified, not assumed.** The batched path is **byte-identical** to the old per-deme path at
thin=1 on both years the old code had finished writing (2015, 2019) — same seed, same rng
consumption order, `diff` clean. The per-deme biallelic filter stays per-deme on the sliced
columns: a site can be biallelic in one deme and not another, so folding it into the year-wide
matrix would have quietly changed the statistic.

**And thinning is unbiased on the simulated side too:** `ld_loss` is **0.01117 at thin=1 vs
0.01126 at thin=25**, a 0.8% difference over the 13 fitted bins. That matches the empirical-side
check (§7.5.1 pt 2) and is what licenses running production thinned.

> **Consequence for the §7.5.1 cut.** Cost is no longer the binding argument for `bin_lo ≥ 562` —
> at 5.4 min even full resolution is affordable. **The time-depth argument is now doing all the
> work**, and it is the better argument anyway: bins below ~561 bp are set by the fixed ancestral
> phase and cannot carry POPMULT signal at any price. Keep the cut; drop "we can't afford it" as
> the reason.
>
> **Generalisable lesson:** a micro-benchmark of the inner loop is not a cost model of the call.
> The 2% that `thin` actually bought was hiding behind a step nobody timed.

#### 7.5.3 The gate PASSES — LD discriminates POPMULT, but it disagrees with F_st [VERIFIED 2026-09-07]

Two points, `diagnostics/ld_probe.py`, thin=25, seed 1, `bin_lo ≥ 562`, mask applied:

| POPMULT | `ld_loss` | 2015 | 2019 | 2023 | sim halfway bin (2015/19/23) |
|---|---|---|---|---|---|
| 2000 | **0.00631** | 0.00575 | 0.00586 | 0.00734 | 56–100 / 562–1000 / 56–100 |
| 5000 | 0.01126 | 0.00940 | 0.01110 | 0.01327 | 32–56 / 32–56 / 32–56 |
| *observed* | — | — | — | — | 56–100 / 32–56 / 56–100 |

**1. The gate passes.** Movement in the fitted bins between the two POPMULTs, against the sim-obs
gap it would be fitting: **0.98 / 1.15 / 0.82** by year. The fitted region moves by about as much
as the distance it has to close, so **LD does discriminate POPMULT** — unlike π, which §7.2.1
showed is pinned at a flat floor. Steps 4–6 are worth building.

**2. But it points the opposite way from F_st, and that is the headline.** `ld_loss` is
**monotonically decreasing toward the bottom of the prior**: 0.0113 at POPMULT=5000 → 0.0063 at
2000. At POPMULT=5000 the simulated curve already decays *faster* than observed (halfway 32–56 bp
against 56–100), so the fit wants **smaller** N. Against F_st's implied **12,469** (§6.7) and
batch 1's median 12,300–12,600, that is a **≥6× disagreement**, and the LD optimum may well sit
*below* the prior floor — two points cannot locate it, and 2000 is the floor.

**3. This is not necessarily a contradiction — it is plausibly a MEASUREMENT of the §6.8 `r`
error. [INFERRED, and the most useful thing here]** LD constrains `4·N_e·r`; F_st constrains
`4·N_e·m`. If `r` is inflated by a factor k, the LD-preferred N is deflated by exactly k while
F_st's is untouched. So the ratio **F_st's N / LD's N ≈ 6+ is an estimate of k** — an independent,
internal handle on how wrong `DEFAULT_RECOMBINATION_RATE` is, from the project's own data rather
than from Cohen et al.'s arithmetic. §6.8 argued on external grounds that `r` is ~100× too high;
this says ≳6× from the inside. **The two are the same sign, and that agreement is the finding.**

**4. What must NOT be done with this.** Do not fit `ld_loss` and `fst_loss` jointly and report the
compromise N. With `r` and `m` both free and both confounded with N, that compromise is set by
the *relative weights*, not by the data — the §7.4.2 failure in a new costume. LD earns its place
only if `r` is pinned on defensible external grounds, exactly as §7.4.3 says of `m`.

**Caveats, stated plainly.** Two POPMULTs, one seed, no replicates. 2019's simulated halfway
(562–1000 bp at POPMULT=2000) is out of line with 2015/2023 and with its own observed bin;
unexplained.

#### 7.5.4 The sweep — `ld_loss` HAS an interior minimum, at POPMULT ≈ 1500 [VERIFIED 2026-09-07]

Six points, thin=25, seed 1, `bin_lo ≥ 562`, mask applied. This is the measurement §7.5.3 said had
to come next, and it answers the question that gates everything else.

| POPMULT | `ld_loss` | 2015 | 2019 | 2023 |
|---|---|---|---|---|
| 500 | 0.02409 | 0.02468 | 0.02481 | 0.02278 |
| 1000 | 0.00691 | 0.00742 | 0.00576 | 0.00754 |
| **1500** | **0.00299** | **0.00271** | **0.00301** | **0.00326** |
| 2000 | 0.00631 | 0.00575 | 0.00586 | 0.00734 |
| 3500 | 0.00908 | 0.00779 | 0.00859 | 0.01086 |
| 5000 | 0.01126 | 0.00940 | 0.01110 | 0.01327 |

**A clean V over 10× in POPMULT, monotone on both arms, minimum at ~1500 in all three years
independently.** So **LD is an identifying statistic, not a one-sided bound** — this is exactly
what π failed to be (§7.2.1: monotone, saturating, no optimum). It is the second equation §7.5
was after.

**But the minimum is below the prior floor (2000) and 8.3× below F_st's 12,469 (§6.7).** That
does not resolve §7.5.3's conflict, it sharpens it. Read through the `r` confound: if `r` is
inflated by factor k, LD's preferred N is deflated by k, so **k ≈ 8.3 → `r` ≈ 3.3e-7** — between
the current 2.75e-6 and §6.8's externally-argued 2.75e-8, and the same direction as §6.8. Two
independent routes to "`r` is too high" now agree in sign and differ in magnitude by ~12×.

**RETRACTED: floor subtraction.** On the two points of §7.5.3 it looked as though 30–55% of
`ld_loss` was a constant `1/n` pedestal, and subtracting each curve's asymptote looked like the
fix. **Six points do not support it.** The offset **flips sign** between POPMULT 500 (+0.0033) and
1000 (−0.0020), so it is not a fixed pedestal; and subtracting it *flattens the minimum badly* —
1500 vs 2000 goes from 111% apart (raw) to 18%. **Keep the raw L1 on levels.** The lesson is
§7.4.1's again: a decomposition fitted on two points is not a decomposition.

**Cost of the sweep, and how to run one.** SLiM's output path is hardcoded in
`CPBSampleSim*.slim`, so forward runs must be **serial** — but they are cheap (4.9 s at POPMULT
500 to 85.8 s at 8000; 2.6 min for five points). Recapitation is the expensive half and
**parallelises** once each point has its own tree: `ld_probe.py --slim-only --raw-ts` writes them,
then N processes run concurrently. 500/1000/1500/3500 took 305/398/489/1007 s and the batch cost
about what its slowest member cost alone. **Cores are not the limit; memory is** — analysis peak
goes as POPMULT^1.10 (§3.1), so ~9 GB for those four on a 15.2 GB box.

> **POPMULT=8000 was OOM-killed locally [2026-09-07]**, at a projected ~12.8 GB. §3.1 recorded
> 12000 failing on this box; **the local ceiling is actually between 5000 and 8000.** Not a
> problem for the sweep — the right arm is established by 2000/3500/5000 — but it tightens the
> §3.1 note, and it means any local work above POPMULT 5000 needs CHTC.

#### 7.5.5 The `ld_loss` noise floor — 2.8% of the signal range. It clears [VERIFIED 2026-09-07]

§7.3's test, applied to LD. Four seeds at **POPMULT=1500** (the minimum), each re-rolling all
three dice — SLiM's forward mating, the recapitation genealogy, and the mutation overlay:

| seed | `ld_loss` |
|---|---|
| 1 | 0.00299 |
| 1001 | 0.00348 |
| 1002 | 0.00411 |
| 1003 | 0.00426 |
| **n=4** | **mean 0.00371, sd 0.00058, CV 15.7%** |

| ratio | value | read against |
|---|---|---|
| noise sd / across-sweep range (0.0211) | **2.8%** | `fst_loss` 2% (fine), `pi_loss` 30% (marginal) |
| noise sd / depth of the minimum vs POPMULT=2000 | 17.6% | ~5.7σ separation |
| noise sd / depth vs POPMULT=1000 | 14.9% | ~6.7σ separation |

**Verdict: `ld_loss` sits with `fst_loss`, not with `pi_loss`.** Its run-to-run spread is 2.8% of
the range the statistic moves across the sweep — inside the band §7.3 set out in advance as "fine"
— and the minimum is separated from both its neighbours by ~6 noise sd. **The V of §7.5.4 is
signal, and the location of its minimum is resolved at the 500-POPMULT grid spacing.**

**Caveats.** Four replicates at ONE point (POPMULT=1500), as with §7.3's three at 5000 — the floor
is not known to be constant across the prior, and `ld_loss`'s level varies 7× across the sweep, so
its variance plausibly does too. Note also that **seed 1 (0.00299) was the low draw**: the sweep in
§7.5.4 is all seed 1, so the minimum's true depth is nearer 0.0037 than 0.0030. That does not move
the minimum's *location*, which is what the inference uses.

**Still open before `ld_loss` can be fitted:** weights from a pilot batch (§7.4.1, never guessed).

#### The matrix-size problem, and why it is not one

LD is pairwise over **SNPs**, so the empirical matrix is ~6.5M² on chr1 and the simulated one is a
few thousand², with no correspondence between their SNPs. **You never compare the matrices.**
Reduce both sides to `mean r² per physical-distance bin` — a short vector indexed by *distance*,
not by SNP identity, so it is dimension-free and works at any SNP count. This is the same move
`ibd_slope` already makes: an n×n F_st matrix plus an n×n distance matrix collapse to one scalar.

#### What already exists

`ToUseOnBeagles/CalculateLD.py` **computes the empirical side, and has now been run** (§7.5.1):
per-year, per-subpop mean r² binned by physical distance, pooled across chromosomes, phased 0/1
haplotype r², verified against `tskit.ld_matrix(stat="r2")`. Binning, MAF and stages all live in
`ld_common.py` (24 log bins, 1 bp → 1e6 bp, `MIN_MAF = 0.05`). Output is bin lo/hi/mid plus an
`r2_` and an `n_` column per population, in popfile order.

> **The original `MAX_DIST = 100_000` / `BIN_SIZE = 1_000` linear bins implied `r ≈ 1e-8` and were
> wrong under every reading of §6.8** — re-binned to log spacing 2026-09-06, which is what let the
> run *measure* the scale instead of assuming it. Settled: the decay is at ~60 bp (§7.5.1).

#### The build — steps 1-4 DONE, 5-6 remain

**The bin scale could not be derived** (§6.8: ρ's units are ambiguous over 1000×), so log-spaced
bins spanning the ambiguity were used and **one run measured it**: the decay is at ~60 bp (§7.5.1),
and neither reading of `ρ_HAN` was right. `BIN_EDGES` is 24 log edges, 1 bp → 1e6 bp, `d = 0`
excluded (msprime can stack mutations on one integer position). The 1e6 bp simulated sequence is
ample. Linear 1-kb bins, the original design, were wrong under every reading.

| # | where | status |
|---|---|---|
| 1 | `ToUseOnBeagles/CalculateLD.py` | **DONE 2026-09-06** — log edges, parallel across chromosomes (`LD_WORKERS`, default `min(16, cpu_count)`). Parallel output **byte-identical** to serial (18/18) because chromosomes pool in chromosome order, not completion order — float addition is not associative. **Memory, not CPU, is the limit:** ~2 GB int8 per worker for chr1, so 16 workers can want 50+ GB; 9 workers is nearly the same wall time at half the memory. |
| 2 | `Python_Code/ld_common.py` | **DONE 2026-09-06.** Shared binning; both sides import it, so it must be **copied** next to `ToUseOnBeagles/`. Hash `4d1d1d92b25b`. |
| 3 | `AnalyzeTreeSeq.py::calculate_ld_decay` | **DONE 2026-09-06.** Per-deme binned r² at real `n_i` → `data/Output_Data/ld_{year}.csv`. |
| — | **the empirical run** | **DONE 2026-09-07** (§7.5.1). 17/17 chromosomes, 7.7 h, spec matched. |
| 3.5 | `diagnostics/ld_probe.py` | **DONE 2026-09-07.** The gate; it passed (§7.5.3). |
| 4 | `ABCAnalysisNoRedis.py` | **DONE 2026-09-07.** `ld_loss` in `calculate_losses`; pools demes `Σsum/Σcnt` (invariant 4), applies `get_keep_mask` (§7.5.1 pt 1), restricts to `bin_lo ≥ 562` (pt 5). |
| 5 | `abc_standardize.py` | **WILL NOT BE DONE.** Weights exist (§7.6.1) but §7.6.2 shows adding `ld_loss` to `D` makes F_st stop contributing, and §7.8 shows why pinning `r` does not rescue it. **`ld_loss` is retired as a fitted statistic** — this step is closed, not pending. |
| 6 | `diagnostics/` | Verified against tskit for the empirical binning and the simulated extractor (§7.5a). A committed regression test is still outstanding. |

#### Simulation side — the four things that are not obvious

**(a) LD is the first statistic here that REQUIRES subsampling — and the effect is ~5–10×, not
~1%. [VERIFIED 2026-09-06]** All four existing statistics share one `pop_samples` list built from
`ts.samples(population=idx, time=time)` — whole demes, 301–714 diploids (§6.6). Fine for π
(unbiased at any n) and ~1% on F_st. **Not remotely fine for r².** Measured on a known-truth
3-deme msprime model (Ne=400, r=1e-6, ρ=1.6e-3), far-field mean r²:

| sample size | haplotypes | far-field r² | `1/n_hap` |
|---|---|---|---|
| **7 diploids** (the real n) | 14 | **0.0799** | 0.0714 |
| 50 diploids | 100 | 0.0165 | 0.0100 |
| 60 diploids (deme cap) | 120 | 0.0156 | 0.0083 |

The far-field level is **almost entirely the `1/n` floor** (Hill 1981, §11) — and it moves **5×**
between n=7 and n=50. Whole-deme simulated demes would sit near 0.008–0.016 against an empirical
0.080: a 5–10× level mismatch, dwarfing any real signal. Compare F_st's whole-deme bias of ~2%
(§6.6). **This is the single biggest trap in the LD work.**

So `calculate_ld_decay()` takes its own sample sets, cut to each site's real `n_i` from
`popFile{year}` (same source and order as `get_keep_mask`). Do **not** repoint the other four
statistics at them: §6.6 measured that subsampling makes `pi_loss` worse.

**Verification passed [2026-09-06]:** subsampling returns whole individuals (7 diploids → 14
nodes, exactly 2 per individual) and caps at deme size; the recovered decay crosses halfway in
bin **562–1000 bp** against a theoretical `d_half` of **625 bp**; and both sides report the same
`ld_common.spec_hash()`. The empirical binning was separately checked against a brute-force pair
loop — identical counts, max sum-r² difference 1.7e-13.

**Expect a short-distance plateau well below r²=1** (~0.45 measured, not ~0.99). That is the
`MIN_MAF` ceiling — two SNPs at different allele frequencies cannot reach r²=1 — and it is a
second reason the MAF filter must be identical on both sides.

**(b) Do not call `ts.genotype_matrix()` on the full tree.** 255,690 sites × 70,078 samples is
1.8e10 entries. Extract per deme instead — 14 haplotypes × 255k sites is 3.6M, trivial. Then
dedupe stacked positions (msprime places mutations on integer positions, so several can share one)
and drop non-biallelic sites, matching the empirical filter.

**(c) Cost — the real model is §7.5.2, not the original analysis.** That timed `decay_one_pop`
in isolation and named `thin` as the cost knob. **In the real call `ts.genotype_matrix()` is 98%
of the cost** (51.3 s vs 1.2 s per deme) and thinning does not touch it; the old code paid that
extraction once per deme, 61 times over the same tree. Batching it to one call per year takes
3 years from 53 min → 20.6 s at `LD_THIN = 25`. What survives from the old table is only the
*shape* of the pair-enumeration cost — quadratic in site density — not any per-trial figure.

`calculate_ld_decay(..., thin=)` is an input-level site thin, distinct from `ld_common.STAGES`
(which thins only beyond 10 kb). It is the right place for the knob for two reasons: **it is
unbiased**, measured at −2.3/−2.5/−2.7% across the 10 kb boundary where the pair count drops 625×,
so the two sides may thin differently and stay comparable in expectation; and **it is
simulated-side-only and NOT in `spec_hash()`**, so it cannot invalidate the 7.7 h empirical run.
At the fitted cut (`bin_lo ≥ 562`) thin=25 gives ~98 bp spacing — ample. `ld_loss` moves 0.8%
between thin=1 and thin=25.

> **Do NOT change `STAGES`.** It is inside `spec_hash()`, so any edit invalidates the finished
> empirical run and buys another 7.7 h on the Beagle machine. A cheaper variant exists (an
> intermediate stage at 1 kb) but would only buy resolution *below* the fitted region.
> `ld_probe.py` hard-stops on a spec mismatch for exactly this reason.

> **Generalisable lesson:** a micro-benchmark of the inner loop is not a cost model of the call.
> The 2% that `thin` actually bought was hiding behind a step nobody timed.

**(d) Write per-deme, decide pooling later.** Per-deme curves at n=5–8 are noisy, and the POPMULT
signal lives in the *decay position*, which is common across demes — so `calculate_losses` will
probably pool. Write the per-deme matrix anyway and pool at scoring time; that keeps the choice
reversible and matches the empirical file's layout.

**(e) Column ordering — CHECKED 2026-09-06, no remap needed.** `CalculateLD.py` writes columns in
**popfile order**, while every other matrix in the project is in **specifier-matrix order** (§4) —
and getting that wrong was a real bug in `CalcGenRel.py`. Verified directly: popfile order and
specifier order are **identical in all three years** (24/17/20 rows, same names, same sequence), so
the empirical LD columns line up with everything else as written. **This is a checked fact, not a
guarantee** — it would break silently if `GeneratePopulationsFile.py` were ever re-run with a
different site ordering, so `calculate_losses` should assert it rather than assume it.

**Cost trap to design around:** `ts.ld_matrix()` is O(L²). The POPMULT=5000 tree carried **255,690
sites** → 6.5e10 entries, which will not allocate. Only pairs within `MAX_DIST` are wanted, so use
the same **block scheme `CalculateLD.py` already uses** — every pair within `MAX_DIST` lies in the
same or an adjacent `MAX_DIST`-wide block, which turns it into a sequence of small BLAS matmuls.
Port that approach rather than calling `ld_matrix` on the full site set.

**Must match on both sides or the curves are incomparable:** distance bins, `MIN_MAF`,
biallelic-only, `MAX_DIST` (capped by the simulated sequence — 1e6 bp against ~55 Mb empirical
chromosomes, so compare only over the overlap), and the per-site sample sizes `n_i`.

**Known systematic that does NOT cancel:** Beagle imputation inflates empirical LD relative to
tskit's perfect phase. Unlike the r² small-sample bias, subsampling does not fix this. The honest
fix is to push simulated genotypes through the same imputation (mask to the empirical missingness
pattern first) — expensive, and worth deciding deliberately. Until then, treat the LD *level* as
suspect and lean on the *shape* of the decay curve.

---

### 7.6 Batch 2 (the LD pilot) — LANDED 2026-09-08. Weights derived; LD must NOT be fitted yet

**200 jobs x 3 trials = 600 trials**, submitted 2026-09-07 on the `ld_loss` wiring. Its only
purpose was producing `WEIGHTS` for `abc_standardize.py` (§7.4.1). It did that, and it also
produced the reason not to use them yet (§7.6.2).

**Batch health: clean.** 200/200 jobs, **600/600 trials, zero lost**, `pop` coverage uniform
across all ten deciles (max |z| = 1.8). Pooled to `out/batch2/abc_results_pilot.csv`.

#### 7.6.0 Two delivery bugs, both silent, both caught by a guard rather than by inspection

**1. Every returned CSV carried batch 1 stapled to the front.** Each of the 200 files was 2,498
rows: batch 1's *entire pooled result* (2,495 rows) followed by that job's 3 real trials. Verified
byte-identical (mod CRLF) over the first 2,496 lines in all 200 files.

**Cause: `out/abc_results.csv` is TRACKED IN GIT** (committed in `03f4e68`). `run_code.sh` clones
`origin/main`, so every job started with batch 1's pooled file already sitting at its output path;
`needs_header` correctly saw a non-empty file and skipped the header, the job appended its 3 rows,
and the whole thing came back. **Untrack the pooled result files or this recurs on every batch.**

`collect_batch.py`'s exact-header guard rejected all 200 rather than pooling batch 1 into batch 2 —
that guard earned its keep. Recovery: originals preserved in `out/pilot_raw_asreturned/`, real rows
(lines 2497+, job id from the filename) rebuilt into `out/pilot_raw/`.

> **New standing rule:** never commit a file the CHTC wrapper writes to. The wrapper's output path
> and the repo's tracked contents share a namespace, and a collision is invisible in the output —
> the CSV looks fine, it is just 2,495 rows of the wrong batch.

**2. `collect_batch.py` silently dropped `ld_loss` from every section.** It was in
`EXPECTED_FIELDS` (so the files parsed) but missing from `LOSSES` and `FITTED`, so the first run
reported weights for pi and F_st only — i.e. nothing for the batch's sole purpose. **This is
§10.2's bug class in a fourth file.** Fixed, along with an `ld_loss` entry in `NOISE_FLOOR`
(0.00074, the §7.5.5 replicates under the same mean-pairwise-|diff| convention as §7.3's entries).

#### 7.6.1 The weights — LD carries the most `pop` information of anything measured

Rank-space variance decomposition, n=600, unique R² (`collect_batch.py` section 7):

| loss | R2_tot | **pop** | total_migr | m | numClust | mu | demographic | mu-nuisance |
|---|---|---|---|---|---|---|---|---|
| `pi_loss` | 0.292 | 0.080 | 0.003 | 0.048 | 0.016 | **0.149** | 0.147 | 0.149 |
| `fst_loss` | 0.614 | 0.079 | 0.089 | **0.403** | 0.014 | 0.000 | 0.585 | 0.000 |
| **`ld_loss`** | 0.651 | **0.258** | 0.186 | 0.115 | 0.050 | 0.000 | **0.610** | 0.000 |
| `ibd_loss` | 0.839 | 0.064 | 0.010 | **0.741** | 0.019 | 0.000 | 0.833 | 0.000 |
| `genrel_loss` | 0.612 | 0.035 | 0.000 | 0.579 | 0.002 | 0.001 | 0.616 | 0.001 |

```
WEIGHTS = {"pi_loss": 0.110, "fst_loss": 0.436, "ld_loss": 0.455}
```
Frozen sigma: `pi_loss` 0.0098757, `fst_loss` 0.0035568, `ld_loss` 0.0014379.

**Three things this measured.**

1. **`ld_loss`'s unique R² on `pop` is 0.258 — 3.3x `fst_loss`'s 0.079, and the largest `pop`
   loading of any statistic in the project.** It takes 0.000 from mu. This is the LD statistic
   doing exactly what §7.5 was built for.
2. **§7.4.1's finding replicates on an independent batch.** pi still loads more on `mutation_rate`
   (0.149) than on `pop` (0.080), so the noise-floor rule would over-weight it 3x (0.353 vs 0.110).
3. **LD is NOT migration-free, correcting §7.5's premise.** It loads 0.186 on `total_migration` and
   0.115 on `m` — together larger than its `pop` loading. Migration changes deme structure, which
   changes LD. So LD is a *better* handle on N than F_st, **not a clean second equation.**

> **Why the demographic-R² rule and not the noise floor.** The floor rule counts all reproducible
> spread as signal, and mu-driven spread is perfectly reproducible — it cannot see the problem.
> Demographic R² is a share of *total* variance, so it discounts replicate noise automatically:
> `ld_loss`'s (0.515)² = 0.265 noise share sits inside its 0.349 unexplained. The rule's real
> weakness is that it **cannot tell useful signal from useless** — `ibd_loss` scores the table's
> highest 0.833 and is ~89% `m`, against a target §7.1 shows is statistically zero. It only
> apportions among statistics already judged admissible on other evidence.
> **Keep the noise floor as a VETO, not a weight**: if a fitted statistic's floor/sigma ever
> approaches 1, drop it rather than down-weight it. That is the job §7.3 actually assigned it.
> Note `ld_loss` has the worst floor/sigma of the three fitted (0.515 vs F_st's 0.086) and its
> floor was measured at ONE POPMULT — 1500, *below* this prior's floor and at the sweep minimum
> where the level is 3-4x lower than anywhere in the prior. That is the report's least reliable row.

#### 7.6.2 Conditioning: adding LD to `D` is a SURRENDER, not a compromise [VERIFIED 2026-09-08]

`diagnostics/condition_slices.py` — §7.4.3's analysis, now a script, run on the pilot under two
weightings. No simulation; runs in seconds.

Adding LD tightens `pop` dramatically, and makes the median far less sensitive to what you pin:

| slice | pop IQR ratio, pi+F_st | pop IQR ratio, pi+F_st+LD |
|---|---|---|
| *(unconditioned)* | 0.93 | **0.78** |
| `m` in [3e-5, 5e-5] | 0.40 | **0.23** |
| `m`[3e-5,8e-5] x `tm`[.10,.20] | 0.71 | **0.13** |
| `m`[3e-5,8e-5] x `tm`[.02,.20] x nC=1/2/3 | 0.41 / 1.17 / 0.69 | **0.13 / 0.17 / 0.14** |

Median `pop` swings 9,686 -> 20,348 (2.1x) across slices under pi+F_st, but sits at 4,320-6,026
under pi+F_st+LD. On its face that is §7.4.3's problem solved.

**It is not. The decisive number:**

```
D built from   n_acc   pi_loss  fst_loss   ld_loss  pop med
pi+fst           120   0.02213   0.00580   0.01413    15200
pi+fst+ld        120   0.02488   0.00800   0.01225     8025
batch median     600   0.02843   0.00806   0.01343    14010
```

**With LD in `D`, the accepted set's median `fst_loss` is 0.00800 against a whole-batch median of
0.00806 — the top 20% fit F_st no better than a random draw.** LD at weight 0.455 has taken over
the ranking entirely. Under pi+F_st alone, F_st improves 0.00806 -> 0.00580, most of the way to
the best available 0.00541.

Corroborating: **corr(`ld_loss`, `fst_loss`) = -0.686** across the prior, and within this batch
their individual optima sit at `pop` ~ 9,686 (F_st) and ~ 3,408 (LD) — 2.8x apart, understating
the gap since §7.5.4's controlled sweep put LD's minimum at 1,500, below the prior floor.

**So the `pop` posterior is not tightening because two statistics agree. It is piling against the
low edge of the prior because LD is monotone across all of it, with F_st providing weak opposing
pressure.** A tight posterior from two conflicting statistics is worse than a broad one from a
single statistic, because it looks like a result. Batch 1 at least failed *visibly* (IQR ratio
1.0); this would fail silently.

**Caveat:** slices hold 24-77 rows, so 5-15 accepted points. Read the pattern, not any single cell.

**Verdict: do not run a large batch, and do not add `ld_loss` to `FITTED_STATS`.** A pi+F_st batch
reproduces batch 1 (unconditioned `pop` IQR 0.93 on 600 trials — non-identification confirmed
independently); a pi+F_st+LD batch produces a sharp answer set by the weight ratio. Neither is
worth the compute. **The blocker is `r` — see §7.7.**

### 7.7 Should `r` be a free ABC parameter? NO — it removes the only thing LD was added for

Asked and settled 2026-09-08. Count the equations:

| statistic | constrains |
|---|---|
| pi | `4*Ne_anc*mu` — and §6.2 shows ~7/8 of it is the *fixed* ancestral phase, so nearly nothing about N |
| F_st | `4*N*m` |
| LD | `4*N*r` |

With `r` **fixed**, `4*N*r` is an equation for N — that is the entire reason LD was added. Free
`r` and it becomes one equation in two unknowns: three unknowns (N, m, r) against two products.
**Adding `r` adds an unknown without adding information**, and strictly worsens the §7.4.2
confound rather than relieving it.

**The curve shape cannot rescue it.** Sved's `E[r²] ~ 1/(1+4*N*r*d) + 1/n` — N and `r` appear only
as a product, at *every* distance. Doubling N is indistinguishable from doubling `r`; there is
nothing in the shape to separate them. (Hayes et al.'s time-depth argument would help if forward N
varied over time; it does not in this model.)

Two lesser objections: `r` reaches SLiM since §6.3, so a wide prior would swing trial cost with the
expensive draws being the ones §6.8 says are wrong; and it would be inferring a parameter nothing
in the data can validate — §7.1's objection to `m`, repeated.

**Pin `r` from a LINKAGE MAP instead.** A genetic map gives cM/Mb from crosses with **no
population-size assumption anywhere in it** — which is the whole problem with every current
estimate, all of which divide by an `N_e` §6.1 shows is ~52x low. §6.8's corrected figure of
~2.8 cM/Mb is an ordinary insect value and a plausible target to confirm. **This single external
number unblocks the LD statistic, and with it N.**

Worth taking to the professor alongside it: the project now has **two independent estimates of how
wrong `r` is, agreeing in sign** — ~100x from Cohen et al.'s own arithmetic (§6.8) and 3-8x from
the internal LD-vs-F_st disagreement (§7.5.4, §7.6.2). That is a finding, not just a question.

---

### 7.8 LD RETIRED as a fitted statistic — the target is a mixture the model cannot be [2026-09-08]

**Decision: `ld_loss` stays computed and stays out of `FITTED_STATS`, permanently unless the
simulation's structure changes.** This supersedes §7.7's "pin `r` and LD is unblocked". `r` was
pinned (§6.8.1) and the rescaled value was tested end to end. It did not rescue LD, and the
measurements below say why.

#### 7.8.1 The rescaled `r` was tested and the prediction FAILED

Four points at **r = 8.0e-7** (the §6.1.3 value), thin=25, seed 1, against §7.5.4's sweep at
2.75e-6. Prediction from `4·N·r`: dividing `r` by 3.44 should multiply LD's preferred POPMULT by
3.44, moving the minimum from 1500 to ≈5,150.

| POPMULT | `ld_loss` @ r=2.75e-6 | `ld_loss` @ r=8.0e-7 |
|---|---|---|
| 1500 | **0.00299** (min) | 0.00649 |
| 2000 | 0.00631 | **0.00322** (min) |
| 3500 | 0.00908 | 0.00557 |
| 5000 | 0.01126 | 0.00615 |

**The minimum moved 1500 → 2000, a factor of 1.33 against a predicted 3.44.** The implied scaling
exponent is ≈**0.23**, where `4Nr` theory says 1.

Checked and **not** the explanation: the 562 bp cut is itself derived from the old `r`, and at
8.0e-7 the time-depth threshold is 1,929 bp. Rescoring at the correct cut gives the same answer,
minimum still 2000.

**Consequence — §7.5.3's central hypothesis is dead.** It read the LD-vs-F_st gap as "a
measurement of the `r` error". At an exponent of 0.23, closing the 8.3× gap would need `r` wrong
by ~10⁴×, which the linkage map excludes outright. **The conflict is not an `r` error.**

#### 7.8.2 The cut should be tuned, not derived — and the answer is invariant to it

Rescoring the stored curves at cuts from 100 bp to 10,000 bp (free; `ld_probe.jsonl` stores the
per-bin curves):

- **Range lost to the `r` correction:** 12.8% of the observed decay's dynamic range survives a
  562 bp cut, only **4.0%** survives 3162 bp. The correction pushed the window into flatter curve.
- **But the formal cut is not the best one.** V-depth (how far the minimum sits below its
  neighbours) at r=8.0e-7 is **36% at the formal 3162 bp cut and 111% at 316 bp**. At r=2.75e-6 it
  peaks at 121% at 562. **So each `r` at its own best cut gives a comparable V** — correcting `r`
  costs little once the cut is re-tuned. `1/(2·G·r)` marks where forward control becomes
  *marginal*, not where it stops.
- **The reassuring part, and the one to quote:** across a **100× range of cuts** the minimum's
  LOCATION never moves — 2000 at the new `r`, 1500 at the old, every time. The cut changes
  contrast, never the answer. State that whenever a cut is reported, because choosing a cut to
  maximise a V is otherwise a researcher degree of freedom.

#### 7.8.3 Why `ld_loss` barely responds to `r`: it is fitting the plateau

The fitted region (≥562 bp) sits ~10× past the ~60 bp decay, so most of it is plateau. Measured on
the per-deme curves (`diagnostics/ld_chrom_check.py` and the plateau check, POPMULT=2000):

- **Both sides obey Hill 1981.** The far-field plateau tracks `1/n_hap` at **+0.969** observed and
  **+0.959** simulated. The subsampling machinery works; the floors are built the same way.
- **The plateau height is mostly that floor** (~0.075) and carries no information — both sides
  reproduce it identically.
- **The signal is the small excess above it:** ~**0.013 observed against ~0.009 simulated**, a 45%
  gap. That difference *is* what `ld_loss` measures — at POPMULT=2000 `ld_loss` is 0.00322 and the
  pooled plateau offset is 0.0026–0.0032 by year, the same number.

Decay *position* is what scales as `4Nr`; plateau *height* does not. That explains the 0.23
exponent and the shallow V together.

> **Correction to an earlier reading.** "`ld_loss` is fitting the wrong part of the curve" was too
> strong. The ≥562 bp region is the **right** region — it is precisely the forward-controlled part,
> the only part that can respond to POPMULT. The problem is narrower: in that region the signal is
> a ~0.004 difference sitting on a floor eighteen times larger.

#### 7.8.4 The killer: the empirical target is a 17-chromosome MIXTURE [VERIFIED 2026-09-08]

`diagnostics/ld_chrom_check.py`, on `data/ld_per_chr/` (17 chr × 3 years). No simulation.

**Test 1 — is the 45% gap significant?** Leave-one-chromosome-out jackknife: the gap is
**+0.00392 / +0.00333 / +0.00268** by year, at **5.2 / 5.1 / 4.2 jackknife SE**. Real, not noise.

**Test 2 — is the excess homogeneous across chromosomes?** No. Chromosomes differ **5.3–5.9×**
more than the jackknife allows, and the ranking is identical in all three years — **chr6, chr5,
chr15 highest; chr17, chr16, chr13 lowest** — with chr6 at ~0.024 against ~0.009, a 2.5× spread.
Three independent sample sets giving the same chromosome order means this is a property of the
chromosomes, not of the beetles. `corr(excess, density proxy)` is negative (−0.33/−0.43/−0.42),
and chr6 is the chromosome §5.4 flags as lowest SNP density.

**That pattern is also what Beagle imputation would produce**, so Test 2 alone is ambiguous.

**Test 3 — distance shift or level shift? This is the discriminator.** Low recombination rescales
*distance* (r² depends on d only through `4Ne·r·d`), so it must slow the decay AND raise the far
field together. Imputation adds a *level* without changing decay speed. Measured: per-chromosome
decay distance against far-field excess,

```
corr(log decay distance, far-field excess) = +0.923 / +0.935 / +0.938
decay distance spans 1,985 -> 26,063 bp (2015), i.e. a 10-19x range
```

**Near-perfect correlation. It is recombination-rate variation, not imputation.** chr6 having both
the slowest decay and the lowest SNP density fits ordinary linked selection.

**And that dissolves the sim-obs gap without invoking imputation.** The simulation applies **one
uniform `r`** to its 1e6 bp locus. The real genome mixes ≥10× in effective `r`, and because the LD
curve is convex, a mixture's far-field average exceeds any single-`r` curve at the mean rate. The
numbers line up: the simulated excess (0.0075 / 0.0079 / 0.0099) sits **at or just below the
FASTEST-decaying real chromosome** (0.0088 / 0.0086 / 0.0100 — nearly exact in 2023). The
simulation reproduces the high-recombination end of the genome; the pooled empirical curve is
higher because it includes everything slower.

#### 7.8.5 Why this retires the statistic rather than prompting another fix

**The decisive arithmetic: the mixture gap (0.0027–0.0039) is the same size as `ld_loss` at its
minimum (0.0032).** So essentially the entire quantity being minimised is explainable by a known
mis-specification, which means the minimum's *location* — POPMULT ≈ 2000 — is set by that
artifact rather than by population size.

Three fixes were considered on 2026-09-08 and **all three rejected** (Sohan's call, and he is
right on each):

| option | why not |
|---|---|
| fit LD **per chromosome** | needs 17 unknown `r_k`, estimated from the very curves being fitted — circular, since LD only constrains `4N·r` |
| fit only chromosomes **near the simulated `r`** | selecting data by how well it matches the model, using the quantity being fitted. Cherry-picking, and it discards most of the genome |
| accept `r` as an **effective average** | **fails in principle, not merely in elegance.** A mixture of rates changes the curve's SHAPE, not just its position — it is a superposition, broader than any single-rate curve. No single `r` reproduces both the far-field level and the decay position. An effective average does not exist |

**What would actually fix it** is a per-chromosome recombination map derived independently of LD.
Hawthorne's 18 linkage groups and the assembly's 18 chromosomes could in principle give
per-chromosome cM/Mb from crosses, with no circularity — but that needs the 2001 AFLP markers
anchored to the 2023 assembly, which almost certainly was never done. **Worth one look; do not
count on it.**

**Imputation is NOT ruled out** — §7.5.1 pt 7 stays open, and it cannot be settled at all with the
files in hand: the Beagle inputs *are* the base data (confirmed 2026-09-08), and
`ConvertBeagleToVCF.py` shows they are hard-called alleles coded 0–3 with **no genotype
probabilities, no DR2, no missing-data marker**, so nothing records what was imputed. The only
route is asking the data provider for the pre-Beagle VCFs. What §7.8.4 does is remove the *need*
to invoke imputation, which moves it well down the list.

#### 7.8.6 What this means for the project

**The `r` line of attack is exhausted.** §6.8 → §7.5.3 → §7.5.4 → §7.7 → §6.8.1 → §7.8 is one long
thread that ends here: `r` was wrong, it was found, it was corrected, and correcting it did not
buy the identifiability LD was recruited for. **Sohan's read 2026-09-08, and it is the right one:
recombination-based approaches are probably not the productive direction.** Next session should
look for a different route to the §7.4.2 N/m confound, not another recombination refinement.

**What survives, and it is not nothing:**
- **`r` is genuinely pinned** for the first time, from a source with no `N_e` in it (§6.8.1).
- **A measured CPB mutation rate exists** (§6.1.2), retiring the midge stand-in.
- **The model's scale is understood**: Q ≈ 78.5, ancestral phase correctly rescaled, forward phase
  not, time not compressible (§6.1.3).
- **Four independent routes agree on `N_e` ≈ 2–5e5**, ~50× above Cohen's 6700.
- **The per-chromosome recombination heterogeneity is measured** — a real property of the CPB
  genome that was not known here before, and worth reporting on its own.

**And the honest bottom line is unchanged from §7.4.3: Nm is identified (≈78–103, stable across
every well-pinned slice); N is whatever the dispersal assumption makes it.** High gene flow across
the Wisconsin landscape is a real, defensible result. An N posterior is not, and nothing measured
on 2026-09-07/08 changes that.


---

## The LD reference blocks, from the old §11

### LD as an ABC summary statistic (§7.5)

**Boitard, S., Rodríguez, W., Jay, F., Mona, S., & Austerlitz, F. (2016).** Inferring population
size history from large samples of genome-wide molecular data — an approximate Bayesian
computation approach. *PLoS Genetics* **12**(3):e1005877.
[doi:10.1371/journal.pgen.1005877](https://doi.org/10.1371/journal.pgen.1005877) ·
[PMC4778914](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4778914/) ·
[code](https://forge-dga.jouy.inra.fr/projects/popsizeabc)
— **The direct methodological precedent, and the one to follow.** "PopSizeABC" infers population
size history by ABC using the allele frequency spectrum plus **average LD in bins of physical
distance between SNPs** as its summary statistics. That is exactly the construction §7.5 needs,
and it is what resolves the "the two LD matrices are different sizes" problem: you never compare
matrices, you compare binned decay curves, which are dimension-free.

**Hayes, B. J., Visscher, P. M., McPartlan, H. C., & Goddard, M. E. (2003).** Novel multilocus
measure of linkage disequilibrium to estimate past effective population size. *Genome Research*
**13**(4):635–643. [link](https://genome.cshlp.org/content/13/4/635)
— **Why LD can see what π cannot.** LD at different distances reflects `N_e` at different times
in the past: long-distance LD reflects *recent* `N_e`, short-distance LD the more distant past.
The decay curve is therefore a time series of `N_e`, not one number. This is the principled reason
to expect LD to carry POPMULT information where π does not — §6.2 shows ~7/8 of π is set by the
fixed ancestral phase, whereas long-range LD is dominated by exactly the recent window the forward
simulation controls.

**Tenesa, A., et al. (2007).** Recent human effective population size estimated from linkage
disequilibrium. *Genome Research* **17**:520–526.
[pdf](https://genome.cshlp.org/content/early/2007/03/09/gr.6023607.full.pdf)
— A worked application of the Hayes binning/curve-fitting approach; useful as an implementation
reference.

### The LD–N_e estimator and its small-sample bias (§7.5, invariant 1)

**Sved, J. A. (1971).** Linkage disequilibrium and homozygosity of chromosome segments in finite
populations. *Theoretical Population Biology* **2**:125–141. — `E[r²] ≈ 1/(1+4N_e·c)`.

**Hill, W. G. (1981).** Estimation of effective population size from data on linkage
disequilibrium. *Genetical Research* **38**:209–216.
[link](https://www.semanticscholar.org/paper/Estimation-of-effective-population-size-from-data-Hill/bb2561666f63910571f9b5dc388b5d2bec6f1445)
— The estimator. Note the form is `E[r²] ≈ 1/(1+4N_e·c) + 1/n`: **the small-sample bias is in the
original equation, not an afterthought.** At our n = 5–8 diploids that `1/n` floor is ~0.06–0.10
and may exceed the real signal.

**Waples, R. S. (2006).** A bias correction for estimates of effective population size based on
linkage disequilibrium at unlinked gene loci. *Conservation Genetics* **7**:167–184.
[link](https://www.researchgate.net/publication/227119323_A_bias_correction_for_estimates_of_effective_population_size_based_on_linkage_disequilibrium_at_unlinked_gene_loci)
— The standard correction for that `1/n` term, with Weir's `[S/(S−1)]²` sample-size adjustment.
**Background reading rather than something to implement**: in an ABC we subsample the *simulated*
demes to each site's real `n_i` so the identical bias appears on both sides (invariant 1), which
is more robust than correcting analytically at n = 5–8.

**Waples, R. K., Larson, W. A., & Waples, R. S. (2016).** Estimating contemporary effective
population size in non-model species using linkage disequilibrium across thousands of loci.
*Heredity* **117**:233–240. [PMC5026751](https://pmc.ncbi.nlm.nih.gov/articles/PMC5026751/)
— Confidence intervals, and the correlated-pairs problem: SNP pairs are not independent, so naive
degrees of freedom badly overstate precision. Same class of error as the Mantel argument in §7.1.

**Ragsdale, A. P., & Gravel, S. (2020).** Unbiased estimation of linkage disequilibrium from
unphased data. *Molecular Biology and Evolution* **37**(3):923–932.
[link](https://academic.oup.com/mbe/article/37/3/923/5614437)
— Relevant only if the phasing assumption is ever revisited. Our data is Beagle-phased, so we use
phased haplotype r² on both sides; the concern here is imputation bias, not phase uncertainty.


---

# §B — resolved defects [ARCHIVED 2026-09-10]

All of these are fixed and verified in code. `CLAUDE.md` §6.3–6.7 keeps the short list of what must
not be undone; the write-ups follow.

### 6.3 `recombination_rate` now reaches the forward simulation [FIXED 2026-07-29]

Was: `CPBSampleSim{Linux,Win}.slim:12` hardcoded `initializeRecombinationRate(1e-8)`, so the ABC
parameter reached **only** `pyslim.recapitate()` and could not affect forward dynamics. Both files
now take `initializeRecombinationRate(RECOMB)` and `Main.py` passes `-d RECOMB=`. Verified against
SLiM 5.1: `-d RECOMB=2.75e-06` (the format Python's `!r` emits) parses; omitting it errors with
`undefined identifier RECOMB` and exit 1, which `check=True` on the subprocess turns into a raise.

**This is a 275× increase in forward recombination** (1e-8 → 2.75e-6), so expect many more edges
per generation and materially higher SLiM memory/runtime than any historical measurement. It is
the main reason §6.4 had to be fixed at the same time.

Still fixed, not inferred, at `DEFAULT_RECOMBINATION_RATE = 2.75e-6`: recombination has **no
signal at all** in π/d_xy/F_st — it shows up only in linkage disequilibrium, so inferring it needs
an LD summary statistic first.

### 6.4 `simplificationRatio=INF` removed [FIXED 2026-07-29, unrun]

Was: `CPBSampleSim{Linux,Win}.slim:4` set `initializeTreeSeq(simplificationRatio=INF)`, telling
SLiM to **never simplify during the forward run**, so the edge table grew unbounded for all 324
generations — almost certainly the real cause of the OOM, not the mutation rate. Now plain
`initializeTreeSeq(timeUnit="generations")`, i.e. SLiM's default ratio of 10.

**Safe — checked against the docs, not assumed.** Simplification is lossless for the genealogy of
retained samples; SLiM retains all living individuals plus everything permanently Remembered
(`treeSeqRememberIndividuals(..., T)` at gens 308/316 — `permanent=T` marks them as real samples),
and future generations descend only from living individuals, so nothing needed later is discarded.
Recombination is unaffected: breakpoints are recorded at reproduction, and simplification runs
afterwards on the already-written tables. The manual frames `simplificationRatio` purely as a
speed/memory tradeoff.

> The recapitation hazard is real but belongs to the **other** mechanism. Recapitation needs the
> input roots, and pyslim is explicit that a *Python-side* `ts.simplify()` before recapitating must
> pass `keep_input_roots=True` — that is why §3 declined simplify-before-recapitate. SLiM's own
> runtime simplification already uses `keep_input_roots`, which is why ordinary SLiM output (SLiM
> simplifies every ~20 ticks by default) is routinely recapitable. Do not conflate the two.

**[OPEN] Not yet run end-to-end.** Output will not be bit-identical (nodes are renumbered,
redundant edges merged); compare **branch-mode diversity under a fixed recapitation seed**, not
file hashes.

**Memory floor this cannot touch:** gens 308/316 permanently Remember *every individual in every
subpop*, pinning all of their ancestry. Downstream only 2–19 individuals per site are ever
sampled, so Remembering a bounded subset (say 50/subpop) would cut retained ancestry a lot. That
changes what is available to the sampler, so it is a design decision, not a free win.

### 6.5 Resolved [VERIFIED]

- **pixy's `count_comparisons` overflowed int32** (fixed 2026-07-28). The field saturates to
  `INT32_MIN` once `comparisons_per_site × no_sites > 2^31`, i.e. at **≥14 diploid individuals**.
  Only `H53-2015` (19) ever crossed it, corrupting **244 of 2015's rows** — but just 18 showed a
  visible sign flip (π = −0.025); the rest, including two `Alsum25` d_xy pairs that saturated on a
  single chromosome, were merely **~18% high and looked entirely plausible**. Fixed by computing
  denominators analytically (§5.1). **Still live in two ways:** any future year with ≥14 sampled
  individuals re-triggers it (2023's max of 308 comparisons/site sits just under the 327.9
  threshold), and **any `empiricalStats` output produced before 2026-07-28 is suspect.**
- **Migration saturation.** The old code normalized each row of `exp(-d·modifier)` to sum to 1,
  making ~76% of each subpop's offspring immigrants every generation — effectively panmixia, which
  is why simulated F_st was tiny and d_xy ≈ π. Fixed by the `total_migration`/`scale` split (§2).
- **KMeans re-randomization.** `random_state=random.randint(0,1000)` gave a different cluster
  layout every run. Seed pinned at 42.
- **`csv.DictReader` ate subpop 0.** The π files have no header, so DictReader consumed the first
  value as a column name. Replaced with `_read_vector`/`_read_matrix`. (The remaining DictReader in
  `read_parameters_from_csv` is correct — that input CSV genuinely has a header.)
- **Missing header on `abc_results.csv`.** `csv_exists = Path(output_csv).exists()` skipped the
  header whenever the file existed — and CHTC's `run_code.sh` **pre-creates** it (so an evicted job
  fails with a real exit code instead of an errno-2 transfer hold), so every job emitted data rows
  with no column names. Now `needs_header` also treats a zero-byte file as needing one.
- **Coalescent rescaling cannot be applied as the model stands.** Standard rescaling (N→N/Q, μ→μQ,
  r→rQ, m→mQ) needs `m → m·Q`, but m was ≈0.76 — even Q=2 exceeds a probability; and the dominant
  term is fixed at `ancestral_Ne`, so the largest contribution to π is unaffected. Both blockers
  are parameter problems, not coalescent subtleties. Worth re-asking once §6.1/§6.4 are settled.

### 6.6 Whole-deme sampling — measured, ~1% on F_st, nothing on π [RESOLVED 2026-08-14]

Numbered last only to keep §6.5's cross-references stable.

`AnalyzeTreeSeq.py:29-31` (and `mu_calibrate.py`, which mirrors it) builds each sample set as
`ts.samples(population=i, time=t)` — *every* node in the deme at that timepoint, because gens
308/316 permanently Remember every individual in every subpop (§6.4). The same `pop_samples` list
then feeds **all four** statistics: π, d_xy, F_st, relatedness. Measured on the POPMULT=5000 tree
(17,372 diploids over 33 demes):

| year | simulated diploids/deme | observed n |
|---|---|---|
| 2015 | 301–714 | 4–19 |
| 2019 | 301–714 | 5–7 |
| 2023 | 301–714 | 5–11 |

**This section previously read that mismatch as a violation of §8 invariant 1 and prescribed
"subsample the simulated demes and recompute." That prescription was wrong for π and unproven for
F_st.** The two statistics differ in kind, and conflating them is what produced the bad advice.

**π — not a problem. Do NOT subsample.** Pairwise diversity is unbiased at any n ≥ 2, so the
whole-deme and n=5 estimates target the *same number*; only the precision differs. `pi_loss` is a
plain L1 distance with no variance-matching term, so degrading the precise side would add a
near-constant penalty to every draw while inflating each trial's Monte-Carlo variance — i.e. it
would raise the noise floor the pass has to clear, in exchange for nothing. §7.2.1 already
measured this exact geometry: uncorrelated scatter scored **worse** than no scatter (shuffled
0.0224 vs flat 0.0209). **There is no requirement to replicate sampling noise on the simulated
side.** §8 invariant 1's "match the sampling design" is aimed at estimators that are *biased* at
small n; π is not one, and the invariant should be read that way.

**F_st — the real question, and it is about bias, not noise.** F_st is a *ratio*, so its estimator
carries an O(1/n) bias whose size depends on the true level. That is systematic: it does **not**
average out over replicates, and it shifts the very quantity the fit reads — simulated 0.0079
against observed 0.00645 is what puts POPMULT ≈ 6000 (§6.2). Genome-wide averaging over many
independent trees suppresses ratio bias substantially, so it may well be negligible; it had simply
never been measured.

**Measured [VERIFIED 2026-08-14] — `diagnostics/fst_subsample.py`, POPMULT=5000, numClusters=33,
seed 1, μ=4.646e-7, 100 replicates.** Recapitate + simplify + mutate once (1405 s, peak 4.7 GB,
70,078 samples, 255,690 sites), then draw each deme's *real* n_i diploid **individuals** (not
random nodes — the empirical unit is a diploid) and recompute π and F_st on those sample sets from
scratch (§8 invariant 2). The whole-deme baseline reproduces production **exactly**
(`fst_loss` 0.00830 against §6.1.1's 0.00830; `pi_loss` 0.02110 against 0.02132), which is what
licenses reading the rest.

**F_st — the bias is real, positive, and ~1%. Not enough to matter.**

| year | whole-deme | subsampled | bias | vs obs |
|---|---|---|---|---|
| 2015 | 0.007579 | 0.007727 | **+1.95%** | obs 0.009149 |
| 2019 | 0.008182 | 0.008174 | −0.09% | obs 0.003199 |
| 2023 | 0.007853 | 0.007933 | **+1.03%** | obs 0.007007 |

`fst_loss` moves **0.00830 → 0.00847 ± 0.00022**, a **+2.0%** shift. The sign is as theory
predicts (smaller within-deme samples slightly deflate within-population diversity, inflating the
ratio), so this is a genuine small-sample bias, not Monte-Carlo wobble. But since F_st tracks
`1/(1+4Nm)` (§6.2), a +1–2% level shift moves the implied POPMULT by +1–2% — **POPMULT ≈ 6000
becomes ≈ 6060**, against a prior of `U(2000, 12000)`. It is far inside the per-replicate spread
(per-pair sd across replicates 0.0016–0.0018, i.e. ~20% of the level).

**Verdict: whole-deme sampling stays. §6.6 requires no pipeline change.** Note the direction if it
is ever wanted: correcting it would raise simulated F_st slightly and therefore push inferred
POPMULT slightly *up*.

**π — subsampling makes the fit WORSE, exactly as predicted.** `pi_loss` **0.02110 → 0.02161 ±
0.00068**, i.e. injecting sampling noise moves it *away* from the observed vector and further above
the flat-simulation floor of 0.02094. Log bias is −0.0003 to −0.00002, confirming unbiasedness at
every n. This is the direct measurement that the old prescription would have degraded the
distance.

**And it kills the §7.2.1 alternative explanation.** Sampling sd at the real n_i is 0.0050–0.0053
in log units against an observed between-site log sd of 0.0423/0.0161/0.0477 — so observed-side
sampling noise is only **1.2% / 9.6% / 1.2% of the observed between-site variance**, attenuating a
site-level correlation by a factor of **0.992 / 0.955 / 0.994**. §7.2.1's Pearsons
(−0.091/+0.538/−0.021) de-attenuate to −0.092/+0.563/−0.021. **Nothing changes: the missing
site-by-site covariation is genuinely absent, not buried in observed-side noise.** By the same
arithmetic the observed spread is *not* mostly sampling noise, so "the sim is too flat" survives
as a real (if smaller, per §7.2 impl. 3) statement.

**Independent of sample size, and untouched by any of this — ~~[OPEN]~~ RESOLVED 2026-08-26, see
§6.7:** the empirical F_st comes from pixy's Weir–Cockerham (SNP-weighted mean, §5.2) while the
simulated side used `ts.Fst`. That estimator mismatch survived this test and was indeed **the
larger of the two concerns** — it was worth a factor of ~2, against ~1% here. Note the old wording
of this paragraph called `ts.Fst` "tskit's Hudson-style ratio", which was itself wrong: it is
Nei/Slatkin, and that error is exactly what let the mismatch sit unexamined.

---

### 6.7 The F_st estimator mismatch — `ts.Fst` is Nei, not Hudson [FIXED 2026-08-26]

**This was worth a factor of ~2 in the inferred POPMULT.** It is the largest defect found since
the §5.1 denominator bug, and like that one it was invisible in the output: both sides emitted
plausible small numbers on the same scale.

**What each side computed. [VERIFIED]**

| side | estimator | formula |
|---|---|---|
| simulated (`AnalyzeTreeSeq.py`, `ts.Fst`) | **Nei (1973) / Slatkin (1991)** | `1 − 2(π_X+π_Y)/(π_X + 2d_xy + π_Y)` = `(d_xy − Hw)/(d_xy + Hw)` |
| empirical (pixy, `--stats fst`, WC default) | **Weir–Cockerham** | variance components `a/(a+b+c)` |

writing `Hw = (π_X + π_Y)/2`. Hudson is `(d_xy − Hw)/d_xy`, so **Nei and Hudson differ by the
factor `(1 + Hw/d_xy) ≈ 2`** at low differentiation — the denominators differ, the numerators do
not. Nei is not a "Hudson-style ratio"; it is roughly half of one.

Verified from the project's own output, not from the docstring: rebuilding the Nei identity from
`diversities_*.csv` and `divergences_*.csv` reproduced `fst_*.csv` to **2.2e-16**, and the implied
`fst_loss` reproduced noise-floor rep1's 0.008715 exactly.

**WC and Hudson estimate the same parameter; Nei does not. [VERIFIED — known-truth msprime
2-deme model, §10 convention.]** WC implemented directly from Weir & Cockerham (1984) —
scikit-allel is not in `cpb-env`, it lives on the Beagle machine.

| n diploids | `ts.Fst` (Nei) | Hudson | WC (ratio-of-sums) | **WC/Hudson** |
|---|---|---|---|---|
| 12 | 0.02008 | 0.03936 | 0.03951 | **1.004** |
| 25 | 0.00645 | 0.01281 | 0.01346 | 1.051 |
| 100 | 0.00538 | 0.01070 | 0.01082 | **1.011** |

with **Nei/Hudson = 0.5027** on the same data. So the simulated side was reporting ~half the
statistic it was being fitted to.

**Size of the error.** Year-averaged, fitted mask applied, the sim-vs-obs estimator gap is
**0.008026** — against an across-prior `fst_loss` range of 0.01026 (**78%**) and a run-to-run
noise floor of 0.00017 (**47×**). This is not a few percent, and §7.3 put essentially the whole
inference on F_st.

**Consequence — it halved the inferred POPMULT.** Using `F_st ≈ 1/(1+k·POPMULT)` fitted to the
§6.1.1 sweep (k = 0.0247, reproduces all three sweep points to <3%):

| convention | sim F_st @ POPMULT=5000 | obs target | implied POPMULT |
|---|---|---|---|
| before (Nei vs WC) | 0.00833 | 0.00645 | **6,234** |
| after (Hudson vs WC) | 0.01647 | 0.00645 | **12,469** |

**Confirmed end-to-end [VERIFIED 2026-08-26]** — the real `analyze_tree_sequence()` re-run on the
POPMULT=5000 `out/simTreeSeq.trees` (1738 s), scored by the real `calculate_losses()`:

| year | new (Hudson) | old (Nei, rep1) | ratio |
|---|---|---|---|
| 2015 | 0.016042 | 0.008102 | **1.980** |
| 2019 | 0.017087 | 0.008641 | **1.977** |
| 2023 | 0.016272 | 0.008254 | **1.971** |

The ratio is `1 + Hw/d_xy`, predicted ≈1.98 and measured 1.971–1.980 in all three years
independently — the cleanest possible confirmation that this was a pure estimator swap and not a
change in what the simulation produces. **`fst_loss` 0.008715 → 0.015714.** The other four losses
behave exactly as they should: `pi_loss` 0.0209→0.0222, `dxy_loss` 0.000152→0.000188 and
`genrel_loss` 0.000115→0.000118 all move only by the mutation-overlay redraw (this run re-drew
mutations), while **`ibd_loss` roughly doubles, 0.004825 → 0.009761 — expected**, since the IBD
slope is a regression on `F_st/(1−F_st)`, so doubling F_st doubles the slope. IBD is diagnostic
only (§7.1), but any stored `ibd_loss` is on the old scale too.

Per year, so it is not an averaging artifact: 2015 → 8,769; 2023 → 11,474; 2019 → 25,229.
**Two of three years sat at or past the old `POPMULT_MAX = 12000`.** **Resolved 2026-08-26: the
ceiling was raised to 25000** on 64 GB machines (§3.1 for the sizing, `ABCAnalysisNoRedis.py` for
the rationale). Total N at the new ceiling is ~83k, beyond both Cohen et al. figures — deliberate,
since §6.1 shows those fail their own paper's internal check by ~52× with every named bias pointing
up. Cost: **2.3× the prior volume**, so a fixed trial budget puts 2.3× fewer draws near the mode.

**The fix.** `AnalyzeTreeSeq.py` now computes Hudson from π and d_xy, which are already in hand:

```python
fsts[i, j] = 1 - 0.5 * (diversities[i] + diversities[j]) / dv
```

It **drops** the `ts.Fst` traversal rather than adding one, so it is free. Verified against tskit
on a fresh msprime simulation (§10): exact against Hudson-from-π/d_xy, and equal to the algebraic
bridge `2·Nei/(1+Nei)` to 2.15e-16. `dv == 0` returns 0.0 rather than a NaN that would silently
poison `fst_loss`.

#### 6.7.1 Two things ruled out on the way [VERIFIED]

**Average-of-ratios pooling in pixy — ruled out.** AoR pooling deflates WC 2–4× and is strongly
n-dependent (measured on the msprime model: AoR/RoS = 0.236 at n=6, 0.657 at n=12, 0.803 at
n=100 — a 3.4× trend). 2015's n runs 4–19, so that trend would be unmissable; measured
`corr(WC/Hudson, harmonic n)` = **−0.010** (2015), +0.310 (2019), +0.072 (2023). **pixy pools
properly, and its WC is the trustworthy empirical number.** This is what licenses keeping the
empirical targets untouched.

**Reconstructing empirical Hudson from `averaged_pi` and `averaged_dxy` — ruled out as a target,
and this matters independently. [OPEN as a data question]** It is tempting (free, no VCFs) but
contaminated: F_st at this level is a ~2–3% difference between two ~0.0125 numbers, so a **1%**
relative offset between π and d_xy moves it by **0.010** — larger than the whole signal. And such
an offset is present. The derived Hudson has **zero negative pairs out of 1,114** with a hard
floor at +0.011/+0.016/+0.014 by year, while pixy's WC properly scatters through zero (24/24/64
negative pairs). Scaling d_xy down **1.13% / 1.60% / 1.39%** puts the floor at zero in all three
years independently — a single systematic offset, not noise.

Leading candidate [INFERRED]: **Beagle imputation depresses within-sample π more than
between-sample d_xy**, since imputation borrows haplotypes from the sample's own panel. That
inflates any π/d_xy-derived F_st while leaving pixy's WC (computed from the same genotypes, but
not as a π/d_xy ratio) comparatively intact. Not blocking — the fix above never uses the derived
quantity — but it is a real property of the empirical data and worth raising with the professor.

**Generalisable lesson, now §8 invariant 9:** "F_st" names a family, not a statistic. Two
estimators can differ by 2× while both look like ordinary small F_st values.


---

# §C — the 99-deme cap [ARCHIVED 2026-09-10]

Lifted by `Python_Code/recapitate_util.py`. `CLAUDE.md` §7.9.5 keeps the stub and the open
modelling question about deme resolution.

#### 7.9.5 The 99-deme cap IS recapitation — an msprime API limit, and it is now lifted [FIXED 2026-09-09]

**The cap was real and it was where Sohan said it was.** `numClusters ∈ {1,2,3}` (×33 = 99 demes)
has always sat exactly one under an msprime hard limit:

```
InputError: Input error in population split: Cannot have more than 100 populations in one event.
```

`pyslim.recapitate(ts, ancestral_Ne=...)` builds a demography in which **every** SLiM
subpopulation splits from one ancestral population in a **single** `population_split` event, and
msprime caps that event at 100 derived populations.

**Measured 2026-09-09, full pipeline, POPMULT=2000, `data/` backed up and restored:**

| numClusters | outcome |
|---|---|
| 99 | **OK** — 147 s, 954 MB peak, 20 deme rows written |
| 200 | SLiM ran fine (112 s), then **died in recapitation** on the limit above |
| 400 | died earlier, in SLiM (see the deme-size wall below) |

**There are TWO walls, at different places, and they are easy to confuse.**

**Wall 1 — deme size rounds to zero, in SLiM.** `CPBSampleSim*.slim:30` does
`sim.addSubpop("p"+i, asInteger(Average Count[i] * POPMULT / numSubpops))`, `asInteger`
**truncates**, and SLiM refuses an empty subpopulation:
`ERROR (Population::AddSubpopulation): subpopulation p38 empty.` At numClusters=200, POPMULT=500
puts 2 of 200 demes under 1.0 individual. The same 200 demes at POPMULT=2000 clears it easily
(smallest Average Count 0.272 → 2.7 individuals), so **the ABC never meets this wall**; it only
bites the interactive `python Main.py` route, whose prompt defaults to POPMULT=500.
**Guarded 2026-09-09 in `Main.main`**, before SLiM runs, with a message naming the minimum POPMULT
for the layout — SLiM's own message names a subpop id and says nothing about clusters.

**Wall 2 — the 100-population event, in recapitation. This is the actual cap.**
**Fixed by `Python_Code/recapitate_util.py`:** merge in **stages**. Split the demes into groups of
≤99, merge each group into a throwaway intermediate population at the recapitation time, then
merge the intermediates into the real ancestral population at `np.nextafter` of that time.

**It is exact, not an approximation.** The second merge is ~1e-13 generations after the first and
the intermediates have size 1.0 (the same value pyslim assigns the SLiM populations, for the same
reason), so the probability of any coalescence inside an intermediate is ~1e-13. Every lineage
reaches the ancestral population at the instant it would have under a single split, and the
ancestral coalescent that follows is identical. Groups nest, so 99 groups of 99 covers 9,801 demes
at two levels and the loop adds levels beyond that.

**Verified 2026-09-09** on the 200-deme POPMULT=2000 tree: `pyslim.recapitate` raises the
InputError; `recapitate_util.recapitate` completes in **127 s** and every tree has
**`num_roots == 1`** — fully coalesced. At ≤99 populations it **delegates to `pyslim.recapitate`
with the same arguments**, so it is identity by construction at the current prior and cannot
perturb any existing result. Wired into `AnalyzeTreeSeq` and all six diagnostics that recapitate.

**WHY YOU WOULD WANT MORE DEMES — the resolution argument, and it corrects an earlier claim here.**
§7.4.2 says "no dispersal structure below ~6.5 km is representable" and pairs it with §7.1's
"sites span 1.7–160 km". **That pairing is wrong**: 1.7–160 km is the *pairwise* range, and the
quantity that matters is nearest-neighbour spacing. Measured over the **43 unique sequenced
sites**: **minimum NN 0.39 km, median NN 1.28 km.** That is exactly the scale CPB dispersal
operates at (rotation studies: 0.3–0.9 km halves insecticide need; 0.5 km can cut postdiapause
abundance ~100×; long-distance flight reaches a few km). **The data has the resolution; the model
does not.**

KMeans centroid spacing against cluster count (ceiling is 1,790 unique field coordinates):

| demes | min NN | median NN | median pairwise |
|---|---|---|---|
| 33 | 6.05 km | 10.93 km | 45.1 km |
| 99 | 2.54 km | 4.74 km | 41.5 km |
| 200 | 1.22 km | 3.02 km | 38.3 km |
| 400 | 0.69 km | 1.68 km | 36.7 km |
| **sequenced sites** | **0.39 km** | **1.28 km** | 29.7 km |

At 33 demes two fields 0.39 km apart are forced into clusters ≥6 km apart, so the short-range part
of the kernel — the only part with real biological content — cannot exist in the model at all.
Raising the count to ~200–400 puts the model's resolution at the biological scale for the first
time, and it also relieves §7.9.2: a good one-to-one site→deme matching is easy once clusters
comfortably outnumber sites.

**THE COST, and it is not the crash.** Deme size is `Average Count × POPMULT / numSubpops`, so
**total simulated N stays ≈3.33 × POPMULT no matter how many demes there are** — more demes means
the same beetles cut into finer pieces. To hold deme size at today's 505 you would need
POPMULT ≈ 30,000 at 200 demes, just past the 25,000 ceiling (§3.1). At POPMULT=25,000 and 200
demes the deme size is 416, comparable to what runs now and inside the memory envelope; **400
demes at a realistic deme size is not affordable.**

**So `numClusters` and POPMULT are entangled by construction**, which is exactly why §7.4.3
measured `numClusters` swinging median `pop` 1.7× (9,650 → 16,740 across {1,2,3}). Raising the
deme count is **not** a free improvement — it changes what POPMULT means, and any comparison
across `numClusters` is a comparison at different deme sizes.

> **[OPEN] The prior has NOT been widened.** `numClusters` is still `randint(1, 4)`. Lifting it is
> now a modelling decision rather than a technical one, and it needs deciding together with the
> POPMULT floor: at 200 demes the prior floor of 2000 gives demes of 33 individuals, which drift
> hard. A defensible pairing would be **numClusters=200 fixed** with POPMULT ≥ ~10,000, reported as
> a sensitivity axis against the current 33 — not a wider prior over both.


---

# §D — the per-file code table [ARCHIVED 2026-09-10]

The old §9. Read this when you need to know what a given script does before touching it. Note it
predates the 2026-09-10 trim, so its `ld_common` / `calculate_ld_decay` / `ld_probe.py` /
`ld_chrom_check.py` rows describe machinery that is still wired and still runs, but whose output is
not fitted.

## 9. Existing code

| File | What it does |
|---|---|
| `Python_Code/Main.py` | Pipeline entrypoint (§2). Threads `total_migration`, `ancestral_Ne`; `KMEANS_SEED = 42`. **`mutation_rate` AND `recombination_rate` have no defaults and raise if omitted** — the latter carried a hardcoded 2.75e-6 until 2026-09-09, which §10.1 claimed had been removed and had not (§6.8.1). `ancestral_Ne` defaults to `scale_constants.ANCESTRAL_NE`. **Also guards deme size before invoking SLiM** (§7.9.5): SLiM refuses a subpop that rounds to zero individuals, and its own message names a subpop id rather than the cluster count. |
| `Python_Code/AnalyzeTreeSeq.py` | Recapitate → simplify → mutate → π, d_xy, F_st, relatedness (batched `indexes=`). F_st is **Hudson, computed from π and d_xy — `ts.Fst` is deliberately not used** (§6.7, invariant 9). **`mutation_rate` and `recombination_rate` have no defaults and raise if omitted** (§10.1); `ancestral_Ne` defaults to `scale_constants.ANCESTRAL_NE`. Recapitates through **`recapitate_util`, not `pyslim` directly** (§7.9.5) — a no-op at ≤99 demes. |
| `Python_Code/ABCAnalysisNoRedis.py` | ABC driver, `calculate_losses`, IBD helpers (`get_site_geo_distances`, `ibd_slope`), `get_keep_mask` + `EXCLUDE_SMALL_SUBPOPS`/`MIN_SUBPOP_N` (§7.0). **LD as of 2026-09-07:** `_read_ld_decay`, `_ld_pooled_mean_abs_diff`, `LD_MIN_BIN = 562`, `LD_EMPIRICAL_SPEC` (a hard stop if `ld_common` drifts from the spec the empirical targets used), and a column-order assertion against the specifier matrix in `getObservedData`. CHTC entrypoint: `python ABCAnalysisNoRedis.py <job_id> [num_trials]` → `../out/abc_results.csv`. **`job_id` is a LABEL ONLY as of 2026-09-07** — it names the output file and nothing else. The RNG is deliberately unseeded: `np.random.seed(job_id)` was removed because the reproducibility it implied was never real (SLiM, `recapitate()` and `sim_mutations()` are all unseeded, so a re-run reproduced the *parameters* but never the losses), while it did create a live footgun — a job id reused across batches silently re-draws that batch's exact parameters. Parameters are recorded per row anyway. numpy seeds from OS entropy at import, not the clock, so simultaneous jobs do not collide [VERIFIED]. |
| `Python_Code/abc_standardize.py` | **Offline** post-processing: σ=1.4826·MAD per fitted stat → standardized `D` → ranked CSV + frozen σ JSON. Not run in the pass loop. **`FITTED_STATS` here is what actually decides what enters `D`.** `WEIGHTS` is **no longer `None`/equal** — set to `{"pi_loss": 0.125, "fst_loss": 0.875}` from batch 1 (§7.4.1), with the decomposition recorded in the config comment. **Still the batch-1 two-statistic set as of 2026-09-09:** the pilot's three-statistic weights (0.110/0.436/0.455) are derived and recorded in §7.6.1 but deliberately NOT installed — see §7.6.2/§7.7. **BOTH weight sets are now STALE**: they were fitted with μ free, and μ left the prior on 2026-09-09 (§7.9.3), so a large share of `pi_loss`'s variance no longer exists. **Re-derive from batch 3; do not reuse either.** Note it has **no CLI**: `RESULTS_CSV`/`RANKED_CSV`/`SIGMAS_JSON` are module constants pointing at batch 1's files, so running it unedited reads batch 1 and overwrites its ranked CSV and frozen sigmas. |
| `Python_Code/ld_common.py` | **Shared LD-decay binning — imported by BOTH sides (§7.5).** `BIN_EDGES` (24 log bins, 1 bp→1e6 bp), `MIN_MAF`, `STAGES`, `accumulate`, `decay_one_pop`, `write_decay`, `report_halfway`. The two sides run on different machines, so this file is **copied** next to `ToUseOnBeagles/` — `spec_hash()` (currently `4d1d1d92b25b`) is printed by both and a mismatch means the curves are not comparable. **Never reimplement any of this on one side only.** |
| `Python_Code/AnalyzeTreeSeq.py::calculate_ld_decay` | Simulated LD decay, gated behind `COMPUTE_LD` (default **ON since 2026-09-07**, `LD_THIN = 25`; was off until `ld_probe.py` shows the fitted region moves with POPMULT; §6.8's bin-scale half is now settled, §7.5.1). The one statistic here that **must** be subsampled to real `n_i` — r²'s `1/n` bias is ~5–10×, against F_st's ~2% (§7.5a), and the empirical plateau tracks `1/n_hap` at r=0.995 (§7.5.1). |
| `Python_Code/GenerateSimulationParams.py` | `determine_migration_rates(distances, total_migration, scale, ...)`. |
| `Python_Code/GenerateClusterData.py` | KMeans clustering, distance matrix, genome→cluster assignment. **`assign_genomes_to_clusters_idv_year` uses the OPTIMAL one-to-one matching (Hungarian, `scipy.optimize.linear_sum_assignment`) as of 2026-09-09** — it was greedy in specifier-row order, which put 60–71% of sites off their nearest cluster (§7.9.2). |
| `SLiM_Code/CPBSampleSim{Linux,Win}.slim` | Forward sim. Neutral, `mutationRate(0)`. Takes `-d POPMULT` and `-d RECOMB` (§6.3), default simplification (§6.4). The two files are identical apart from path separators — **fix both or neither.** |
| `diagnostics/ridge_sweep.py` | §6.2.1 harness. `--setup` builds a `.trees` via the production path; each subsequent call runs one `4·Ne·μ = const` ridge point and appends JSON to `out/ridge_*.jsonl` (one point per process, so a failure costs only that point). Reports branch-mode diversity mean/sd/**CV**, site π, F_st, plus wall time and peak RSS. |
| `diagnostics/mu_calibrate.py` | §6.1.1 harness. Recapitates **once** (μ-free), then sweeps μ over cheap mutation overlays to solve `pi_loss` exactly (weighted median in log space, §7.0 mask applied). One POPMULT per process; appends JSON to `out/mu_calibration.jsonl`. `--skip-slim` reuses the `.trees` on disk. Also reports `fst_loss` at the calibrated μ, so both fitted statistics land in one run. **Since 2026-08-12 it also stores the per-subpop VECTORS** (`pi_sim`, `pi_obs`, `branch_div`, `deme_rel_size`, site names, both sample sizes) under `vectors`, which is what `pi_covary.py` consumes — summaries alone cannot answer a site-by-site question. |
| `diagnostics/pi_covary.py` | §7.2.1 harness. Reads the per-subpop vectors stored by `mu_calibrate.py` — **no simulation, no recapitation, runs in seconds.** Three tests: site-label permutation correlation, the flat-simulation `pi_loss` floor (a property of the *observed* vectors alone, so it applies at every POPMULT without re-running), and a shuffle null in the units of the objective. `--dump` prints the per-site table. |
| `diagnostics/fst_subsample.py` | §6.6 harness. Recapitate + simplify + mutate **once** (`--save-ts`/`--load-ts` make re-runs seconds), then R replicates drawing each deme's real n_i diploid **individuals** and recomputing π and F_st from scratch. Reports **bias** (mean over replicates − whole-deme) separately from **spread** (sd over replicates), in `fst_loss`/`pi_loss` units. Answers whether whole-deme sampling shifts the F_st *level*; π is reported only to size the noise, not to justify changing it. |
| `diagnostics/noise_floor.py` | §7.3 harness, and the end-to-end integration test. Holds every parameter fixed and re-runs the pipeline R times with different seeds, through the **real** `analyze_tree_sequence` + `calculate_losses`. Varies the three dice production rolls (SLiM `-s`, recapitation, mutation overlay); reuses `cluster_data.csv`/`migration_rates.csv` on disk, since KMeans is deterministic in production anyway. **Appends each replicate to `out/noise_floor.jsonl` as it completes** (a 3-rep run is ~1.6 h and OSPool evicts), and `--summarize` pools replicate records across processes — refusing to pool if they disagree on landscape/μ/Ne, or if seeds repeat. Guards the inputs on entry: deme count vs matrix dimension, and constant off-diagonal row sums. `--fixed-tree` gives a strict lower bound (recapitation+mutation only) for when SLiM is unavailable — good coverage for π, poor for F_st. |
| `diagnostics/mu_calibrate_summary.py` | Tabulates `out/mu_calibration.jsonl`: μ vs POPMULT, `branch_div` against its hard ceiling, the Monte-Carlo spread of the μ iterates, per-year sim-vs-obs log spread, and the saturation extrapolation. |
| `diagnostics/collect_batch.py` | §7.4 harness. Concatenates a CHTC batch's per-job CSVs into one results file and reports the landing checks. **No simulation — runs in seconds.** Dedupes the per-job header rows (a plain `cat` interleaves 500 of them) and **recovers `job_id` from the filename**, which is otherwise lost: `job_id` is not a CSV column and `iteration` restarts at 0 per job, so without it you cannot tell *which* jobs came up short. Reports: inventory (absent/short/over-long jobs, trials lost), whether failures track `pop` (the §3.1 memory signature), prior coverage by decile, per-statistic σ against the §7.3 noise floor, **rank-space variance decomposition onto the parameters** (§7.4.1 — this is what sets the weights), the loss gradient in `pop` with SE-aware flattening, and a suggested `WEIGHTS` block. `--no-write` reports without touching anything. **`EXPECTED_FIELDS` is the mechanical anti-pooling guard** — matched exactly, so a batch-1 file (12 cols) and a post-`ld_loss` file (13 cols) cannot be read by the same invocation; it names batch-1 files explicitly and exits cleanly when all files are rejected. `write_concat` refuses to overwrite a pooled CSV with a different layout, so forgetting `--out` cannot destroy batch 1 — **both guards fired for real on the pilot** (§7.6.0). **`LOSSES`/`FITTED`/`NOISE_FLOOR` gained `ld_loss` on 2026-09-08**; they had been missed when `EXPECTED_FIELDS` was updated, which silently dropped LD from every analysis section (§10.2). |
| `diagnostics/ld_probe.py` | §7.5.1 harness, and the gate on steps 4–6. Recapitate + simplify + mutate once (`--skip-slim` reuses the `.trees` on disk; `--save-ts`/`--load-ts` make re-runs seconds), then run the **real** `AnalyzeTreeSeq.calculate_ld_decay` per year and compare the pooled curve against `averaged_ldDecay_{year}.csv` over bins `bin_lo ≥ --min-bin` (default 562). One POPMULT per process, appending to `out/ld_probe.jsonl`; `--summarize` reports whether the fitted region **moves** with POPMULT against the sim-obs gap — that ratio is the decision number. **It writes per-deme curves to `out/ld_probe/ld_{year}_pop{P}_thin{T}_s{seed}.csv`, NOT to `data/Output_Data/ld_*.csv`** — reading the latter after an `ld_probe` run gets you whatever the last full-pipeline run left there, which cost a wrong conclusion on 2026-09-08. Note the filename does **not** encode `r`, so two runs at the same popmult/thin/seed and different `r` overwrite each other. The jsonl **stores the full per-bin curves**, so any cut can be rescored for free without re-simulating (§7.8.2). Hard-stops if `ld_common.spec_hash()` ≠ the spec the empirical run used, and asserts the LD column order equals specifier-matrix order (§7.5e). Also times the LD step alone, which is the per-trial cost §7.5c wants. |
| `diagnostics/condition_slices.py` | §7.4.3/§7.6.2 harness. **No simulation — runs in seconds.** Reads a pooled batch CSV, freezes sigma exactly as `abc_standardize.py` does, then ranks draws **within** narrow slices of `m` / `total_migration` / `numClusters` — which emulates a batch actually run at those pinned values. Reports posterior/prior IQR ratios for `pop` and `Nm` (`Nm` from the DEME size, `3.33*POPMULT/(33*numClusters)`). Its point is the `WEIGHTINGS` dict at the top: it scores the same batch under two fitted sets at once, which is what exposes a statistic that tightens `pop` by taking over the ranking rather than by agreeing with the others. Read the `pop` column DOWN — a median that swings with the slice is the answer changing with the assumption, not an error bar. |
| `diagnostics/ld_chrom_check.py` | §7.8.4 harness. **No simulation — seconds.** Runs three tests on `data/ld_per_chr/` (17 chr × 3 years): a leave-one-chromosome-out **jackknife** on the far-field excess (is the sim-obs gap significant?), a **heterogeneity** check against that SE, and — the discriminator — **decay distance vs far-field excess**. That last one separates the two candidate causes: low recombination rescales *distance* so it must slow the decay AND raise the far field together, while imputation adds a *level* and leaves decay alone. Measured +0.92 to +0.94, i.e. recombination. Also the tool that showed the empirical target is a ≥10× `r` mixture the single-locus model cannot be. |
| `diagnostics/qpost.py` | Post-process an existing `.trees` → branch diversity, site π, F_st. Fast; no SLiM needed. |
| `diagnostics/temporal_fc.py` | **§7.9.4/§7.9.6 harness — the temporal statistic.** Waples (1989) F_c over the **19 field-pairs resampled at identical coordinates** across years, matched on COORDINATES not names (the same field is `BrilowskiHome-2015` and `BrikalskiPats-2023`). Two arms: `n_real` (each timepoint cut to the field's real n_i — what is measurable) and `n_big` (both cut to `--big-n`, so the ~0.14 sampling pedestal shrinks and the N signal is readable). Also reports Hudson F_st on the same sample sets, which is what makes the §7.9.6 kernel test a single-pass comparison. **`--deme-choice` is NOT cosmetic** — the production site→deme mapping puts the same field in different demes in different years (§7.9.2), so it defaults to `nearest`. `--max-bytes` chunks the genotype extraction (§7.9.7); verified bit-exact against the unchunked path. Runs on an already-mutated tree in minutes. |
| `Python_Code/fc_common.py` | **The shared temporal-F_c spec, imported by BOTH sides (§7.9.8).** MIN_MAF, the year/generation map, the coordinate-matched field-pair list, the n<4 mask, `waples_fc`, and `spec_hash()`. Must be **COPIED** next to `ToUseOnBeagles/` exactly as `ld_common.py` is. The MAF cut here is an ASCERTAINMENT, not cosmetic -- at n=7 it conditions on the quantity forming F_c's denominator, which is tolerable only because both sides do it identically. Never reimplement any of this on one side. |
| `ToUseOnBeagles/CalcTemporalFc.py` | Empirical F_c from the Beagle VCFs -- **written 2026-09-09, NOT YET RUN.** Per-chromosome outputs are kept deliberately (`fc_out/chr*_temporalFc.csv`): they are the only route to a between-chromosome jackknife or a re-pool without re-reading every VCF, which is precisely what §7.8.4 needed for LD. Pools sums over total loci, never a mean of per-chromosome F_c (invariant 4). |
| `diagnostics/test_fc_roundtrip.py` | **Round-trip regression test for the empirical path. Seconds, no VCFs, no SLiM.** Writes a synthetic VCF + popfiles + specifier matrices, runs `CalcTemporalFc.main()` over them, and compares against F_c straight from the genotype matrix. **It caught the multi-allelic mismatch between the two sides on its first run** (§7.9.8C); tolerance is 1e-8 because the two paths accumulate in different orders. Run it after ANY change to `fc_common` or either side. |
| `Python_Code/scale_constants.py` | **The three scale constants: μ, r, `ancestral_Ne`. ONE home — import, never re-declare** (§7.9.3). They were literals in seven files and §6.8.1 recorded a planned `r` change being silently ignored in most of them. Q = 100 is a **declared** choice and `ancestral_Ne` is derived from it; `MU_TRUE`/`R_TRUE`/`PI_OBS` are the measured inputs, each with its citation. Also carries `branch_div_ceiling()` and `ld_time_depth_bp()`, both of which MOVE when a constant does. Dependency-free so anything can import it. |
| `Python_Code/recapitate_util.py` | **Recapitation past msprime's 100-populations-per-event limit** (§7.9.5) — the real reason `numClusters` was capped at 3. Staged merge: groups of ≤99 into throwaway intermediates, then into the ancestral population `np.nextafter` later. Exact, not approximate (intermediates live ~1e-13 generations at size 1.0). **Delegates to `pyslim.recapitate` unchanged at ≤99 populations**, so it is identity at the current prior. Verified at 200 demes: 127 s, every tree `num_roots == 1`. Used by `AnalyzeTreeSeq` and all six diagnostics that recapitate — `grep -rn "pyslim.recapitate" --include=*.py .` should only hit this file. |

**Every F_st in `diagnostics/` is Hudson as of 2026-08-26 (§6.7).** `grep -rn "\.Fst(" --include=*.py .`
must return **nothing** — a hit means something drifted back to Nei and its F_st is ~2× low. Treat
that grep as standing, like `.simplify(` in §2. Note the *recorded* numbers in `out/*.jsonl` and in
the §6.1.1/§6.2.1/§6.6 tables predate the fix and are still Nei-scale.
| `ToUseOnBeagles/*` | Empirical-side pipeline (§5). Runs on the Beagle machine, not here. |

Both hand-rolled estimators (`CalcGenRel.py`, `CalculateLD.py`) were **verified bit-identical to
tskit** (`genetic_relatedness`, `ld_matrix(stat="r2")`, ~1e-16/1e-18). **Keep it that way** —
re-run that verification if you touch them. Note this is an *estimator* check on identical input;
it says nothing about which set of sites each side is computed over in production, which is
exactly where the §5.1 denominator mismatch lived.

`run_code.sh` is the CHTC wrapper and is **not in the repo**. It clones from `origin/main`, so
**any code or empirical-target change must be committed and pushed before a CHTC submission.**

> **And the converse, learned the hard way on batch 2 (§7.6.0): never COMMIT a file the wrapper
> writes to.** `out/abc_results.csv` was tracked, so every job cloned batch 1's pooled results into
> its own output path and appended to them — 200 files each returning 2,495 rows of the wrong batch
> plus 3 real ones. Nothing in the job's own output looks wrong.
>
> **RESOLVED 2026-09-09 by a layout change, not just a gitignore.** Each job still writes to
> `../out/abc_results.csv` (`ABCAnalysisNoRedis.py:470,583,705` — unchanged, that is the job's own
> path). What moved is the **pooled** result: batch 1 now lives in **`out/batch1/`** and the pilot
> in **`out/batch2/`**, which are tracked and safe because the wrapper never writes there. The
> `out/` root is gitignored at `abc_results.csv`, `abc_results_ranked.csv` and `abc_sigmas.json`,
> and all raw per-job data (`out/*_raw/`) is gitignored too — `out/pilot_raw_asreturned/` alone is
> 144 MB. Constants updated in `abc_standardize.py`, `collect_batch.py --out` and
> `condition_slices.py --results`; both tools verified against the new paths.
> **Keep pooled results under `out/<batch>/` and never at the `out/` root.**


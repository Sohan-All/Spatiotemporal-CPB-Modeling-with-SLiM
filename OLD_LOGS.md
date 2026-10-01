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

Second archive pass, 2026-09-27 (CLAUDE.md cut from 3,809 lines to a session-pickup brief). Each
section below is the full text of the old CLAUDE.md section as it stood that day; where CLAUDE.md
now carries a short version, the long one is here:

- **§E — the prediction reframe** (old §1.1).
- **§F — pipeline, environment, per-run cost** (old §2, §3, §3.1).
- **§G — the empirical denominators** (old §5, §5.1–5.4).
- **§H — the scale constants and how they were reached** (old §6.1–6.2.1, §6.8, §6.8.1).
- **§I — ABC distance, kinship, noise floor, batch 1, LD** (old §7–§7.5.6).
- **§J — information budget, temporal F_c, batch 3, re-founding** (old §7.9–§7.9.12).
- **§K — conventions history** (old §10–§10.2).

References (old §11) moved to `REFERENCES.md`.

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


---

# §E — the prediction reframe, full text (old §1.1) [ARCHIVED 2026-09-27]

CLAUDE.md §1.1 keeps the decision and the acceptance criterion.

### 1.1 The goal is PREDICTION, and it reframes every open problem [DECIDED 2026-09-17, Sohan]

**The deliverable is a calibrated, validated model of how genetic material moves among Wisconsin
CPB populations across years, which can then project any allele's spread given a selection
coefficient. It is NOT a posterior on N.** Everything below follows from that, and several
conclusions reached under the old framing change sign.

**The end use is insecticide-resistance spread, and the model is neutral. That is a
FACTORIZATION, not a contradiction — but only if it is stated as one.** Spread speed of a
selected allele goes roughly as `σ·sqrt(2·s)`: the neutral data calibrate the dispersal term `σ`,
which cannot be had from first principles, and `s` is supplied separately. "Neutral model
predicting resistance" is indefensible; "neutral model calibrating the gene-flow term in a
resistance forecast" is standard practice. **Use the second phrasing.** No selection is added to
the simulation — SLiM is a pure neutral recorder (§2) and stays one.

> **Have this ready, because the professor will raise it.** Cohen et al. 2022 — §11, and very
> likely his own lab — is titled *"Evidence of hard selective sweeps suggests **independent**
> adaptation to insecticides"*: resistance arose SEPARATELY in WI and NY rather than spreading
> between them. If CPB resistance is dominated by repeated independent origins, a spread model
> may be predicting the wrong process. **The reframe that survives it:** a gene-flow calibration
> tells you **how large a region a single origin can supply**, which is exactly what decides
> whether to expect spread or independent origins. That is a better question than "how fast does
> it spread", and these data can speak to it.

#### What this DEMOTES

**The N/m identifiability failure (§7.4.2, §7.4.3) stops being the blocker.** A model can carry
unidentifiable parameters and still have well-identified predictions, provided the forecast
depends on them through the combination that *is* identified. What is identified here:

| quantity | status | what it governs |
|---|---|---|
| **Nm ≈ 78–103** | tight — IQR ratio 0.14–0.32 across every pinned slice (§7.4.3) | migrants per deme per generation, i.e. the gene-flow rate |
| **θ = 4Nμ** | pinned by the π calibration (§6.1.1, §7.9.3) | standing variation available to selection |
| **the temporal drift rate** | measured directly over 19 field pairs (§7.9.8) | how far a field's frequencies move per generation |

**`N` alone enters a resistance forecast mainly through `Ns` and through establishment
probability (≈ `2s`, which has no N in it). For a strongly selected allele — and resistance is —
`Ns >> 1` everywhere in the prior, so the dynamics are deterministic and insensitive to which N
was picked.** [INFERRED, textbook] So the ridge may simply not reach the deliverable.

**§7.9.1's information budget is also demoted.** It bounds how much the *genealogy* can tell you
about N, i.e. it constrains FITTING. It does not constrain forward simulation. It remains the
explanation for why every coalescent statistic failed; it is not a limit on forecasting.

#### What this PROMOTES — and both are now critical

**1. The temporal-dynamics misfit (§7.9.11D/H) goes from nuisance to disqualifying.** Under
parameter inference it merely pushed the posterior around. Under prediction it is an error in the
forecast quantity itself: simulated demes accumulate drift over 16 generations and real fields do
not. Run forward, that model predicts allele frequencies in individual fields wandering — and
alleles drifting to local loss or fixation — substantially faster than reality. For resistance
that is wrong in the most consequential direction. **You cannot extrapolate from a model that
misses the in-sample version of your own target.** §7.9.12's re-founding probe is therefore
clearly worth its ~2 h, and this is the reason.

> **THE PROBE RAN 2026-09-19/20 AND IT WORKS — §7.9.12F/G.** Re-founding each field annually from
> 40 kernel-drawn colonists gives drift growth **−0.00043 ± 0.00044** against an observed
> **−0.00151**, where a persistent deme gives **+0.0101 ± 0.0037** and §7.9.11H showed no prior
> value reaches the data. A bottleneck *without* outside founders makes it 2.5× worse, so the
> mechanism is the immigration, not the small size. **The disqualifying misfit now has a
> demonstrated fix.**
>
> **It costs the N-identifiability story, and that lands squarely on this section.** Under
> re-founding **neither fitted statistic reads deme size** — F_c's elasticity in N goes −0.996 →
> −0.08/+0.04 and F_st's −0.95 → −0.05 — both instead reading founder number `K` and the kernel.
> The table above lists **Nm** and **the temporal drift rate** as identified quantities; both were
> measured under the persistent-deme model and **neither carries over unexamined.** This section's
> argument survives in form — unidentifiable parameters, identified predictions — but its specific
> contents have to be re-derived under whichever structure is adopted. **The hold-out (§10.2) is
> still the test, and it matters more now, not less.**

**2. Spatial resolution now has to match the prediction scale, and it does not.** Previously
underweighted because it only blurred the fit (§7.9.2). It is worse than that for a forecast: the
model's demes are 6–10 km clusters and everything inside one is instantaneously mixed, while real
fields are **0.39 km apart at closest, 1.28 km median**, and CPB dispersal is sub-kilometre
(§11). **So the model as it stands cannot predict field-to-field spread at all — it assumes that
spread has already happened.** If the deliverable is spatial, deme count must rise toward field
resolution. The 99-deme cap is lifted (§7.9.5), so this is now a deliberate choice rather than a
constraint — and it is entangled with N, since deme size is `Average Count × POPMULT /
numSubpops`.

#### The data already give the MARGINAL forecast, with no model at all

The excess of F_c over its pedestal *is* "how much does a field's allele frequency change over t
generations", measured over 19 pairs: **≈0.0283 at t=8 and ≈0.0268 at t=16** (§7.9.11D). So a
field at frequency `p` today is, 8 generations out, a draw with variance ≈ `0.028·p(1−p)` — and
16 generations out, the same. That is a forecast, and it needed no simulation.

**What the model must therefore earn is the CONDITIONAL forecasts the data cannot measure:**
where an allele introduced at one field shows up and when, counterfactual management, horizons
beyond 16 generations, and the selection overlay. That is a real contribution — but it sharpens
the point above, because the model currently fails to reproduce the one predictive quantity the
data measure directly.

#### The validation: hold out 2023 [pair counts VERIFIED 2026-09-17 from `averaged_temporalFc.csv`]

Out-of-sample forecast accuracy replaces parameter recovery as the acceptance criterion — which
this project has never had, and whose absence is part of why every result has felt inconclusive.
"2023 falls inside the 95% predictive interval" is a criterion.

The 19 temporal pairs split:

| span | t | pairs |
|---|---|---|
| 2015 → 2019 | 8 | **2** (`Arlington`, `GarrisonNE × OkrayGarrisonNE`) |
| 2019 → 2023 | 8 | 8 |
| 2015 → 2023 | 16 | 9 |

**Fit on** the 2015 and 2019 π and F_st matrices plus those 2 early pairs; **hold out** the 2023
π and F_st and **17 of the 19 temporal pairs**. The held-out set is the big half, which is what
makes it a strong test. **The thin part, and say it plainly: only 2 pairs inform the temporal
fit.** The two full spatial years carry the demographic fit, and the test is whether that fit
predicts temporal behaviour it never saw.

> **And those 2 pairs are the least representative in the set [VERIFIED 2026-09-19 from
> `averaged_temporalFc.csv`].** They are `Arlington-2015 × Arlington-2019` at **n=4/7** — the
> smallest retained sample anywhere in the target, and the UW research station, i.e. the one site
> §7.9.12A finds is plausibly **not** rotated — and `GarrisonNE-2015 × OkrayGarrisonNE-2019`,
> which is a **rename of one site**, not two fields. So if re-founding is the mechanism, the
> hold-out's temporal fit rests on the one pair that may not be re-founded plus one relabel.
> **Note also the thinness is an artifact of the hold-out choice, not of the collection:** the
> data hold **10 pairs at t=8** and 9 at t=16; dropping 2023 deletes 17 of 19 because every other
> pair touches it.
>
> **A second, sharper caveat if re-founding turns out to be right.** Its whole signature is that
> F_c is nearly **flat in t** — so all 19 pairs would be measuring one level rather than a rate,
> and "fit 2, predict 17" becomes close to trivially satisfiable. In that world the temporal half
> of the hold-out tests almost nothing and **the real out-of-sample content is the 2023 π and
> F_st matrices.** Say which half of the hold-out carried the test; do not report the pair count
> as if all 17 were independent evidence.

Mechanically cheap — mask the 2023 columns through the existing `get_keep_mask` path and restrict
the F_c pairs. No new machinery.

**This is what makes the project defensible regardless of whether N and m ever separate:** *the
parameters are individually unidentifiable; the out-of-sample predictions are well constrained.*
That converts the largest apparent failure into the methodological result.

#### Where the old blockers land now

- **N vs m — carry it, do not invert it.** For spread you want the migration *fraction*
  (`total_migration` × kernel weight). F_c constrains N well (elasticity −0.996, §7.9.4); F_st
  constrains Nm; **dividing them is ill-conditioned — the gradients are 5.7° apart,
  `det(J) = 0.098` (§7.9.4)** — so the migration fraction stays poorly determined. Under a
  predictive goal the answer is a **sensitivity surface**, not an inversion, and §7.4.3 is
  already its prototype.
- **The dispersal kernel.** Under inference its lack of a biological referent at 6 km resolution
  was fatal (§7.4.2). Under prediction it is ordinary: **specify it from the dispersal ecology
  (§11), state it as an assumption, and report sensitivity.** That is routine in predictive
  ecology. It does NOT license reading `1/m` as a beetle dispersal distance — that rule stands.
- **`total_migration`'s ceiling** is still the remaining N confound for F_c (§7.9.10F) and still
  wants an external bound — but as a sensitivity axis, not as a blocker.
- **Selection is NOT added.** The factorization above is what avoids it.

#### What this does NOT license

- **It does not rescue a `pop` posterior.** §7.9.11G stands: no weights are installed and none
  should be reported. The reframe says the posterior may not be *needed*, not that it is now
  valid.
- **It does not excuse the drift misfit.** The opposite — see above.
- **It does not mean "fit whatever, the predictions will be fine."** That is an empirical claim,
  and the hold-out is exactly the test of it. **If 2023 falls outside the predictive interval,
  the reframe has failed and must be reported as having failed.**

---


---

# §F — pipeline, environment and per-run cost, full text (old §2, §3, §3.1) [ARCHIVED 2026-09-27]

CLAUDE.md §2–3 keep the pipeline, the env and the current CHTC sizing.

## 2. How the simulation pipeline works [VERIFIED]

`Main.py::main(num_clusters, migration_rates_modifier, population_modifier, total_migration=0.05,
mutation_rate, recombination_rate, ancestral_Ne=6700)`:

1. `GenerateClusterData.cluster_coordinates(..., random_state=KMEANS_SEED)` — KMeans over field
   coordinates from `data/final_data_for_modeling.csv`. Seed pinned at 42.
2. Writes `data/cluster_data.csv` and `data/cluster_distances.csv` — **overwriting in place.**
3. `GenerateSimulationParams.determine_migration_rates(distances, total_migration, scale, ...)`
   → `data/migration_rates.csv`.
4. `slim -d POPMULT=<pop> -d RECOMB=<r> SLiM_Code/CPBSampleSim{Win,Linux}.slim` →
   `out/simTreeSeq.trees`. Neither constant has a default inside the `.slim` files — an absent
   `-d` is a loud `undefined identifier` error, not a silent fallback to a different scale.
5. `AnalyzeTreeSeq.analyze_tree_sequence(..., ancestral_Ne=6700)` — recapitate → simplify →
   overlay mutations → write, per year, `diversities`, `divergences` (d_xy), `fst`, and
   `relatedness` matrices under `data/Output_Data/`. Pairwise stats use **single batched
   `indexes=` traversals** (verified bit-identical to per-pair loops, ~24× faster at k=8).

**Key architectural fact:** SLiM runs with `initializeMutationRate(0)`. It is a pure neutral
tree-sequence recorder. **All mutations are overlaid afterwards** by `msprime.sim_mutations()`,
*after* `simplify()`.

Consequences:
- Mutation rate has **zero** effect on SLiM's memory or runtime.
- Site-level π is **exactly linear in μ**: `π = μ × (branch-mode diversity)`. Branch-mode
  diversity (`ts.diversity(mode="branch")`) is `2·E[T_pair]` in generations and is the
  mutation-free view of the demography. Use it when diagnosing.

**Migration parameterization.** `total_migration` is the total immigration fraction per subpop per
generation (a real bounded rate); `scale` (the old `migration_rates_modifier`) is only the
dispersal-kernel decay, reshaping *where* migrants come from. The kernel is built over sources
(self excluded), normalized to 1, then scaled by `total_migration`, so each row's off-diagonals
sum to exactly `total_migration`; SLiM fills retention as `1 − total_migration`. [VERIFIED]

**`ts.simplify()` renumbers populations.** It defaults to `filter_populations=True`, which drops
unreferenced populations and renumbers survivors contiguously — while
`ts.samples(population=idx, ...)` uses the *original* cluster-row index. This produced the CHTC
"2-of-5 rows" bug (`Sample sets must contain at least one element`) at `numClusters=3` (×33 = 99
demes), where many demes go unreferenced. Fixed with `filter_populations=False`
(`AnalyzeTreeSeq.py:146`, commit `c5963ae`). **This was a correctness bug, not just a crash:**
when filtering occurred but the index stayed in range, the old code silently returned the *wrong
deme's* statistics. Treat any pre-`c5963ae` results as suspect. [VERIFIED]

> **The same bug survived in the diagnostics harness until 2026-08-04** — `qdriver.py:154` (since
> deleted) and `qpost.py:52` both simplified with the default `filter_populations=True` and then
> queried `ts.samples(population=i, ...)` with original cluster-row indices. `qpost.py` now passes
> `filter_populations=False`. Exposure was probably small at numClusters=33 (most demes are
> referenced), but **any §6.2-style sweep number that came from these two scripts rather than the
> full pipeline should be re-measured**, and this is a prerequisite for using them on the §6.1
> `--anc-ne` sweep. Grep for `.simplify(` before trusting any new script.

---

## 3. Environment [VERIFIED]

**`cpb-env` (from `environment.yml`, renamed from `environment3.yml` in `929f1c2`): SLiM 5.1,
pyslim 1.1.1, tskit 1.0.2, msprime 1.4.1,
python 3.12.** Internally consistent SLiM-5 stack; the pipeline runs end to end. Invoke via
`conda run -n cpb-env python ...` from `Python_Code/` (paths are relative to that dir).

> An older note claimed "SLiM 4.3 tree sequences need pyslim 1.0.x; pyslim 1.1+ fails." That
> described the *old* env. Do **not** downgrade pyslim on this env.

> **Windows Smart App Control blocks binaries INTERMITTENTLY, `slim.exe` included. Do not plan
> around it in either direction.**
> `HKLM:\SYSTEM\CurrentControlSet\Control\CI\Policy\VerifiedAndReputablePolicyState = 1` and has
> **never** been turned off. What changes is the *verdict*, which is reputation-based and therefore
> time-varying: `slim.exe` was refused 08-12 and 08-23, then ran normally on 08-24 minutes after a
> refusal, at an unchanged policy state. `sklearn.cluster` and `scipy.stats` behaved the same way
> and now import cleanly. **A refusal does not predict the next invocation.**
>
> **Probe at the moment you need it** — `slim -v`, `conda run -n cpb-env python -c "import
> sklearn.cluster"`. Two successive notes in this file were wrong in opposite directions because
> both generalised one probe into a standing fact. **Record what you observed and when, never what
> is "blocked".**
>
> What has *always* worked locally: everything starting from an existing `.trees` — recapitation,
> simplification, mutation overlay, every statistic. That is why `--skip-slim` / `--fixed-tree`
> exist in the diagnostics. A full noise floor (§7.3) needs the forward genealogy re-drawn and
> **that has run here** (~28 min/replicate), so CHTC is the fallback, not the only route. Turning
> SAC off is a one-way change via Windows Security → App & browser control; given the block is
> intermittent, not obviously worth doing.
>
> This is why `import Main` in `ABCAnalysisNoRedis.py` is **lazy, inside `model()`** — it keeps the
> loss/statistics half (`_read_vector`, `_read_matrix`, `get_keep_mask`) importable without the
> simulation stack, so diagnostics reuse the **real** mask and readers.

### 3.1 Per-run cost [VERIFIED 2026-07-31, measured post-§6.3/§6.4, numClusters=33]

All on one 15.2 GB Windows box. Peak memory is the OS high-water mark (`peak_wset`), tracked
separately for the SLiM child and the Python process — **they do not overlap**, SLiM has exited
before recapitation starts.

| POPMULT | SLiM time | SLiM peak | `.trees` | analysis time | **analysis peak** | total |
|---|---|---|---|---|---|---|
| 500 | 4.0 s | 173 MB | 21.1 MB | ~274 s | 611 MB | 4.5 min |
| 5000 | 50.9 s | 1672 MB | 216.5 MB | 1764 s | **7632 MB** | **30.3 min** |
| 12000 | 169.6 s | 3938 MB | 519.8 MB | — **OOM** — | ≈20 GB (est.) | — |

**Scaling, measured, not assumed:**
- **SLiM memory and `.trees` are linear in POPMULT** (9.7× and 2.36× against 10× and 2.4×). This
  is what confirms §6.4's fix works — with `simplificationRatio=INF` the edge table grew with
  generations×individuals instead.
- **Analysis-phase memory is SUPERLINEAR**, exponent ≈1.10 (12.5× for 10×). A linear
  extrapolation *under*-predicts: it gave 6.1 GB at POPMULT=5000 against 7632 MB measured.
- **Analysis time is SUBlinear** (6.4× for 10×), projecting **~60 min/trial at POPMULT=12000**.

**Recapitation is still the bottleneck** — ~97% of wall time, and it is where the memory goes.

> **POPMULT=12000 OOMs on a 16 GB machine.** Killed inside `pyslim.recapitate`, before the
> `Simplifying tree sequence...` print. This is **not** a §6.4 regression: the forward phase
> completed fine at 3.9 GB. It is the memory floor §6.4 says simplification cannot touch — gens
> 308/316 permanently Remember *every* individual in *every* subpop, so recapitation faces ~240k
> sample nodes.

> **The table above is at the OLD constants (μ 4.646e-7, r 2.75e-6, `ancestral_Ne` 6700).** The
> CHTC table below is **measured** at the current Q=100 constants and supersedes the 2026-08-26
> projection (which put 25000 at 43.5 GB and 1.93 h/trial).

**CHTC sizing [MEASURED 2026-09-15 from batch 3's HTCondor logs, `out/batch3_log/`, Q=100].** All
500 jobs terminated normally; no evictions, no holds. `MemoryUsage` is a **job-level** figure —
the maximum over that job's 5 trials, sampled periodically — so it is binned by the largest
POPMULT in the job.

| largest POPMULT in job | jobs | peak memory, median | peak memory, max |
|---|---|---|---|
| 5000–10000 | 3 | 9.4 GB | 10.0 GB |
| 10000–15000 | 24 | 14.4 GB | 18.5 GB |
| 15000–20000 | 110 | 21.9 GB | 28.9 GB |
| 20000–23000 | 173 | 27.6 GB | 36.9 GB |
| **23000–25000** | 190 | **32.2 GB** | **40.9 GB** |

- **Time per trial** (fit of job `TimeExecute` to Σ c·POPMULT^k over its 5 trials): ≈
  5.24e-4 · POPMULT^1.559 s → **5 min at 5000, 20 min at 12000, 44 min at 20000, 63 min at 25000**;
  27 min averaged over the prior. A 5-trial job's median is **1.8 h**. One job (413) took 8.6 h at
  ordinary POPMULTs — probably a slow node. [INFERRED]
- **Disk** peaked at 1.6 GB against 10 GB requested.
- **`request_memory = 48 GB` is enough for the next batch** — 17% over the highest logged peak.
  Because `MemoryUsage` is sampled, the true instantaneous peak can sit a little above the log;
  that is the headroom's job. Undersized memory holds jobs mid-pass and biases the pooled pass
  toward small draws, so do not cut further than that without a reason.

**Simplify-before-recapitate** was considered and **declined** (§9 `AnalyzeTreeSeq.py:126,146`
recapitates the full tree sequence and only then simplifies to the ~2.4k genomes actually
sampled). Two-thirds of that reasoning is now known to be wrong: it is *safe* with
`keep_input_roots=True` (pyslim documents exactly this), and the benefit is **not** local-only —
it is the difference between running and OOMing at POPMULT ≳ 9000, on CHTC nodes too. Left
unchanged deliberately (2026-07-31): nothing is broken, and CHTC can request more memory.
**Revisit if memory or throughput becomes binding.**

Total simulated N is **not** POPMULT. Subpop size is `Average Count × POPMULT / numSubpops`, so
`total N ≈ POPMULT × mean(Average Count) ≈ 3.33 × POPMULT`. The older OOM (POPMULT≈40000 on a
128 GB machine) predates the §6.4 fix; the forward phase is now linear and cheap, so if that
recurs it will be recapitation again, not the edge table.

---


---

# §G — the empirical denominators, full text (old §5 incl. §5.1–5.4) [ARCHIVED 2026-09-27]

CLAUDE.md §5 keeps the script list and the rules.

## 5. The empirical statistics pipeline (`ToUseOnBeagles/`)

Runs on the machine holding the Beagle VCFs and pixy output, **not** in this repo. All paths are
relative to that working directory. Its outputs are copied into `data/empiricalStats/`.

```
ConvertBeagleToVCF.py  Beagle -> VCF (one record per line)
PixyTheFiles.py        runs pixy per chromosome per year -> statsChr{i}_{year}/
CallableSites.py       per-chromosome callable-site denominators (shared config)
AverageData.py         pools pixy output -> averaged_{pi,dxy,fst}_{year}.csv
                       denominators are ANALYTIC, not read from pixy (5.1)
CalcGenRel.py          per-year relatedness from VCFs -> averaged_genRel_{year}.csv
CalculateLD.py         per-year, per-subpop LD decay -- log bins, PARALLEL across chromosomes
                       (one process per chr, all 3 years per VCF read; LD_WORKERS to tune).
                       Imports ld_common.py, which must be COPIED next to it. Renamed from
                       "CalculateLD,py" on 2026-09-06 -- the comma was a typo, never a .py
CalcTemporalFc.py      Waples temporal F_c over the 19 matched field pairs -> averaged_temporalFc.csv
                       + chr*_temporalFc.csv. Imports fc_common.py (COPIED next to it). Drops
                       fc_common.EXCLUDED_SAMPLES. Current spec 49afd0877028 (7.2.2H); the run
                       on file was made under it on 2026-09-12
CalcKinship.py         individual KING kinship within and between sites -> kinship_out/ (7.2.2).
                       Its raw values and printed reading guide are WRONG for this data until the
                       per-individual het-deficit correction is applied
```

One command runs the recalculation:

```bash
python AverageData.py && python CalcGenRel.py
```

### 5.1 π and d_xy are per-SNP unless corrected [VERIFIED]

`PixyTheFiles.py:28` runs pixy with `--bypass_invariant_check` on a **variant-sites-only** VCF, so
`count_comparisons` counts comparisons at SNPs only and `avg_pi` is per-SNP heterozygosity, not
per-site nucleotide diversity. Confirmed directly from pixy's own output: chr1 `no_sites =
6,549,657`, exactly the VCF record count, with zero invariant sites.

`AverageData.py` corrects this by extending the denominator over callable sites:

```
comparisons_per_site = C(2n, 2)        for pi     <- ANALYTIC, from sample size
                     = (2n_i)(2n_k)    for d_xy
denominator          = comparisons_per_site * callable_sites
```

The numerator needs no adjustment — invariant sites contribute zero differences.

**The missingness assumption is exactly satisfied.** [VERIFIED] `count_missing = 0` and
`count_comparisons / no_sites` divides to exact integers: 91 = C(14,2) → 7 diploid individuals,
231 = C(22,2) → 11. Beagle imputes everything, so comparisons-per-site is constant and the
extrapolation is exact, not approximate.

That last fact is *why* the denominator is computed from sample size rather than read out of pixy's
`count_comparisons`, which overflows int32 for large subpops (§6.5). `AverageData.py` still reads
the field to cross-check and **raises** if it disagrees where pixy could have been right — do not
soften that guard into a warning.

**Callable sites are provisional.** True callable counts are unrecoverable — the Beagle files hold
only variant sites, so the upstream filtering was never recorded. `CallableSites.py` uses pixy's
`window_pos_2` (position of the last SNP per chromosome), a *lower* bound on length and therefore
an *upper* bound on π. Assembly chromosome lengths would give the other bound; the gap is order
10–20%, against the ~2400× scale error being corrected, so no downstream conclusion turns on it.

Genome-wide: **81,141,632 SNPs over 930,522,190 bp = 8.72% density.** Corrected mean π ≈ **0.0122**
(from 0.1403 uncorrected) — an ordinary insect value.

### 5.2 Pool across chromosomes, never average [VERIFIED]

These statistics are ratios, so genome-wide = `Σ numerators / Σ denominators`. Averaging 17
per-chromosome ratios over-weights sparse chromosomes. `AverageData.py` pools; `CalcGenRel.py`
always did.

**F_st is the exception, by necessity.** Proper pooling (Bhatia et al. 2013) is `Σa / Σ(a+b)` over
the Weir–Cockerham variance components, but pixy's `fst.txt` emits only `avg_wc_fst` and `no_snps`
— the components are not in the file. `AverageData.py` uses a SNP-count-weighted mean, strictly
better than a plain mean. To pool properly, recompute from the VCFs with
`allel.weir_cockerham_fst`, which returns the components. Error is a few percent, and F_st is
scale-invariant anyway.

### 5.3 Which statistics the denominator affects

**The dividing line is whether the statistic is a ratio.** [VERIFIED by reasoning]

| statistic | fitted? | affected by per-SNP denominator? | why |
|---|---|---|---|
| π | **yes** | **yes** | level; scales with diversity |
| F_st | **yes** | no | ratio of variance components; scaling cancels |
| IBD slope | **yes** | no | regression built on F_st |
| d_xy | diagnostic | **yes** | level; same denominator as π |
| genetic relatedness | diagnostic | **yes** | sum of centred products; not μ-invariant |

Invariant sites contribute zero to both numerator and denominator of F_st, so computing it over
variant sites only is correct. **The structure/migration side of the inference was never
compromised by the scale issue.**

Relatedness is a *level*: tskit applies `span_normalise=True`, making the simulated side per base
pair over the full 1e6 bp. `CalcGenRel.py` therefore normalizes by callable sites
(`NORMALISE_BY = "callable"`); `"segregating"` preserves the old behaviour for the tskit estimator
check only.

### 5.4 Cosmetic quirks that are not bugs [VERIFIED]

- **`ConvertBeagleToVCF.py:33` hardcodes `CHROM = 9`** for every file. Chromosome identity lives
  only in filenames. Harmless — every script keys off the filename consistently — but it means the
  VCF contents cannot tell you which assembly chromosome a file is.
- **`UserWarning: 'GT' FORMAT header not found`** from scikit-allel. `generate_VCF_header`
  (`ConvertBeagleToVCF.py:82-85`) omits the `##FORMAT=<ID=GT,...>` declaration, so allel falls
  back to its default GT spec (diploid), which matches the data. If allel had actually failed to
  parse GT, `chunk[0]["calldata/GT"]` would raise `KeyError`, not warn.
- **chr6 has the lowest SNP density** (6.26% vs 8.72% genome-wide). All chromosomes were processed
  identically, so there is no differential artifact to correct. Accepted as-is.

---


---

# §H — the scale constants and how they were reached (old §6.1, 6.1.1–6.1.3, 6.2, 6.2.1, 6.3–6.7 stub, 6.8, 6.8.1) [ARCHIVED 2026-09-27]

CLAUDE.md §6 states the current constants and the traps. Several values below are SUPERSEDED (ancestral_Ne 6700, mu 4.646e-7, r 2.75e-6): the Q=100 reset of 2026-09-09 replaced them.

## 6. Known defects — read before changing the ABC

Ordered by how much they distort the inference.

### 6.1 `ancestral_Ne = 6700` — provenance found; the value fails its own source's internal check

**Source [VERIFIED 2026-08-04]: Cohen et al. 2022, *Evolutionary Applications* 15:1691–1705
(doi:10.1111/eva.13498), Figure 3a.** PDF is in the repo root. 6700 is `N_a`, the **dadi**-inferred
ancestral effective size of the common ancestor of the Hancock, WI and Long Island, NY pest
populations, at a split **325 generations (160 yr ± 1.6)** ago. The same figure gives
`N_WI = 15,000 (±163)` and `N_NY = 40,000 (±760)`.

**This one paper is the provenance of nearly every fixed constant in the project:**

| project constant | source in Cohen et al. |
|---|---|
| `ancestral_Ne = 6700` | Fig. 3a `N_a = 6,700 (±10)` |
| `DEFAULT_RECOMBINATION_RATE = 2.75e-6` | §3.2 `r_HAN = 2.75e-6` — the **Wisconsin** population, i.e. the right one for us |
| SLiM run length 324 generations (`CPBSampleSim*.slim:43,48,53`) | the 325-generation divergence |
| "biological μ ≈ 2.1e-9" | the *Chironomus riparius* midge rate Cohen used — **SUPERSEDED 2026-09-08, see §6.1.2: a measured CPB rate now exists** |
| `pop` prior (was `U(2000, 12000)` → total N 6.7k–40k) | bracketed `N_WI = 15,000` and `N_NY = 40,000` — **superseded: now `U(2000, 25000)` → 6.7k–83k, §6.7** |

**The old guess in this section was half wrong.** 6700 is *not* a contemporary or LD-based
estimate. In dadi's `no_mig` model the ancestral population is **constant-size extending
infinitely into the past**, which is conceptually exactly what `recapitate(ancestral_Ne=)` wants.
**The role is right.** The value is not.

**The ~217× conflict is internal to Cohen et al., not a disagreement between their data and ours.
[VERIFIED — arithmetic]** Watterson's θ from *their own reported* 11.8M polymorphic sites over an
~870 Mb genome at n=28 diploids (a_56 = 4.594):

```
theta_W = 11.8e6 / (870e6 * 4.594)   = 2.95e-3 per site
Ne      = theta_W / (4 * 2.1e-9)     ~ 3.5e5      <- from Cohen's own SNP count
```

That is **52× their own reported `N_a` = 6700**, using their own mutation rate. Our π-derived
`Ne = 0.0122/(4·2.1e-9) ≈ 1.45e6` is only **4×** from *that* — same order of magnitude.
**Our π is not the outlier; 6700 is.**

**The paper says so itself, twice.** §4.1: sequencing was low-coverage (>5× average, ≥3× per
individual), so "possible heterozygous sites [are] mistaken as homozygous… reducing singletons and
causing a bias that results in **underestimating demography**… our results might lead to an
**underestimate in effective population size**." And their two methods disagree with each other —
"the dadi estimates… were ~4-fold larger than estimates from the Stairway plot," whose Figure 1
sits at 1k–20k. The authors treat this as a modest caveat; the arithmetic above says it is ~50×.

**Leading mechanistic candidate — and it is the mirror image of the §5.1 bug we just fixed on our
own side. [INFERRED]** `N_a` is derived, not measured: `N_a = θ_dadi / (4μL)` with
**L = 840 Mb, "all intergenic sequence data."** But the 2D-SFS was built from intergenic regions
"with **stringent quality thresholds for coverage and likelihood**." If θ was fit to an SFS drawn
from a heavily filtered subset while L was set to the full intergenic span, `N_a` is deflated by
exactly that fraction — the same denominator/numerator mismatch as §5.1, pointing the other way.
Reconciling with our π would need L ≈ 3.9 Mb, which is too small to be the whole story on its own,
so this is likely **compounded with** the low-coverage singleton loss the authors name. Confirming
it would need their supplement. Either way both named biases push the same direction: *up*.

How μ has been absorbing the error:

| scenario | π target | required μ at Ne=6700 | × biological |
|---|---|---|---|
| both errors | 0.140 | 5.2e-6 | ~2400× |
| denominator fixed (now) | 0.0122 | 4.6e-7 | ~217× |
| both fixed (Ne ≈ 1.5e6) | 0.0122 | 2.1e-9 | 1× |

`ancestral_Ne` is a threaded parameter, **not** an ABC free parameter — it is confounded with μ in
π (they enter only as `4·Ne·μ`), so inferring it would build a second ridge. Exposed for
sensitivity analysis only. Scaling it with `population_modifier` was considered and **rejected**
(fabricates N-identifiability, conflates two demographic epochs).

**Resolution [VERIFIED 2026-08-04 by the ridge sweep, §6.2.1 — this REVERSES the earlier
recommendation to set `ancestral_Ne ≈ 1.4e6`]:**

**Keep `ancestral_Ne = 6700`. Do not raise it.** Two measured facts force this:

1. **It would not change the inference.** π is invariant along the `4·Ne·μ = const` ridge in
   *both* its level (±2.6%) and its between-subpop relative spread (CV flat to <1%). So moving
   along the ridge is very nearly a no-op for the ABC.
2. **It is computationally impossible.** Recapitation cost scales as **Ne^2.34** (measured).
   Ne = 1.452e6 extrapolates to **~600 days per trial.**

**What must change is the language, not the constant.** μ = 5e-6 is **not a mutation rate** — it
is half of a calibration constant whose only meaningful content is the product `4·Ne_anc·μ`, set
to match observed π. Never report μ on its own, never present it as biological, and do not infer
it: report **θ = 4Nμ** (§6.2 already said this). The critique of Cohen's 6700 above stands as a
statement about CPB biology; it just does not license a code change, because the biologically
"correct" value cannot be simulated.

> **Calibration is POPMULT-dependent.** Setting `4·Ne·μ = 0.0122` does *not* yield π = 0.0122 —
> forward-phase coalescence pulls `branch_div` below `4·Ne_anc`, so the realised π is lower, and
> by a POPMULT-dependent factor (81% of target at POPMULT=500, ~50% at POPMULT=150). μ and POPMULT
> are therefore mildly coupled through the π level. Calibrate μ at the POPMULT you intend to run.
> **Done — see §6.1.1.**

### 6.1.2 A CPB-SPECIFIC mutation rate now exists — 5.8e-9 [FOUND 2026-09-08]

**Xu et al. (2026), *Genome Biology and Evolution* 18(2):evag027** — trio sequencing,
**16 parent-offspring trios** from two families, ~32.8× coverage, 92 de novo mutations,
~491.8 Mb callable per trio, all 14 Sanger-checked variants confirmed.

```
mu = 5.8e-9 per site per generation   (95% CI 4.7e-9 - 7.2e-9)
```

**This retires the midge stand-in.** §6.1's table above says "there is no CPB-specific rate" and
the project has been quoting *Chironomus riparius*' 2.1e-9. The measured CPB rate is **2.76×
higher**, and the paper notes it is ~2× the median for insects generally.

**It does NOT change `DEFAULT_MUTATION_RATE`.** §6.1's resolution stands: μ in this model is half
of a calibration constant, set so `4·Ne_anc·μ = π_obs`. What the trio rate changes is the
*interpretation* — see §6.1.3.

**Two independent corroborations fall out of it.** The paper reports **π = 0.005** against our
corrected 0.0122 (§5.1) — same order, where the pre-fix per-SNP value was 0.1403, so this is
outside confirmation that the denominator fix was right. And their own implied ancestral
`N_e ≈ 2.4e5` sits with the other estimates in §6.1.3.

> **Caveat, and it matters for POPMULT.** The same paper puts *present-day* US `N_e` at around
> 1e4, which is a different quantity from the ancestral figure and ~25× smaller. **5.26e5 is the
> right comparison for `ancestral_Ne` and the recapitation phase; it is NOT a target for POPMULT**,
> which represents contemporary deme sizes.

### 6.1.3 What the model's scale actually is: Q ≈ 78.5, and only `r` is out of step [VERIFIED 2026-09-08]

§6.5 records that coalescent rescaling "cannot be applied as the model stands". That was written
when `m ≈ 0.76`, before the `total_migration` refactor, and it is now out of date in a more
interesting way than expected: **two-thirds of the rescaling had already happened by accident.**

Standard rescaling simulates `N/Q` and sets `μ→μQ`, `r→rQ`, `m→mQ`, `G→G/Q`. Asking what `Q` each
knob implies, with the trio μ and the §6.8 linkage-map `r`:

| knob | implied Q |
|---|---|
| ancestral Ne (`π/(4μ_true)` = 5.26e5, vs simulated 6700) | **78.5** |
| μ (4.646e-7 / 5.8e-9) | **80.1** |
| `r` (2.75e-6 / 1.02e-8) | **269.6** |
| G (325 / 324) | **1.00 — time was never compressed** |

**The μ calibration IS a correct rescaling of the ancestral phase, by algebra rather than by
luck.** μ was calibrated so `4·Ne_anc_sim·μ_sim = π_obs`, and `Ne_anc_true` is *defined* as
`π_obs/(4·μ_true)`. Substituting gives `μ_sim = μ_true·Q` identically. Measured agreement: 4.55e-7
predicted against 4.646e-7 calibrated, **−2.0%**.

**So `r` is the one knob out of step, and by 3.4×, not 269×.** At Q = 78.5 it should be
**8.0e-7**. Setting it to the real-world 1.02e-8 would be **78× too small for this model** — that
is the unrescaled value, correct for a beetle and wrong for a simulation running at 1/78 scale.

**Three independent routes now agree on the true ancestral `N_e`, and the LD one is genuinely
independent of the other two:**

| route | `N_e` |
|---|---|
| our π = 0.0122 with the trio μ | 5.3e5 |
| the trio paper's own π = 0.005 with the same μ | 2.4e5 |
| **LD**: `ρ_obs ≈ 0.017/bp` (§7.5.1) with the linkage-map `r` | **4.2e5** |
| Watterson from Cohen's own SNP count (§6.1) | 3.5e5 |

All within ~2× of each other and all ~50× above Cohen's 6700, which is §6.1's conclusion reached
from four directions instead of one.

**What rescaling still CANNOT fix: time.** A consistent Q = 78.5 needs `G → 325/78.5 ≈ 4.1`
generations. The model runs **324**, anchored to the real WI/NY split, with gens 308/316/324 being
the 2015/2019/2023 collections at 2 generations/year. Compressing would put those three sampling
points 0.1 generations apart. So the forward phase accumulates far more drift than reality —
`G/N` per deme runs **0.13 to 4.8** across the prior against **~0.02** for a real metapopulation of
`N_e ≈ 5e5` over 33 demes, i.e. **6× to 200× too much**. That is structural, not a parameter, and
it is a candidate explanation for residuals the other knobs cannot absorb.

> **[OPEN] The one idea worth costing.** The LD cut is `1/(2·G·r)`, so it moves with G, not with
> `ancestral_Ne`. Running the forward phase longer would push forward control down toward the
> ~60 bp decay where the real LD signal is (§7.5). It needs ~10,400 generations against 324 —
> 32× — and SLiM's forward phase is cheap while recapitation would get *cheaper*. The obstacle is
> biological: 325 generations is the invasion timescale, so a much longer structured phase needs
> justifying. **A professor question, not a switch to flip.**
>
> **And scaling `ancestral_Ne` with POPMULT is NOT the alternative** — asked and rejected
> 2026-09-08, for the reasons §6.1 already gives plus a decisive one: recapitation costs
> `Ne^2.34` (§6.2.1), so tying it to POPMULT makes the prior's upper half unrunnable (39× cost
> from POPMULT 5000 → 25000). It would also wire the constant to the unknown, manufacturing a
> sharp N posterior out of an assumption.

### 6.1.1 μ recalibrated: 5e-6 → 4.646e-7 [VERIFIED 2026-08-11]

`diagnostics/mu_calibrate.py`, `numClusters=33`, `total_migration=0.05`, `ancestral_Ne=6700`,
seed 1. Raw records in `out/mu_calibration.jsonl`; `mu_calibrate_summary.py` tabulates them.

**Method.** Recapitate + simplify **once** (~97% of cost, and entirely μ-free), then sweep μ over
cheap mutation overlays on the ~2.4k-sample simplified tree. Branch-mode diversity `b_i` gives an
analytic first μ with no re-running; the loop then corrects the multiple-hit deficit. The
objective is the *actual* fitted loss, not mean-matching: since `log π_sim,i = log μ + log b_i`,
minimising `pi_loss` over `log μ` is a weighted-L1 problem whose exact minimiser is the **weighted
median** of `log π_obs,i − log b_i` with weights `1/(3·n_year)` — the same year-normalisation
`calculate_losses` uses, with the §7.0 mask applied to both sides.

| POPMULT | subpop N | `branch_div` | /ceiling | **μ_calib** | 4Ne·μ | `pi_loss` | `fst_loss` | sim F_st |
|---|---|---|---|---|---|---|---|---|
| 500 | 50 | 21647 | 0.789 | 5.564e-7 | 0.01491 | 0.0672 | 0.07145 | 0.0766 |
| 2000 | 202 | 25629 | 0.934 | 4.819e-7 | 0.01291 | 0.0290 | 0.01856 | 0.0197 |
| **5000** | 505 | 26847 | 0.978 | **4.646e-7** | 0.01245 | 0.0213 | **0.00830** | 0.0079 |


> **Nei-scale banner (§6.7): every `fst_loss` / `sim F_st` figure in this section predates the
> 2026-08-26 Hudson fix and is ~2× low. Orderings hold, LEVELS do not — never compare one
> against a post-fix number. Not re-measured, by decision. `pi_loss` and `branch_div` are
> unaffected.**

(observed F_st = 0.00645 in every row.) **`DEFAULT_MUTATION_RATE` is now the POPMULT=5000 value,
4.646e-7 — the old 5e-6 was 10.8× too large**, being calibrated against the pre-2026-07-28
per-SNP π target (§5.1).

**Three independent cross-checks passed**, which is what licenses trusting the harness: implied
`branch_div` at POPMULT=500 reproduces `ridge_sweep.jsonl`'s 21985; `fst_loss` at POPMULT=500
reproduces §6.2's 0.0719; and `fst_loss` at POPMULT=5000 reproduces §6.2's 0.00829 **exactly**.

**`site_pi = μ · branch_div` is not exact — it runs 1–2% low**, and the deficit grows with `μ·b`
(multiple hits at a site). Purely analytic calibration lands ~1.8% low, so the iteration matters.
Monte-Carlo precision of the result is **0.44–0.99%** (each iterate re-draws mutations, so the
loop oscillates rather than converging to a fixed point; `mu_calibrated` is the last iterate).

**`branch_div` saturates, which is why ONE fixed μ covers the whole prior.** Recapitation
coalesces any surviving pair at rate `1/(2·Ne_anc)`, so
`branch_div ≤ 2·(324 + 2·6700) = 27448` — a hard ceiling, giving μ a floor of ~4.45e-7. Fitting
the deficit `1 − b/ceiling ≈ 94.2·POPMULT^−0.973` (a **fit, not a derivation** — an `exp(−P/τ)`
form was tried and rejected for being unable to fit both ends) extrapolates to μ = 4.51e-7 at
POPMULT=8000 and 4.49e-7 at 12000. Holding μ at 4.646e-7 across the prior `U(2000, 12000)`
therefore drifts π by only **−3.2% to +2.1%** — and that drift is *signal*: it is how the π level
carries POPMULT information now that μ is effectively pinned.

> **Still true after the ceiling moved to 25000 [VERIFIED 2026-08-26 by the same saturation fit].**
> `branch_div` at POPMULT=25000 is 27312 against 26797 at 5000 — **+1.9%**, versus +1.4% at 12000.
> Across the widened prior `U(2000, 25000)` the drift is **−3.5% to +1.9%**, barely worse than the
> −3.2%/+2.1% over the old one, because `branch_div` is already at 99.5% of its hard ceiling.
> **μ does NOT need recalibrating for the wider prior** — saturation is exactly what buys this.

**Consequence for §6.2 — the blocker is gone.** With μ calibrated per POPMULT, π and F_st no
longer pull in opposite directions; both improve monotonically together (π 0.067→0.029→0.021,
F_st 0.071→0.019→0.0083). The conflict was an artifact of the miscalibrated μ, not a feature of
the data.

### 6.2 Is N identifiable? — must be re-derived [OPEN]

The earlier conclusion was that population size is nearly unidentifiable from π: at
`ancestral_Ne = 6700` and μ = 5e-6, recapitation alone predicts `4·Ne·μ = 5.63e-5` against a
pipeline output of `6.38e-5`, i.e. **~7/8 of simulated diversity was set by a constant that
`population_modifier` does not touch.**

**That arithmetic is conditional on values now known to be wrong.** If ancestral Ne rises ~217×
and μ drops to biological, the forward/ancestral balance shifts and N may become identifiable.
**Re-derive, do not assume.** Use branch-mode diversity to diagnose. A POPMULT sweep
is the way to settle it; a previous attempt OOM'd at POPMULT=40000. (`qdriver.py`, the old
harness for this, was deleted in `929f1c2` — `ridge_sweep.py`/`mu_calibrate.py` supersede it.)

Regardless: π depends on N and μ only through **θ = 4Nμ** — they are confounded. If both stay
free, the posterior is a *ridge*, not a peak, and a peaked-looking N marginal is an artefact.
Report **θ=4Nμ** and **Nm**.

**Partial answer [VERIFIED 2026-07-31]: N *is* identifiable — through F_st, not through π.**
A two-point POPMULT sweep at μ=5e-6, `total_migration=0.05`, numClusters=33:

> **[NARROWED 2026-09-06 by §7.4.2 — read this before quoting the claim below.]** The sweep that
> follows holds the **dispersal kernel and `total_migration` FIXED** and varies only POPMULT. Under
> that condition N is identifiable through F_st, exactly as stated. **It is not identifiable in the
> actual pass**, where `m` is free over four orders of magnitude: batch 1 returned a `pop` posterior
> whose IQR is 82–102% of the prior's at every acceptance level, while `m` collapsed 20×. F_st
> constrains **Nm**; this sweep separates N from m only because it fixed m by construction. The
> conclusion is conditional on that, and the condition does not hold in production.

| POPMULT | subpop size | Nm | sim F_st | **fst_loss** | sim π | **pi_loss** |
|---|---|---|---|---|---|---|
| 500 | ~50 | 2.5 | 0.0754–0.0778 | 0.0719 | 9.5e−2 | 2.043 |
| 5000 | ~505 | 25 | 0.0076–0.0082 | **0.00829** | 1.17e−1 | **2.262** |
| *observed* | — | 31–83 | 0.0032–0.0083 | — | 1.22e−2 | — |


> **Nei-scale banner (§6.7): every `fst_loss` / `sim F_st` figure in this section predates the
> 2026-08-26 Hudson fix and is ~2× low. Orderings hold, LEVELS do not — never compare one
> against a post-fix number. Not re-measured, by decision. `pi_loss` and `branch_div` are
> unaffected.**

F_st tracks `1/(1+4Nm)` almost exactly and `fst_loss` improves **8.7×**. So the structural side
carries real information about N. Do not conclude "N is unidentifiable" from the π argument alone.

**~~But the two fitted statistics currently pull in opposite directions.~~ RESOLVED 2026-08-11 by
§6.1.1.** The old reading was: the POPMULT that fits F_st drives π *further* from target (9.7× too
high at POPMULT=5000), so **no POPMULT satisfies both**. That was true only at the miscalibrated
μ = 5e-6. With μ recalibrated to 4.646e-7, π and F_st improve **together** with POPMULT
(π 0.067→0.029→0.021, F_st 0.071→0.019→0.0083 over POPMULT 500→2000→5000). The tradeoff was an
artifact of μ, not a feature of the data, exactly as this section suspected.

This section also **predicted the recalibration, and the arithmetic route won**: matching π at
POPMULT=5000 and Ne=6700 was estimated at μ ≈ **5.2e-7** by sweep extrapolation, against
**4.6e-7** from §6.1's π arithmetic. The measured value is **4.646e-7** — the arithmetic route was
right to within 1%, the extrapolation 12% high. Both remain ~220× the biological rate, which is
the §6.1 problem restated and is *not* resolved by this (nor can it be — §6.2.1).

### 6.2.1 The ridge sweep — π survives, but Ne_anc cannot be raised [VERIFIED 2026-08-04]

`diagnostics/ridge_sweep.py`, POPMULT=150, numClusters=33, seed 1, two points on
`4·Ne·μ = 0.0122`: **(Ne=6700, μ=4.55e-7)** and **(Ne=20000, μ=1.53e-7)**, a 2.985× step.

| quantity | Ne=6700 | Ne=20000 | ratio |
|---|---|---|---|
| `recap_s` | 183.8 | 2369.7 | **×12.89 → exponent 2.34** |
| `recap_peak_mb` | 292 | 291 | ×1.00 |
| site π (2015/19/23) | .00582/.00616/.00631 | .00570/.00600/.00615 | **−2.1/−2.6/−2.5%** |
| `branch_div` mean | 12918/13649/13992 | 37665/39635/40616 | ×2.92/2.90/2.90 |
| `branch_div` **sd** | 4010/4286/3816 | 11687/12508/11170 | ×2.92/2.92/2.93 |
| `branch_div` **CV** | .3104/.3140/.2727 | .3103/.3156/.2750 | **×1.000/1.005/1.008** |
| mean F_st | .2503/.2277/.2045 | .2518/.2291/.2063 | +0.6/+0.6/+0.9% |

**Test 1 — π is invariant along the ridge. PASSES.** π moves only −2.1 to −2.6% for a 3× change in
Ne_anc. The small residual drift is the forward-phase contribution shrinking in relative terms,
exactly as theory says it should.

**Test 2 — the between-subpop spread does NOT collapse. The prediction above was WRONG.** Mean and
SD of `branch_div` both scale ×2.9 with Ne_anc, so the **CV is flat to under 1%.** The mechanism:
the ancestral contribution to a subpop's coalescence time is weighted by the probability its pairs
did *not* already coalesce in the forward phase, and that probability varies by subpop. So the
ancestral phase **multiplies** the forward structure rather than **adding** a constant to it, and
relative structure is preserved no matter how deep the ancestral phase gets.

Consequences, all favourable:
- **π stays a fitted statistic.** It keeps its full relative information about POPMULT and
  migration. Do not demote it.
- **Log space is exactly right** (§7): `pi_loss` measures relative differences, which is precisely
  the quantity shown to be ridge-invariant.
- **The circularity worry is defused.** The element-wise variation is genuine independent signal
  and it does not shrink, so calibrating the *level* does not hollow out the statistic.
- **F_st is confirmed insensitive to `ancestral_Ne`** (<1%), as a forward-phase ratio should be.

**Test 3 — cost. This is the new blocker.** `recap_s` scales as **Ne^2.34** while memory stays
flat at ~291 MB. So the constraint is **wall time, not RAM** — the opposite of the §3.1 OOM
problem. Extrapolated from the Ne=6700 baseline at POPMULT=150:

| target Ne_anc | factor | projected recapitation |
|---|---|---|
| 2e5 | ×2,805 | ~143 h (6 days) |
| 1.452e6 (π-implied) | ×288,691 | **~14,700 h (614 days)** |

**Raising `ancestral_Ne` to the biologically-implied value is computationally out of reach**, by
about four orders of magnitude, and no amount of CHTC memory fixes a wall-time wall. Since §6.2.1
also shows raising it would barely move the inference, the correct move is to **keep 6700** — see
§6.1's resolution.

**Caveats, stated honestly.** Measured across 3× in Ne (6700→20000), not the full 217×; the
mechanism (mean and sd both ∝ Ne_anc) is clear and the extrapolation is principled, but the top of
the ridge is unverified. POPMULT=150 is far below the prior range, chosen because cost is
dominated by Ne_anc — dropping POPMULT 3.3× (500→150) cut recapitation only 19%, so POPMULT is
**not** a useful cost lever. The absolute CV is strongly POPMULT-dependent (0.27–0.31 at
POPMULT=150 vs 0.081–0.093 at POPMULT=500); what was shown invariant is its *insensitivity to
Ne_anc*, not its value.

### 6.3-6.7 Resolved defects — the standing rules only [full write-ups in `OLD_LOGS.md` §B]

Six defects, all fixed and verified in code. What must not be undone:

- **`r` reaches SLiM.** `CPBSampleSim*.slim` takes `-d RECOMB`; before 2026-07-29 forward
  recombination was hardcoded at 1e-8 and the ABC parameter reached only `pyslim.recapitate()`.
  An absent `-d` is a loud `undefined identifier` error, not a silent fallback.
- **`simplificationRatio=INF` is gone.** SLiM simplifies at its default ratio, so the edge table no
  longer grows unbounded across all 324 generations. Simplification is lossless for retained
  samples; SLiM's own runtime simplification uses `keep_input_roots`, which is why its output is
  routinely recapitable. Do not confuse that with a *Python-side* `ts.simplify()` before
  recapitating, which must pass `keep_input_roots=True` (§3).
- **pixy's `count_comparisons` overflowed int32** at ≥14 diploid individuals, corrupting 244 of
  2015's rows — only 18 of them visibly. Denominators are now analytic (§5.1). **Any
  `empiricalStats` output produced before 2026-07-28 is suspect**, and any future year with ≥14
  sampled individuals re-triggers it.
- **Migration saturation** (row-normalising `exp(-d·modifier)` to 1 made ~76% of each subpop's
  offspring immigrants — effectively panmixia), **KMeans re-randomization** (seed now pinned at
  42), **`csv.DictReader` eating subpop 0** (the π files are headerless; use
  `_read_vector`/`_read_matrix`), and **the missing header on `abc_results.csv`** (a zero-byte
  pre-created file now counts as needing one) — all fixed. §2 has the current migration
  parameterization.
- **Whole-deme sampling stays.** Measured at POPMULT=5000 over 100 replicates: it biases `fst_loss`
  by **+2.0%** (implied POPMULT ≈ 6000 → ≈ 6060, far inside the prior) and makes `pi_loss` **worse**
  if "corrected", because π is unbiased at any n ≥ 2 and injecting sampling noise only raises the
  floor the pass must clear. **Read invariant 1 as a rule about bias, not about noise.**
- **F_st is HUDSON, computed from π and d_xy. `ts.Fst` is Nei/Slatkin and is ~2× low.** This was
  worth a factor of ~2 in the inferred POPMULT (6,234 → 12,469) and is why `POPMULT_MAX` was raised
  to 25000. The ratio `1 + Hw/d_xy` was predicted ≈1.98 and measured 1.971–1.980 in all three years
  independently, which is what confirms it was a pure estimator swap. **Standing check: a grep for
  `.Fst(` across the repo must return nothing.** Recorded numbers in `out/*.jsonl` and in the
  §6.1.1 / §6.2.1 tables predate the fix and are Nei-scale — never compare one against a post-fix
  number.
- **pixy pools F_st properly** (average-of-ratios was ruled out — the strong n-dependence it would
  produce is absent), so the empirical WC targets are trustworthy and stay untouched. **Do not
  reconstruct empirical Hudson from `averaged_pi` and `averaged_dxy`** — a ~1% relative π/d_xy
  offset moves F_st by more than the whole signal, and such an offset is present (1.13/1.60/1.39%
  by year). **Mechanism found 2026-09-11 (§7.2.2E):** a ~24% per-individual heterozygote deficit
  in the genotype calls removes within-individual differences from π and not from d_xy; correcting
  π for it removes the Hudson floor in all three years and closes ~91% of the gap to WC. **The
  rule stands** — the residual is still a large fraction of the signal.

---

### 6.8 `DEFAULT_RECOMBINATION_RATE = 2.75e-6` inherits the §6.1 Ne error [OPEN 2026-09-06]

Found while scoping the LD statistic (§7.5). **Not a transcription bug — the code faithfully uses
what Cohen et al. report.** The problem is one level up, and it is §6.1 again.

**What the paper actually says** (p6, verified from the PDF in the repo root):

> "The per-generation per base recombination rate (r) differed only slightly yet was also
> significantly different (**r_HAN = 2.75e-6** and r_LI = 1.95e-6; p-value <2.2e-16)."

**How they got it (p3, §2.4).** pyrho v0.1.6 estimates the *population* recombination rate; the
paper then converts with **`ρ = 4·N_e·r`**, taking `N_e` from dadi. So `r` is **derived, not
measured** — the measured quantity is `ρ_HAN = 0.163`. The arithmetic is self-consistent:
`4 × 14,777 × 2.75e-6 = 0.1625 ≈ 0.163`. ✓

**Therefore `r` inherits the dadi `N_e` error exactly, and §6.1 shows that error is ~52×:**

| `N_e` used | value | implied `r` | in cM/Mb |
|---|---|---|---|
| Cohen dadi `N_WI` (what the paper used, and what we use) | 14,777 | **2.758e-6** | **276** |
| Watterson from Cohen's **own** SNP count (§6.1) | 3.5e5 | 1.164e-7 | 11.6 |
| our corrected π (§6.1) | 1.452e6 | **2.806e-8** | **2.81** |

**276 cM/Mb is not a possible recombination rate** — it is ~275× a typical eukaryotic value, and
the paper itself notes its ρ estimates are "approximately an order of magnitude lower than
*Drosophila*", which is only coherent if `r` is ordinary and `N_e` is not. Correcting `N_e` gives
**2.81 cM/Mb, an entirely ordinary insect rate.**

**Third, independent line of evidence.** `ToUseOnBeagles/CalculateLD.py` was written with
`MAX_DIST = 100_000`, `BIN_SIZE = 1_000`. LD half-decays where `4·N_e·c = 1`, i.e. at
`1/(4·N_e·r)` bp. At a deme `N_e` ≈ 600 that is **152 bp at r=2.75e-6** (everything collapses into
bin 1 — the 100-bin structure would be flat and useless) versus **15.2 kb at r=2.75e-8** (bin 16,
squarely in the design's sweet spot). **Whoever chose that binning was implicitly assuming
r ≈ 1e-8.**

**Consequences, in order of how much they bite:**

- **For LD (§7.5) this is decisive.** LD decay is governed by **`ρ = 4·N_e·r`, not by `r` alone**,
  and `ρ = 0.163` is what pyrho actually measured, so it does not depend on the dadi `N_e`.
  **But ρ's UNITS are themselves ambiguous, so it is not the clean anchor it first appears
  [CORRECTED 2026-09-06]:**

  | reading of `ρ_HAN = 0.163` | ρ per bp | LD half-decay `1/ρ` | verdict |
  |---|---|---|---|
  | per **bp** (what the paper's own arithmetic implies) | 1.63e-1 | **6.1 bp** | **not a real organism** |
  | per **kb** | 1.63e-4 | **6.1 kb** | plausible — human ρ≈4e-4/bp gives ~2.5 kb |

  The paper is internally consistent in per-bp units (`4 × 14,777 × 2.75e-6 = 0.163`), but per-bp
  implies LD vanishing within ~6 bp, which no organism does. So **either the dadi `N_e` is wrong
  (§6.1: ~52× low) or pyrho's output units were misread on the way in.** Both roads end at the
  same place: **the derived `r` is unreliable, and ρ cannot be used as a numeric anchor until its
  units are confirmed.** Settling it needs pyrho's documented output units, or — far easier —
  **one run of `CalculateLD.py` with log bins, which measures the empirical decay position
  directly and makes the whole arithmetic unnecessary (§7.5 step 1).**

  > **MEASURED 2026-09-07 — NEITHER reading is right, and ρ is now moot as an anchor.** The
  > empirical run happened (§7.5.1): half-decay of the excess-over-plateau is **~60 bp** in all
  > three years, i.e. `ρ_obs ≈ 1/60 = 0.017 per bp`. That is **10× below** the per-bp reading of
  > `ρ_HAN = 0.163` and **~100× above** the per-kb reading. **Do not use `ρ_HAN` as a numeric
  > anchor in either unit** — the direct measurement supersedes it, and the residual 10× gap is
  > unexplained (candidates: Cohen's `ρ` is Hancock-only while ours pools 17 sites and 3 years;
  > and Beagle imputation, which inflates LD and would push our `1/ρ` *longer*, i.e. the wrong
  > way to close the gap). What survives from this table is the *conclusion*, not the arithmetic:
  > `r` derived as `ρ/(4N_e)` inherits whichever `N_e` you accept, and that question is still
  > §6.1's.
- **For π and F_st, nothing changes.** Neither has any recombination signal (§6.3); `r` only sets
  how many independent trees the genome contains, which affects variance, not expectation. **No
  existing result is invalidated.** This is why it went unnoticed.
- **For SLiM cost, it is pure waste.** §6.3 raised forward recombination 1e-8 → 2.75e-6 ("a 275×
  increase"), and that 275× is *exactly* the §6.1 error factor. The old hardcoded 1e-8 was closer
  to right than the value that replaced it — for the wrong reason, since it was a placeholder, not
  an estimate. Extra recombination inflates edge counts and forward-phase memory for no signal.

#### 6.8.1 The linkage map — `r_true ≈ 1.0e-8`, so the model wants 8.0e-7 [FOUND 2026-09-08]

**The map exists, and it is the `N_e`-free anchor this section asked for.**

**Hawthorne, D. J. (2001), *Genetics* 158(2):695–700** — AFLP map of *L. decemlineata* from 74
backcross individuals (Long Island × Freeville, NY):

- **total map length 1,032 cM**, **18 linkage groups** (= haploid chromosome number)
- 172 AFLP + 10 codominant markers; 86 well-spaced AFLPs in the published map
- mean intermarker distance 11.1 cM (±6.72)

**Yan et al. (2023), *Scientific Data* 10:36** — chromosome-level assembly, **~1,008 Mb anchored
to 18 chromosomes (17 + XO)**, scaffold N50 58.32 Mb.

```
1,032 cM / 1,008 Mb = 1.02 cM/Mb   ->   r_true ~ 1.0e-8 per bp per generation
```

| source | `r` | vs the map |
|---|---|---|
| **linkage map (this)** | **1.02e-8** | — |
| `DEFAULT_RECOMBINATION_RATE` (from Cohen) | 2.75e-6 | **269× high** |
| §7.5.4's internal LD estimate | 3.3e-7 | 32× high |
| §6.8's external correction (ρ / 4·π-Ne) | 2.81e-8 | 2.7× high |
| **the old hardcoded placeholder (§6.3 "bug value")** | **1e-8** | **0.98× — right by accident** |

This section's suspicion that "the old hardcoded 1e-8 was closer to right than the value that
replaced it" is confirmed to ~2%, and it vindicates whoever chose the original `MAX_DIST`/
`BIN_SIZE` LD bins, which this section deduced were "implicitly assuming r ≈ 1e-8."

**Caveats.** 1,032 cM is a **lower bound** — 86 placed markers over 18 groups leaves chromosome
ends uncovered; adding half the mean intermarker spacing per end gives ~1,230 cM → 1.22 cM/Mb, so
realistically **r_true ≈ 1.0–1.5e-8**. Genome-size choice moves it <20%. **Unchecked from the
abstract:** whether the map is sex-averaged or sex-specific (the BC1 was 33 males, 46 females and
the paper is partly about sex chromosomes), and whether the X linkage group should be dropped to
match our 17-chromosome data. And **§7.5 records a ≥10× spread in effective `r` BETWEEN
chromosomes**, so 1.02e-8 is a genome-wide average over a wide distribution, not a precise
constant.

**THE VALUE THE MODEL WANTS IS NOT 1.02e-8.** §6.1.3 shows the simulation runs at Q ≈ 78.5, so a
consistently rescaled `r` is `r_true × Q` = **8.0e-7**. Setting `DEFAULT_RECOMBINATION_RATE` to
the real-world 1.02e-8 would be **78× too small**. The current 2.75e-6 is **3.4× too high**, not
269×.

**[RECOMMENDED, not yet applied] Set `DEFAULT_RECOMBINATION_RATE = 8.0e-7.** Justified
independently of any LD fit — a linkage map plus §6.1.3's algebra — and it made the four-point
four-point LD sweep run 3–4× faster, since less recombination means fewer edges. It breaks comparability
with batches 1 and 2 and with §3.1's timings; Sohan's call 2026-09-08 was that re-running batches
is cheap enough that this is acceptable.

> **Three places would silently ignore the change. Fix them in the same commit.**
> - **`Main.py:24` still has `recombination_rate=2.75e-6` as a hardcoded default** — §10.1 claims
>   all three scale-setting parameters became `None` + raise, and `AnalyzeTreeSeq` did, but
>   `Main.main` was missed. It is a literal, not the constant, so it would not follow. The ABC path
>   always passes `r` explicitly, so no trial was affected — the §10.1 story exactly.
> - **Five diagnostics carry their own copies**: `fst_subsample.py:64`, `mu_calibrate.py:62`,
>   `noise_floor.py:308`, `qpost.py:32`, `ridge_sweep.py:172`. Only `ld_probe.py` reads the
>   constant.
> - **`LD_MIN_BIN = 562` is DERIVED from `r`** (`ABCAnalysisNoRedis.py:44`): the time-depth cut is
>   `1/(2·G·r)`, which at 8.0e-7 becomes **1,929 bp** (first bin edge 3162), not 562.

**The `r` question is now settled enough to act on, and it did NOT rescue LD — see §7.5.**

> **~~STRENGTHENED 2026-09-07~~ — that note was WRONG, and the simulation falsified it the same
> day (§7.5.3). Retained because the error is instructive.**
>
> It argued analytically that the observed ~60 bp half-decay implies `N_e,deme ≈ 1515` via
> `1/(4·N_e,deme·r)`, i.e. **POPMULT ≈ 15,000** at the current `r` — pleasingly close to §6.7's
> F_st-implied 12,469. **Measurement says the opposite.** At POPMULT=5000 the *simulated* curve
> already crosses halfway at **32–56 bp**, i.e. it decays FASTER than the observed 56–100 bp, so
> the fit wants POPMULT **lower**, not 3× higher. At POPMULT=2000 the simulated halfway lands at
> 56–100 bp and matches the observed bin in 2015 and 2023.
>
> **Why the arithmetic failed:** `report_halfway` is the midpoint between the curve's max and min,
> and the min is the `1/n_hap` floor (§7.5.1 pt 1) — it is *not* the Sved `1/ρ`. So `1/(4·N_e·r)`
> was never the quantity being compared, and the analytic map from ρ to POPMULT is invalid on
> both sides. **The simulation is the only valid map.** This is the second time in this section a
> confident derivation from ρ has been overturned by a measurement; treat ρ arithmetic as a
> sanity check, never as a result.
>
> **What survives, and it is the part that matters:** `r` and POPMULT are confounded in LD exactly
> as μ and POPMULT are in π — LD constrains `4·N_e·r`, so an `r` inflated by factor k deflates the
> LD-preferred POPMULT by k. Treat `r` as **half of a calibration constant** whose meaningful
> content is `ρ = 4·N_e·r`, never as a biological recombination rate — the same language rule §6.1
> imposed on μ. §7.5.3 turns this into a quantitative handle on the `r` error.

---


---

# §I — the ABC distance, kinship, noise floor, batch 1 and LD (old §7 through §7.5.6) [ARCHIVED 2026-09-27]

CLAUDE.md §7 keeps the fitted set, the losses, the priors and one-line conclusions.

## 7. The ABC distance — design

In ABC there is no gradient-descent "loss." This is the **distance** `d(S_sim, S_obs)` used for
rejection/weighting; it only has to *rank* parameter draws sensibly.

**Fitted:** element-wise **log-π**, off-diagonal **F_st**.
**Computed but NOT fitted — RETIRED 2026-09-08:** binned **LD decay**. `ld_loss` is wired
end-to-end and is still the most `pop`-informative statistic ever measured here (unique R² 0.258).
**It is not fitted, and the decision is closed — see §7.5.** `r` was pinned (§6.8.1) and it did
not help: the simulation is a single-locus, single-`r` model and the empirical target pools 17
chromosomes spanning ≥10× in effective `r`, which are structurally different objects.
**Diagnostic:** **IBD slope**, **d_xy**, **genetic relatedness**.
**COMPUTED, NOT YET FITTED — temporal F_c** (§7.9.4, §7.9.6, §7.9.10). It is the first statistic
measured here that reads N nearly independently of the dispersal assumption — **but that is a
property of the PERSISTENT-DEME model, not of the statistic [narrowed 2026-09-20, §7.9.12F].**
Under the re-founding structure F_c reads **founder number `K`**, with an elasticity in N of
−0.08/+0.04 against −0.996 here, and its kernel sensitivity rises 5.2×. Everything this
paragraph and §7.9.6/§7.9.8B/§7.9.11B claim about F_c is conditional on demes being persistent. **`fc_loss` is ON by
default since 2026-09-12** (`COMPUTE_FC` defaults to 1 in both modules, commit `3818aa7`), so
every trial computes it. **Batch 3 LANDED 2026-09-15 (§7.9.11); the corner replicates LANDED
2026-09-16 (§7.9.11H).** `abc_standardize.py` already lists it in `FITTED_STATS`, but `WEIGHTS` is still `None`, so nothing has been scored: **no weights
are installed and no `pop` posterior may be reported.** F_c tightens `pop` (IQR ratio 0.49 at
top 1% against batch 1's 0.82), but every version of `fc_loss` pushes `pop` toward the ceiling,
because simulated F_c sits above observed at both gaps across the whole prior and gains drift
from t=8 to t=16 while the data gains none (4.8 sd at fixed parameters). **The corner replicates
settle that this is MODEL STRUCTURE, not a prior problem** — at the prior's most favourable corner
the drift growth falls 4× and then stops at +0.0011, still above observed in 10 of 10 replicates
(§7.9.11H). Read §7.9.11 before installing weights.

> **IBD was demoted from fitted on 2026-07-29** (§7.1). The switch that actually matters is
> `abc_standardize.py::FITTED_STATS`; `calculate_losses` still returns `ibd_loss`.

> **START AT §7.9 IF YOU ARE LOOKING FOR A NEW STATISTIC.** §7.9.1 measures the ceiling on all of
> them — only 0.5–7% of the genealogy is POPMULT-dependent, and that fraction is itself an `Nm`
> quantity — which is why π, d_xy, relatedness, IBD and LD all failed the same way. It also gives
> the filter a replacement must pass.

### 7.0 Small-subpop exclusion [VERIFIED]

`ABCAnalysisNoRedis.py` has `EXCLUDE_SMALL_SUBPOPS = True`, `MIN_SUBPOP_N = 4`. Subpops with
fewer than 4 diploid individuals are dropped from the **fitted** statistics via
`get_keep_mask(year)`, a boolean mask in specifier-matrix row order applied identically to the
observed and simulated sides.

Why: at n ≤ 3, **46.7% of 2015's pairs return a negative F_st** (vs 5.3% at 4≤n≤7 and 0% at n≥8) —
those entries scatter around a noise floor rather than measuring differentiation.

Only 2015 is affected. It drops exactly two sites — **`Arlington2015` (n=2) and `H67-2015` (n=2)**
— removing exactly the 45 noise-floor pairs. **Side effect worth knowing:** `Arlington2015` is the
typo duplicate (§4), so for the fitted statistics the duplicated-`Arlington` ambiguity is now moot
— only `Arlington-2015` (n=4) survives. The underlying data question is still open.

**Relatedness is deliberately NOT masked** — it is centred on the populations present when it was
computed, so slicing rows/cols is not the same as recomputing on the subset (invariant 2). Both
sides stay full-size. π, F_st and d_xy are all safe to slice.

`calculate_losses` returns per-statistic, un-standardized, count-normalized distances (mean over
entries within a year, then mean over years):

```
pi_loss     = mean_years( mean_i | log pi_sim,i - log pi_obs,i | )   <- fitted, log-space
fst_loss    = mean_years( mean_pairs | Fst_sim - Fst_obs | )         <- fitted
ibd_loss    = mean_years( | slope_sim - slope_obs | )                <- diagnostic (was fitted)
dxy_loss    = mean_years( mean_pairs | dxy_sim - dxy_obs | )         <- diagnostic
genrel_loss = mean_years( mean_pairs | R_sim - R_obs | )             <- diagnostic
ld_loss     = mean_years( mean_bins | r2_sim,b - r2_obs,b | )        <- fitted
```

`ld_loss` uses only bins with `bin_lo >= LD_MIN_BIN` (562 bp), and pools demes as
`sum(sum_r2)/sum(cnt)` over the kept demes — never a mean of means (§5.2, invariant 4).

> **`ld_loss` is in `calculate_losses` and in `CSV_FIELDNAMES`, but is NOT in
> `abc_standardize.py::FITTED_STATS`, and that is a CLOSED DECISION, not a missing step (§7.5).**
> Its weight was derived from the pilot batch and is unusable; adding it to `D` makes `fst_loss`
> stop contributing. Do not wire it in.

### 7.1 Why IBD is not fitted [VERIFIED 2026-07-29]

**The observed IBD slope is indistinguishable from zero in all three years.** Mantel test, 9999
permutations of site labels, using the project's own `ibd_slope`/`get_site_geo_distances`:

| year | n | pairs | slope | Mantel r | **p (two-sided)** | \|slope\|/null_sd |
|---|---|---|---|---|---|---|
| 2015 | 24 | 550 | −1.42e−03 | −0.069 | **0.639** | 0.46 |
| 2019 | 17 | 272 | +5.03e−04 | +0.214 | **0.148** | 1.43 |
| 2023 | 20 | 380 | +5.53e−05 | +0.004 | **0.985** | 0.02 |

The sign flips across years — noise, not a weak real effect. **Identical conclusion on the old
per-SNP targets**, so this does not depend on the denominator fix.

A Mantel test is required here, not an OLS p-value: `ibd_slope` fits over all n(n−1) ordered
off-diagonal pairs (552 for 2015, from 24 sites), which are massively non-independent — every site
appears in 23 of them. An OLS p-value would treat correlated pairs as independent observations and
report significance almost regardless of signal. Permuting **site labels** is the exchangeable unit.

Consequences:
- **`scale` (dispersal-kernel decay) is unidentifiable from these data.** Fix it or report it as
  unidentified. Do not present it as inferred.
- Fitting a target that is noise adds a pure-noise term to the standardized `D`, costing
  acceptance efficiency and blurring the parameters that *are* identifiable.
- **`total_migration` is probably still identifiable, via the F_st level rather than IBD:** at
  mean F_st ≈ 0.003–0.008, **Nm ≈ 31–83**.

**Geographic scale is not the explanation** — sites span 1.7–160 km, median pairwise ~34 km, in all
three years. That is ample range to detect IBD in an insect. And **Rousset's slope assumes
drift–dispersal equilibrium**, which is questionable for a recent fast-spreading invader
(Whitlock & McCauley 1999), so a null result does not cleanly mean "no IBD" — it can equally mean
"not at equilibrium yet."

### 7.2 F_st carries real signal, but it is site-coherent, not distance-structured [VERIFIED]

F_st was checked for sample-size noise-domination before being left as the sole structural fitted
statistic. **It is not noise-dominated.** Spearman rho between `min(n_i,n_j)` and `|F_st|`, with
site-label permutation: 2015 −0.082 (p=0.63), 2019 +0.242 (p=0.23), 2023 +0.082 (p=0.71). No
association. (2019 has almost no leverage — sizes are 5–7.)

Instead, **each year's extremes converge on a single site**, which is what real differentiation
looks like and noise does not:

- **`Mortensen9-2015`** (n=5): mean F_st 0.0686, **5.7× the next-highest site.**
- **`H41-2023`** (n=7): mean F_st 0.0472, **3× the next-highest.**
- 2019: no isolate at all; max pairwise F_st 0.0092.

**These are almost certainly a within-site sampling artifact, not landscape structure. [INFERRED,
strong]** Both isolates are simultaneously the **lowest-π** and (near-)**highest within-site
relatedness** site in their year:

| year | corr(mean F_st, π) | corr(mean F_st, self-relatedness) | isolate |
|---|---|---|---|
| 2015 | **−0.721** | +0.339 | `Mortensen9-2015`: π rank 1 (lowest), selfRel rank 3 |
| 2019 | −0.380 | −0.089 | — |
| 2023 | **−0.915** | **+0.748** | `H41-2023`: π rank 1 (lowest), selfRel rank 1 |

Low within-site diversity + high within-site relatedness + high F_st against everything is the
signature of **a sample of close relatives** (e.g. beetles taken off one plant or one egg mass —
entirely plausible for CPB) or a recent founder/bottleneck event at that field. Either way it is a
*local* phenomenon the landscape migration model cannot and should not reproduce.

> **RESOLVED 2026-09-11 by individual kinship (§7.2.2) — one cause per isolate, and neither is a
> founder event.** `Mortensen9-2015` is a family: all 10 pairs are third-degree or closer after
> correction. `H41-2023` is not related at all: three of its seven samples failed genotyping
> (observed heterozygosity 0.004–0.106 against ~0.27), which fakes the same
> low-π / high-selfRel / high-F_st signature. Still local, still not simulable — the conclusion
> above holds, the mechanism changed.

**Geography rules out the innocent explanation.** The isolates are geographically *ordinary*:
`Mortensen9-2015` ranks 19/24 on mean distance to other sites, `H41-2023` ranks 17/20. Meanwhile
the genuinely remote sites are undifferentiated — `Alsum59-2023` is the most remote (114.6 km) and
ranks 11/20 on F_st. corr(mean F_st, mean distance) = −0.065 / +0.251 / +0.014.

Two implications:
1. **The dispersal kernel may be misspecified.** The data says "one site is an isolate, the rest
   are near-panmictic"; an exponential distance-decay kernel with one global `total_migration`
   produces smooth structure and cannot make a single deme an isolate. Element-wise F_st fitting
   will be dominated by pairs the model structurally cannot match.
2. **π and F_st are not independent statistics here** (r = −0.72 / −0.92 in 2015/2023). Fitting
   both with equal weight after MAD-standardization partly double-counts one signal.
3. **The two isolates also dominate the observed π *spread*** [VERIFIED 2026-08-11]. Dropping the
   single site cuts the fitted between-site log-sd of π from **0.0422 → 0.0212** in 2015
   (`Mortensen9-2015`, 50% of the spread) and **0.0477 → 0.0155** in 2023 (`H41-2023`, 67%).
   2019, which has no isolate, sits at 0.0160 already. So the *genuine* between-site π spread is
   **0.014–0.021** in all three years, and roughly half of the raw spread is the artifact.

   This matters for reading §6.1.1: simulated π spread falls with POPMULT (log-sd 0.10–0.13 at
   POPMULT=500, 0.028 at 2000, 0.010 at 5000). Against the *raw* observed spread the simulation
   looks 4× too flat at POPMULT=5000; against the **cleaned** spread it is only ~1.5×, and the
   observed value is bracketed inside the prior (matching around POPMULT ≈ 3000–4000). The
   simulation is not failing to produce structure — it is failing to produce an artifact, which
   is correct behaviour.

   **~~Caveat, unresolved:~~ ANSWERED 2026-08-12 — they do NOT covary. See §7.2.1.** The
   suspicion was right: `pi_loss` is scoring flatness, not site-level fit.

**There is deliberately no `total_loss` during the pass.** The combined standardized distance
`D = sqrt(Σ_j w_j (loss_j/σ_j)²)` with `σ_j = 1.4826·MAD` is built **offline** by
`Python_Code/abc_standardize.py`, using the run set as its own pilot batch, with equal weights to
start. This commits us to **rejection ABC**, which is what CHTC wants.

Rationale for the choices:

- **F_st, not raw d_xy.** d_xy ≈ π + differentiation, so its entries are nearly redundant with π
  and its differences are dominated by the diversity-level mismatch. Fitting d_xy would
  reintroduce the μ-degeneracy.
- **π element-wise, not mean/SD.** The specifier-matrix mapping makes simulated-subpop ↔ real-site
  a genuine within-year correspondence (§4), so element-wise is meaningful and dimensions always
  match.
- **π in log space; F_st and IBD not.** F_st is bounded near 0 and the IBD slope can be negative
  (log undefined). Log-π gives relative-error semantics and linearises the θ=4Nμ ridge.
- **IBD slope:** regress `F_st/(1 − F_st)` on **ln**(geographic distance) over all off-diagonal
  pairs (Rousset), take the OLS slope. One scalar regardless of k. Real-site coordinates from the
  specifier matrix are used for **both** the observed and simulated regressions.
- **Relatedness is not fitted** — it is linearly dependent by construction (centring makes rows sum
  to ~0) and measures the same signal as F_st. Kept as a **posterior-predictive check**: a
  statistic you didn't fit is far better validation than one you did.
- **Normalize out year entry counts** so 24- vs 17- vs 20-subpop years contribute comparably.
- **IBD is retained as a diagnostic** for the same reason — with the slope no longer fitted, a
  simulated-vs-observed slope comparison becomes an honest posterior-predictive check.

Verified behaviour [VERIFIED]: `calculate_losses(obs, obs)` → all five losses exactly 0.0.
Perturbed sim (π×1.1, F_st+0.02) → `pi_loss = 0.0953 = ln(1.1)` exactly, `fst_loss = 0.02`.
Real-data IBD slopes finite: 2015 −1.42e-3, 2019 +8.04e-4, 2023 +1.44e-4 (weak/mixed).

**Current priors** (`ABCAnalysisNoRedis.py`). **The three scale constants moved to `scale_constants.py` on 2026-09-09 and μ left the prior entirely — §7.9.3.**

| parameter | prior | note |
|---|---|---|
| `m` (kernel decay) | `lognorm(s=1.5, scale=1e-4)` | **MUST BE FIXED BEFORE BATCH 2 (§7.4.2).** Not merely "unidentifiable from IBD" (§7.1) — it is *actively* absorbing the F_st signal that should constrain N. Prior spans 6.1e-7–1.4e-2, i.e. `exp(−m·d)` from ~0.98 to `e^−476` at 34 km: uniform kernel to total isolation. Unique R² inside `fst_loss` is **0.429 against 0.091 for `pop`**. Batch 1 inferred it sharply (IQR ratio 0.05 at top 1%) and left N at the prior. **Never report an `m` posterior — IBD is null (§7.1), so nothing validates it.** §7.9.6 measured the damage directly: across this prior's informative band, F_st's implied N swings **5.3×** with the kernel alone — and temporal F_c's swings only 1.23×, which is why F_c may let `m` stay free |
| `total_migration` | `U(0.001, 0.301)` | **placeholder ceiling** — needs a biological bound. Also trades against N through Nm: spearman(log `pop`, `total_migration`) = **−0.336** *within* batch 1's accepted set (§7.4.2). Fixing `m` alone will not free N while this stays 300× wide |
| `pop` (POPMULT) | `U(2000, 25000)` | ≈ 6.7k–83k individuals. **Raised from 12000 on 2026-08-26** — the §6.7 F_st fix moved the implied value onto the old ceiling. ~44 GB / ~1.9 h per trial at the top (§3.1) |
| `numClusters` | {1,2,3} | **CSV records the raw draw; actual count is ×33.** The ×3 ceiling (99 demes) was NOT a modelling choice — it sat one under msprime's 100-populations-per-event limit, **lifted 2026-09-09** (§7.9.5). Still {1,2,3}: raising it is now a deliberate decision, because deme size is `Average Count × POPMULT / numSubpops` and more demes means *smaller* ones |
| `mutation_rate` | **FIXED 5.8e-7** — no longer drawn | **Left the prior 2026-09-09 (§7.9.3).** Was `lognorm(s=0.02, scale=4.646e-7)`, tightened 0.5→0.05→0.02 and still the largest driver of `pi_loss` (unique R² 0.149 on μ against 0.080 on `pop`, pilot batch). There was never an unknown to draw: μ is *defined* by `4·Ne_anc·μ = π_obs` with both other terms fixed. Value is `MU_TRUE × Q` = 5.8e-9 × 100 (§7.9.3); **report θ=4Nμ, never μ** |
| `recombination_rate` | **fixed 1.02e-6** | `R_TRUE × Q` = 1.02e-8 × 100 (§6.8.1, §7.9.3). Was 2.75e-6. Reaches SLiM (§6.3-6.7); not inferred, and **must not be** (§7.5) |
| `ancestral_Ne` | **fixed 5259** | Not a prior entry — listed so the scale is in one place. Derived from the declared Q=100 (§7.9.3); was 6700 |

### 7.2.1 π does NOT covary site-by-site — `pi_loss` scores flatness [VERIFIED 2026-08-12]

`diagnostics/pi_covary.py`, on the per-subpop vectors `mu_calibrate.py` now stores. POPMULT=5000,
numClusters=33, seed 1 (the re-run reproduced μ=4.646e-7, `pi_loss`=0.02132, `fst_loss`=0.00830
exactly, so this sits on the same tree as §6.1.1).

**Test A — correlation of log π_sim against log π_obs**, site-label permutation null, 9999 perms
(a permutation test for the same reason as §7.1: few sites, coupled by migration):

| year | n | Pearson | p | Spearman | p | vs `branch_div` | p |
|---|---|---|---|---|---|---|---|
| 2015 | 22 | −0.091 | 0.61 | −0.068 | 0.77 | −0.092 | 0.60 |
| 2019 | 17 | +0.538 | **0.027** | +0.194 | 0.46 | +0.550 | 0.020 |
| 2023 | 20 | −0.021 | 0.91 | +0.146 | 0.54 | −0.034 | 0.85 |
| pooled | 59 | +0.025 | 0.86 | | | | |

**2019's nominal hit is not real.** Spearman is only 0.194, so the Pearson is leveraged by a few
points; one p<0.05 in three tests is ordinary; and it **does not replicate**. The simulation's
low-diversity demes are *the same three cluster rows every year* (`Alsum18`/`Alsum59`/`Arlington`,
`branch_div` ≈ 26.1–26.4k against ≈27.0k elsewhere) — fixed simulated structure against a shuffling
observed rank. They are low-π in 2019, but `Alsum25`/`Alsum140` are among the **highest** observed
π in 2015, and 2023 splits.

**Nor is the null an artifact of noise on the observed side [VERIFIED 2026-08-14, §6.6].** The
obvious objection — observed π is estimated from 4–19 diploids, so a real correlation would be
attenuated — was measured with `diagnostics/fst_subsample.py`: sampling sd at the real n_i is
0.0050–0.0053 in log units, only 1.2%/9.6%/1.2% of the observed between-site variance, giving
attenuation factors of 0.992/0.955/0.994. De-attenuated, the Pearsons above become
−0.092/+0.563/−0.021. **The correlation is absent, not hidden.**

**Test B — the flat-simulation floor. This is the decisive number.** `pi_loss` for a
*structureless* simulation, levelled by the same weighted-median rule, is **0.02094**. It depends
only on the *observed* vectors, so it applies at every POPMULT already measured — no re-simulation:

| POPMULT | `pi_loss` | vs flat floor | sim log-sd |
|---|---|---|---|
| 500 | 0.06717 | **3.21× worse** | 0.1187 |
| 2000 | 0.02898 | 1.38× worse | 0.0295 |
| 5000 | 0.02132 | **1.02× — at the floor** | 0.0105 |

**No POPMULT tested beats a simulation with no between-site structure at all.** The whole
0.067→0.021 improvement is excess *uncorrelated* scatter decaying to that floor, and it saturates.

**Test C — shuffle null in the units of the objective** (permute the sim within year, re-calibrate
μ each time): real 0.02052, null mean 0.02238 ± 0.00113, **percentile 5.5, p = 0.055**. The
ordering `real (0.0205) < flat (0.0209) < shuffled (0.0224)` is exactly the L1 geometry — scatter
in the wrong place is worse than no scatter. Correctly-placed structure buys **+2.0%** of
`pi_loss`. (Verified on synthetic data first, per §10: a flat sim beat an uncorrelated
equal-spread sim in **400/400** draws.)

**Consequence — π's POPMULT information is one-sided and saturating, not peaked.** It genuinely
excludes small POPMULT (500 is 3.2× worse than saying nothing), but above ≈5000 it is already at
the floor and cannot discriminate further. L1 **never penalises the simulation for being too
flat**, so the preference runs monotonically upward and identifies no optimum. **Do not read a
large-N π signal as corroborating F_st** — F_st carries the actual peak (§6.2, wanting POPMULT
≈6000), and the two are not independent anyway (§7.2: r = −0.72/−0.92). Keep π fitted — the lower
bound is real — but treat it as a bound plus a posterior-predictive check, not as a second vote.

**Mechanism probe, unexplained:** corr(log `branch_div`, log deme size) is **negative** in all
three years (−0.217/−0.284/−0.156). Larger demes show *slightly lower* diversity, which is
backwards for drift. The effect is tiny (sim log-sd 0.010) and the three low-`branch_div` demes are
also among the largest, so isolation is probably confounded with size here. [OPEN] — but note
this means the simulation's per-site π variation is **not** deme-size-driven, i.e. its only
plausible route to matching observed π site-by-site is not the one operating.

### 7.2.2 Individual kinship — Mortensen9 is a family, H41-2023 is failed samples, and every sample carries a ~24% heterozygote deficit [VERIFIED 2026-09-11]

`ToUseOnBeagles/CalcKinship.py` (KING, genome-wide pooled counts, MAF ≥ 0.05 over the whole
year) was run on the Beagle machine; output is in `out/kinship_out/`. The known-truth test
`diagnostics/test_kinship.py` passes all four checks (re-run 2026-09-11), so **the estimator is
right. What is wrong is the raw output's scale.**

> The correction is `diagnostics/kinship_correct.py` (added later the same day, §7.2.2H). It
> reproduces every number in A–D in seconds, with no VCFs.

#### A. Unrelated pairs read −0.30, not 0

Within-site site medians run −0.19 to −0.39 (H41-2023 −0.65), and between-site pairs sit in the
same place. That is not relatedness. It is a **per-individual heterozygote deficit**.

Model: individual `i` has a rate `e_i` at which a true heterozygote is called homozygous for a
random allele. Derived, exact in expectation under HWE within a year:

```
KING_obs = (phi_true - ebar) / (1 - ebar),     ebar = (e_i + e_j) / 2
het_obs,i = (1 - e_i) * H
```

so an unrelated pair reads `−ē/(1−ē)`, and −0.30 means `ē ≈ 0.23`. **Fit:** `e_i` by least
squares on **between-site pairs only** (`phi_true ≈ 0`), then correct every pair with
`phi_true = phi_obs·(1−ē) + ē`.

| year | individuals | between-site pairs | fit R² | median `e` | `e` range | corr(obs het, 1−e) | CV of het: raw → het/(1−e) | implied true het |
|---|---|---|---|---|---|---|---|---|
| 2015 | 167 | 13,247 | 0.960 | 0.242 | −0.097 – 0.367 | 0.955 | 0.063 → 0.019 | 0.2699 |
| 2019 | 115 | 6,220 | 0.969 | 0.242 | −0.094 – 0.456 | 0.980 | 0.071 → 0.018 | 0.2713 |
| 2023 | 138 | 9,032 | 0.983 | 0.240 | 0.051 – 0.954 | 0.994 | 0.141 → 0.063 | 0.2672 |

The het column is the check: raw heterozygote counts are not what the fit minimises, yet they
fall into line with the fitted `e_i` (CV shrinks 3–4×), and all three years imply the same true
heterozygosity at MAF ≥ 0.05. Corrected between-site baseline: median −0.0004/−0.0004/−0.0005,
1st–99th percentile about ±0.02.

**What the counts CANNOT tell apart: miscalling versus genuine inbreeding.** Per individual, a
random-allele miscall at rate `e` and inbreeding at `F = e` give identical genotype expectations.
**Miscalling is the reading [INFERRED]**, because `F_IS ≈ 0.24` uniformly across every field in
every year is not plausible for an outbreeding pest in demes of hundreds; because Cohen et al.
§4.1 warns of exactly this bias in low-coverage CPB data; and because S221's observed
heterozygosity of 0.004 is not a living outbred insect.

#### B. CalcKinship.py's reading guide is backwards for this data

Its docstring says imputation biases kinship **up**, so a null is strong and a positive is
merely suggestive. **The measured bias is strongly DOWN.** The raw `kinship_site_summary.csv`
gives `Mortensen9-2015` **zero** pairs at second degree or closer, and the script's printed
decision rule would have sent the conclusion to "founder/bottleneck event" — which is wrong
(C below). **Do not read the raw `degree`/`ge_*` columns or the printed rule** until the script
reports corrected values (TODO §8.3). Separately, it calls its formula KING's "between-family
form"; it is the homogeneous (within-family) form. That does not matter here, since both give
`−ē/(1−ē)` for unrelated pairs of equal true heterozygosity.

#### C. The two isolates have DIFFERENT causes, and neither is landscape structure

**`Mortensen9-2015` is a family.** All 10 pairs clear the third-degree cut, 9 are second-degree
or closer and 2 first-degree; median **+0.125** (half-sib level), max +0.288. Its `e_i` are
ordinary (0.25–0.37), so the correction is well-conditioned — S176 × S179 (`e` 0.249/0.261)
reads +0.214 on its own. **§7.2's "sample of close relatives" is confirmed; no founder event is
needed.**

**`H41-2023` is three failed samples, not relatives.**

| sample | observed het (expect ≈0.27) | `e` |
|---|---|---|
| S221 | 0.004 | 0.954 |
| S222 | 0.040 | 0.816 |
| S223 | 0.106 | 0.558 |
| S219, S220, S224, S225 | normal | 0.21–0.27 |

Their corrected kinship with each other reads 0.56–0.82 ("duplicate"), but the correction is
ill-conditioned as `e → 1`, and **they read +0.48 to +0.65 against S212 of `H15-2023` (`e` =
0.631), a different field.** Failed samples in different fields resembling one another is a
shared artifact, not kin — e.g. poorly covered sites being called toward one allele, which the
random-allele model does not cover. [INFERRED] **Any kinship value involving `e > 0.4` is
meaningless.**

That also explains why H41 fakes an isolate. A random-allele miscall removes only
within-individual differences, predicting a π deficit of `ē/(2n−1)` ≈ **3.6%** at H41. Observed
π is 0.01002 against a 2023 typical 0.0124, **≈19% low** — so the H41 π and F_st are consistent
with the one-allele artifact rather than with random miscalling. [INFERRED]

Site-level corr(mean `e`, log π) is −0.53/−0.08/**−0.90** and corr(mean `e`, mean F_st)
+0.56/−0.23/**+0.80**, both driven by the two isolates.

#### D. Everywhere else, relatedness is local

Every within-site pair corrected above the second-degree cut, excluding failed samples:

| site | pair | observed | `e` | corrected | note |
|---|---|---|---|---|---|
| `Alsum25-2015` | S15 × S16 | +0.306 | 0.268/0.237 | **+0.481** | duplicate level with ordinary `e` — probably one beetle sequenced twice or a relabel [INFERRED] |
| `Alsum25-2015` | S9 × S10 | +0.143 | 0.164/**−0.097** | +0.172 | hinges on a negative-`e` sample |
| `Alsum18-2023` | S26 × S27 | +0.256 | 0.197/0.058 | +0.351 | first degree |
| `Alsum18-2023` | S21 × S22 | +0.006 | 0.205/0.051 | +0.133 | second degree |
| `OkrayGrosheks40-2019` | S306 × S307 | +0.097 | **−0.094**/0.266 | +0.175 | hinges on a negative-`e` sample |

Plus one marginal third-degree pair in `H66-2023` (S254 × S255, +0.046). **Negative `e` means more
heterozygotes than HWE allows**, whose textbook cause is a mixed-DNA sample — and a mixture looks
related to whoever contributed. Those two flags are less secure than they look. S22 and S27 are
the two lowest-`e` samples in 2023, which is the same concern in milder form.

**Every other site-sample has a corrected median ≈ 0.00 and nothing above the third-degree cut.**
That answers the question the script was written for: **a handful of samples, not the
collection.** Dropping them is a correction, not a caveat on every estimate.

Between sites, the only pair beyond second degree outside the failed samples is S85 (`C5-2015`,
`e` 0.312) × S177 (`Mortensen9-2015`, `e` 0.367 — the highest in 2015) at +0.270. S85 also reads
~+0.045 against several unrelated individuals, i.e. it is mildly off-model, so this one is
unresolved.

#### E. It explains ~91% of the §6.7.1 π/d_xy offset

A per-individual deficit removes within-individual comparisons from π — one in every `2n−1`
haplotype pairs — and d_xy has none. Predicted per-site π deficit `ē_site/(2n−1)`: mean
**1.96/1.96/2.02%**, range 0.64–3.60%. Hudson rebuilt from `averaged_pi` and `averaged_dxy`, n ≥ 4
sites, unordered pairs, before and after dividing π by `1 − deficit`:

| year | as-is: min, negatives | π corrected: min, negatives | pixy WC: min, negatives | mean(WC − Hudson) as-is → corrected | corr with WC as-is → corrected |
|---|---|---|---|---|---|
| 2015 | +0.01132, 0/231 | −0.00411, 58/231 | −0.00190, 12/231 | −0.01786 → +0.00152 | 0.982 → 0.998 |
| 2019 | +0.01601, 0/136 | −0.00297, 42/136 | −0.00148, 12/136 | −0.01802 → +0.00159 | 0.697 → 0.986 |
| 2023 | +0.01394, 0/190 | −0.00357, 80/190 | −0.00206, 32/190 | −0.01829 → +0.00167 | 0.985 → 0.999 |

The floor is gone in all three years. It slightly overshoots, plausibly because `e` is fitted at
MAF ≥ 0.05 while π pools every SNP, and because between-site KING runs slightly negative from real
differentiation, nudging `e` up. **The §6.3–6.7 rule stands — do not reconstruct empirical Hudson
from π and d_xy.** The residual 0.0016 is still a large fraction of mean WC F_st. What changed is
that the offset now has a mechanism. (Genuine inbreeding would produce it identically; see A.)

**What it does NOT touch:**
- **Between-site π spread.** With the two isolates out, log-sd of π moves 0.0212 → 0.0205,
  0.0160 → 0.0164 and 0.0155 → 0.0154. §7.2.1 stands.
- **The π level.** The ~2% offset is absorbed by μ calibration, whose target is `π_obs`.
- **WC F_st, the fitted target.** WC's θ carries a separate within-individual component, so it
  is built to be consistent under within-individual correlation, and the corrected-Hudson ≈ WC
  agreement above fits that. [INFERRED]
- **`averaged_genRel` self-values** are inflated by the deficit. They are diagnostic only.

#### F. Temporal F_c: the flagged samples ARE the outliers [VERIFIED 2026-09-11]

`data/empiricalStats/averaged_temporalFc.csv` is present (§7.9.8E). Flag a matched pair if either
sample contains a within-site pair corrected to second degree or closer (both `e < 0.4`), or any
individual with `e > 0.4`. Among the **15 pairs with n ≥ 7 on both sides**:

| pair | t | F_c − pedestal | flag |
|---|---|---|---|
| `H41-2015 × H41-2023` | 16 | **0.0885** | failed samples S221–S223 |
| `Alsum18-2019 × Alsum18-2023` | 8 | 0.0411 | relatives |
| `OkrayGrosheks40-2019 × -2023` | 8 | 0.0365 | relatives |
| `Alsum59-2019 × -2023` | 8 | 0.0336 | failed sample S32 (`e` 0.456) |
| `H15-2015 × H15-2023` | 16 | 0.0314 | failed sample S212 |
| the other 10 | 8/16 | 0.0216 – 0.0289 | — |

**Perfect separation; exact one-sided rank-sum p = 1/3003.** Caveat, stated plainly: the `e` cut
was chosen with these numbers in view. Separation holds for cuts of 0.35–0.45 and fails at 0.30
(pulls in `C5-2015`'s S85 at 0.312 and `H10-2015`'s S182 at 0.340) and at 0.50 (drops `Alsum59`).
The kinship cut is the script's own pre-set second-degree boundary.

**Mechanism.** Relatives and failed samples both cut the number of effectively independent
haplotypes, which inflates the sampling variance of a site frequency and so raises F_c. The
simulation draws unrelated individuals with perfect genotype calls, so it reproduces neither.

**Three consequences, all binding before a pilot batch:**

1. **These pairs (or their samples) must come out of the target before `fc_loss` sees it.**
   **→ DONE by sample exclusion on both sides; half-sibs deliberately kept — §7.2.2H.**
   `H41` alone sits ~0.065 above the typical excess; pooled within the t=16 gap over nine pairs,
   that is roughly **0.007** — about 40% of F_c's whole-prior drift range (~0.017, §7.9.8).
   Pooling does not protect against it. Any sample removal must update `n_i` on **both** sides
   (invariant 1).
2. **The uniform ~24% deficit is NOT reproduced on the simulated side. [OPEN — measure before
   fitting F_c.]** **→ MEASURED AND FIXED: +0.034–0.035, 5.7× the signal; production now imitates
   the miscalls — §7.2.2H.** A per-individual deficit inflates the variance of a sample frequency by about
   `(1+F)`, so the naive n=7 pedestal would go 0.143 → ~0.177, twice F_c's whole-prior drift
   range. **Do not trust that figure:** §7.9.4 showed MAF ascertainment at n=7 already distorts
   the pedestal decomposition, so the real size has to be measured. The direction is known:
   empirical F_c reads as a smaller N. Measure it by applying per-individual random-allele
   miscalls — `e_i` matched to each field's measured values — to the simulated genotypes in
   `temporal_fc.py`. One cheap alternative is one random allele per individual on both sides,
   which is exactly immune to random-allele miscalls at the cost of half the haplotypes. It fails
   wherever the random-allele model fails, which the failed samples show it does.
3. **No drift is visible between the gaps. [OPEN]** With flagged pairs out, mean excess is
   **0.0248** at t=8 (n=5, sd 0.0023) against **0.0231** at t=16 (n=5, sd 0.0017). Drift predicts
   t=16 > t=8. Possible causes are five pairs a side, n=7 ascertainment compressing the drift term
   (§7.9.4), or genuinely large N. Get the simulated t=8 and t=16 at real n before reading
   anything into it.

#### G. An unverified number to re-derive before acting

CalcKinship.py's docstring says removing both isolates moves the F_st-implied POPMULT **12,500 →
26,400** (2.12×) and takes the three years from 2.88× to 1.35× agreement. **It is recorded nowhere
else and was not re-derived here.** 26,400 is **above `POPMULT_MAX = 25000`**, so if it holds,
dropping the isolates from the fitted statistics re-opens the §6.7 ceiling problem. Verify before
dropping sites.

#### H. F_c made usable: exclusions on both sides, miscalls reproduced in simulation [VERIFIED 2026-09-11]

Two decisions (Sohan, 2026-09-11), both resting on measurements from `diagnostics/fc_miscall.py`.
It ran on `out/ld_probe_p2000.trees` and `ld_probe_p5000.trees`, 16 reps each, at the real
post-exclusion n, with all conditions computed on the same draws so the shifts are paired. Records
are in `out/fc_miscall.jsonl`; `--summarize` pools only records under the current spec.

**1. Relatives: drop only what the simulation cannot produce.** SLiM records each individual's
parents (`pedigree_p1/p2`), so relatives in a random simulated sample are counted exactly:

| POPMULT | field samples holding a shared-parent pair | full-sib pairs | F_c shift from drawing unrelated individuals |
|---|---|---|---|
| 2000 | 100/512 (19.5%) | 0 of 9,392 | −0.0021 ± 0.0002 |
| 5000 | 41/512 (8.0%) | 0 of 9,392 | −0.0009 ± 0.0001 |

Half-sibs are ordinary in a simulated deme, and they carry ~20% of the POPMULT 2000→5000 F_c
signal (0.0015 of 0.0061). Purging them from the empirical side alone would bias F_c low
(Waples & Anderson 2017). **So S22 and S306, both half-sib level, are KEPT, and the full sib S27 is
dropped**, along with the five failed samples. The rule in `fc_common.EXCLUDED_SAMPLES` is now
`e > 0.4` or corrected kinship above the **first-degree** cut (0.177). `kinship_correct.py`
re-derives the list from the kinship files and exits non-zero if the two disagree.

Final list: S32 (`Alsum59-2019`), S27 (`Alsum18-2023`), S212 (`H15-2023`), and S221/S222/S223
(`H41-2023`). n becomes 6 for `Alsum59-2019`, `Alsum18-2023` and `H15-2023`, and 4 for `H41-2023`;
all 19 field pairs survive. **Spec `ee863fff3bbf` → `49afd0877028`.**

**2. Miscalls: reproduced on the simulated side.**

| POPMULT | clean pooled F_c | shift from miscalls | shift from miscalls + unrelated draws |
|---|---|---|---|
| 2000 | 0.16273 | +0.03387 ± 0.00008 | +0.03187 |
| 5000 | 0.15666 | +0.03506 ± 0.00005 | +0.03425 |

| POPMULT 2000→5000 signal | pooled | t=8 | t=16 |
|---|---|---|---|
| clean genotypes | 0.00607 | 0.00312 | 0.00934 |
| with miscalls (**the production rule**) | **0.00488 (80%)** | 0.00213 | 0.00794 |

The shift is **5.7× the clean signal**, so F_c computed on clean simulated genotypes would read far
too small an N. It is nearly constant in N; its small growth (+0.0012 from 2000 to 5000) is the 20%
of signal lost. Note that **t=8 carries little signal in either condition**, so t=16 does most of
the work.

**Code changes:**
- **`fc_common.py`:**
  - `EXCLUDED_SAMPLES` is part of the spec hash.
  - `popfile_members` and `popfile_counts` apply it. `check_exclusions` raises on an absent id or
    a wrong site.
  - `read_miscall_rates`, `field_miscall_rates` and `apply_miscall` are **not in the hash**. The
    empirical side computes nothing differently, so hashing them would force a multi-hour empirical
    re-run for a change that cannot move the empirical number.
- **`CalcTemporalFc.py`:** takes its VCF columns through `popfile_members` and prints the exclusions.
- **`AnalyzeTreeSeq.calculate_temporal_fc`:**
  - One draw per **field-year**, via the new `_subsample_individuals`, which keeps each diploid's
    node pair.
  - Each simulated individual is given one of its field's real retained rates, in random order, and
    `apply_miscall` runs before frequencies are taken.
  - It raises if `data/empiricalStats/miscall_rates.csv` is missing.
- **`ABCAnalysisNoRedis.FC_EMPIRICAL_SPEC = "ee863fff3bbf"`:** the hash the current target was
  computed under, so `COMPUTE_FC=1` raises until that target is re-run and the constant updated.
- **`diagnostics/kinship_correct.py`:** the §7.2.2A fit as committed code; `--write-rates` writes
  the rates file (420 samples).
- **`diagnostics/fc_miscall.py`:** the measurement above.

**Verified:**
- **`test_fc_roundtrip.py`, extended:** exclusions shrink n identically in the script's columns, its
  n_a/n_b and the direct computation (4 pairs at n=6); agreement 4.5e-11; a typo or wrong site in
  the list raises.
- **`apply_miscall` on known genotypes:** e=0 is the identity, e=1 leaves no heterozygotes,
  homozygotes are untouched, the realised rate is within 0.004 of e, and allele-1 frequency moves
  by +0.0003.
- **Production `calculate_temporal_fc`**, run with and without miscalls on `ld_probe_p2000.trees` (6
  draws each): pooled shift **+0.0338**, matching the diagnostic. All 19 rows are finite, and the
  reduced n appear.

**Caveats:**
- **Old constants.** The trees are at μ 4.646e-7, r 2.75e-6, `ancestral_Ne` 6700, with one tree per
  POPMULT. The shift is a sampling effect and should not depend on those, but it has not been
  re-measured at Q = 100.
- **The miscall model is random-allele.** It fits ordinary samples well (§7.2.2A) and does not fit
  the failed ones (§7.2.2C), which is why those are dropped rather than imitated. Rates are
  estimated at MAF ≥ 0.05, the band F_c uses.
- **`diagnostics/temporal_fc.py` does NOT apply miscalls.** It probes the statistic itself, so the
  §7.9.4 / §7.9.6 / §7.9.8 numbers stand as clean-genotype results. **Never compare its F_c levels
  with an empirical target.**
- **No floor yet.** `fc_loss`'s floor has not been measured with miscalls on.

**Before `COMPUTE_FC=1`:**
- ~~(a) Copy the new `fc_common.py` and `CalcTemporalFc.py` to the Beagle machine and re-run; it must
  print `49afd0877028`.~~ **DONE 2026-09-12.**
- ~~(b) Copy the new target in and set `FC_EMPIRICAL_SPEC` to that hash.~~ **DONE 2026-09-12.**
  Exactly the four pairs containing an excluded sample changed (Alsum18 and Alsum59 at t=8, H15
  and H41 at t=16); the other 15 rows are identical to the `ee863fff3bbf` run. H41's excess over
  its pedestal fell 0.0885 → 0.0352. The pooled excess is 0.0283 at t=8 and 0.0268 at t=16, with
  all 19 pairs included.
- (c) Run one live POPMULT=500 trial with `COMPUTE_FC=1`, asserting no empty cells (§10.2).
- (d) Measure the `fc_loss` floor through the production path.
- (e) Run a pilot batch for weights.

F point 3, the missing t=8 vs t=16 difference, is still open.

---

### 7.3 The noise floor — what it is, and why it binds harder here [VERIFIED 2026-08-24]

**The question.** Rejection ABC keeps the draws with the smallest distance `D`. That infers
something only if a draw scores well *because its parameters are right*. But every trial rolls
dice three times — SLiM's forward mating, the recapitation genealogy, and the msprime mutation
overlay — so in truth `D = signal(θ) + luck`. If the luck term is comparable to how much `D`
moves across the parameter range, then "smallest `D`" selects lucky seeds, and the pass returns a
tight, confident-looking posterior built out of coincidence. **Nothing in the output would reveal
that had happened**, which is why this is measured before the pass, not after.

**The floor** is the spread of `D` over replicates at *identical* parameters. It is the resolution
limit of the whole apparatus, and **ε must sit above it** — below that you are slicing finer than
the simulator can resolve.

**Why it binds harder here than in a generic ABC.** The fitted set is only two statistics, and
§7.2.1 showed π is pinned at its flat floor: it excludes small POPMULT but cannot peak. So
**F_st is carrying essentially the whole inference**, and the floor on `fst_loss` specifically is
the number that decides whether the pass means anything. Its signal is `fst_loss` 0.0186 → 0.0083
over POPMULT 2000 → 5000 (§6.1.1), a range of ~0.0103. Roughly:

| run-to-run `fst_loss` spread | reading |
|---|---|
| ~0.0003 | signal ≈ 35× noise — fine |
| ~0.002 | ~20% of the range — usable, ε must be generous |
| ~0.005 | half the POPMULT range is noise — F_st cannot identify N |

**Harness:** `diagnostics/noise_floor.py`, which deliberately routes through the **real**
`AnalyzeTreeSeq.analyze_tree_sequence()` and `ABCAnalysisNoRedis.calculate_losses()`, so it
doubles as the end-to-end integration test §1 wants before any cluster spend.

**MEASURED 2026-08-24 — 3 replicates. `fst_loss` clears decisively; `pi_loss` is marginal.**
Seeds 1001/1002/1003, POPMULT=5000, numClusters=33, total_migration=0.05, μ=4.646e-7,
ancestral_Ne=6700, run locally (20.3–27.6 min each; SLiM 50.1/50.0/62.7 s, `.trees` 206.2–206.3 MB).
Records in `out/noise_floor.jsonl`; pooled with `noise_floor.py --summarize --popmult 5000`.

> **Nei-scale banner (§6.7): every `fst_loss` / `sim F_st` figure in this section predates the
> 2026-08-26 Hudson fix and is ~2× low. Orderings hold, LEVELS do not — never compare one
> against a post-fix number. Not re-measured, by decision. `pi_loss` and `branch_div` are
> unaffected.**
> **The `fst_loss` floor was NOT re-measured before batch 1** (decided 2026-08-26). The working
> assumption is that level and range both roughly double, leaving the 2% verdict intact — one
> measured point is consistent (×1.80, §6.7) but one point is not a floor. **Superseded in
> practice:** §7.4.1 sets the weights from the variance decomposition, which needs no floor at
> all. `pi_loss`'s 30% stands as measured. `ibd_loss` also doubles and is on the old scale.


| statistic | 1001 | 1002 | 1003 | mean | sd | mean\|diff\| | CV% |
|---|---|---|---|---|---|---|---|
| **`pi_loss`** | 0.020863 | 0.024467 | 0.020963 | 0.02210 | 0.00205 | **0.00240** | 9.29 |
| **`fst_loss`** | 0.008556 | 0.008456 | 0.008715 | 0.00858 | 0.00013 | **0.00017** | 1.52 |
| `ibd_loss` | 0.004825 | 0.004691 | 0.004902 | 0.00481 | 0.00011 | 0.00014 | 2.23 |
| `dxy_loss` | 0.000152 | 0.000233 | 0.000153 | 0.00018 | 0.00005 | 0.00005 | 25.84 |
| `genrel_loss` | 0.000115 | 0.000110 | 0.000118 | 0.00011 | 0.00000 | 0.00001 | 3.64 |

**`fst_loss`: run-to-run 0.00017 against an across-prior range of 0.01026 — 2%. The pass is
viable.** Signal sits ~60x above noise, inside the "fine" band this section set out in advance,
and nowhere near the ~0.005 that would have meant F_st cannot identify N. ε has room to be tight.
This is the number that mattered: §7.2.1 established that F_st carries essentially the whole
inference, and all three replicates re-drew the forward genealogy through SLiM, which is exactly
where F_st's run-to-run variance lives and precisely what `--fixed-tree` cannot capture (§3).

**`pi_loss`: run-to-run 0.00240 against a range of 0.00810 — 30%. Marginal, and it is §7.2.1
restated from the other direction.** π saturates at its flat floor above POPMULT ≈ 5000, so its
usable range is small and coalescent noise eats a large fraction of it. Note **two of the three
replicates (0.020863, 0.020963) sit within 0.5% of the §7.2.1 flat-simulation floor of 0.02094**,
and 1002's 0.024467 is an excursion *upward* — the expected geometry, since L1 never penalises a
too-flat simulation, so π can only leave the floor by bad luck.

**Consequence for `abc_standardize.py` — now measured, not predicted.** Equal weights after
MAD-standardization would put a statistic that is 30% noise at parity with one that is 2% noise.
TODO §4 anticipated that π would "warrant *less* than half"; this floor is the direct evidence.
**Set unequal `WEIGHTS` once the pilot batch exists, and inspect the per-statistic MADs before
accepting the default.** Keep π fitted — its lower bound on POPMULT is real (§7.2.1) — but do not
let it vote equally with F_st.

**Caveat, stated honestly.** Three replicates at ONE point in parameter space (POPMULT=5000). The
floor is not guaranteed constant across the prior; F_st's variance in particular could differ at
POPMULT=2000, where the level is 2.2x higher. Nothing here has been measured at the prior's edges.

---

### 7.4 Batch 1 landed — the pass ran, and it identified the KERNEL, not N [VERIFIED 2026-09-06]

500 jobs × 5 trials submitted 2026-08-26 on the code as of `929f1c2` (Hudson F_st, §6.7;
`POPMULT_MAX = 25000`). Collected with `diagnostics/collect_batch.py`; pooled results in
`out/batch1/abc_results.csv`, ranked in `out/batch1/abc_results_ranked.csv`, frozen σ in
`out/batch1/abc_sigmas.json`. (Moved out of `out/` root 2026-09-09 — see §9's layout note.)

**Batch health: clean. [VERIFIED]** 2,495 of 2,500 rows. All 499 present files hold exactly 5 rows
— no short files, so no per-trial failures at all. **Job 254 is absent entirely** (no file): a
whole job died or its transfer failed, and because a lost job records no rows it records no
parameters either, so an OOM at high `pop` cannot be distinguished from a random eviction without
the HTCondor log. At 1/500 nothing turns on it.

**`pop` coverage is uniform across all ten deciles** (max |z| = 1.8; top decile z = +0.54).
**This closes TODO §4's "validate memory at POPMULT=25000"** — the §3.1 extrapolation held and
jobs ran at the ceiling. There is no truncated upper tail, so the pool is not biased toward the
small draws that fit.

#### 7.4.1 The weights: π's spread is majority μ-nuisance, so it gets 0.125

The §7.3 rule (σ against the replicate noise floor) was applied first and **gives near-equal
weights, 0.488/0.512** — contradicting §7.3's own prediction. That rule is wrong here, and why is
the actual finding. A rank-space variance decomposition of each loss onto all five parameters
(`collect_batch.py` §7, unique R², n=2,495):

| loss | R²_total | **demographic** | **μ-nuisance** | unexplained |
|---|---|---|---|---|
| `pi_loss` | 0.240 | **0.089** | **0.153** | 0.760 |
| `fst_loss` | 0.622 | **0.624** | 0.000 | 0.378 |
| `ibd_loss` | 0.846 | 0.851 | 0.000 | 0.154 |
| `dxy_loss` | 0.043 | 0.033 | 0.010 | 0.957 |
| `genrel_loss` | 0.640 | 0.645 | 0.000 | 0.360 |

("demographic" = unique R² of `pop` + `total_migration` + `m` + `numClusters`.)

**`pi_loss`'s batch spread is majority μ-draw** — 63% of the variance the parameters explain comes
from `mutation_rate`, only 10% from `pop`. μ is a *deliberate* nuisance dimension (§7.2 priors
table; only θ=4Nμ is ever reported), so weighting on that variance puts the μ draw into `D`. The
noise-floor rule cannot see this, because μ-driven spread is not replicate noise. This is TODO
§1's `s=0.02` concern surviving into the batch: even at the tightened prior, **μ is still the
single largest driver of `pi_loss`.**

`fst_loss` takes **0.000** from μ (ρ = +0.011, n.s.) — the empirical confirmation, in the real
pipeline, of §5.3's claim that F_st is denominator-invariant.

Set in `abc_standardize.py`: **`WEIGHTS = {"pi_loss": 0.125, "fst_loss": 0.875}`**, i.e. weight by
the demographic share. Frozen σ: `pi_loss` 0.0101535, `fst_loss` 0.00316778. These weights need no
noise floor, which also sidesteps the deferred Nei→Hudson floor caveat (§7.3 banner, TODO §4).

> **This rule replicated on an independent batch** (the LD pilot, 600 trials): `pi_loss` still
> loaded more on `mutation_rate` (0.149) than on `pop` (0.080), so the noise-floor rule would have
> over-weighted it 3×. **Keep the noise floor as a VETO, not a weight** — if a fitted statistic's
> floor/σ ever approaches 1, drop the statistic rather than down-weight it. Demographic R² only
> apportions among statistics already judged admissible on other evidence: it cannot tell useful
> signal from useless, which is why `ibd_loss` tops the table below and must still stay out.

> **Do NOT promote `ibd_loss` on the strength of its 0.851.** It is the highest demographic R² in
> the table and it is ~90% `m`. The observed IBD slope is statistically zero in all three years
> (§7.1), so fitting it would measure distance from a noise target. Same for `genrel_loss` (~92%
> `m`). Both stay diagnostic — but they are now well-characterised posterior-predictive checks
> **for `m` specifically**, which is more than was known before.

#### 7.4.2 The result that matters: N is unidentified, `m` is sharply inferred

Posterior IQR as a fraction of prior IQR — 1.0 means the data said nothing:

| acceptance | n | `pop` median | **`pop` IQR ratio** | `m` median | **`m` IQR ratio** | `total_migration` IQR ratio |
|---|---|---|---|---|---|---|
| top 20% | 499 | 14299 | **0.97** | 4.26e-05 | 0.23 | 0.76 |
| top 10% | 249 | 13910 | **1.02** | 4.79e-05 | 0.24 | 0.74 |
| top 5% | 124 | 12572 | **1.01** | 5.23e-05 | 0.22 | 0.67 |
| top 2% | 49 | 12305 | 0.82 | 6.92e-05 | 0.12 | 0.65 |
| top 1% | 24 | 12647 | **0.82** | 6.99e-05 | **0.05** | 0.73 |

**`pop` never tightens** — at the top 1% its spread is still 82% of the prior's. **`m` collapses
20×.** The cleanest single view is the top 20 runs by `D`: **`fst_loss` sd across them is
0.000033** — they are hitting the *identical* F_st — while `pop` spans **2,959 to 24,938**
(essentially the whole prior) and `m` sits in a narrow band 5.7e-5–8.9e-5.

**The pass pinned F_st with the dispersal kernel and left population size free.** It inferred the
one parameter §7.1 says must not be reported as inferred, and learned nothing about the one the
project exists to estimate.

**This is a PARAMETERIZATION failure, not a distance-weighting failure.** F_st ≈ `1/(1+4Nm)`, so
with `m` free over four orders of magnitude, N and the kernel are confounded. The compensation is
visible inside the accepted set: **spearman(log `pop`, `total_migration`) = −0.336** — bigger N,
less migration, same F_st. The 0.125/0.875 weights are still correct on the evidence (equal
weights would be worse); they simply cannot rescue a confounded parameterization.

**Why `m` wins.** Its prior `lognorm(s=1.5, scale=1e-4)` spans 6.1e-7 to 1.4e-2. The kernel runs
on **cluster centroids** (`data/cluster_distances.csv`, 33 clusters: min 6.5 km, **median 44.1 km**,
max 173.7 km — not the 1.7–160 km *site* spread quoted in §7.1), so `1/m` is an e-folding distance
in metres. Across the prior that runs from 1,639 km (kernel flat — a pure island model) to 0.07 km
(all migration to the single nearest cluster). Against that, `pop`'s 12.5× range has almost no
leverage on F_st.

Measured kernel behaviour, per `1/m` (mean over the 33 demes; max possible eff. sources = 32):

| `m` | `1/m` | eff. source demes | weight on nearest | regime |
|---|---|---|---|---|
| 6.1e-7 | 1639 km | 32.0 | 0.032 | island model — **`m` has no effect** |
| 5e-6 | 200 km | 31.4 | 0.038 | island model — **`m` has no effect** |
| 2e-5 | 50 km | 26.2 | 0.065 | informative |
| **4.3e-5** | **23 km** | 18.5 | 0.121 | informative — *batch-1 posterior* |
| **7e-5** | **14 km** | 12.9 | 0.186 | informative — *batch-1 posterior* |
| 1e-4 | 10 km | 9.2 | 0.250 | informative |
| 1e-3 | 1.0 km | 1.4 | 0.854 | informative |
| 1.4e-2 | 0.07 km | 1.0 | 0.992 | nearest-neighbour only — **`m` has no effect** |

**The prior is not badly placed** — 85.4% of its mass lands in the informative band ≈[6e-6, 6e-4],
with 3.0% below and 11.6% above. The problem is not a wasteful prior; it is that inside that band
`m` simply has more leverage on F_st than `pop` does.

> **Note the discretisation limit this exposes.** The nearest two cluster centroids are 6.5 km
> apart, so **no dispersal structure below ~6.5 km is representable at all** — any `1/m` shorter
> than that collapses to "all migrants from the one nearest cluster" and the model stops
> distinguishing them. Real CPB dispersal is largely sub-kilometre. The kernel is therefore
> describing between-cluster exchange, not beetle movement, and `1/m` must not be read as a
> dispersal distance in the biological sense. (`numClusters` ∈ {1,2,3}×33 changes the centroid
> spacing, so this floor moves with it.) Inside `fst_loss`, unique R² for `m` is **0.429** against **0.091**
for `pop`: the kernel is ~4.7× more influential than population size in the statistic carrying the
whole inference.

**Two things that did go right.** The posterior *centre* is consistent with §6.7 — median `pop` ≈
12,300–12,600 at tight acceptance against §6.7's predicted 12,469 — so raising the ceiling was
right; there is simply no peak under it. And the `fst_loss` gradient in `pop` flattens from decile
~8 onward with deciles 8/9/10 all within 1 SE of the minimum, so **`POPMULT_MAX = 25000` is not
truncating** — the §6.7 failure at 12000 does not repeat.

**Residual misfit worth recording.** Best achievable `fst_loss` is **0.0054**, ~18× the (Nei-scale)
noise floor. The model cannot match observed F_st even at its best, consistent with §7.2's
suspicion that an exponential distance-decay kernel structurally cannot produce "one isolate, rest
near-panmictic."

**What this licenses saying, and what it does not.** This batch constrains **Nm**, not N. Reporting
an N posterior from it would be reporting the prior. Separating N from m needs an external
constraint on dispersal, which these data cannot supply — that is the professor conversation in
TODO §3, not an ABC-tuning problem.

> **Consequence, and it is the new blocker: `m` must be fixed before batch 2.** TODO §3's "decide
> whether `scale` should stay a free parameter" is no longer a prior-volume economy — it is the
> difference between inferring N and not. §6.2 measured the identifiability that is missing here
> (POPMULT 500→5000 improving `fst_loss` 8.7×) **at a fixed kernel and fixed `total_migration`**;
> freeing `m` is what destroyed it. Fixing `m` alone is not sufficient either — `total_migration`
> is free over a 300× range and trades against N through the same Nm. More draws under this
> parameterization buy resolution on `m` and none on N.

#### 7.4.3 Pinning the kernel does NOT rescue N — it converts it into an assumption [VERIFIED 2026-09-06]

Batch 1 answers this for free, with no new simulation: condition on narrow slices of `m` and
`total_migration` and see what `pop` does. Scored with the real weights (0.125/0.875) and frozen σ.

| slice | n | top-20% `pop` median | `pop` IQR ratio | Nm median | Nm IQR ratio |
|---|---|---|---|---|---|
| *(unconditioned)* | 2495 | 14299 | **0.97** | 106 | 0.77 |
| `m` ∈ [3e-5, 5e-5] | 285 | 12451 | 0.77 | 85.9 | 0.26 |
| `m` ∈ [5e-5, 8e-5] | 303 | 13194 | 0.85 | 91.2 | 0.31 |
| `m` ∈ [8e-5, 1e-4] | 439 | 17965 | 0.67 | 172.5 | 1.00 |
| `m` ∈ [1e-4, 3e-4] | 411 | 21051 | 0.50 | 285.0 | 1.63 |
| `m`[3e-5,8e-5] × `tm`[0.02,0.10] | 135 | **19680** | 0.39 | 78.3 | 0.18 |
| `m`[3e-5,8e-5] × `tm`[0.10,0.20] | 214 | **12846** | 0.68 | 86.1 | 0.17 |
| `m`[3e-5,8e-5] × `tm`[0.20,0.31] | 197 | **7852** | 0.59 | 95.5 | 0.42 |
| `m`[3e-5,8e-5] × `tm`[0.02,0.20] × nC=1 | 132 | 9650 | 0.80 | 103.4 | 0.32 |
| `m`[3e-5,8e-5] × `tm`[0.02,0.20] × nC=2 | 115 | 13902 | 0.73 | 81.2 | 0.22 |
| `m`[3e-5,8e-5] × `tm`[0.02,0.20] × nC=3 | 102 | 16740 | 0.79 | 81.3 | 0.14 |

**Two things, and the second is the important one.**

1. **Pinning does help.** `pop`'s IQR ratio falls from 0.97 to 0.39–0.85. So conditioning on the
   nuisance parameters does recover *some* resolution on N — fixing `m` is worth doing.

2. **But the inferred N is a direct function of what you pin it to.** Holding `m` fixed and moving
   only `total_migration` across its prior swings median `pop` **19,680 → 12,846 → 7,852, a 2.5×
   range**. Moving `m` alone swings it 12,451 → 21,051 (1.7×). `numClusters` adds another 1.7×
   (9,650 → 16,740). **These are not error bars; they are the answer changing with the assumption.**

Meanwhile **Nm is stable at ≈78–103 across every well-pinned slice**, and its IQR ratio tightens to
0.14–0.32 — far tighter than `pop` ever gets. (The high-`m` slices are the exception, running
Nm 172–285; those are the slices where the kernel is pulling toward nearest-neighbour exchange.)
For reference §7.1's anchor from the observed F_st level is **Nm ≈ 31–83**, so ≈80–100 sits just
above it — consistent, not contradictory.

**This is the honest statement of what the data supports: Nm is identified; N is whatever the
dispersal assumption makes it.** Fixing `m` and `total_migration` does not turn an unidentified
posterior into an inferred N — it turns it into a point estimate whose value is set by an
assumption the genetics cannot check. Two consequences bind on batch 2:

- **The pinned values must be defensible on external grounds, not chosen.** Picking them to land
  on a plausible N would be circular, and given the 2.5× lever it is very easy to do accidentally.
- **A sensitivity analysis over the assumption is mandatory, not optional.** Report N *as a
  function of* assumed (`m`, `total_migration`) — the surface above is the prototype. A single
  fixed pin would report one column of that surface and hide the other four.

This is what makes the grid design in TODO §4 the right shape: it costs the same as one fixed-pin
batch and it produces the sensitivity surface as its primary output rather than as an afterthought.

---

### 7.5 LD decay — BUILT, MEASURED, and RETIRED. Do not restart it [closed 2026-09-08]

`ld_loss` is wired end to end — `ld_common.py` imported by both sides, `calculate_ld_decay` in
`AnalyzeTreeSeq` (`COMPUTE_LD` default ON, `LD_THIN = 25`), `_ld_pooled_mean_abs_diff` and
`LD_MIN_BIN` in `ABCAnalysisNoRedis` — and it is the most `pop`-informative statistic ever measured
in this project: unique R² **0.258** on `pop`, against **0.079** for `fst_loss`. **It is still not
in `abc_standardize.py::FITTED_STATS`, and that is a closed decision, not a pending step.**

**The reason is structural, and it is that every chromosome has a different recombination rate.**
The simulation applies one uniform `r` to one 1e6 bp locus. The empirical target pools 17
chromosomes whose *effective* `r` spans **≥10×** — per-chromosome decay distance runs 1,985 →
26,063 bp, with the same chromosome ranking in all three years (chr6, chr5, chr15 slowest), which
is what proves it is a property of the genome and not of the beetles. Because the LD curve is
convex, a mixture's far field sits **above** any single-`r` curve, so no effective average `r`
exists: a mixture changes the curve's *shape*, not merely its position. Measured, the sim–obs gap
this creates (0.0027–0.0039) is the **same size as `ld_loss` at its own minimum** (0.0032) — so
essentially the entire quantity being minimised is a known mis-specification, and the minimum's
location (POPMULT ≈ 2000) is set by that artifact rather than by population size.

> **CONTRADICTED 2026-09-17 by §7.5.6 below.** Scoring against a single chromosome — any of the
> 17, including the fastest-decaying — puts the minimum in the same place. The mixture is real,
> but it is NOT what sets the minimum's location. The retirement stands for a different reason.

Everything tried against it failed, and each failure is recorded in `OLD_LOGS.md` §A:

- **Pinning `r`** (§6.8.1's linkage-map value) was tested end to end. `4Nr` predicts the LD minimum
  moves 3.44×; it moved **1.33×**. Closing the LD-vs-F_st gap that way would need `r` wrong by
  ~10⁴×, which the linkage map excludes.
- **Freeing `r`** cannot work in principle: LD constrains `4Nr`, so a free `r` adds an unknown
  without adding an equation, and Sved's form puts N and `r` in a product at *every* distance —
  there is nothing in the shape to separate them.
- **Fitting per chromosome** needs 17 unknown `r_k` estimated from the very curves being fitted.
- **Selecting chromosomes near the simulated `r`** is choosing data by how well it matches the
  model, using the quantity being fitted.
- **Adding `ld_loss` to `D` anyway** is a surrender, not a compromise: at its derived weight the
  accepted set's median `fst_loss` equals the whole-batch median — F_st stops contributing — and
  `pop` tightens because LD is monotone across the prior, not because two statistics agree.

**What survives and still matters.** `r` is genuinely pinned from a linkage map with no `N_e` in it
(§6.8.1). The per-chromosome recombination heterogeneity is a real measured property of the CPB
genome that was not known here before, and is worth reporting on its own. The empirical run
(`data/empiricalStats/averaged_ldDecay_*.csv`, `data/ld_per_chr/`, 7.7 h on the Beagle machine) and
its per-chromosome files are kept — **do not delete them**; they are the only route to a
between-chromosome jackknife or a re-pool without another 7.7 h run. And `ld_common.spec_hash()`
(`4d1d1d92b25b`) still guards both sides, so nothing may edit the binning spec casually.

**Beagle imputation inflating empirical LD is still [OPEN]**, and cannot be settled with the files
in hand — the Beagle inputs *are* the base data, hard-called with no genotype probabilities and no
missingness marker, so nothing records what was imputed. ~~The recombination finding above removes
the *need* to invoke it, which moves it well down the list.~~ **Back on the list — §7.5.6.**

#### 7.5.6 One chromosome does not rescue LD — the misfit is a LEVEL shared by every chromosome [VERIFIED 2026-09-17]

`diagnostics/ld_one_chrom.py`, no simulation: re-scores all 2,500 batch-3 trials (their saved
`ld_{year}.csv`) against one chromosome's curve from `data/ld_per_chr/`, and the ten
`out/fc_loss_floor/` replicates for a floor. Record in `out/ld_one_chrom.jsonl`.

**The chromosome is chosen blind to the simulation.** All chromosomes share one Ne, so decay
distance ∝ 1/r_k; the chromosome at the MEDIAN decay distance of the 17 has the most typical r,
which is what the model's linkage-map genome-average r stands for. That is **chr2** (~2,900 bp,
geometric mean over years). This is not the rejected "choose chromosomes near the simulated r"
above: the rule compares chromosomes only with each other. `--all` prints every chromosome as a
sensitivity display and must never be used to choose one. Caveat: the median is unweighted by
chromosome length (not in the repo), while the map r is length-weighted.

**Cross-checks passed.** Re-scoring against the genome-wide target reproduces the batch's own
`ld_loss` column to **7e-18**; the 2015 decay-distance range reproduces 1,985 → 26,063 bp.

Single-statistic rejection, median accepted `pop`:

| statistic | top 5% | top 1% | IQR ratio at top 1% | floor/signal |
|---|---|---|---|---|
| LD, 17-chromosome pool | 4,710 | 3,136 | 0.21 | 0.14 |
| **LD, chr2 only** | **4,633** | **3,136** | 0.23 | 0.14 |
| F_c alone | 18,381 | 16,718 | 0.68 | — |
| F_st alone | 13,184 | 10,297 | 0.76 | — |

- **The chromosome barely matters.** Over all 17 the top-5% median runs **4,600–5,900**, including
  chr6, which decays 13× slower than chr17. spearman(decay distance, preferred `pop`) = −0.63 is
  the r–N trade-off in the expected direction, but tiny against a 13× spread in decay.
- **Noise is not the problem** (floor/signal 0.14). **The problem is the answer:** every version's
  lowest-loss `pop` decile is the FIRST (2,000–4,400), i.e. LD is pressed against the prior floor,
  3–5× below F_c and F_st.
- **So the mixture does not set the minimum.** The fitted bins (≥ `LD_MIN_BIN` = 1,778 bp) sit
  higher on the observed side than the simulated side on EVERY chromosome, and smaller N adds
  forward-phase LD, so the loss rewards smaller N regardless of chromosome.

**Reading [INFERRED]:** a raised long-range LD level common to all chromosomes is what imputation
would produce, and the old between-chromosome Test 3 (`ld_chrom_check.py`) is blind to a uniform
shift by construction — it only compares chromosomes with each other. So **imputation is a
candidate again**, and still cannot be settled with the files in hand. A genuinely small recent N
is the other reading, and both F_c and F_st argue against it.

**LD stays retired, for a sharper reason:** a level mismatch shared by every chromosome, not the
recombination mixture. Do not re-open it via chromosome choice — this closes that route.

---


---

# §J — the information budget, temporal F_c, batch 3 and the re-founding model (old §7.9–7.9.12) [ARCHIVED 2026-09-27]

CLAUDE.md §7.9 keeps the current state. Every N-identifiability claim about F_c below is conditional on the persistent-deme model (see §7.9.12F).

### 7.9 The information budget — WHY N is unidentifiable, mechanically [VERIFIED 2026-09-09]

§7.4.2 measured that N is unidentified and §7.5 exhausted the recombination route. This section
answers the prior question — *why* — and it turns out to be derivable rather than empirical. The
answer changes what a replacement statistic has to look like, so **read this before proposing
one.**

#### 7.9.1 Only 0.5–7% of the genealogy is POPMULT-dependent, and that sliver is itself `Nm`

Recapitation hands the deep genealogy to a **fixed** `ancestral_Ne = 6700`. So POPMULT can only
act on pairs that coalesce inside the 324-generation forward window. That fraction is directly
readable from `branch_div`, whose hard ceiling is `2·(324 + 2·6700) = 27448` (§6.1.1): a pair
coalescing forward saves ~`2·Ne_anc` of pair-time, so

```
frac coalescing forward  ~=  (ceiling - branch_div) / (2 * 2 * Ne_anc)
```

From `out/mu_calibration.jsonl`, plus §6.1.1's saturation fit at the prior ceiling:

| POPMULT | deme N | `branch_div` | deficit | **frac coalescing in the forward window** |
|---|---|---|---|---|
| 500 | 50 | 21647 | 5801 | 21.7% |
| **2000** (prior floor) | 202 | 25629 | 1819 | **6.8%** |
| 5000 | 505 | 26847 | 601 | 2.2% |
| **25000** (prior ceiling) | 2525 | 27312 (fitted) | 136 | **0.5%** |

**Across the production prior, between 6.8% and 0.5% of the genealogy is POPMULT-dependent. That
is the total information budget for N, for ANY statistic computed from the tree sequence.** It is
not a property of π or F_st; it is a ceiling on all of them, and it explains why every statistic
tried so far has behaved the same way.

**And the sliver is not even an N quantity — it is an `Nm` quantity.** The strict scattering-phase
bound is `1/(1 + 4·N_deme·m)`, the probability two lineages in one deme coalesce before *either*
migrates out (coalescence rate `1/2N` against total migration rate `2m`), at
`m = total_migration = 0.05`:

| POPMULT | deme N | `4Nm` | `1/(1+4Nm)` bound | measured | ratio |
|---|---|---|---|---|---|
| 500 | 50 | 10.0 | 9.1% | 21.7% | 2.4 |
| 2000 | 202 | 40.4 | 2.4% | 6.8% | 2.8 |
| 5000 | 505 | 101.0 | 1.0% | 2.2% | 2.2 |

The measured fraction sits at a **roughly constant 2.2–2.8× the bound** across a 10× range in
POPMULT. The multiplier is pairs that separate and later *re-meet* inside the 324-generation
window, which the strict bound drops by construction. **The point is the proportionality, not the
constant: the forward window's coalescence fraction tracks `1/(1+4Nm)`, so it carries `Nm`, not
N.**

**This is Wakeley's separation of timescales** (§11), and stating it in those terms is what makes
§7.4.2 predictable rather than surprising. A structured genealogy has a fast **scattering phase**,
whose only parameter is `Nm`, and a slow **collecting phase** — a Kingman coalescent on the whole
metapopulation — which is where N appears *on its own*. The collecting-phase timescale here is
`~2·D·N ≈ 33,000` generations at POPMULT=5000. **The forward run is 324.** So the forward phase is
essentially pure scattering, and the entire collecting phase has been replaced by a constant.

**Consequences, and they are the useful part:**

- **N is not hiding in a statistic nobody tried.** The simulator does not express it. Any summary
  of that tree inherits the confound, which is why d_xy, relatedness, IBD and LD all failed in the
  same direction.
- **The filter on any candidate statistic is now specific**, not "is it informative": (a) does its
  signal live inside the forward window rather than being diluted by the ancestral phase, and
  (b) does it read an **absolute rate** rather than a coalescence *probability*? Probabilities of
  two lineages meeting are always `Nm`. A rate measured over a known number of generations is N.
- **The upper half of the prior is nearly informationless.** At POPMULT 25000 only 0.5% of
  coalescences are POPMULT-dependent. **Do not "fix" this by lowering `POPMULT_MAX`** — that would
  be choosing the prior to sharpen the answer, exactly the circularity §7.4.3 warns about. Report
  the budget curve as a diagnostic instead.
- **It is a MODEL-STRUCTURE result, and that makes it publishable.** "Only 0.5–7% of the genealogy
  is sensitive to the parameter we are inferring, and that fraction is itself a function of `Nm`"
  is a far stronger thing to report than a list of statistics that did not work.

> **What would actually change it** is bringing the collecting phase inside the forward run —
> `G ≳ 2·D·N`, i.e. ~33,000 generations against 324. That is the same open idea §6.1.3 costed for
> LD, now with a second and independent motivation. Still a professor question, because 325
> generations is the invasion timescale.

#### 7.9.2 A resampled field is a DIFFERENT DEME in different years [FOUND 2026-09-09]

Found while building §7.9.4's probe, and independent of it.
`GenerateClusterData.assign_genomes_to_clusters_idv_year` is **greedy in specifier-row order and
one-to-one per year** — `if cluster.genome_assignments[yearIdx] is not None: continue`, so each
site takes its nearest cluster *that no earlier site of that year has already claimed*.

**Measured** (33 clusters, `data/cluster_data.csv` as of this writing):

| year | sites NOT on their nearest cluster | median site→deme distance | if assigned to NEAREST | optimal one-to-one (Hungarian) |
|---|---|---|---|---|
| 2015 | 16/24 (67%) | 10.74 km | 3.23 km | 7.41 km |
| 2019 | 12/17 (71%) | 11.31 km | 4.22 km | 8.75 km |
| 2023 | 12/20 (60%) | 8.01 km | 3.31 km | 7.48 km |

Worst single case: `Paramount652-2015` is represented by a deme **26.7 km** away when a cluster
sits **4.8 km** from it.

**Decompose the cost, because it decides what to do about it.** Most of the displacement is the
*one-to-one constraint itself* (median 3.2 → 7.4 km), not the greedy ordering (7.4 → 10.7 km).
One-to-one is defensible: without it two sites could share a deme and would then have *identical*
simulated statistics, which is worse for element-wise fitting. **Greedy ordering is not** —
`scipy.optimize.linear_sum_assignment` gives the optimal one-to-one matching at no cost and cuts
total displacement 12–18%.

**Two separate consequences:**

1. **Geography is blurred on the simulated side.** The kernel runs on centroid distances and the
   IBD regression uses real site coordinates for *both* sides (§7.2), so a median 10.7 km
   site→deme error is injected noise in exactly the quantity `m` is supposed to explain. A
   candidate partial explanation for weak simulated IBD, and not covered by any existing check.
   **[OPEN]**
2. **There is no persistent field across years.** Arlington is cluster 10 in 2015 and cluster 21
   in 2019. Anything comparing the same field across years — §7.9.4 above all — must choose its
   own deme mapping and cannot use the production one.

**FIXED 2026-09-09 (Sohan's call): `assign_genomes_to_clusters_idv_year` now uses the optimal
one-to-one matching** (`scipy.optimize.linear_sum_assignment`, the Hungarian algorithm). Same
one-site-per-cluster constraint, globally minimal total displacement, no dependence on row order.
Verified against the old greedy path on all three years: every site still assigned exactly once,
no cluster takes two, and total displacement falls **18.0% / 14.9% / 12.3%** (2015/2019/2023).
It breaks comparability with batches 1 and 2, which is accepted — it ships alongside the Q = 100
rescaling (§7.9.3), which breaks it anyway.

**What it does NOT fix.** Each year is still matched independently, so a field can still land in
different demes in different years when its neighbours differ; anything comparing a field across
years must pick its own deme (`diagnostics/temporal_fc.py --deme-choice`). And the one-to-one
constraint itself still costs more than the greediness did — median 3.2 → 7.4 km, against 7.4 →
10.7 km for the ordering. **The real cure is more clusters**, which makes a good matching easy;
see §7.9.5.

#### 7.9.3 `mutation_rate` is now FIXED, not drawn [DONE 2026-09-09]

`prior_distributions` no longer contains `mutation_rate`; `sample_prior()` returns
`DEFAULT_MUTATION_RATE = 4.646e-7` as a constant. The CSV column stays, so the output layout is
unchanged and every row still records the scale it ran at.

**The justification is that there was never an unknown to draw.** μ is *defined* by
`4·Ne_anc·μ = π_obs` with both other terms fixed — `Ne_anc = 6700` by fiat (§6.1), `π_obs = 0.0122`
measured (§5.1) — so it is a deterministic function of an observation and a chosen constant. A
constant like `ancestral_Ne`, not a parameter.

**§6.1.2's trio rate is what makes this external rather than self-referential.** Xu et al.
(2026) measured `μ_true = 5.8e-9`; §6.1.3 puts the model at Q ≈ 78.5; a consistent rescaling needs
`μ_sim = μ_true·Q = 4.55e-7` against the calibrated **4.646e-7**, i.e. **−2.0%**. The trio rate
does not replace the constant, it **confirms** it — the first support that value has had beyond
"it makes π come out right."

**Measured payoff, from the pilot batch's own decomposition:** `pi_loss`'s single largest driver was the μ
draw (unique R² **0.149**) against **0.080** for `pop`. Deleting the draw deletes that nuisance
variance outright. π becomes a demographically clean statistic for the first time, and its weight
can rise on demographic signal instead of being held to 0.110 to keep the μ draw out of `D`.
**Re-derive `WEIGHTS` from batch 3 rather than reusing §7.4.1's or the pilot's** — both were fitted
with μ free and are not valid under a fixed μ.

**Two traps, stated because both are the obvious next move:**
- **Do NOT set `DEFAULT_MUTATION_RATE = 5.8e-9.** That is the *unrescaled* rate — correct for a
  beetle, 78× too small for a model running at 1/78 scale. Same error class as §6.8.1's `r`.
- **Do NOT widen the prior to the trio CI** (4.7–7.2e-9, log-sd ≈ 0.11, **5× the old s=0.02**).
  That uncertainty belongs to `Ne_true = π_obs/(4·μ_true)`, which is *reported*, not simulated.
  Propagating it into μ_sim would make `pi_loss` noisier for nothing.

> **`collect_batch.py` was made robust to this in the same commit.** `_zrank` divides by the rank
> standard deviation, so a fixed parameter makes the design matrix singular; the decomposition now
> drops zero-variance parameters and **prints which**, so a *silently* fixed parameter (a bug)
> stays as visible as a deliberately fixed one. Verified behaviour-preserving — re-run on batch 2
> it reproduces the pilot's table exactly.

> **Where the trio rate DOES earn its keep is reporting.** It converts the fixed constants into
> biology — `Ne_anc,true = π_obs/(4·μ_true) ≈ 5.3e5`, with a real CI — and it pins Q. Note the
> scale question this raises and settles: Hudson F_st is `(E[T_b] − E[T_w])/E[T_b]`, so `Ne_anc`
> cancels top and bottom, which is *why* §6.2.1 measured F_st moving <1% for a 3× change in
> `ancestral_Ne`. **Therefore `Nm ≈ 78–103` is a REAL-WORLD number and no factor of Q should ever
> be applied to it.** Easy to get wrong at write-up time.

#### 7.9.4 Temporal F_c — reads N exactly [MEASURED 2026-09-09]

> **The "does it break the confound" half of this section was answered later the same day: YES, on the axis that matters. See §7.9.6.** What follows tested `total_migration`, which §7.4.2 shows is NOT the knob that absorbed F_st in batch 1. The elasticities below are correct and worth keeping — they are what shows F_st is not a pure `Nm` statistic — but the identifiability verdict in point 4 is superseded.

The statistic §7.9.1's filter asks for: signal entirely inside the forward window, and an
**absolute rate** rather than a coalescence probability.

**What it is.** 19 field-pairs are resampled at *identical coordinates* across years — 9 at t=16
generations (2015 vs 2023), 10 at t=8 — nearly all n=7 vs n=7, after the §7.0 n<4 mask drops the
`Arlington2015` typo duplicate. Waples (1989) `F_c = mean_loci (x−y)²/((x+y)/2 − xy)` on
common variants (pooled MAF ≥ 0.05). Harness: `diagnostics/temporal_fc.py`.

**Why it should escape the budget.** The ancestral phase is *shared history* between two
timepoints, so it cancels exactly from a frequency **change**: F_c is 100% forward-phase, against
0.5–7% for π and F_st. It is also μ-free and denominator-free (§5.1's problem does not arise), and
r-free in expectation — r sets only the number of effectively independent loci, i.e. precision.

**Why §11's rejection of the temporal method does not apply.** That rejection was that at n=5–8
the sampling correction `1/(2S_a)+1/(2S_b) ≈ 0.14` is 11–21× the drift signal. Decisive if you are
*correcting* it to estimate Ne. In ABC you **reproduce** it — subsample the simulated demes to the
same n_i (invariant 1) and the identical offset appears on both sides and cancels. Only precision
survives as an objection, and precision is bought with loci (Waples: doubling S, t or L buys about
the same).

**Measured**, three points, seed 1, 5 subsample draws each, all on trees that already existed or
cost one SLiM run. `n_big` = every deme cut to 50 diploids so the pedestal is 0.02 instead of 0.14:

| pt | POPMULT | `total_migration` | deme N | Nm | F_st | F_c(t=16) | pedestal | **F_c − pedestal** |
|---|---|---|---|---|---|---|---|---|
| A | 2000 | 0.05 | 202 | 10.1 | 0.04659 | 0.03848 | 0.02 | **0.01848** |
| B | 5000 | 0.05 | 505 | 25.2 | 0.01952 | 0.02742 | 0.02 | **0.00742** |
| C | 5000 | 0.02 | 505 | 10.1 | 0.02190 | 0.02920 | 0.02 | **0.00920** |

**1. F_c reads deme size essentially exactly. PASSES, and this is the real result.** A→B is a
2.50× change in deme N and the drift term moves 0.01848 → 0.00742, a ratio of **2.49**. As an
elasticity, `d ln(F_c−ped)/d ln N = −0.996` against a theoretical **−1.000**. At POPMULT=2000 the
raw value 0.0388 sits 2% from the textbook `t/(2N) = 0.0396`. **No other statistic in this project
tracks N that cleanly** — and it does it while sitting outside the §7.9.1 budget.

**2. At the REAL sample sizes it is marginal, as predicted.** With n=7 the pedestal is ~0.143 and
F_c moves only 0.16197 → 0.15480 across A→B, i.e. **4.4% of the level, 1.9 sd of subsample noise**.
Across the full prior (POPMULT 2000→25000, 12.5×) the drift term would run 0.0185 → 0.0015, a
range of ~0.017 against a per-draw sd of ~0.0037 — roughly 5σ, and usable. Note the pedestal
cancels exactly in `|F_c,sim − F_c,obs|`, so it costs nothing but the noise it carries.

**3. F_st is NOT a pure `Nm` statistic in this model, and that CORRECTS a standing claim.**
Elasticities from the same three points:

```
d ln F_st        / d ln N = -0.949      d ln F_st        / d ln m = -0.126
d ln (F_c - ped) / d ln N = -0.996      d ln (F_c - ped) / d ln m = -0.235
```

so the iso-lines are `N·m^0.132 = const` for F_st and `N·m^0.236 = const` for F_c. **A pure `Nm`
confound would put the exponent at 1.000 for F_st.** It is 0.13 — F_st barely responds to
`total_migration` here (+12% for a 2.5× *reduction*, where `1/(1+4Nm)` predicts +147%).

§6.2's "F_st tracks `1/(1+4Nm)` almost exactly" is not contradicted so much as **shown to have
been untestable**: that sweep varied only POPMULT at fixed migration, and over that path
`1/(1+4Nm)` and `1/N` are indistinguishable. This is the first measurement in the project that
moves migration and population size separately.

**4. THE DECISIVE TEST WAS OUTSTANDING WHEN THIS WAS WRITTEN — it has since been run (§7.9.6).** The design varied
`total_migration`. **That is the wrong knob.** §7.4.2 measured that what absorbs `fst_loss` is the
**kernel decay `m`** — unique R² **0.429**, against **0.089** for `total_migration`. So the axis
that destroyed batch 1's identifiability is the one not tested here. The gradients above are only
5.7° apart (`det(J) = +0.098`, non-zero but ill-conditioned), which on this axis is not enough to
claim the confound is broken.

> **NEXT MEASUREMENT — DONE 2026-09-09, see §7.9.6. F_st tripled; F_c did not move.** Hold POPMULT and `total_migration` fixed and move
> the kernel `scale` from 5e-5 (1/m = 20 km, ~18 effective source demes) to ~1e-3 (1/m = 1 km,
> ~1.4 effective source demes, §7.4.2's table). F_st should move a lot. **The question is whether
> F_c moves at all** — it should not, because where an immigrant comes from is second-order for
> drift within a deme over 8–16 generations, whereas it is first-order for between-deme
> differentiation. If F_c is flat while F_st swings, F_c pins N and F_st is free to pin the
> kernel, and that *is* the break. Reuse `iso_nm_point`-style setup: regenerate
> `migration_rates.csv` at the new scale, one SLiM run, recapitate, then `temporal_fc.py`.

**One more caveat, found while re-running: the pedestal decomposition is NOT reliable at n=7.**
The `F_c − pedestal` column is trustworthy in the `n_big` arm (drift ratio 2.49 against a deme-size
ratio of 2.50) but not in `n_real`, where it collapses to 0.0078 / 0.0006 — far below `t/(2N)` and
not in the right ratio. Cause: the pooled-MAF ≥ 0.05 filter is an *ascertainment* on the sample,
and at 14 haplotypes the observable frequencies are multiples of 1/14, so conditioning on the very
quantity that forms F_c's denominator biases it. **Consequence: never read `ne_hat` at real n.**
It does not affect the statistic's use in ABC — the identical filter and identical n_i apply to
both sides, so `|F_c,sim − F_c,obs|` stays valid — but it does mean the *raw* F_c is the statistic
and the decomposition is a diagnostic only.

**Caveats, stated plainly.** One seed per point; no §7.3-style noise floor (the sd quoted is
subsample-draw spread on a FIXED tree, so it omits SLiM, recapitation and the mutation overlay).
F_st here is Hudson over the ~15 matched-field demes at the subsample size, not the production
whole-deme value over all demes, so its LEVEL is not comparable to §6.7's — only the elasticities
are, and those are ratios of like against like. `total_migration` moved only 2.5×, against a prior
spanning 300×.

**Which deme stands for a field** is not a free choice here — see §7.9.2. The probe defaults to
`--deme-choice nearest` (the cluster nearest the field's real coordinates, identical at both
timepoints) precisely because the production mapping puts the same field in different demes in
different years.

#### 7.9.5 The 99-deme cap is lifted [FIXED 2026-09-09 — full account in `OLD_LOGS.md` §C]

`numClusters ∈ {1,2,3}` (×33 = 99 demes) always sat exactly one under an msprime hard limit:
`pyslim.recapitate` splits **every** SLiM subpopulation from the ancestral population in a
**single** `population_split` event, and msprime caps that event at 100 derived populations.
`Python_Code/recapitate_util.py` merges in stages of ≤99 into throwaway intermediates and is
**exact, not approximate** (the intermediates live ~1e-13 generations at size 1.0). It **delegates
to `pyslim.recapitate` unchanged at ≤99 populations**, so it is identity at the current prior and
cannot perturb any existing result. Verified at 200 demes: 127 s, every tree `num_roots == 1`.
A grep for `pyslim.recapitate` should only hit that file.

A **second, unrelated wall** lives in SLiM and is easy to confuse with it: deme size is
`Average Count × POPMULT / numSubpops`, `asInteger` truncates, and SLiM refuses an empty
subpopulation. The ABC never meets it; the interactive `python Main.py` route (POPMULT=500) does.
`Main.main` guards it before SLiM runs, naming the minimum POPMULT for the layout — SLiM's own
message names a subpop id and says nothing about clusters.

> **[OPEN] The prior has NOT been widened, and lifting it is now a modelling decision.** The 43
> unique sequenced sites are **0.39 km apart at closest, 1.28 km median**, which is the scale CPB
> dispersal actually operates at — while 33 demes have a minimum centroid spacing of **6.05 km**,
> so the short-range part of the kernel, the only part with real biological content, cannot exist
> in the model at all. **But total simulated N is ≈3.33 × POPMULT no matter how many demes there
> are**, so more demes means smaller ones: holding today's deme size (505) at 200 demes needs
> POPMULT ≈ 30,000, just past the 25,000 ceiling. `numClusters` and POPMULT are entangled by
> construction, which is why §7.4.3 measured `numClusters` swinging median `pop` 1.7×. A defensible
> pairing is **numClusters=200 fixed with POPMULT ≥ ~10,000**, reported as a sensitivity axis
> against the current 33 — not a wider prior over both.

#### 7.9.6 The kernel test — F_st swings 3×, F_c does not move. THE CONFOUND BREAKS [MEASURED 2026-09-09]

§7.9.4 tested the wrong knob. This is the right one.

**Why this axis and not `total_migration`.** §7.4.2 measured what actually absorbs `fst_loss` in a
real batch: the **kernel decay `m`** at unique R² **0.429**, against **0.089** for
`total_migration`. Batch 1 failed because `m` free over four orders of magnitude soaked up F_st and
left `pop` at the prior. So the question is not "is F_c a different functional of (N, m)" in
general — it is specifically **does F_c survive the kernel swinging.**

**The swing, quantified before running** (from `cluster_distances.csv`, 33 demes):

| kernel `scale` | 1/m | effective source demes | weight on nearest | median source distance |
|---|---|---|---|---|
| 5e-5 *(baseline, the trees already on disk)* | 20 km | **16.7** | 0.139 | 28.2 km |
| 1e-3 *(this run)* | 1 km | **1.4** | 0.854 | 12.1 km |

A 12× change in effective source count, spanning most of the informative band §7.4.2 identifies
(≈[6e-6, 6e-4]) and running just past its top. If F_st does not move on this axis, nothing moves it.

**Built at the OLD constants on purpose** — μ 4.646e-7, r 2.75e-6, `ancestral_Ne` 6700 — because
the baseline trees were, and mixing in the Q = 100 rescaling (§7.9.3) would confound the kernel
effect with the scale change. POPMULT and `total_migration` (0.05) held fixed. Seed 1, 5 subsample
draws. `data/migration_rates.csv` backed up and verified restored.

| arm | POPMULT | F_st (5e-5 → 1e-3) | **ΔF_st** | F_c(t=16) (5e-5 → 1e-3) | **ΔF_c** | ΔF_c / subsample sd |
|---|---|---|---|---|---|---|
| n_big | 2000 | 0.04659 → 0.13105 | **+181%** | 0.03848 → 0.04105 | **+6.7%** | **1.5** |
| n_big | 5000 | 0.01952 → 0.05611 | **+187%** | 0.02742 → 0.02780 | **+1.4%** | **0.8** |
| n_real | 2000 | 0.05019 → 0.13568 | **+170%** | 0.16197 → 0.16649 | **+2.8%** | **1.1** |
| n_real | 5000 | 0.01947 → 0.05850 | **+201%** | 0.15480 → 0.15390 | **−0.6%** | **0.2** |

**F_st roughly triples in ALL FOUR cells (+170% to +201%). F_c moves −0.6% to +6.7%, every one of
them within 1.5 subsample standard deviations of its own noise — i.e. statistically
indistinguishable from not moving at all.**

**The elasticities, measured at BOTH POPMULTs, and this is the number to quote:**

```
                       d ln / d ln N          d ln / d ln kernel
                    P=2000    P=5000        P=2000    P=5000
F_st                -0.949    -0.926        +0.345    +0.352
F_c (drift term)    -0.996    -1.084        +0.043    +0.017
```

Both elasticities are **stable across a 2.5x change in POPMULT**, which is what licenses reading
them as properties of the statistics rather than of one point. `det(J) = +0.303` at POPMULT=2000
against **+0.098** on the `total_migration` axis (§7.9.4); the two gradient directions are
**17.5°** apart rather than 5.7°.

**What that means in the units that matter.** F_st's implied N scales as `kernel^0.364`; F_c's as
`kernel^0.043`. Across the kernel prior's informative band (~100×):

| statistic | swing in the N it implies, across the whole kernel range |
|---|---|
| F_st | **5.3–5.8×** — this IS the batch-1 failure |
| F_c | **1.07–1.22×** |

(the two figures are POPMULT=2000 and 5000; F_st's exponent on the kernel is 0.364/0.380, F_c's
0.043/0.016.)

**F_c reads N almost independently of the dispersal assumption.** That is exactly the property
§7.4.3 said was missing, and it is the first statistic measured here that has it. Paired with
F_st — which is now free to do what it is good at, pinning the kernel — this is a genuine second
equation rather than a second vote.

**Mechanism, and it is the reason to believe the result rather than just the numbers.** F_c
measures drift *inside* a deme over 8–16 generations. Where an immigrant came from is second-order
for that; it is first-order for between-deme differentiation. The kernel reshapes *which* demes
exchange migrants without changing how fast a deme's own allele frequencies wander, so F_st feels
it and F_c does not. Predicted before the run, and measured.

**Caveats, stated plainly.**
- **The 2×2 is complete** (added later on 2026-09-09 with `--simplify-first`, §7.9.7). The
  re-runs reproduced the already-measured cells **exactly** — n_real at POPMULT=5000 returned
  0.15390 and F_st 0.05850, digit for digit — which is what licenses mixing them.
- **One seed per point, and no §7.3-style noise floor.** The sd quoted is subsample-draw spread on
  a FIXED tree; it omits SLiM, recapitation and the mutation overlay. A real floor for `fc_loss`
  is a prerequisite before fitting, exactly as §7.5.5 was for LD.
- **Two points on the kernel axis**, so the elasticity is a two-point slope, not a curve. §7.4.2's
  table says the kernel stops mattering below ~6e-6 and above ~6e-4, so the response is certainly
  not log-linear across the whole prior.
- **The pedestal decomposition remains unusable at n=7** (§7.9.4's caveat): the n_real drift term
  goes *negative* at POPMULT=5000. Use raw F_c as the statistic; the decomposition is a diagnostic
  for the n_big arm only.
- F_st here is Hudson over the ~15 matched-field demes at the subsample size, not the production
  whole-deme value, so its LEVEL is not comparable to §6.7's. Only the elasticities are, and those
  are ratios of like against like.

> **What this does NOT settle.** F_c is not yet a fitted statistic. Before it can be:
> **(1)** the empirical side does not exist — it needs a script on the Beagle machine computing
> per-field per-year allele frequencies and Waples F_c over the 19 matched pairs on common
> variants, with the *identical* MAF filter (the filter is an ascertainment, §7.9.4);
> **(2)** a replicate noise floor (§7.3's test, applied to `fc_loss`);
> **(3)** weights from a pilot batch — never guessed (§7.4.1).
> And §7.9.1's ceiling is untouched by any of this: F_c escapes the budget because it measures a
> *change*, but it is still reading the same 324-generation window.

#### 7.9.7 Probe memory: chunked genotype extraction [ADDED 2026-09-09]

`temporal_fc.py --max-bytes N` caps the genotype matrix at N int8 entries and chunks the
extraction. At POPMULT=5000 with `--big-n 50` the matrix is ~243k sites × 3.3k nodes = **800 MB**,
which was enough to get the probe **killed by the OS** with ~2.3 GB free on the 15.2 GB box.

Chunking costs **one extra `genotype_matrix` traversal per chunk** — extraction tracks sites and
trees, not sample count (§7.5.2) — so it is a real time-for-memory trade, not free. One chunk
reproduces the unchunked path exactly.

**Verified bit-exact**, not assumed: on a toy msprime tree at three chunk sizes (≈11, 28 and 82
chunks) the maximum absolute difference against the unchunked path is **0.0**, including with a
deliberately duplicated group to exercise the dedup path.

**`--simplify-first` is the better lever, and it is the one to use.** Chunking trades time for
memory; simplifying to the rep's own nodes before extracting reduces BOTH. Measured on the
POPMULT=5000 kernel tree: **70,078 samples / 5.55M edges / 242,685 sites → 3,300 samples / 1.02M
edges / 57,816 sites.** Wall time for a full point went 64 s + 93 s → **28 s + 60 s**.

It is exact for this statistic: allele frequencies among the retained samples are untouched by
`simplify`, and the sites `filter_sites=True` drops are monomorphic in the retained set, which the
MAF filter removes anyway. **Confirmed empirically** — re-running the POPMULT=2000 kernel point
with `--simplify-first` reproduced F_c 0.16649 / 0.04105 and F_st 0.13568 / 0.13105 **exactly**.

> The tree itself is NOT the memory problem — measured at **0.32 GB** loaded. Both OOM kills came
> from `genotype_matrix` on top of a box with ~2.5 GB free (browser + Slack, no single hog), and
> the harness kills background tasks on SYSTEM pressure, so a modest job can be collateral damage.
> Check free memory before launching a long point.

#### 7.9.8 F_c validation — noise floor measured, empirical side written and round-trip tested [2026-09-09]

Two of the three things standing between §7.9.6 and fitting F_c are now done. The third is blocked
on a run that has to happen on the Beagle machine.

**A. The replicate noise floor. §7.3's test, applied to F_c.**

4 seeds at POPMULT=2000 and 2 at POPMULT=5000, at the **new Q=100 constants**, each re-rolling all
three dice production rolls — SLiM's forward mating (`-s`), the recapitation genealogy and the
mutation overlay (both unseeded). Every replicate recorded `fc_spec = ee863fff3bbf`.

| arm | POPMULT | n | mean F_c(t=16) | **floor sd** | CV% |
|---|---|---|---|---|---|
| n_real | 2000 | 4 | 0.16549 | **0.00124** | 0.75 |
| n_real | 5000 | 2 | 0.15430 | 0.00200 | 1.30 |
| n_big | 2000 | 4 | 0.03928 | **0.00100** | 2.53 |
| n_big | 5000 | 2 | 0.02742 | 0.00016 | 0.57 |

Against the POPMULT signal:

| arm | signal, POPMULT 2000→5000 | floor / signal | separation | floor / full-prior range |
|---|---|---|---|---|
| n_real | 0.01119 | **11.1%** | 9.0 sd | 11.9% |
| n_big | 0.01187 | **8.4%** | 11.9 sd | 5.6% |

**Verdict: usable, and it sits between the good statistics and the marginal one.** §7.3's bands
were `fst_loss` 2% (fine), `ld_loss` 2.8% (fine), `pi_loss` 30% (marginal). F_c at 8–12% is
noisier than F_st in relative terms but nowhere near π, and POPMULT 2000 and 5000 are separated by
**9 replicate standard deviations** in the arm that is actually measurable.

**B. The number that matters is not the floor — it is SELECTIVITY.**

A floor only says whether a statistic can resolve anything. The batch-1 failure was not a floor
problem: `fst_loss`'s floor was excellent (2%) and the pass still learned nothing about N, because
F_st moves *more* with the nuisance than with the parameter. Measured here on the same replicates,
at POPMULT=2000, each statistic in units of **its own** floor:

| change | F_c | F_st |
|---|---|---|
| POPMULT 2000 → 5000 (2.5×) | **11.1 sd** | 13.2 sd |
| kernel 5e-5 → 1e-3 (12× in source demes) | **2.6 sd** | **41.0 sd** |
| **selectivity = POPMULT response / kernel response** | **4.28** | **0.32** |

**F_st responds to the dispersal kernel three times more strongly than to population size. F_c
responds to population size four times more strongly than to the kernel. That 13× difference in
selectivity is §7.4.2's failure and its remedy, in one table.**

Note **selectivity is floor-independent** — the standard deviations cancel, so it is the ratio of
raw responses — which makes it robust to the fact that the F_st computed here is a noisier
estimator than production's (see caveats).

**C. The empirical side exists, and building it found a real bug.**

`ToUseOnBeagles/CalcTemporalFc.py` computes F_c over the 19 matched field pairs from the Beagle
VCFs; the spec lives in `Python_Code/fc_common.py`, imported by both sides and **copied** next to
the Beagle scripts exactly as `ld_common.py` is, with a `spec_hash()` both sides print.

`diagnostics/test_fc_roundtrip.py` simulates a tree, writes it as a VCF with popfiles and
specifier matrices laid out as the Beagle machine has them, runs the empirical script's real
`main()`, and compares against F_c computed straight from the genotype matrix. **The first run
disagreed by up to 2.6e-3** — too small to notice in a results table, far too large to be
rounding.

**Cause: multi-allelic sites, handled differently on each side.** The simulated side folded every
derived allele together as `G > 0`; the empirical side passed `max_allele=1` to scikit-allel,
which does not *exclude* a third allele but silently shrinks that site's denominator.
**7.1% of sites in a real simulated tree carry more than two alleles** (recurrent mutation), so
this was not a corner case. Both sides now drop non-biallelic sites **per pair** (invariant 7),
through `fc_common` — per pair for the same reason the LD code judges it per deme (§7.5.2).
Agreement is now **3.8e-11**, which is float64 summation order.

> **The effect on the simulated numbers alone is negligible — measured, not assumed.** Re-running
> a known cell under the fix moved F_c 0.04105 → 0.04107, **0.05%**. So every number in §7.9.4 and
> §7.9.6 stands as written. What the bug would have corrupted is the *comparison between sides*,
> which is the only thing the fit actually uses.

**D. The spec hash caught a second problem within the hour of being written.**

The first floor run was started before `fc_common` existed, and each replicate spawns a fresh
process — so replicates were picking up whatever `temporal_fc.py` was on disk at that moment.
Three completed replicates recorded either no spec or the pre-biallelic `bf930becbbb4` against a
current `ee863fff3bbf`: **two different statistics pooled into one floor.** All six stale records
were discarded and the run restarted; the driver now pins and prints the spec at the start.

**Generalisable lesson, and it is not the one `ld_common` already taught.** A spec hash is usually
framed as protecting against the *other machine* drifting. Here it caught **this** machine drifting
mid-run, because a long sweep re-reads its own code between replicates. **Any diagnostic that
spawns per-replicate subprocesses must record the spec in every record, and the summariser must
refuse to pool records that disagree** — `fc_floor_summary` asserts exactly that.

**Caveats, stated plainly.**
- **n=4 at one POPMULT and n=2 at the other.** The floor is not known to be constant across the
  prior, and the n_big arm's sd at POPMULT=5000 (0.00016) is 6× smaller than at 2000 on two
  replicates, which is too few to believe.
- **This is a floor on the STATISTIC, not on `fc_loss`.** The loss is `|F_c,sim − F_c,obs|` with a
  constant observed value, so the two floors coincide away from zero — but `fc_loss` does not
  exist yet and has not been measured.
- **The F_st column above is NOT production `fst_loss`.** It is Hudson over the ~15 matched-field
  demes at the subsample size, computed on the same sample sets so the comparison is like for
  like. Its floor (CV 9.6%) is much worse than production's 1.52% (§7.3) precisely because it uses
  far fewer demes and smaller samples. Selectivity is unaffected, being floor-free.
- One kernel point and one POPMULT step feed the selectivity numbers; both are two-point slopes.

**E. What remains, and why it cannot be done here.**

**The pilot batch is blocked on the empirical target file.** `fc_loss` is deliberately NOT wired
into `calculate_losses` yet: §10.2's rule is that the only check catching that class of bug is a
live trial asserting no empty cells, and without `averaged_temporalFc.csv` a live trial cannot
exercise `fc_loss`. Wiring it blind is the exact mistake §10.2 documents twice.

**The unblocking sequence, on the Beagle machine:**
1. ~~Copy `Python_Code/fc_common.py` next to `ToUseOnBeagles/`.~~ **DONE.**
2. ~~`python CalcTemporalFc.py`.~~ **STARTED 2026-09-09 evening — RUNNING, expect hours.**
   (The LD equivalent took 7.7 h over 17 chromosomes; F_c reads the same VCFs but does far less
   per site — no pair enumeration — so it should be well under that. Not yet timed.)
3. ~~Confirm the spec hash.~~ **CONFIRMED 2026-09-09: it printed `ee863fff3bbf`, matching the
   simulated side exactly.** The two sides are computing the same statistic. This was the gate
   that mattered, and it passed on the first try because the spec was shared rather than
   reimplemented.
4. ~~Copy the target and keep the per-chromosome files.~~ **DONE — present by 2026-09-11:**
   `data/empiricalStats/averaged_temporalFc.csv`, and the 17 per-chromosome files in
   `data/Fc_per_chr/` — keep them, they are the only route to a between-chromosome jackknife or a
   re-pool without re-reading every VCF (§7.5). **Before using the target, read §7.2.2F: its five
   largest pairs among those with n ≥ 7 are exactly the five containing relatives or failed
   samples.** ~~That target is now STALE~~ **Re-run under `49afd0877028` on 2026-09-12
   (§7.2.2H).**

Then: one live POPMULT=500 trial asserting no empty cells, and a pilot batch for weights.
`fc_loss` itself is wired — see §7.9.9.

#### 7.9.9 `fc_loss` wired, the ten lists collapsed, and the live trial that found the real bug [2026-09-09]

**A. Wired, and OFF by default.** `AnalyzeTreeSeq.calculate_temporal_fc` computes the simulated
side (subsampled to each field's real n_i, nearest-centroid deme so a field does not move between
years, spec from `fc_common`); `ABCAnalysisNoRedis` gained `_read_temporal_fc`, `_fc_loss`, and
loud failures if either side's file is missing. Both modules read the same `COMPUTE_FC`
environment variable, default **0**. With it off the CSV layout is byte-identical to before, so
nothing existing is disturbed. Turn it on only once the empirical target exists (§7.9.8E).

**B. The ten hand-written lists are GONE.** §10.2 asked for exactly this "the next time the
statistic set changes", and adding `fc_loss` was that time. `PARAM_NAMES`, `LOSS_NAMES` and
`CSV_FIELDNAMES` are defined **once**; `_build_row`, `_format_losses` and `_copy_raw_features`
derive from them, replacing four hand-written blocks duplicated across the two runner functions.
`collect_batch` now imports the same constants, closing the fourth file in that bug class.

Two guards make the failure mode impossible rather than merely unlikely:
- `_build_row` **raises** on a missing loss. `csv.DictWriter` fills a missing key with an empty
  string, which is precisely how a batch completed looking healthy with a blank `ld_loss` column.
  A crash costs one trial; a blank column costs a batch.
- `calculate_losses` asserts its own return keys equal `LOSS_NAMES` before returning.

**C. The live trial (§10.2's rule) — PASSED, and it earned its keep.**

Two POPMULT=500 runs, `COMPUTE_FC=1`, the first used as a stand-in observed target held in scratch
and **never** written into `data/empiricalStats/` (a fabricated target sitting where the real one
belongs is how someone ends up fitting to noise). `data/` backed up and restored.

- **obs-vs-obs identity: all seven losses exactly 0.0**, `fc_loss` included.
- A real trial produced a complete 14-column row: **no missing keys, no NaN, no blank cells in the
  written CSV.**
- Both spec hashes printed (`ld_common 4d1d1d92b25b`, `fc_common ee863fff3bbf`).

> **The first attempt failed, and the failure was the good kind.** `Main.main` takes the ACTUAL
> deme count while `model()` multiplies the raw draw by 33; passing `num_clusters=1` tripped the
> §7.9.5 deme-size guard **before SLiM ran**, with a message naming the real problem ("20
> sequenced sites but only 1 clusters"). Without that guard it would have been a SLiM error naming
> a subpopulation id.

**D. THE FINDING: `fc_loss` must POOL BEFORE DIFFERENCING. [MEASURED]**

The live trial ran two simulations at **identical parameters**, which makes it an accidental
one-replicate noise measurement — and it exposed that the obvious construction is wrong.

| construction | value between two IDENTICAL-parameter runs |
|---|---|
| per-pair: mean over pairs of \|F_c,sim − F_c,obs\| | **0.0280** |
| per-gap: mean over gaps of \|pooled F_c,sim − pooled F_c,obs\| | **0.0118** |
| all 19 pooled into one number | 0.0005 |

**Per-pair is 55× the fully-pooled difference.** Individual field pairs at n=7 scatter enormously;
the POPMULT signal is a level shift common to all of them, which survives pooling while the
scatter cancels. Fitting per-pair differences would be fitting mostly sampling noise.

**This is §7.2.1's geometry again**, and that section *measured* the consequence rather than
arguing it: under L1, uncorrelated scatter scores **worse** than no scatter, so a per-element loss
over noise-dominated elements penalises the simulation for having the right amount of variance.
It is also what LD already does — pool the demes, then compare (§7.5d).

**But do NOT pool all 19 into one number either.** In the live trial the two gaps moved in
*opposite* directions (t16 rose, t8 fell), so a single pool cancels real signal along with the
noise. The gaps carry different amounts of drift (8 vs 16 generations) and must stay separate —
the same argument that normalises out per-year entry counts (§7). **`_fc_loss` pools within a gap,
differences once per gap, then averages the gaps.**

Unit-checked: identity → 0.0; a uniform +0.01 shift → 0.01; swapping which *pair* holds which
value inside a gap → 0.0 (pooling is meant to be blind to that); a +0.10 level shift in one gap
only → 0.05.

The per-pair values are still written to `data/Output_Data/temporal_fc.csv` and copied into the
raw-feature store, so a per-pair variant can be scored after the fact without re-simulating.

**E. Consequence for the noise floor — §7.9.8's number does NOT transfer.**

§7.9.8 measured the floor on the **pooled statistic** (0.00124 sd at POPMULT=2000). `fc_loss` is a
*difference* of two pooled statistics, so its floor is larger — and the only figure in hand is
**≈0.0118 at POPMULT=500**, from the two identical-parameter runs above. That is one replicate, at
a POPMULT **below the prior floor** where demes hold ~50 individuals and drift noise is ~4× what it
is at POPMULT=2000.

**So the `fc_loss` floor is NOT yet measured at a prior-relevant POPMULT, and 0.0118 should be
read as a caution rather than an estimate.** ~~Re-run §7.9.8's replicates through the real
`calculate_losses` once the empirical target exists.~~ **Done 2026-09-12 — §7.9.10.**

#### 7.9.10 F_c goes live: target re-run, live trial, `fc_loss` floor, batch 3 submitted [2026-09-12]

**A. Empirical target re-run under `49afd0877028`.** `CalcTemporalFc.py` was re-run on the Beagle
machine with the new `fc_common.py`, and `FC_EMPIRICAL_SPEC` was updated (commit `aca8845`).
**Exactly the four pairs holding an excluded sample changed**, and the other 15 rows are identical to
the `ee863fff3bbf` run, which is what a pure exclusion change must produce:

| pair | n before → after | excess over pedestal |
|---|---|---|
| `Alsum18-2019 × -2023` (t=8) | 7/7 → 7/6 | 0.0411 → 0.0324 |
| `Alsum59-2019 × -2023` (t=8) | 7/7 → 6/7 | 0.0336 → 0.0327 |
| `H15-2015 × -2023` (t=16) | 8/7 → 8/6 | 0.0314 → 0.0324 |
| `H41-2015 × -2023` (t=16) | 7/7 → 7/4 | **0.0885 → 0.0352** |

Pooled observed F_c: **0.18463 at t=8, 0.18825 at t=16.** (The raw levels rose slightly because
smaller n raises the pedestal. The excess fell.)

**B. Live trial (§10.2) — PASSED.** POPMULT=500, m=5e-5, total_migration=0.05, 33 demes, Q=100
constants, through the real `model()` → `calculate_losses()` → `_build_row` → DictWriter.
- Both spec hashes matched, and obs-vs-obs gave exactly 0.0 on all seven losses.
- The 14-column row had no empty or non-finite cells, and all 19 simulated pairs were finite.
- `model()` took 63 s.

**C. The `fc_loss` floor — 30% of signal, marginal. [MEASURED]** Harness
`diagnostics/fc_loss_floor.py`.
- **Method.** Every replicate is a real `ABC.model()` → `calculate_losses()` call, so all five dice
  are re-rolled: SLiM, recapitation, the mutation overlay, the per-field-year subsample and the
  miscall assignment. m=5e-5, total_migration=0.05, 33 demes, Q=100 constants, 5 replicates at each
  POPMULT.
- **Safety.** It backs up and restores the `data/` files the pipeline overwrites.
- **Records.** `out/fc_loss_floor.jsonl`, with per-replicate outputs in `out/fc_loss_floor/`.
  `--summarize` refuses to pool mismatched specs, constants or parameters.

| quantity | POPMULT 2000 | POPMULT 5000 | change | pooled sd | separation | floor/signal |
|---|---|---|---|---|---|---|
| **`fc_loss`** | 0.01132 ± 0.00171 | 0.00550 ± 0.00180 | −0.00581 | 0.00175 | **3.3 sd** | **30%** |
| pooled sim F_c, t=8 | 0.19015 ± 0.00146 | 0.18720 ± 0.00120 | −0.00295 | 0.00134 | 2.2 sd | 45% |
| pooled sim F_c, t=16 | 0.20536 ± 0.00323 | 0.19669 ± 0.00240 | −0.00868 | 0.00285 | 3.0 sd | 33% |
| `fst_loss` | 0.01973 ± 0.00034 | 0.00888 ± 0.00022 | −0.01086 | 0.00029 | 37.9 sd | 3% |
| `pi_loss` | 0.03360 ± 0.00915 | 0.02270 ± 0.00142 | −0.01091 | 0.00655 | 1.7 sd | 60% |

- **Worse than §7.9.8's 11%, for two reasons.**
  - **The signal is smaller.** The POPMULT signal averaged over gaps is ~0.0058, against 0.0112 on
    clean genotypes. That is consistent with §7.2.2H's prediction that the miscalls cost ~20%.
  - **The noise is larger.** 0.00175 against 0.00124. Production re-draws the subsample and the
    miscall rates in each trial, and H41-2023 is down to n=4.
- **Sim F_c is still ABOVE observed at POPMULT=5000 in both gaps**, so F_c alone points higher.
- **`fc_loss` cannot reach zero.** Simulated t16 − t8 is 0.0095 at 5000, against an observed 0.0036.
  A loss that averages |per-gap difference| therefore has a minimum set by that mismatch. This is
  §7.2.2F point 3 surfacing in the loss, and it is why `fc_loss`'s CV climbs to 33% at 5000.
- **The Hudson-scale `fst_loss` floor at Q=100 (deferred by the §7.3 banner) is 1.7–2.5% CV** —
  §7.3's 2% verdict stands after the §6.7 fix.
- **The floor is a VETO, not a weight (§7.4.1).** At 30% F_c stays admissible; its case is
  kernel-insensitivity (§7.9.6), not precision. What the floor does set is a lower limit on ε: below
  it, rejection keeps lucky draws. The model is misspecified (`fst_loss` minimum ~18× its floor, the
  gap mismatch above), and luck is more available where a statistic is noisier, so a too-tight ε
  would bias toward the noisy regions of the prior. That is a bias, not just extra width.
- **Caveats.** 5 replicates each, so each sd is ±~35%. One kernel point. Nothing near POPMULT 25000.
- **Cost on this box** (other load running):

  | POPMULT | min/rep | peak Python | peak SLiM |
  |---|---|---|---|
  | 2000 | 2.7–3.2 (one contended rep 14.6) | ≤0.9 GB | 0.35 GB |
  | 5000 | 12–17 | 1.5–3.7 GB | 0.86 GB |

  That is about half the §3.1 old-constant cost at 5000 (30 min, 7.6 GB).

**D. μ does NOT need recalibrating at Q=100. [MEASURED, same replicates]** Median simulated π over
median observed π, masked, averaged over years: **0.972 at POPMULT 2000, 0.993 at 5000.** This
closes TODO §7.5b's "recalibrate μ at the new `ancestral_Ne`" — the Q algebra (§7.9.3) put the
level within 1–3%.

**E. `COMPUTE_FC` defaults ON (commit `3818aa7`).** Both modules read the variable with default `"1"`.
It is ON by default rather than set by the CHTC wrapper, because an unset variable would produce a
batch that runs cleanly with no `fc_loss` column — the §10.2 failure class. Verified: unset gives 14
columns ending `fc_loss`, and `COMPUTE_FC=0` gives 13, with obs-vs-obs exactly 0 in both.
**Consequence:** `collect_batch.py`'s header match is exact, so **re-analysing batch 1 or batch 2
now needs `COMPUTE_FC=0`.**

**F. Batch 3 SUBMITTED 2026-09-12 — the F_c pilot.**
- **Design.** 500 jobs × 5 trials = 2,500, on `3818aa7`. The prior is **identical to batch 1's**:
  `m` lognorm(1.5, 1e-4), `total_migration` U(0.001, 0.301), `pop` U(2000, 25000), `numClusters`
  {1,2,3}. μ and `r` are fixed at the Q=100 values.
- **Why that shape.** It is a clean before/after, differing from batch 1 only in F_c and the
  constants. It asks one question, whether `pop` tightens once `fc_loss` enters `D` with the kernel
  free, and it supplies the weights.
- **Memory.** `request_memory` was left at 64 GB. Extrapolating the measured peaks to 25000 gives
  **22–44 GB**, depending on the scaling exponent (1.10 from §3.1, 1.28 from median replicates, 1.54
  from the worst ones). That is too uncertain to cut, and the HTCondor logs will give the real
  number.
- **Expect partial N resolution at best.** F_c breaks the KERNEL confound, not the
  `total_migration` one: its elasticity on `total_migration` is −0.24, larger than F_st's −0.13
  (§7.9.4). Across the 300× prior that alone swings F_c's implied N ~3.8×. Until `total_migration`
  has a defensible ceiling, **read any `pop` tightening against that.**
- **Landed 2026-09-15 — see §7.9.11.** Memory came in at ≤40.9 GB (§3.1).

#### 7.9.11 Batch 3 landed — F_c tightens `pop`, but pushes it toward the ceiling, because simulated drift runs high at both gaps [VERIFIED 2026-09-15]

**Nothing here installs weights or reports a `pop` posterior.** The analyses were one-off scripts,
**not in the repo**. Everything is reproducible from `out/batch3/abc_results.csv`,
`out/batch3_raw/detailed_sim_results_<job>/run<k>/temporal_fc.csv` and `out/fc_loss_floor/`.

**A. Health: clean.**
- **Delivery.** 500 jobs × 5 = **2,500 of 2,500 trials**; every job wrote exactly 5 rows and 14
  columns ending `fc_loss`. HTCondor logs (`out/batch3_log/`): all normal termination, no
  evictions or holds, stderr is git clone progress only. The parameters printed in every `.out`
  match the CSV in all 500 jobs.
- **Coverage.** `pop` is uniform across deciles; the top decile is slightly OVER-full (z = +1.39),
  so no high-POPMULT trials died.
- **Files.** Pooled with plain `python collect_batch.py` → `out/batch3/abc_results.csv`. Raw
  per-job files are in `out/batch3_raw/`. **The first return had all 500 `detailed_sim_results_*`
  folders empty** — the wrapper was transferring the wrong path. Sohan fixed it 2026-09-15; each
  trial now has its 16 files, and the per-trial F_c tables reproduce the CSV's `fc_loss` to 1e-16.
- **Memory and time** are now measured; see §3.1.

**B. Decomposition (`collect_batch.py` §7, μ fixed so no nuisance column).** Unique rank-space R²:

| loss | pop | total_migration | m | numClusters | demographic |
|---|---|---|---|---|---|
| `pi_loss` | 0.113 | 0.023 | 0.177 | 0.019 | 0.331 |
| `fst_loss` | 0.073 | 0.085 | 0.422 | 0.012 | 0.593 |
| **`fc_loss`** | **0.260** | 0.108 | **0.004** | **0.152** | 0.524 |

- **§7.9.6 holds in production:** `fc_loss` carries essentially nothing from the kernel `m`
  (0.004), while `fst_loss` is still mostly `m` (0.422).
- **A new confound: `numClusters` (0.152).** F_c reads deme size, and deme size is
  `Average Count × POPMULT / numSubpops`, so the deme count moves F_c exactly as population size
  does. Like `total_migration` (0.108), it trades against `pop`.
- **Suggested weights** (§7.4.1 rule): **pi 0.229 / fst 0.409 / fc 0.362. NOT installed.**
- **The veto flag is up.** Section 4's floor/σ is **1.23 for `fc_loss`** (π 0.48, F_st 0.075),
  i.e. its replicate floor (§7.9.10C) exceeds its spread across the batch. §7.4.1 says a fitted
  statistic whose floor/σ nears 1 should be dropped. The reason is visible in the gradient: the
  median `fc_loss` falls only 0.0034 → 0.0031 between POPMULT 15–20k and 20–25k, far below its
  ~0.0018 replicate sd. Above ~15,000 it cannot rank draws.

**C. The read-out.** Population-size IQR as a fraction of the prior's. The read-out code reproduces
batch 1's documented 0.97/1.02/1.01/0.82/0.82 exactly.

| accepted | batch 1 (π+F_st) | batch 3, π+F_st (no F_c) | batch 3, π+F_st+F_c |
|---|---|---|---|
| top 20% | 0.97 | 0.88 | 0.77 |
| top 10% | 1.02 | 0.77 | 0.73 |
| top 5% | 1.01 | 0.81 | 0.68 |
| top 2% | 0.82 | 0.86 | 0.67 |
| top 1% | 0.82 | 0.74 | **0.49** |

- **What went right.**
  - F_c is taking over from the kernel: log-`m` IQR ratio stays 0.50–0.74, where batch 1
    collapsed to 0.07.
  - The migration trade-off weakened: spearman(log `pop`, `total_migration`) in the top 20% is
    −0.18, against batch 1's −0.34.
  - F_st still contributes: accepted median `fst_loss` 0.0060 against 0.0080 for the whole
    batch, unlike the LD failure (§7.5).
- **What went wrong.**
  - **The accepted median climbs toward the ceiling:** 18,100 at top 20% → 19,600 at top 1%.
  - **The tight levels select luck.** Accepted `fc_loss` median is 0.0021 at top 2% and 0.0018 at
    top 1%, i.e. at the replicate floor — below it rejection picks lucky draws (§7.9.10C).

**D. Why it pushes up: simulated F_c is above observed at BOTH gaps, and gains drift from t=8 to
t=16 while the data gains none.**

Observed pooled F_c is 0.18463 at t=8 and 0.18825 at t=16. Across batch 3:
- simulated t=16 is above observed in **99.9%** of trials;
- simulated t=8 is above observed in **85.5%** (median +0.0014 even above POPMULT 20,000).

**Most of the raw t16 − t8 difference is PEDESTAL COMPOSITION, not drift.** The t=16 pool holds
more small samples (H41-2023 n=4, Arlington-2015 n=4, LaszSE 5/5), so its pooled pedestal is
**0.0051 higher** than t=8's. That is identical on both sides and cancels in the loss. Subtracting
the pooled pedestal leaves the drift term:

| | excess, t=8 | excess, t=16 | growth 8→16 |
|---|---|---|---|
| observed | 0.0283 | 0.0268 | **−0.0015** |
| batch 3, POPMULT > 20,000 (median) | 0.0297 | 0.0315 | +0.0018 |
| fixed parameters, POPMULT 5000 (5 reps) | | | **+0.0044 ± 0.0012** |
| fixed parameters, POPMULT 2000 (5 reps) | | | +0.0101 ± 0.0037 |

- **The growth mismatch is NOT noise.** At fixed parameters (`out/fc_loss_floor/`: m=5e-5,
  `total_migration` 0.05, 33 demes) all five POPMULT-5000 replicates show growth ≥ +0.0029, and
  observed is **4.8 sd** below their mean (3.2 sd at 2000). The sd across batch trials (0.0021) is
  inflated by parameter variation and must not be used as the noise on this quantity.
- **Nothing in the prior reaches it.** Growth shrinks with deme size (spearman −0.36) and with
  `total_migration` (−0.33), reaching only ~+0.0015 at the prior's favourable corner. Just 3–4%
  of trials are as low as observed. **CONFIRMED DIRECTLY 2026-09-16 by the corner replicates
  (§7.9.11H): +0.00111 ± 0.00033 sem at that corner, 7.9 sem above observed, 0 of 10 replicates
  as low as the data** — so this extrapolation held, slightly optimistically.
- **Do not reason "two 8-generation steps add to one 16-generation change."** Most t=8 pairs are
  2019→2023 fields absent from the t=16 set; only Arlington and GarrisonNE have all three
  comparisons. And at n=7 the excess over pedestal is mostly not drift (§7.9.4), on both sides:
  GarrisonNE's per-pair excess is 0.0292/0.0298/0.0297 simulated and 0.0289/0.0244/0.0241
  observed for 15→19 / 19→23 / 15→23. The mismatch lives in the pooled gaps, not in these
  triangles.
- **Reading [INFERRED, and a professor question]:** real fields accumulate less drift over 16
  generations than any persistent simulated deme in the prior. A field sampled at the same
  coordinates in 2015 and 2023 may not be one persistent population: potato fields rotate yearly
  and are recolonised from overwintering beetles nearby (§11 dispersal refs). That is a model
  structure question. **Do not widen `total_migration` to chase it** — that is choosing the prior
  to fit.

**E. The loss's form does not matter — tested, not assumed.** Batch 3 re-scored under four
versions, each with its own §7.4.1 weights. The fixed-parameter floor is pooled over the ten
`fc_loss_floor` replicates.

| `fc_loss` version | accepted `pop` median, 20% → 1% | `pop` IQR ratio, 20% → 1% | demographic R² | floor/signal |
|---|---|---|---|---|
| **current**: mean over gaps of \|pooled diff\| | 18,100 → 19,600 | 0.77 → 0.49 | 0.52 | 0.30 |
| level: \|mean over gaps of pooled diff\| | 17,900 → 20,600 | 0.76 → 0.38 | 0.48 | 0.30 |
| t=8 only | 16,500 → 20,200 | 0.85 → 0.60 | 0.23 | 0.45 |
| t=16 only | 17,800 → 19,000 | 0.77 → 0.83 | 0.53 | 0.33 |

- **The push is not caused by the growth mismatch.** Even t=8 alone, which never sees t=16, pulls
  `pop` to 16–20k, because simulated t=8 is itself above observed across the prior.
- **"Level" only looks tighter.** Its top 1% (0.38) comes from trading t=8 below observed
  (−0.0006) against t=16 above (+0.0035), with a loss median (0.0014) under its own floor. That is
  compromise plus luck.
- **Keep the current form.** No variant fixes anything.

**F. Unequal-sample-size pairs read observed > simulated — not a bug; tracks leftover relatedness
[pattern VERIFIED, mechanism INFERRED].**
- **Code checked.** All 2,500 trial tables match `averaged_temporalFc.csv` exactly on pair list,
  t, n_a, n_b and pedestal. All 32 field-years get a miscall-rate pool of exactly n with no
  excluded sample. `test_fc_roundtrip.py` re-run passes (4.5e-11).
- **Pattern.** Over the 1,063 trials above POPMULT 15,000, gap-centred residual (obs − sim, each
  gap's mean offset removed):
  - the **6 unequal-n pairs** (Alsum18, Alsum59, H15, H41, both Arlington-2015) average **+0.0032,
    all positive**;
  - the 13 equal-n pairs average −0.0015, 3 positive;
  - one-sided Mann–Whitney **p = 0.0016**.
- **The better predictor is relatedness among RETAINED samples.** spearman(residual, fields' mean
  corrected within-field kinship) = **+0.76** (p < 0.001), from `kinship_correct.fit_year`.
  Unequal-n fields average 0.0048 against 0.0018; unequal n is mostly a proxy, since those fields
  are the ones that had problem samples.
  - The largest residual is OkrayGrosheks40 (equal n, +0.0119, z 2.1), which keeps half-sib S306
    (kinship 0.175).
  - Alsum18 keeps a 0.133 pair.
  - **H15 (+0.0077, z 2.7) has no kinship signal and is unexplained.**
- **The keep-half-sibs rule STILL HOLDS where the fit lands [MEASURED 2026-09-15].** The worry was
  that §7.2.2H's 19.5%/8.0% (POPMULT 2000/5000) would collapse at the larger demes the fit
  prefers, leaving real half-sibs uncancelled.
  - **Method.** A count-only census (one-off script, not in the repo) used fc_miscall.py's
    definition: a random field sample at real post-exclusion n, nearest deme, shares a parent via
    SLiM's `pedigree_p1/p2`. It needs no recapitation. 64 draws × 32 field-years per tree.
  - **The trees.** The existing `out/ld_probe_p{2000,5000}.trees`, plus one fresh SLiM-only run at
    POPMULT 20,000: 33 demes, m=5e-5, `total_migration` 0.05, r=1.02e-6, seed 20260915. It ran
    from a scratch mirror of `data/`, so the repo was untouched; 101 s, ≤2.3 GB.
  - **Check.** At 16 draws the census reproduced fc_miscall.py (21.9%/7.8% vs 19.5%/8.0%).

  | POPMULT | field deme (median) | field samples holding a shared-parent pair | pairs sharing a parent | × deme size |
  |---|---|---|---|---|
  | 2000 | 265 | 22.6% | 1.371% | 3.63 |
  | 5000 | 662 | 9.5% | 0.548% | 3.63 |
  | 20000 | 2651 | **2.6%** | 0.144% | 3.81 |

  - **It is a clean 1/N law:** P(pair shares a parent) ≈ 3.7 / deme size at all three points.
    Full sibs stay essentially absent (1 in 37,568 pairs at 2000 and 5000, none at 20,000).
  - **Where the fit lands.** Batch 3's top 20% implies a field deme of median ~1,470 (IQR
    ~940–2,380, from 0.1325 · POPMULT / numClusters). The law then predicts **~4.4%** of field
    samples holding a pair, 2.8–7.0% over the IQR.
  - **Observed: 2 of 32 real field samples (6.3%)** hold a retained pair at second degree or closer
    — OkrayGrosheks40-2019 (S306×S307, 0.175) and Alsum18-2023 (S21×S22, 0.133). P(≥2 of 32) is
    0.41 at 4.4% and 0.20 at 2.6%, so the data are consistent with the simulation.
  - **The census is a lower bound.** It counts only shared parents, while KING's second degree
    also includes grandparent and avuncular pairs.
- **So relatedness does NOT explain D's level misfit.** The simulation already produces about as
  many half-sib samples as the data hold. The net F_c the kept pairs add beyond what the
  simulation matches is of order 0.0003–0.0006 per gap (rough arithmetic from §7.2.2H's shift per
  unit prevalence) — far below the t=16 offset (+0.0047). The per-pair +0.76 correlation is real,
  but it is about which pairs scatter high, not about the pooled level. **H15 is still
  unexplained.**
- **Arlington2015 (n=2) shares Arlington-2015's coordinates exactly.** It is masked on both sides,
  so it biases nothing, but two beetles of that field are unused — §4's open duplicate question.

**G. What this licenses.**
- **Do not report a `pop` posterior from batch 3**, and do not install the suggested weights as a
  result:
  - the veto flag is up (B);
  - the tight levels sit at the floor (C);
  - the push toward the ceiling comes from a misfit no prior value closes (D, E).
- **Do not raise `POPMULT_MAX` because F_c points up.**
- **What batch 3 DOES establish:**
  - F_c removes the kernel confound in production (B, C);
  - the measured CHTC cost (§3.1);
  - a specific, quantified misfit — too much drift, and drift that grows where the data's does
    not — which is a model-structure finding to take to the professor, not an ABC-tuning problem.

**H. The corner replicates LANDED 2026-09-16 — growth falls 4× at the prior's edge and stops at
+0.0011, still positive in every replicate. The misfit is MOSTLY STRUCTURAL [VERIFIED]**

The question was whether the 8→16 drift growth (D) can reach zero anywhere in the prior. It cannot.

**Health and provenance.** 10 of 10 jobs returned, one row each, 14 columns ending `fc_loss`, no
empty cells. Every row carries the corner spec exactly — m 5e-5, `total_migration` 0.30 (the prior
ceiling), POPMULT 25,000 (the ceiling), numClusters 1 = 33 demes (the largest demes), μ 5.8e-7,
r 1.02e-6 — from `Python_Code/corner_replicates.csv` (`5c8baa4`), run through the same
`model()` → `calculate_losses()` path as a batch trial.

> **Two departures from the plan above, neither affecting the result.** The files came back to
> **`out/batch4_raw/`**, not `out/corner_raw/`, and are named `abc_results_{0..9}.csv` with
> `detailed_sim_results_{0..9}/run1/`. **No HTCondor logs were returned**, so the memory peak was
> NOT read. §3.1's ceiling row stands on batch 3's logs (median 32.2 GB, max 40.9 GB for jobs
> reaching 23–25k); this run would have been a free independent check of it at numClusters=1 and
> `total_migration` 0.30, and that check was lost, not failed.
> The analysis was a one-off script, **not in the repo**; everything below is reproducible from
> `out/batch4_raw/` and `data/empiricalStats/averaged_temporalFc.csv`.

**Cross-checks passed before reading anything.** All 19 pairs present and finite in every
replicate; pair list, `t`, `n_a`, `n_b` and pedestals identical to the empirical target; and
recomputing `fc_loss` from each per-pair table reproduces the CSV column to **1e-9**.

**THE RESULT.** growth = (pooled sim F_c t16 − pooled pedestal t16) − (same at t8):

| | excess t=8 | excess t=16 | **growth** |
|---|---|---|---|
| **observed** | 0.02832 | 0.02681 | **−0.00151** |
| **corner, n=10** | 0.02940 ± 0.00063 | 0.03052 ± 0.00081 | **+0.00111 ± 0.00105** |
| fixed params, POPMULT 5000 / `tm` 0.05 (§7.9.11D) | | | +0.0044 ± 0.0012 |

(± is the replicate sd. sem = 0.00033, so the mean is pinned to ±0.0004 exactly as planned;
95% CI **[+0.00036, +0.00186]**; per-replicate range −0.00029 to +0.00255.)

**Verdict, against the rule written down before the run: it clears the "above ~+0.001 → structural"
threshold, but only just, so the honest reading is MOSTLY structural with a real prior-edge
component.**

- **The prior edge does a lot.** Growth falls **4×**, from +0.0044 to +0.00111 — that is **56% of
  the way** from the reference to observed. The parameters are not irrelevant to this misfit.
- **And then it stops.** The mean sits **3.4 sem above zero** and **7.9 sem above observed**, and
  **not one of the ten replicates is as low as the data.** There is no value in the prior that
  reproduces the data's *flat* drift, which was the actual question.
- **Batch 3's extrapolation held, slightly optimistically:** §7.9.11D projected ~+0.0015 at the
  favourable corner from the batch's own gradients; measured +0.0011.

**The residual is concentrated in the LONG gap.** Simulated t=8 very nearly closes at the corner
(**+0.00108** above observed, against batch 3's +0.0014 median above POPMULT 20,000), while t=16
stays **+0.00371**. So both the level misfit and the growth misfit live in the 16-generation
comparison, and the t=8 offset §7.9.11E used to argue "the push is not caused by the growth
mismatch" is itself a prior-edge effect that this corner mostly removes.

> **DO NOT read this corner as a candidate parameter point. F_st rejects it flatly.** Simulated
> mean off-diagonal F_st is **0.000196 / 0.000193 / 0.000199** by year against observed
> **0.00915 / 0.00320 / 0.00701** — about **35× too low**, i.e. effectively panmictic
> (deme N ≈ 2.5k at `total_migration` 0.30 gives 4Nm ≈ 3000). That is why `fst_loss` is
> **0.00643 in all ten replicates with a spread of 1e-5**: it is not a fit, it is the
> no-structure floor, ≈ mean |F_st_obs|. A constant `fst_loss` across replicates is the
> signature of that floor, not of a quiet estimator. **Even had growth reached zero here, the
> corner could not have been the answer** — it only says where the F_c misfit goes when the prior
> is pushed as hard as it can be pushed.

**The other losses, for the record.** `fc_loss` averages **0.00240 ± 0.00050**, above its 0.0018
replicate floor (§7.9.10C) — so the corner is *not* a place `fc_loss` bottoms out either.
`pi_loss` averages **0.0247 ± 0.0050**, with the expected one-sided excursions above §7.2.1's
0.02094 flat floor (minimum 0.0209, maximum 0.0356).

**Caveats.** Ten replicates at ONE point, so this is a statement about the prior's most favourable
corner, not a surface. Growth was measured on pooled gaps, which is the quantity `fc_loss` uses;
the per-pair scatter is far larger and is not what is being tested (§7.9.9D). And the corner is
favourable *for F_c* by §7.9.11D's gradients — it is the worst corner for F_st, which is the point
of the banner above.

**WHAT THIS LICENSES.**
- **The drift finding is now a result, not a suspicion.** In 19 resampled fields, allele
  frequencies drift no further over 16 generations than over 8 (observed growth −0.0015), while a
  persistent simulated deme gives +0.0044 at POPMULT 5000 and **never reaches the data anywhere in
  the prior, including at its most favourable edge.** Take it to the professor as the top item.
- **The mechanism question is unchanged and is now the whole question:** is a field sampled at the
  same coordinates in 2015 and 2023 one persistent population, given yearly potato rotation and
  recolonisation from overwintering beetles nearby (§11 dispersal refs)? If not, the model needs
  field demes that are **re-founded from their neighbours** — a structure change, not a prior.
- **Still do not widen `total_migration` or raise `POPMULT_MAX`.** This run is the reason: the
  ceiling was *already* tested and it does not close the gap, so widening would buy a worse F_st
  fit for a misfit that stays.
- **Batch 3's read-out is unaffected.** No weights are installed and no `pop` posterior is
  reported (§7.9.11G stands).

#### 7.9.12 The re-founding model — the rotation premise CHECKED, and what building it would take [SCOPING 2026-09-17]

§7.9.11H leaves one question: are these fields persistent populations at all? This section checks
the premise against the records in hand and scopes the model change.

> **STATUS 2026-09-22: the toggle is in production and its floor is measured — H. Batch 5 (a
> 5×5 CHTC smoke test of the new prior) is running.**
>
> **STATUS 2026-09-20: the probe HAS run and it worked — jump to F, then G.** A–E below are the
> scoping written before it, kept because the predictions they made are what the measurement
> tested: B predicted that a bottleneck without outside founders would produce nothing (it made
> the misfit 2.5× worse), E predicted the kernel confound might come back (it did, 5.2×), and C's
> "fiddly two-step dance across generations" turned out to be two ordinary `early()` blocks.
> **Where A–E say the probe has not been run, read F.**

**A. The rotation record exists, and it is in a file already in the repo [VERIFIED]**

`data/final_data_for_modeling.csv` is not only the KMeans coordinate source. It carries `grower`,
`farm`, `field`, `field_fvid`, `croptype` and `year` over **2014–2024** — 3,266 field-years,
1,830 unique fields. `croptype` is a **USDA CDL code**: **43 = potato**, and it is 92% of rows
(2,994/3,266); the rest are corn (1), dry bean (42), soy (5), alfalfa (36), sweet corn (12) and a
long tail. So it is a scouting record *of potato fields*, which is what makes its gaps meaningful
and also what limits them (see the caveat).

Matching the sequenced sites to it on coordinates — they are **not** the same points, median
nearest-field distance **0.39 km**, 35 of 61 site-years within 0.5 km:

- **The sequenced sites are COMMERCIAL GROWER FIELDS, not research plots.** The H-series are
  **Heartland Farms** field codes H-10, H-15, H-41, H-53, H-67 — *not* the Hancock research
  station, which is the natural guess from the names and is wrong. Others map to Alsum Farms
  (their numbered fields 018/025/059/140 match to within 0.16 km), Mortenson Brothers, Okray
  Family Farm, Coloma Farms, Plover River Farms, Signature Farms and KA Farms.
- **`Arlington` is the one exception: 32.8 km from any scouted field** — the UW Arlington
  Agricultural Research Station, outside the commercial dataset entirely.
- **Over 11 years, each F_c-pair field matched within 0.5 km appears in potato in only 1–3 of 11
  years, and NO field ever appears in two consecutive years.** H-10: 2015 and 2020. H-41: 2015,
  2020, 2023. C-05: 2023 only. That is the signature of a 3–4 year rotation, which is standard
  Wisconsin practice.

> **The record is DEMONSTRABLY INCOMPLETE, so absence is NOT proof of non-potato.** Of the 35
> sequenced site-years whose field matched within 0.5 km, only **16 appear in the scouting record
> for that same year** — yet all 35 must have been in potato, because beetles were collected off
> them. The table therefore misses at least ~54% of true potato field-years, and "1–3 of 11" is a
> **lower bound on cropping frequency, not a measurement of it.** Read this as *consistent with*
> heavy rotation, not as proof of it. Confirming it properly is a question for the data provider
> (TODO §8.4), not another analysis.

> **AND THE OBVIOUS CONTROL CANNOT BE RUN.** The clean test would be to compare continuously
> cropped sites against rotated ones: if both show flat F_c, re-founding is not the explanation.
> There is exactly **one** candidate continuously cropped site — `Arlington`, the research station
> — and it carries **n=4** in 2015, the smallest retained sample in the set (its typo duplicate
> `Arlington2015`, n=2, is already masked by §7.0). One site at the smallest n cannot carry that
> contrast. **The within-data control for this hypothesis does not exist**, which means the
> hypothesis has to be judged on mechanism and on what the simulation does under it.

**B. Is re-founding a reasonable model? Yes, and it predicts the SPECIFIC signature [INFERRED]**

It is not an invented patch. Extinction–recolonization metapopulations are standard: **Slatkin
(1977)** for the propagule-vs-migrant-pool distinction, **Wade & McCauley (1988)**, **Whitlock &
McCauley (1990)** — and **Whitlock & McCauley (1999)** is already in §11 for a different reason.

The mechanism maps onto §7.9.11D exactly. If a field is re-founded each year by a small number of
colonists drawn from a large regional pool, its allele frequency is a fresh noisy draw each year,
nearly uncorrelated with the year before. The variance of the between-year difference is then
about twice the founder variance plus whatever the regional pool itself drifted — and the regional
pool is large, so it drifts slowly. **That sum is nearly flat in t**, which is what the data show
(excess 0.0283 at t=8 against 0.0268 at t=16) and what a persistent deme cannot produce at any
size, as §7.9.11H measured. So this is a mechanism whose qualitative prediction is already
observed, not extra freedom bolted on until the fit improves.

> **The bottleneck IS the model. Do not implement re-founding without it.** Replacing a deme from
> a migrant pool with no reduction in founder number just fills it with a blend of its neighbours:
> that *reduces* differentiation and creates no year-to-year independent noise. Small `k` is what
> produces the flat signature. An implementation without a founder bottleneck would produce
> nothing and would look like the idea failing.

**C. What building it would take [SCOPING, not measured]**

`SLiM_Code/CPBSampleSim{Win,Linux}.slim` is **57 lines with no scheduling logic at all** — constant
subpop sizes, one fixed migration matrix, and `treeSeqRememberIndividuals` at gens 308/316/324.
This would be the first real dynamics the forward model has ever had.

- **SLiM (~40–60 lines, BOTH files — §9: fix both or neither).** An event firing every **2
  generations** (2 generations/year, §11) that, for the demes being re-founded, sets that deme's
  migration rates to its kernel row normalised to 1 and shrinks it to `k`, plus a companion event
  restoring size the following generation. **In SLiM's WF model a size change takes effect for the
  next generation's offspring, so the bottleneck is a two-step dance across generations** — that is
  the one fiddly part and it should be checked against the SLiM manual, not taken from this note.
  A nonWF rewrite would give more control and is not worth it.
- **Python (cheap, and this is the §7.9.9B refactor paying off).** New `-d` constants passed from
  `Main.main` with **no defaults** (§10.1); new entries in `prior_distributions` and `PARAM_NAMES`.
  Before the ten-list collapse this was the §10.2 bug class; now it is one edit plus the prior.
  `Main.main`'s deme-size guard needs extending so `k ≥ 1` — a bottleneck walks straight into
  §7.9.5's empty-subpopulation wall.
- **Untouched:** recapitation, `AnalyzeTreeSeq`, `fc_common`, the empirical side, every statistic.

Roughly a day of code. **The validation is the real cost** — a new `fc_loss` floor under the new
model, a check that F_st stays fittable, then a pilot batch.

**D. Sequencing — the cheap probe FIRST. [DONE 2026-09-19 — the answer is in F]**

Hardcode one re-founding rule at one plausible `k`, run **SLiM only** at POPMULT=2000, recapitate,
and point `diagnostics/temporal_fc.py` at it for t=8 and t=16. **Does growth go flat or negative?**
If it does not, none of the parameterization work in C is justified, and the cost was one SLiM run.

> **It went flat. F has the numbers.** The design as executed added two things this plan did not
> have, and both earned their place: a **negative control** (the same bottleneck with ordinary
> migration — it made the misfit worse, which is what shows the immigration is the mechanism),
> and a **paired baseline** at the same POPMULT and constants, which reproduced §7.9.11D's
> independently measured +0.0101 and so validated the harness before anything was read from it.

> **Do NOT try to settle this analytically first.** §6.8 records two confident derivations from ρ
> that measurement overturned, and §7.9.4 shows the n=7 pedestal decomposition is itself unreliable,
> so the arithmetic here would be exactly the kind that has misled this project before. **The
> simulation is the only valid map.**

**E. Risks, stated before anyone starts**

- **It adds 2–3 parameters** (founder number, re-founding rate or schedule, propagule vs migrant
  pool) to a model that already cannot separate N from `m`. §7.4.3 binds: more nuisance dimensions
  are more places for the signal to go.
- **It moves F_st too, and the direction is not obvious.** Recolonization raises differentiation
  under a propagule pool and can lower it under a migrant pool (Wade & McCauley's condition is
  roughly whether `k` is below `2Nm + 1/2`). `fst_loss` already sits ~18× its floor at best
  (§7.4.2), so this could help or hurt and must be watched, not assumed.
- **IT MAY HAND BACK THE KERNEL CONFOUND — and that is the risk that matters most, because
  kernel-insensitivity is F_c's whole reason for existing [RAISED 2026-09-19, Sohan].** §7.9.6
  measured F_c moving <1.5 sd for a 12× swing in effective source demes, and §7.9.11B confirmed
  it in production (unique R² **0.004** on `m`). **That was measured under a PERSISTENT-deme
  model**, where an immigrant's origin is second-order for drift inside a deme. Under re-founding
  it is first-order by construction: the founders *are* the immigrants, so where the kernel says
  they come from sets the deme's new frequency directly. A propagule pool drawn from one near
  neighbour and a migrant pool drawn from the whole region give different founder variance.
  **So the probe (D) must measure BOTH things, not just the drift curve**: run it at two kernel
  values and check whether F_c's kernel elasticity is still ≈0.04. If re-founding fixes the drift
  growth but re-couples F_c to the kernel, the model has traded a temporal misfit for the batch-1
  identifiability failure, and that is not obviously a good trade.
- **Possible upside that cuts the other way:** concentrating drift into single generations should
  put **more** coalescence inside the 324-generation forward window, which would *loosen*
  §7.9.1's 0.5–7% information budget rather than only adding nuisance. That is a reason to measure
  the budget curve again under the new model, not a reason to expect it.
- **A 6–10 km deme is not a field, and re-founding one is not what rotation does.** Rotation moves
  the potato *within* a region; a cluster deme contains many fields. So at 33 demes the probe tests
  the **mechanism** (does a founder bottleneck flatten temporal F_c?) and not a faithful model of
  rotation. A faithful one needs field-resolution demes — which is §1.1's promoted spatial-
  resolution item, so the two are the same decision approached from different sides.

**F. THE PROBE RAN — growth flattens, F_st falls, and BOTH fitted statistics stop reading N [MEASURED 2026-09-19, replicated 2026-09-20]**

Step 1 of D, executed. Harness: `diagnostics/refound_probe.slim` (the production
`CPBSampleSimWin.slim` with one added block, file paths as `-d` constants so the repo's `data/`
is only ever read), scored through the real `diagnostics/temporal_fc.py`. Records in
`out/refound_probe.jsonl`. All arms POPMULT=2000, 33 demes, `total_migration` 0.05, Q=100
constants, SLiM seed 20260919, 5 subsample draws, `--simplify-first`.

**The re-founding rule.** Every **2 generations** (= annually, §11) each deme is cut to
`REFOUND_K` diploids and its migration row is renormalised so a fraction `REFOUND_M` of the
founding cohort is drawn from the *other* demes, kernel-weighted; the next generation restores the
deme's normal size and migration row. It fires on **odd** ticks so the three sampling ticks
(308/316/324, all even) are each one generation *after* a founding — a field colonised in spring
and sampled in summer. `REFOUND_M = 0.05` reproduces the production migration row exactly, i.e. a
**pure bottleneck**, and is the negative control; `REFOUND_M = 1.0` is full
extinction–recolonization from the migrant pool.

> **WF timing, verified here rather than taken from the manual** (§7.9.12C flagged this as the
> fiddly part). In SLiM 5.1 `early()` runs **before** offspring generation, so
> `setSubpopulationSize()` and `setMigrationRates()` called in `early()` of tick *g* govern the
> cohort produced in tick *g* itself. The "two-step dance across generations" the scoping note
> feared is just two consecutive `early()` blocks. Confirmed by printing deme sizes each tick.

**The decision number**, growth = (F_c t16 − pedestal t16) − (F_c t8 − pedestal t8):

| arm | growth, n_real | growth, n_big | F_st (n_big) |
|---|---|---|---|
| **baseline** — persistent deme | **+0.01038** | +0.01007 | 0.0208 |
| **bottleneck** — K=40, M=0.05 | **+0.02618** | +0.02751 | 0.0832 |
| **re-founding** — K=40, M=1.0 | **−0.00221** | +0.00344 | 0.0141 |
| re-founding — K=100, M=1.0 | −0.00060 | +0.00041 | 0.0052 |
| **observed** | **−0.00151** | −0.00151 | — |

- **The baseline validates the harness.** +0.01038 against §7.9.11D's independently measured
  **+0.0101 ± 0.0037** for a persistent deme at POPMULT 2000 — measured through the *production*
  `model()` path with miscalls, where this runs through `temporal_fc.py` without them. **That also
  proves growth is miscall-free**, which §7.9.11D assumed: the miscalls add **+0.0337 to the t=8
  excess and +0.0334 to t=16**, so they cancel in the difference. Only *growth* may be compared
  across the two paths; levels may not (§7.2.2H).
- **The negative control fired, and it is what makes the result mean something.** A bottleneck
  *alone* makes the misfit **2.5× worse**, exactly as §7.9.12B predicted — small `k` without
  outside founders is just faster drift inside a persistent deme, and drift still accumulates with
  t. **The flattening requires founders from elsewhere.** Anyone who implements the bottleneck and
  skips the migration renormalisation will conclude the idea failed.
- **Re-founding collapses growth by 3–5× and lands it at or just below zero.** n_real gives
  −0.0022 against an observed −0.0015; n_big gives +0.0034. **Do not read the sign too hard** —
  one tree per arm, and §7.9.4 warns the pedestal decomposition is unreliable at n=7, which is the
  n_real arm. What is solid is the size of the move: the baseline→re-founding swing is 0.0126,
  **3.4× the ±0.0037 replicate sd** §7.9.11D measured at this POPMULT.
- **F_st FALLS, which is the favourable direction — but this is NOT yet "F_st fits better", and
  the distinction matters [CORRECTED 2026-09-20].** Wade & McCauley's condition resolved the
  way E hoped (migrant pool, `k` below `2Nm + 1/2`): within these runs differentiation drops
  0.0208 → 0.0141 → 0.0052 (n_big; baseline → K=40 → K=100). **That across-arm comparison is
  sound — same estimator, same demes, same sample sizes.** What does NOT follow is any
  comparison to the observed 0.0065:
  - `temporal_fc.mean_hudson_fst` is **Hudson over the ~15 matched-field demes at the subsample
    size**; the observed 0.0065 is **pixy Weir–Cockerham over every deme at full n**. §7.9.6's
    own caveat says the LEVELS are not comparable, only elasticities. An earlier draft of this
    bullet compared them anyway; that was wrong and is retracted.
  - **And a discrepancy fell out of noticing it, which is [OPEN].** This baseline arm gives
    F_st **0.0208** at POPMULT 2000, where §7.9.4 records **0.0466** for a persistent deme at the
    same POPMULT, same `total_migration` 0.05 and same kernel, through the same function — a
    **2.2× gap**. The only known difference is the constants (§7.9.4 predates the Q=100 reset).
    But F_st should be nearly insensitive to `ancestral_Ne` (§6.2.1, <1% for a 3× change) and has
    no `r` dependence in expectation, and `ancestral_Ne` *fell* 6700 → 5259, which should raise
    F_st, not halve it. **So the 2.2× is unexplained and one of the two numbers is not what it
    is labelled.** Resolve it before quoting either in a write-up; a persistent-deme arm re-run
    at the old constants would settle it in ~2 min.
  So: re-founding moves F_st **down**, and the persistent model is the one known to overshoot
  (§7.4.2's best `fst_loss` is ~18× its floor). That is promising and it is not a measured fit.

**BUT IT DISSOLVES THE REASON F_c WAS VALUABLE. Both halves of §7.9.12E's warning landed.**

| | persistent deme | re-founding |
|---|---|---|
| `d ln(F_c drift) / d ln kernel` | **+0.050** | **+0.259** |
| `d ln(F_st) / d ln kernel` | +0.615 | +1.017 |
| F_st response ÷ F_c response | **12.4×** | **3.9×** |

(kernel 5e-5 → 1e-3, ×20; n_big. The baseline's +0.050 reproduces §7.9.6's **+0.043**, measured
at the old constants, which is the cross-check that licenses the comparison.)

- **F_c's kernel sensitivity rises 5.2×.** That is precisely the risk raised on 2026-09-19 and
  written into E before the run.
- **Worse, and this was NOT anticipated: F_c stops reading N and starts reading `K`.** At fixed
  POPMULT, raising `K` 40 → 100 (×2.5) drops the drift term ×2.40 (t=8) and ×2.62 (t=16) —
  **elasticity in `K` of −0.95 and −1.05**, which is exactly the −0.996 F_c had *in N* for a
  persistent deme (§7.9.4). The drift source has moved from deme size to founder number, because a
  deme re-founded from 100% immigrants is a fresh binomial draw of `2K` gametes and its subsequent
  size barely enters. **So §7.9.6's headline — "F_c reads N nearly independently of the dispersal
  assumption" — does not survive this model structure.**
- **CONFIRMED 2026-09-20 — F_c has stopped reading N, and so has F_st.** The direct test: hold
  `K = 40` and double POPMULT 2000 → 4000 (deme 202 → 404). If F_c still read deme size the drift
  term would halve. It does not move.

  | quantity, n_big arm | POPMULT 2000 | POPMULT 4000 | elasticity in N | persistent-deme reference |
  |---|---|---|---|---|
  | F_c drift term, t=8 | 0.02496 | 0.02514 | **+0.010** | −0.996 (§7.9.4) |
  | F_c drift term, t=16 | 0.02840 | 0.02596 | **−0.130** | −0.996 |
  | **F_st** | 0.01406 | 0.01350 | **−0.059** | −0.95 (§7.9.6) |

  - **This was the predicted half.** A deme re-founded from 100% immigrants is a fresh binomial
    draw of `2K` gametes; it expands to its full size in one generation and that size never gets
    the chance to matter. F_c reads `K`, with elasticity −1 (measured above), and N with
    elasticity ≈0.
  - **THIS WAS NOT PREDICTED, AND IT IS THE BIGGER FINDING: F_st stops reading N too.** Its
    elasticity falls from −0.95 to **−0.06**. Under re-founding, differentiation is set by founder
    sampling (≈ `1/2K`) rather than by drift inside a deme, so **BOTH fitted statistics now read
    `K` and the kernel, and NEITHER reads deme size.** Under this model structure `N` is not
    weakly identified — it is very nearly absent from the fitted statistics altogether. Whether
    that matters is §1.1's question, not an ABC-tuning one: `K` and the kernel *are* the
    dispersal-side parameters a gene-flow forecast needs. But it must be checked, not assumed,
    that the forecast really is insensitive to N — **do not carry over any N-identifiability
    claim from §7.9.4, §7.9.6, §7.9.8B or §7.9.11B into the re-founding model.**
  - **REPLICATED 2026-09-20 — the n_real/n_big disagreement was NOISE, and both arms now agree
    the statistics have stopped reading N.** 5 independent replicates per POPMULT (SLiM seeds
    20260921–25), each re-rolling SLiM, recapitation, the mutation overlay and the subsample, so
    these are a **real replicate floor**, not the fixed-tree subsample spread the single points
    carried. Records: `out/refound_reps.jsonl`; harness
    `diagnostics/refound_probe.slim` + `temporal_fc.py --simplify-first --max-bytes 100000000`.

    | arm | quantity | POPMULT 2000 | POPMULT 4000 | **elasticity in N [95% CI]** |
    |---|---|---|---|---|
    | **n_big** | F_c drift, t=8 | 0.02580 ± 0.00051 | 0.02438 ± 0.00033 | **−0.082 [−0.149, −0.015]** |
    | **n_big** | F_c drift, t=16 | 0.02537 ± 0.00023 | 0.02608 ± 0.00024 | **+0.040 [+0.004, +0.076]** |
    | **n_big** | F_st | 0.01342 ± 0.00023 | 0.01299 ± 0.00018 | **−0.047 [−0.111, +0.016]** |
    | n_real | F_c drift, t=8 | 0.01429 ± 0.00132 | 0.01174 ± 0.00152 | −0.283 [−0.732, +0.166] |
    | n_real | F_c drift, t=16 | 0.01435 ± 0.00148 | 0.01319 ± 0.00075 | −0.121 [−0.455, +0.213] |
    | n_real | F_st | 0.01164 ± 0.00135 | 0.01311 ± 0.00206 | +0.172 [−0.380, +0.723] |

    (± is the standard error of the mean over 5 replicates.)

    - **Every CI excludes −1 by a wide margin**, against the persistent deme's −0.996. The
      conclusion stands on both arms independently: **F_c and F_st have both stopped reading
      deme size.**
    - **The single point's −0.625 was noise.** It sits inside the replicated n_real CI
      [−0.732, +0.166]. n_real's standard errors are **4–10× n_big's**, which is §7.9.4's
      n=7 ascertainment warning measured rather than asserted — use n_big for any elasticity.
    - **n_big's two gaps have CIs that just exclude zero, in OPPOSITE directions** (−0.082 and
      +0.040). That is residual structure at the edge of this design's precision, not a signal;
      do not read either sign. What the design supports is "indistinguishable from 0, nowhere
      near −1".
  - **The growth number replicated too, and it is the headline (n_big, same replicates):**
    **−0.00043 ± 0.00044 at POPMULT 2000** against an observed **−0.00151** — 2.5 sem apart,
    i.e. the re-founding model reproduces the data's flat drift. The persistent deme at the same
    POPMULT gives **+0.0101 ± 0.0037** (§7.9.11D), which is ~28 sem above observed. At POPMULT
    4000 growth is +0.00170 ± 0.00053, so the match is not uniform across N and this is one
    (K, kernel) point, not a fit.
    > **The single-tree value was +0.00344, about 4 sd above the replicate mean.** One tree was
    > not enough for this quantity, which is exactly what the replicates were run to find out.
    > Quote the replicated figure, never the single point.
  - Cost: 10 replicates, ~2 min each, at 1.1 GB free. Both memory levers were needed (§7.9.7) —
    `--simplify-first` alone is what turned a 2 h thrash into 85 s.

**G. BUILD IT AS A TOGGLE, NOT A FORK — and the OFF state is VERIFIED IDENTICAL to production
[Sohan's design, 2026-09-20; VERIFIED same day]**

`refound_probe.slim` already carries the toggle: **`REFOUND_K < 0` disables re-founding entirely**
and the event block returns before touching anything. The question that decides whether a toggle
is safe is whether OFF is *the current model* or merely something close to it — because if it is
only close, every existing result silently changes meaning.

**It is exact.** Production `CPBSampleSimWin.slim` and `refound_probe.slim` with `REFOUND_K=-1`,
run at **the same SLiM seed** (424242, POPMULT 2000, r 1.02e-6, same input files):

| | production | toggle OFF |
|---|---|---|
| nodes / edges / individuals | 165,933 / 824,944 / 20,811 | **identical** |
| trees | 214,573 | **identical** |
| every node, edge and individual column | — | **identical** |
| individual metadata (SLiM pedigree + parent ids) | — | **identical** |

The **mating record itself matches**, not just the summary counts, so OFF is not an approximation
of the production model — it *is* it. (The files differ by ~5 KB of provenance and top-level
metadata, which record the script path and the `-d` constants and must differ by construction.
Never compare these trees byte-wise.) Harness: `diagnostics/toggle_check.py`.

**What the toggle buys, and why it is the right shape:**
- **No fork.** One model file, so §9's "fix both `.slim` files or neither" does not become "keep
  four files in sync".
- **Existing results stay comparable.** A batch with the toggle off is the same experiment as
  batch 3, which is what makes a before/after honest.
- **A and B can run in ONE batch**, so the comparison is free of the between-batch differences
  (constants, cluster matching) that already broke comparability twice.
- **It is a permanent regression test.** Re-running the check above after any edit to the forward
  model proves the edit did not disturb the persistent-deme path. Nothing else in this project
  tests that.

**Two rules for wiring it, both from §10.1/§10.2's history:**
- **`REFOUND_K` and `REFOUND_M` get NO defaults, in the `.slim` and in `Main.main`.** A missing
  `-d` must be a loud `undefined identifier`, exactly as `POPMULT` and `RECOMB` are. A toggle is
  precisely the kind of switch that could otherwise flip silently, and the output files would
  carry no record of which model produced them — §10.1's whole lesson.
- **Both must join `PARAM_NAMES` and the CSV**, so every row states which model it ran under.
  With the §7.9.9B refactor that is one edit, not ten.

**How to read this, and it depends entirely on §1.1.** Under the *old* framing this would be close
to fatal: the statistic built to identify N no longer identifies N. **Under the prediction framing
it may be the right answer rather than a loss.** Re-founding *is* gene flow — the founders are the
dispersers — so `K` and the kernel are the σ-side parameters the deliverable actually needs, and
F_c becomes a direct reader of the colonization process rather than an indirect reader of deme
size. That is closer to the product than N ever was. **But it must be stated as a change of target,
not smuggled in**: every claim in §7.9.6, §7.9.8B and §7.9.11B about F_c's selectivity *for N* is
conditional on the persistent-deme model and has to be re-derived, not inherited.

**Caveats, stated plainly.**
- **One tree per arm, one seed, 5 subsample draws.** The per-arm sds are subsample spread on a
  fixed tree and omit SLiM, recapitation and the mutation overlay — a **lower bound** on replicate
  noise (§7.9.6). No §7.3-style floor exists for any of this.
- **`K = 40` was not fitted, and its LEVEL agreement is partly circular.** It was chosen from the
  anchor excess ≈ `1/k` against the observed 0.028. The *flatness* is not what it was chosen for,
  and that is the result. On level, K=40 is in fact **~1.9× too high** once the +0.0335 miscall
  shift is added (0.053 against an observed 0.0283) and K=100 lands nearer (0.037) — **an
  indication that the level wants a larger `K`, not a fit.** Read it through §7.9.4's warning that
  the n=7 decomposition is unreliable.
- **Two kernel points**, so the elasticities are two-point slopes (§7.9.6's caveat, unchanged).
- **33 demes, so this is the mechanism, not rotation** — the bullet above.

**H. The toggle is IN PRODUCTION, and `fc_loss` is admissible on K [DONE 2026-09-22, commit `5b68c89`]**

**What shipped.** The re-founding block lives in both `CPBSampleSim{Win,Linux}.slim`, which still
differ only in path separators. Re-verified with `toggle_check.py` at seed 424242, POPMULT 2000:
the new script OFF is identical to the pre-toggle script, and ON (K=40, M=1.0) is identical to
`refound_probe.slim`, in every table column and the pedigree metadata both times.
- **One departure from the probe:** founders are capped at each deme's normal size (Sohan's
  call). A deme smaller than K is still recolonised wholly from immigrants, just not
  bottlenecked. It only bites when K exceeds a deme's size; `Main.main` prints how many demes it
  hits when not silent.
- **No defaults anywhere (§10.1):**
  - the `.slim` reads both constants even when OFF, so a missing `-d` is `undefined identifier`
    in either arm, and it `stop()`s on 0 ≤ K < 1 or M ∉ (0, 1];
  - `Main.main(refound_k=, refound_m=)` raises on None, a non-integer K, K = 0 or a bad M, all
    before any `data/` write;
  - `model()` reads `parameter["refound_k"]` with no `.get()` fallback.
- **The prior (Sohan, 2026-09-22):**
  - each trial is ON with p = 0.75 (`REFOUND_P_ON`);
  - K ~ log-uniform(10, 400), floored to an integer; OFF rows record **K = −1**;
  - `REFOUND_M = 1.0` is fixed and recorded per row;
  - **`pop` is still sampled in both arms**, so the A/B stays symmetric with batch 3, and
    because N-insensitivity was measured only over POPMULT 2000–4000 at one K.
- **CSV is 16 columns** (`refound_k`, `refound_m` after `recombination_rate`).
  `PRE_REFOUND_PARAM_NAMES` keeps the old list. **`collect_batch.py` needs `--pre-refound` to read
  batch 3 or earlier.** Verified: 2,500 rows, and §7.9.11B's decomposition reproduced exactly.
  With both arms in a batch it prints a warning that its sections pool them — per-arm analysis is
  not written yet.
- **Harness knock-ons:**
  - `fc_loss_floor.py` requires `--refound-k`, and reads pre-toggle records as K = −1;
  - `ld_probe`, `mu_calibrate`, `noise_floor` and `ridge_sweep` pass `REFOUND_K=-1`;
  - `corner_replicates.csv` carries the two columns, set to OFF;
  - `ABCAnalysis.py` (legacy) was not updated and now raises.
- **Live trial passed (§10.2):** POPMULT 500, one ON and one OFF trial, through
  `run_sims_from_csv` and the DictWriter — 2 rows × 16 columns, no empty or non-finite cells.

**The floor under re-founding.** `fc_loss_floor.py --refound-k {40,100}`, 5 reps each, POPMULT
2000, m 5e-5, `tm` 0.05, 33 demes. Records in `out/fc_loss_floor_refound_k{40,100}.jsonl`.
The K 40 → 100 step is 2.5×, the same size as the old POPMULT 2000 → 5000 test.

| loss | K=40 | K=100 | separation | floor/signal |
|---|---|---|---|---|
| `fc_loss` | 0.01867 ± 0.00344 | 0.00930 ± 0.00197 | 3.3 sd | 30% |
| `fst_loss` | 0.01184 ± 0.00024 | 0.00629 ± 0.00010 | 30 sd | 3.3% |
| `pi_loss` | 0.0911 ± 0.0131 | 0.0412 ± 0.0048 | 5.1 sd | 20% |

- **`fc_loss` reads K about as well as it read N under persistent demes** (§7.9.10C: also 3.3 sd,
  30%). Admissible, noisy.
- **`fst_loss` brackets the data and is NOT the §7.9.11H no-structure floor.** Production Hudson
  sim mean F_st is **0.0126 at K=40 and 0.0052 at K=100**; observed is 0.0092/0.0032/0.0070.
- **F_c and F_st may want different K.** Even at K=100, sim F_c is above observed at both gaps
  (t8 +0.0075, t16 +0.011), so F_c's level wants K > 100 while F_st's sits in 40–100. It is a
  K-version of the old N tension, at ONE (m, `tm`) point. Check it first in the pilot.
- **`pi_loss` under re-founding is almost pure LEVEL error — do not treat it as independent
  evidence on K.** Median sim/obs π is **0.912 at K=40 and 0.961 at K=100**, and |ln 0.912| =
  0.092 ≈ `pi_loss`. Sim π spread is near-flat (log-sd 0.001–0.003). The bottlenecks lower π,
  and μ was calibrated under persistent demes (0.97–0.99 there, §7.9.10D). So π reads K only
  through a calibration made for the other model.
- **ON trials are cheaper:** 1.1 min/rep at K=40 and 2.1 at K=100, against ~3 for persistent,
  at POPMULT 2000.

---


---

# §K — conventions history (old §10, §10.1, §10.2) [ARCHIVED 2026-09-27]

CLAUDE.md §10 keeps the rules.

## 10. Conventions

- Python; `numpy`, `scikit-allel`, `tskit`, `msprime`, `pyslim`.
- Scripts carry a `CONFIG` block at the top rather than argparse (`diagnostics/*` are the
  exception — they take flags, since they're meant to be swept).
- When implementing any new statistic, **verify it against tskit on a small msprime simulation
  first**, then apply it to real data. That pattern has caught real bugs three times.
- The pipeline **overwrites files in `data/` in place** (§2). Any experiment should run on a copy
  of the tree, not the repo, or redirect those paths.
- Prefer failing loudly over silent fallbacks for anything that changes units or scale — a default
  that reproduces a known-wrong scale is worse than a crash, because the output files carry no
  record of which scale they are on.

### 10.1 The scale-setting parameters have no defaults [VERIFIED 2026-08-23]

The rule above, applied to the constants that set the diversity and linkage scales. Three call
sites carried defaults reproducing known-wrong scales; all are now `None` + `raise`:

| site | old default | what it actually was |
|---|---|---|
| `Main.py::main` | `mutation_rate=5e-6` | the pre-2026-07-28 per-SNP calibration, **~10.8× too large** (§6.1.1) |
| `AnalyzeTreeSeq.py::analyze_tree_sequence` | `mutation_rate=1e-7` | a third scale, matching nothing |
| `AnalyzeTreeSeq.py::analyze_tree_sequence` | `recombination_rate=1e-8` | **the §6.3 bug value** — the rate that never reached SLiM |

The ABC path always passed all of these explicitly, so **no trial result was ever affected.** What
was exposed was the README's own documented interactive route (`python Main.py`), which silently
produced output at the 5e-6 scale — into files carrying no record of which scale they are on.

`Main.py`'s `__main__` block now prompts for μ, defaulting to `DEFAULT_MUTATION_RATE` imported from
`ABCAnalysisNoRedis`, so the constant has exactly one home. Both guards fire **before** any `data/`
write, so a wrong invocation cannot leave a half-overwritten tree behind. `ABCAnalysis.py` (the
legacy pyabc driver) was updated to pass both explicitly; it is otherwise stale and is a deletion
candidate.

`ancestral_Ne=6700` deliberately **keeps** its default — a fixed project constant with stated
provenance (§6.1), not a calibration that drifts.

### 10.2 Adding a statistic means editing TEN hand-written lists [VERIFIED the hard way, twice]

`ABCAnalysisNoRedis.py` enumerates the statistic set by hand in **eight** places — `fieldnames`,
the `row = {...}` dict, the `detailed_sim_results` copy loop, and the progress `print`, each
duplicated across the two near-identical runner functions (`run_sims_from_csv` and
`run_abc_simulation`). Nothing ties them together, and **`csv.DictWriter` fills a missing key with
an empty string rather than raising.**

Wiring `ld_loss` in hit exactly that. It was correctly in `calculate_losses` and in **both**
`fieldnames` lists — every file-level and import-level check passed — and the CSV still came out
with a **blank `ld_loss` column**, because the `row` dicts are built from explicit keys. A batch
would have completed looking perfectly healthy and been useless for its only purpose. The copy
loop had the same hole, silently dropping the LD curves from the raw-feature store.

**It recurred on 2026-09-08, in a different file.** `diagnostics/collect_batch.py` has the same
shape: `ld_loss` was in `EXPECTED_FIELDS` (so the pilot's files parsed and the inventory looked
perfect) but missing from `LOSSES` and `FITTED`, so **every analysis section silently omitted it**
— including the weights, which were the batch's only purpose. Two more lists, ten in total.

**The rule this buys: a live single-trial run is the ONLY check that catches this class of bug.**
Type checks, import checks, and reading the diff do not. Before any batch that changes the
statistic set, run one real trial at a small POPMULT and assert the output row has **no empty
cells** — that is a stronger check than inspecting any particular column.

```bash
# ~4.5 min at POPMULT=500, exercises Main -> analyze_tree_sequence -> model ->
# calculate_losses -> DictWriter, i.e. every link the batch will use
```

**~~Collapse these lists to one shared constant the next time the statistic set changes.~~ DONE 2026-09-09** (§7.9.9B), when `fc_loss` was added. `ABCAnalysisNoRedis` now defines `PARAM_NAMES` / `LOSS_NAMES` / `CSV_FIELDNAMES` **once**, and `_build_row`, `_format_losses` and `_copy_raw_features` derive from them; `collect_batch` imports the same constants. Two guards turn the failure mode from unlikely into impossible: `_build_row` **raises** on a missing loss instead of letting `DictWriter` write an empty string, and `calculate_losses` asserts its own return keys equal `LOSS_NAMES`.

The lesson that stands regardless: **a missing statistic never raises on its own — it just is not there.** It cost two separate silent failures before the lists were collapsed, and the second was caught only because the numbers it should have printed were conspicuously absent. **The live single-trial run is still the check** (§7.9.9C): it is what confirmed the collapse works, and it is what found the `fc_loss` pooling bug the same day.

---


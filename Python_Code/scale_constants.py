"""The three constants that set the model's SCALE. One home, imported everywhere.

WHY THIS FILE EXISTS. mu, r and ancestral_Ne were literals in seven places
(ABCAnalysisNoRedis, Main, AnalyzeTreeSeq, and five diagnostics), and CLAUDE.md 6.8.1 recorded
that a planned change to r "would be silently ignored" in most of them. That is 10.2's bug class
applied to constants instead of to statistic lists, and the fix is the same: collapse them.
IMPORT FROM HERE. Do not re-declare any of these.

Deliberately dependency-free (no numpy, no sim stack) so anything can import it.

-------------------------------------------------------------------------------------------
THE MODEL RUNS AT 1/Q REAL SCALE, AND Q IS NOW A CHOICE (Q = 100, set 2026-09-09)
-------------------------------------------------------------------------------------------
Standard coalescent rescaling simulates N/Q and compensates with mu -> mu*Q, r -> r*Q, m -> m*Q,
G -> G/Q. This model does that for the ANCESTRAL phase and not for time (CLAUDE.md 6.1.3): the
forward phase runs its true 324 generations, anchored to the real WI/NY invasion.

Until 2026-09-09, Q was DERIVED from Cohen et al.'s ancestral_Ne = 6700 and came out at ~78.5.
CLAUDE.md 6.1 establishes from four independent directions that 6700 fails its own paper's
internal check by ~52x, so a Q derived from it inherits that. **The arrow is now reversed.**
Q is a declared modelling decision (chosen round, for tractability) and ancestral_Ne is DERIVED:

    MU_TRUE  = 5.8e-9      measured. Zhang et al. 2026, GBE 18(2):evag027 -- 16 CPB
                           parent-offspring trios, ~32.8x coverage, 92 de novo mutations
                           (95% CI 4.7-7.2e-9). CLAUDE.md 6.1.2.
    R_TRUE   = 1.02e-8      measured, and Ne-FREE. Hawthorne 2001 (Genetics 158:695-700), 1,032 cM
                           over 18 linkage groups, divided by Yan et al. 2023's 1,008 Mb
                           assembly = 1.02 cM/Mb. A linkage map has no population-size term
                           anywhere in its derivation -- the defect in every other estimate this
                           project has used. CLAUDE.md 6.8.1.
    PI_OBS   = 0.0122      measured, after the 5.1 callable-site denominator fix.

    Ne_true  = PI_OBS / (4 * MU_TRUE)  = 5.259e5     <- what the ancestral population really was
    Q        = 100                                    <- CHOSEN
    ANCESTRAL_NE      = Ne_true / Q    = 5259        <- what we simulate
    MUTATION_RATE     = MU_TRUE * Q    = 5.8e-7
    RECOMBINATION_RATE= R_TRUE  * Q    = 1.02e-6

Self-consistent by construction: 4 * 5259 * 5.8e-7 = 0.0122 = PI_OBS.

WHY THIS IS BETTER THAN THE OLD CHAIN, and it is not just tidiness. Every constant now descends
from two MEASURED quantities and one DECLARED Q. Nothing descends from Cohen's 6700. The old
mu = 4.646e-7 was a number fitted to make pi come out right at an ancestral size nobody believes;
the new one is an external measurement times a stated scale factor, and it lands 2% away from the
fitted value (5.8e-7 vs 4.646e-7 * (6700/5259) = 5.92e-7), which is the consistency check.

WHAT THIS BREAKS, deliberately (Sohan's call 2026-09-09: re-running batches is cheap):
  * Batches 1 and 2, and every timing in CLAUDE.md 3.1 -- all ran at mu 4.646e-7, r 2.75e-6,
    ancestral_Ne 6700.
  * WEIGHTS in abc_standardize.py. 7.4.1's and 7.6.1's were both fitted with mu FREE and are
    invalid anyway (7.9.3); re-derive from batch 3.
  * The 6.1.1 saturation table. The branch_div ceiling is 2*(324 + 2*ANCESTRAL_NE), which moves
    27448 -> 21684, so those numbers must be re-measured before being quoted again.
  * LD_MIN_BIN, which is DERIVED from r as 1/(2*G*r). See ABCAnalysisNoRedis.

WHAT IT DOES NOT CHANGE:
  * The information budget (7.9.1). The fraction of pairs coalescing in the forward window is a
    property of the forward phase (N, m, G) alone; ancestral_Ne sets how DEEP the rest is, not
    how many coalesce early. Still 0.5-7% across the prior, still tracking 1/(1+4Nm).
  * F_st. It is (E[T_b]-E[T_w])/E[T_b], so ancestral_Ne cancels top and bottom -- which is why
    6.2.1 measured F_st moving <1% for a 3x change in it. Nm ~= 78-103 stays a REAL-WORLD number
    and must never have a factor of Q applied to it.

COST, expected rather than measured: recapitation scales as Ne^2.34 (6.2.1), so 6700 -> 5259 is
about 0.6x, and r dropping 2.7x means fewer edges in the forward phase (6.8.1 measured 3-4x
faster runs at the lower r). Both changes make trials CHEAPER.

**RECALIBRATE BEFORE TRUSTING pi's LEVEL.** MUTATION_RATE here is the ANALYTIC value. 6.1.1
measured that the analytic value lands ~2% low because forward coalescence pulls branch_div below
4*Ne_anc and multiple hits eat a further 1-2%. Re-run diagnostics/mu_calibrate.py at the new
ANCESTRAL_NE to get the fitted value; until then treat pi's level as good to a few percent.
"""

# ---- measured inputs (never edit without changing the citation) --------------------------
MU_TRUE = 5.8e-9            # Zhang et al. 2026, GBE 18(2):evag027 (95% CI 4.7-7.2e-9)
R_TRUE = 1.02e-8            # Hawthorne 2001 map / Yan et al. 2023 assembly = 1.02 cM/Mb
PI_OBS = 0.0122             # CLAUDE.md 5.1, genome-wide mean after the denominator fix

# ---- the declared scale ------------------------------------------------------------------
Q = 100                     # CHOSEN, not derived. See the header.

# ---- derived, and used everywhere --------------------------------------------------------
ANCESTRAL_NE = 5259         # = round(PI_OBS / (4 * MU_TRUE * Q))
MUTATION_RATE = MU_TRUE * Q         # 5.8e-7. NOT a mutation rate -- report theta = 4*N*mu.
RECOMBINATION_RATE = R_TRUE * Q     # 1.02e-6. NOT a recombination rate -- see the header.

# The forward run length, in generations. 325-generation WI/NY split (Cohen et al. 2022 Fig 3a)
# at 2 generations/year; SLiM Remembers at 308/316/324 for 2015/2019/2023. Hardcoded in
# CPBSampleSim*.slim -- this is a MIRROR for the arithmetic below, not the source of truth.
FORWARD_GENERATIONS = 324


def branch_div_ceiling():
    """Hard upper bound on branch-mode diversity, 2*(G + 2*Ne_anc) (CLAUDE.md 6.1.1).

    A pair that never coalesces forward waits an ancestral 2*Ne_anc. Used by the 7.9.1
    information-budget arithmetic, which reads the forward-coalescence fraction off the gap
    between this ceiling and the measured branch_div.
    """
    return 2 * (FORWARD_GENERATIONS + 2 * ANCESTRAL_NE)


def ld_time_depth_bp():
    """Shortest distance the forward window controls, 1/(2*G*r) (CLAUDE.md 7.5.1 pt 5).

    Below this, LD is set by the FIXED ancestral phase and carries no POPMULT signal. Derived
    from r, so it MUST move whenever RECOMBINATION_RATE does -- CLAUDE.md 6.8.1 lists forgetting
    this as one of the three places an r change is silently ignored.
    """
    return 1.0 / (2 * FORWARD_GENERATIONS * RECOMBINATION_RATE)

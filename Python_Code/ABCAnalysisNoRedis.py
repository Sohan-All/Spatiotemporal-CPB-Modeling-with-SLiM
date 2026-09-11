import csv
import math
import numpy as np
import sys
import shutil

from pathlib import Path
from scipy import stats

# Shared LD binning. Safe to import at module level: ld_common pulls in only numpy/csv/hashlib,
# NOT the simulation stack, so diagnostics keep importing this module without SLiM (CLAUDE.md 3).
import ld_common as ldc
import fc_common as fcc

# Main is imported lazily inside model() so this module stays importable without the simulation
# stack; diagnostics/*.py reuse the readers and get_keep_mask (CLAUDE.md 9).

# THE THREE SCALE-SETTING CONSTANTS LIVE IN scale_constants.py. Import, never re-declare --
# they were literals in seven files and CLAUDE.md 6.8.1 recorded a planned change being silently
# ignored in most of them. Re-exported here under their historical names so callers keep working.
#
# Q = 100 as of 2026-09-09 (CLAUDE.md 7.9.3): mu and r are external measurements times a DECLARED
# scale factor, and ancestral_Ne is derived from them -- nothing descends from Cohen's 6700 any
# more. Neither is a biological rate; only theta = 4*N*mu and rho = 4*N*r are meaningful.
#
# mutation_rate is FIXED, not inferred (it left prior_distributions on 2026-09-09). There is no
# unknown to draw: it is DEFINED by 4*Ne_anc*mu = pi_obs with Ne_anc chosen and pi_obs measured.
# Do NOT set it to MU_TRUE = 5.8e-9 -- that is the unrescaled rate, Q times too small for a model
# running at 1/Q scale. Do NOT widen a prior to the trio CI: that uncertainty belongs to
# Ne_true = pi_obs/(4*mu_true), which is REPORTED, not simulated.
import scale_constants as sc

DEFAULT_MUTATION_RATE = sc.MUTATION_RATE            # 5.8e-7
DEFAULT_RECOMBINATION_RATE = sc.RECOMBINATION_RATE  # 1.02e-6
DEFAULT_ANCESTRAL_NE = sc.ANCESTRAL_NE              # 5259

# Total N ~ 3.33*POPMULT. Raised 12000 -> 25000 on 2026-08-26 (CLAUDE.md 6.7).
# ~44 GB and ~1.9 h per trial at the ceiling -- size CHTC requests from CLAUDE.md 3.1.
POPMULT_MAX = 25000

# recombination_rate is intentionally absent -- fixed at DEFAULT_RECOMBINATION_RATE.
prior_distributions = {
    "m": stats.lognorm(s=1.5, scale=np.exp(np.log(0.0001))),   # kernel decay; unidentifiable (7.1)
    "total_migration": stats.uniform(loc=0.001, scale=0.3),     # U(0.001, 0.301)
    "pop": stats.uniform(loc=2000, scale=POPMULT_MAX - 2000),   # POPMULT ~ U(2000, 25000)
    "numClusters": stats.randint(1, 4),                         # 1, 2 or 3; scaled x33 in model()
    # mutation_rate is intentionally absent -- fixed at DEFAULT_MUTATION_RATE as of 2026-09-09.
    # It was lognorm(s=0.02, scale=DEFAULT_MUTATION_RATE), tightened 0.5 -> 0.05 -> 0.02 because a
    # wider draw let mu rather than POPMULT dominate pi_loss (CLAUDE.md 7.2.1). Even at s=0.02 it
    # still did: the batch-2 decomposition put pi_loss's unique R2 at 0.149 on mu against 0.080 on
    # pop (CLAUDE.md 7.6.1). Deleting the draw deletes that nuisance variance outright, which is
    # why pi_loss's weight can now rise on demographic signal rather than being held down to keep
    # the mu draw out of D. See DEFAULT_MUTATION_RATE above for why there is no unknown to draw.
}

# THE STATISTIC AND PARAMETER SETS, DEFINED ONCE. CLAUDE.md 10.2 records this bug class costing
# two separate silent failures: the set was enumerated by hand in TEN places (fieldnames, the row
# dict, the raw-feature copy loop and the progress print, each duplicated across the two runner
# functions, plus collect_batch's LOSSES/FITTED), nothing tied them together, and csv.DictWriter
# fills a missing key with an EMPTY STRING rather than raising. A batch completed looking healthy
# with a blank ld_loss column. **Derive everything from these lists; never retype a statistic
# name.** 10.2 asked for exactly this collapse "the next time the statistic set changes" -- adding
# fc_loss on 2026-09-09 is that time.
PARAM_NAMES = ["m", "total_migration", "pop", "numClusters", "mutation_rate", "recombination_rate"]
_BASE_LOSSES = ["pi_loss", "fst_loss", "ld_loss", "ibd_loss", "dxy_loss", "genrel_loss"]

# Temporal F_c is OFF until the empirical target exists (7.9.8E). AnalyzeTreeSeq reads the SAME
# variable and the two must agree -- with it on, both sides are required and a missing file raises.
COMPUTE_FC = bool(int(__import__("os").environ.get("COMPUTE_FC", "0")))
LOSS_NAMES = _BASE_LOSSES + (["fc_loss"] if COMPUTE_FC else [])
CSV_FIELDNAMES = ["iteration"] + PARAM_NAMES + LOSS_NAMES

# Raw per-year features copied into the detailed store so offline sigma has values, not just
# losses. Temporal F_c is NOT per-year (it spans years), so it is copied separately.
RAW_FEATURE_STATS = ["diversities", "divergences", "fst", "relatedness", "ld"]

# Empirical target for temporal F_c. Produced by ToUseOnBeagles/CalcTemporalFc.py, which must be
# RUN on the Beagle machine -- writing this path by hand would be fabricating data.
TEMPORAL_FC_OBS = Path("../data/empiricalStats/averaged_temporalFc.csv")
TEMPORAL_FC_SIM = Path("../data/Output_Data/temporal_fc.csv")


# Subpops too thinly sampled for a usable pairwise Fst are dropped from the FITTED statistics --
# at n<=3, 46.7% of 2015's pairs return a negative Fst. Only 2015 is affected (CLAUDE.md 7.0).
EXCLUDE_SMALL_SUBPOPS = True
MIN_SUBPOP_N = 4

# Lowest distance bin entering ld_loss, in bp. DERIVED FROM r, not a free choice: the
# 324-generation forward window only controls d >= 1/(2*G*r), and everything shorter is set by
# the FIXED ancestral phase and carries no POPMULT signal at any price (CLAUDE.md 7.5.1 pt 5).
# It therefore MOVES whenever RECOMBINATION_RATE does -- 6.8.1 lists forgetting that as one of
# three places an r change is silently ignored. At r = 1.02e-6 the cut is 1513 bp, so the first
# admissible ld_common bin edge is 1778 (was 562 at r = 2.75e-6).
# ld_loss is NOT fitted and will not be (7.8); this is kept correct so the diagnostic stays
# meaningful, not because anything downstream depends on it.
LD_MIN_BIN = 1778

# The ld_common spec the empirical targets were computed under (ldCalcOut.txt, 2026-09-07).
# Different spec on the two sides => the curves are binned differently and ld_loss is meaningless.
# Checked once in getObservedData(), where it is cheap and fires before any trial runs.
LD_EMPIRICAL_SPEC = "4d1d1d92b25b"

# The fc_common spec the empirical temporal-F_c target was computed under. The target now in
# data/empiricalStats/ predates fc_common.EXCLUDED_SAMPLES (CLAUDE.md 7.2.2F), so this stays at the
# OLD hash -- and fc_loss raises -- until ToUseOnBeagles/CalcTemporalFc.py is re-run under the
# current spec and this is updated to the hash that run prints.
FC_EMPIRICAL_SPEC = "ee863fff3bbf"


# ---------------------------------------------------------------------------
# Feature I/O + helpers
# ---------------------------------------------------------------------------

def _read_vector(path):
    '''Read a headerless single-column CSV into a 1-D float array. Reads EVERY row (fixes the
    csv.DictReader bug that silently dropped subpop 0 -- CLAUDE.md 5.7).'''
    vals = []
    with open(path, mode='r', newline='', encoding='utf-8') as f:
        for row in csv.reader(f):
            if not row:
                continue
            v = row[0].strip()
            if v != "":
                vals.append(float(v))
    return np.array(vals, dtype=float)


def _read_matrix(path):
    '''Read a headerless square CSV into a 2-D float array; blank cells -> NaN.'''
    matrix = []
    with open(path, mode='r', newline='', encoding='utf-8') as f:
        for row in csv.reader(f):
            if not row:
                continue
            matrix.append([np.nan if v.strip() == "" else float(v) for v in row])
    return np.array(matrix, dtype=float)


def _read_ld_decay(path):
    '''Read an LD-decay table written by ld_common.write_decay -> (sum_r2, cnt, labels).

    Returns SUMS, not means. Pooling demes is sum(numerators)/sum(denominators), never a mean of
    means (CLAUDE.md 5.2, invariant 4); write_decay stores the mean and the pair count, so the
    numerator is recovered as mean*count.

    Handles both sides: the empirical file's labels are site names, the simulated file's are deme
    indices "0".."K-1". Both are in specifier-matrix row order (7.5e), which is what makes the
    column-wise mask valid on both.
    '''
    rows = []
    with open(path, mode='r', newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            rows.append(row)
    if len(rows) != ldc.N_BINS:
        raise ValueError(f"{path}: {len(rows)} bins vs ld_common's {ldc.N_BINS} -- the two sides "
                         f"were computed under different specs, so ld_loss is meaningless")
    labels = [c[3:] for c in rows[0] if c.startswith("r2_")]
    lo = np.array([int(r["bin_lo"]) for r in rows], dtype=np.int64)
    if not np.array_equal(lo, ldc.BIN_EDGES[:-1]):
        raise ValueError(f"{path}: bin edges differ from ld_common's -- specs have drifted")
    cnt = np.array([[float(r[f"n_{p}"]) for p in labels] for r in rows], dtype=float)
    mean = np.array([[np.nan if r[f"r2_{p}"].strip() == "" else float(r[f"r2_{p}"])
                      for p in labels] for r in rows], dtype=float)
    return np.nan_to_num(mean, nan=0.0) * cnt, cnt, labels


def _ld_pooled_mean_abs_diff(sum_sim, cnt_sim, sum_obs, cnt_obs, keep):
    '''Mean |r2_sim - r2_obs| over the fitted distance bins, demes pooled by pair count.

    RAW levels, deliberately. Subtracting each curve's asymptote to fit "shape" instead was tried
    and REJECTED (CLAUDE.md 7.5.4): the sim-obs floor offset flips sign across the POPMULT range
    rather than acting as a fixed pedestal, and removing it flattens the minimum from 111% to 18%.
    '''
    fit = ldc.BIN_EDGES[:-1] >= LD_MIN_BIN

    def _pool(s, c):
        s = np.asarray(s)[:, keep].sum(axis=1)
        c = np.asarray(c)[:, keep].sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(c > 0, s / np.where(c > 0, c, 1), np.nan)

    m_sim, m_obs = _pool(sum_sim, cnt_sim), _pool(sum_obs, cnt_obs)
    ok = fit & np.isfinite(m_sim) & np.isfinite(m_obs)
    if not np.any(ok):
        return np.nan
    return float(np.mean(np.abs(m_sim[ok] - m_obs[ok])))


_GEO_DIST_CACHE = {}


def _haversine_m(lat1, lon1, lat2, lon2):
    '''Great-circle distance in metres (matches GenerateClusterData.distance formula).'''
    r = 6371000.0
    p = math.pi / 180.0
    a = (0.5 - math.cos((lat2 - lat1) * p) / 2
         + math.cos(lat1 * p) * math.cos(lat2 * p) * (1 - math.cos((lon2 - lon1) * p)) / 2)
    return 2 * r * math.asin(math.sqrt(a))


def get_site_geo_distances(year):
    '''Pairwise REAL-SITE geographic distances (metres) for a year's subpops, indexed by
    subpop = specifier-matrix row order -- the same ordering as the pi/Fst/relatedness matrices
    (see GenerateClusterData.assign_genomes_to_clusters_idv_year). Used for BOTH the observed and
    simulated IBD slopes. Real-site coords are cols 1 (lat), 2 (lon) of the specifier.
    Cached (the specifier files are fixed).'''
    if year in _GEO_DIST_CACHE:
        return _GEO_DIST_CACHE[year]
    coords = np.genfromtxt(Path(f"../data/Genetic_Data/specifier_matrix_{year}.csv"),
                           delimiter=",", usecols=(1, 2))
    lats, lons = coords[:, 0], coords[:, 1]
    n = len(lats)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                D[i, j] = _haversine_m(lats[i], lons[i], lats[j], lons[j])
    _GEO_DIST_CACHE[year] = D
    return D


def ibd_slope(fst_matrix, geo_dist):
    '''Rousset isolation-by-distance slope: OLS slope of Fst/(1-Fst) on ln(distance) over all
    off-diagonal pairs. Masks non-finite Fst, Fst>=1, and non-positive distances (e.g. the
    duplicate-site typo -- CLAUDE.md 4). Returns NaN if <2 usable pairs.'''
    fst = np.asarray(fst_matrix, dtype=float)
    n = fst.shape[0]
    xs, ys = [], []
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            d = geo_dist[i, j]
            f = fst[i, j]
            if not np.isfinite(f) or f >= 1.0 or not np.isfinite(d) or d <= 0:
                continue
            xs.append(math.log(d))
            ys.append(f / (1.0 - f))
    if len(xs) < 2:
        return np.nan
    slope, _ = np.polyfit(np.array(xs), np.array(ys), 1)
    return float(slope)


def _offdiag_mean_abs_diff(A, B):
    '''Mean absolute difference over off-diagonal entries (count-normalized within year);
    NaN entries skipped.'''
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    n = A.shape[0]
    mask = ~np.eye(n, dtype=bool)
    return float(np.nanmean(np.abs(A[mask] - B[mask])))


def _pi_log_mean_abs_diff(pi_sim, pi_obs):
    '''Mean absolute difference of log(pi), element-wise (count-normalized); non-finite/non-positive
    entries skipped. Log-space gives relative error and linearises the theta=4Nmu ridge.'''
    pi_sim = np.asarray(pi_sim, dtype=float)
    pi_obs = np.asarray(pi_obs, dtype=float)
    mask = np.isfinite(pi_sim) & np.isfinite(pi_obs) & (pi_sim > 0) & (pi_obs > 0)
    if not np.any(mask):
        return np.nan
    return float(np.nanmean(np.abs(np.log(pi_sim[mask]) - np.log(pi_obs[mask]))))


_KEEP_MASK_CACHE = {}


def _specifier_site_names(year):
    '''Site names in specifier-matrix row order (col 0) -- the canonical subpop ordering (4).'''
    names = []
    with open(Path(f"../data/Genetic_Data/specifier_matrix_{year}.csv"), encoding="utf-8") as f:
        for line in f:
            if line.strip():
                names.append(line.split(",")[0].strip())
    return names


def _build_row(iteration, parameters, losses):
    """One CSV row, built from PARAM_NAMES + LOSS_NAMES rather than by hand.

    RAISES on a missing loss. That is the whole point: csv.DictWriter fills a missing key with an
    EMPTY STRING, so the old hand-written dicts turned "I forgot to add the new statistic" into a
    batch that completed, looked healthy, and was useless for its only purpose -- twice
    (CLAUDE.md 10.2). A crash here costs one trial; a blank column costs a batch.
    """
    row = {"iteration": iteration}
    for k in PARAM_NAMES:
        if k == "total_migration":
            row[k] = parameters.get(k, 0.05)
        elif k == "recombination_rate":
            row[k] = parameters.get(k, DEFAULT_RECOMBINATION_RATE)
        elif k == "mutation_rate":
            row[k] = parameters.get(k, DEFAULT_MUTATION_RATE)
        else:
            row[k] = parameters[k]
    for k in LOSS_NAMES:
        if k not in losses:
            raise KeyError(f"calculate_losses did not return {k!r}; got {sorted(losses)}. "
                           f"LOSS_NAMES and calculate_losses must agree (CLAUDE.md 10.2).")
        row[k] = losses[k]
    return row


def _format_losses(losses):
    """Progress line, derived from LOSS_NAMES so a new statistic appears without an edit."""
    return " ".join(f"{k.replace('_loss', '')}={losses[k]:.4g}" for k in LOSS_NAMES)


def _copy_raw_features(iteration_dir):
    """Keep this trial's raw features so offline sigma has values, not just losses.

    Per-year matrices plus, when it is on, the temporal F_c table -- which is NOT per-year, since
    it spans them. Copying raw features is what makes a batch re-scorable without re-simulating
    (it is how a different LD_MIN_BIN could be tried after the fact, CLAUDE.md 7.8.2).
    """
    for year in ["2015", "2019", "2023"]:
        for stat in RAW_FEATURE_STATS:
            src = Path(f"../data/Output_Data/{stat}_{year}.csv")
            if src.exists():
                shutil.copy2(src, iteration_dir / f"{stat}_{year}.csv")
    if COMPUTE_FC and TEMPORAL_FC_SIM.exists():
        shutil.copy2(TEMPORAL_FC_SIM, iteration_dir / TEMPORAL_FC_SIM.name)


def _read_temporal_fc(path):
    """{(site_a, site_b): (t, fc)} from a temporal-F_c file. Same layout on both sides."""
    out = {}
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            out[(row["site_a"].strip(), row["site_b"].strip())] = (int(row["t"]), float(row["fc"]))
    if not out:
        raise ValueError(f"{path} holds no field pairs")
    return out


def _fc_loss(obs, sim):
    """POOL each generation gap first, THEN difference. Mean over gaps of |pooled difference|.

    THE ORDER MATTERS AND IT IS NOT THE OBVIOUS ONE. Differencing per field pair and then
    averaging -- the natural reading of "mean absolute difference", and what this function did
    first -- fits mostly sampling noise. Measured directly by the 10.2 live trial: two runs at
    IDENTICAL parameters gave pooled F_c of 0.19127 and 0.19076, a difference of 0.0005, while the
    per-pair loss between the same two runs read **0.0280 -- 55x larger**. Individual field pairs
    at n=7 scatter enormously; the POPMULT signal is a level shift common to all of them, and it
    survives pooling while the scatter cancels.

    This is 7.2.1's geometry again, and that section measured the consequence rather than arguing
    it: under L1, uncorrelated scatter scores WORSE than no scatter, so a per-element loss over
    noise-dominated elements penalises the simulation for having the right amount of variance.
    It is also what the LD statistic already does -- pool the demes, then compare (7.5d).

    Pooled by GAP rather than over all 19 pairs at once for the reason calculate_losses normalises
    out per-year entry counts (CLAUDE.md 7): t=8 has 10 field pairs and t=16 has 9, they carry
    different amounts of drift (8 vs 16 generations), and neither should outvote the other.
    Note the two gaps can move in OPPOSITE directions -- in the live trial t16 rose while t8 fell
    -- so pooling all 19 together would cancel real signal, not just noise.

    The sampling pedestal is NOT subtracted from either side -- it is identical by construction
    (the simulated demes are cut to the same n_i) and cancels in the difference. Subtracting it
    would be right for a point estimate of Ne and is wrong here (fc_common).

    The per-pair values are still written to data/Output_Data/temporal_fc.csv and copied into the
    raw-feature store, so a per-pair variant can be scored after the fact without re-simulating.
    """
    missing = set(obs) ^ set(sim)
    if missing:
        raise ValueError(
            f"temporal F_c field pairs differ between the two sides: {sorted(missing)[:4]}. "
            f"Both sides build the list from fc_common.matched_field_pairs, so this means they "
            f"are running different specs or different specifier matrices.")
    by_gap = {}
    for key, (t, fc_obs) in obs.items():
        t_sim, fc_sim = sim[key]
        if t_sim != t:
            raise ValueError(f"{key}: generation gap {t_sim} != {t}")
        if np.isfinite(fc_obs) and np.isfinite(fc_sim):
            by_gap.setdefault(t, []).append((fc_sim, fc_obs))
    if not by_gap:
        raise ValueError("temporal F_c: no usable field pairs on either side")
    # pool WITHIN the gap, then take one absolute difference per gap
    return float(np.mean([abs(np.mean([a for a, _ in v]) - np.mean([b for _, b in v]))
                          for v in by_gap.values()]))


def get_keep_mask(year):
    '''Boolean mask over specifier-matrix rows: True = subpop retained in the FITTED statistics.
    Drops subpops with fewer than MIN_SUBPOP_N diploid individuals when EXCLUDE_SMALL_SUBPOPS is
    set. Indexed by specifier row order, so the SAME mask is valid for the observed and the
    simulated matrices (CLAUDE.md 4) -- that is what keeps the comparison element-wise.
    Raises if a specifier site is missing from the popfile rather than silently keeping it.'''
    year = str(year)
    if year in _KEEP_MASK_CACHE:
        return _KEEP_MASK_CACHE[year]

    names = _specifier_site_names(year)

    if not EXCLUDE_SMALL_SUBPOPS:
        mask = np.ones(len(names), dtype=bool)
    else:
        counts = {}
        with open(Path(f"../data/Genetic_Data/popFile{year}"), encoding="utf-8") as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2:
                    counts[parts[1].strip()] = counts.get(parts[1].strip(), 0) + 1
        missing = [s for s in names if s not in counts]
        if missing:
            raise ValueError(f"{year}: specifier sites absent from popFile{year}: {missing}")
        mask = np.array([counts[s] >= MIN_SUBPOP_N for s in names], dtype=bool)

    _KEEP_MASK_CACHE[year] = mask
    return mask


def model(parameter):
    '''
    The model function that runs the SLiM simulation with the given parameters: 
    1. migration rate (m)
    2. population size (pop)
    3. number of clusters (numClusters)
    4. mutation rate (mutation_rate)
    5. recombination rate (recombination_rate)
    
    :param parameter: This is a dictionary containing the parameters for the simulation.
    '''
    
    #Get the parameters
    m = parameter.get("m", prior_distributions["m"].rvs())
    total_migration = parameter.get("total_migration", prior_distributions["total_migration"].rvs())
    pop = int(np.floor(parameter.get("pop", prior_distributions["pop"].rvs())))
    numClusters = parameter.get("numClusters", prior_distributions["numClusters"].rvs()) * 33  #scale to 33, 66, or 99
    mutation_rate = parameter.get("mutation_rate", DEFAULT_MUTATION_RATE)
    recombination_rate = parameter.get("recombination_rate", DEFAULT_RECOMBINATION_RATE)

    #Run the model - change silent to true for actual runs
    import Main   # lazy: see the note beside the imports at the top of this file
    Main.main(num_clusters=numClusters, migration_rates_modifier=m, population_modifier=pop,
              total_migration=total_migration, mutation_rate=mutation_rate, recombination_rate=recombination_rate, silent=True)
    
    
    #Read in the simulated output data (pi vector; dxy, Fst, relatedness matrices)
    outDict = {}
    for year in ["2015", "2019", "2023"]:
        outDict[f"{year}_diversity"] = _read_vector(Path(f"../data/Output_Data/diversities_{year}.csv"))
        outDict[f"{year}_divergence"] = _read_matrix(Path(f"../data/Output_Data/divergences_{year}.csv"))
        outDict[f"{year}_fst"] = _read_matrix(Path(f"../data/Output_Data/fst_{year}.csv"))
        outDict[f"{year}_relatedness"] = _read_matrix(Path(f"../data/Output_Data/relatedness_{year}.csv"))

        # LD decay (CLAUDE.md 7.5). Fails loudly rather than defaulting: an absent file means
        # AnalyzeTreeSeq.COMPUTE_LD is off, and a trial silently missing a FITTED statistic is
        # exactly the kind of scale/completeness error 10 says must crash instead (10.1).
        ld_path = Path(f"../data/Output_Data/ld_{year}.csv")
        if not ld_path.exists():
            raise FileNotFoundError(
                f"{ld_path} missing -- ld_loss is a fitted statistic but the simulated LD curve "
                f"was not written. Set AnalyzeTreeSeq.COMPUTE_LD (or COMPUTE_LD=1 in the "
                f"environment).")
        outDict[f"{year}_ld_sum"], outDict[f"{year}_ld_cnt"], _ = _read_ld_decay(ld_path)


    if COMPUTE_FC:
        # Written by AnalyzeTreeSeq.calculate_temporal_fc in the same run. If it is absent the
        # simulated side did not compute it, which means COMPUTE_FC disagrees between the two
        # modules -- fail loudly rather than emit a trial missing a fitted statistic.
        if not TEMPORAL_FC_SIM.exists():
            raise FileNotFoundError(
                f"COMPUTE_FC is on but {TEMPORAL_FC_SIM} was not written. AnalyzeTreeSeq reads the "
                f"same COMPUTE_FC environment variable -- set it for the whole process, not just "
                f"this module.")
        outDict["temporal_fc"] = _read_temporal_fc(TEMPORAL_FC_SIM)
    return outDict



def calculate_losses(x, x0):
    '''
    Per-statistic distances between observed (x) and simulated (x0) feature sets.

    Each statistic is averaged over its entries within a year (so the 24/17/20-subpop years
    contribute comparably), then across years. pi is compared in LOG space. There is deliberately
    NO total_loss -- abc_standardize.py builds the combined distance offline (CLAUDE.md 7).

    Returns (all un-standardized):
      - pi_loss     : FITTED  (log-space, element-wise)
      - fst_loss    : FITTED  (off-diagonal)
      - ld_loss     : FITTED  (binned LD decay, demes pooled, bins >= LD_MIN_BIN)
      - ibd_loss    : DIAGNOSTIC only (|IBD slope difference|)
      - dxy_loss    : DIAGNOSTIC only (off-diagonal)
      - genrel_loss : DIAGNOSTIC only (off-diagonal)

    IBD is diagnostic, not fitted: the observed slope is indistinguishable from zero in all three
    years (CLAUDE.md 7.1).
    '''
    pi_terms, fst_terms, ld_terms = [], [], []
    ibd_terms, dxy_terms, genrel_terms = [], [], []

    for year in ["2015", "2019", "2023"]:
        keep = get_keep_mask(year)
        kk = np.ix_(keep, keep)

        # pi (fitted, log-space, element-wise)
        pi_terms.append(_pi_log_mean_abs_diff(x[f"{year}_diversity"][keep],
                                              x0[f"{year}_diversity"][keep]))

        # Fst (fitted, off-diagonal)
        fst_terms.append(_offdiag_mean_abs_diff(x[f"{year}_fst"][kk], x0[f"{year}_fst"][kk]))

        # LD decay (fitted). Unlike the matrices above this is (bins x demes), so the mask
        # selects COLUMNS and the demes are then pooled by pair count inside the helper.
        ld_terms.append(_ld_pooled_mean_abs_diff(
            x0[f"{year}_ld_sum"], x0[f"{year}_ld_cnt"],
            x[f"{year}_ld_sum"], x[f"{year}_ld_cnt"], keep))

        # IBD slope (diagnostic) -- same real-site distances for observed and simulated
        geo = get_site_geo_distances(year)[kk]
        ibd_terms.append(abs(ibd_slope(x[f"{year}_fst"][kk], geo)
                             - ibd_slope(x0[f"{year}_fst"][kk], geo)))

        # dxy (diagnostic, off-diagonal)
        dxy_terms.append(_offdiag_mean_abs_diff(x[f"{year}_divergence"][kk],
                                                x0[f"{year}_divergence"][kk]))

        # Deliberately NOT masked: relatedness is centred on the populations present when it was
        # computed, so slicing it is not the same as recomputing on the subset (CLAUDE.md 8#2).
        genrel_terms.append(_offdiag_mean_abs_diff(x[f"{year}_relatedness"],
                                                   x0[f"{year}_relatedness"]))

    out = {
        "pi_loss": float(np.nanmean(pi_terms)),
        "fst_loss": float(np.nanmean(fst_terms)),
        "ld_loss": float(np.nanmean(ld_terms)),
        "ibd_loss": float(np.nanmean(ibd_terms)),
        "dxy_loss": float(np.nanmean(dxy_terms)),
        "genrel_loss": float(np.nanmean(genrel_terms)),
    }
    # Temporal F_c spans years rather than sitting inside one, so it is computed outside the loop
    # and carries no year mask -- the n<4 cut is applied when the pair list is built (fc_common).
    if COMPUTE_FC:
        out["fc_loss"] = _fc_loss(x0["temporal_fc"], x["temporal_fc"])
    if set(out) != set(LOSS_NAMES):
        raise ValueError(f"calculate_losses returned {sorted(out)} but LOSS_NAMES is "
                         f"{sorted(LOSS_NAMES)} -- the two must agree (CLAUDE.md 10.2)")
    return out

def getObservedData():
    '''Load empirical features: pi vector, plus dxy / Fst / genetic-relatedness matrices per year.
    Uses _read_vector (which reads every row, fixing the csv.DictReader drop of subpop 0 -- 5.7).'''
    outDict = {}
    for year in ["2015", "2019", "2023"]:
        outDict[f"{year}_diversity"] = _read_vector(Path(f"../data/empiricalStats/averaged_pi_{year}.csv"))
        outDict[f"{year}_divergence"] = _read_matrix(Path(f"../data/empiricalStats/averaged_dxy_{year}.csv"))
        outDict[f"{year}_fst"] = _read_matrix(Path(f"../data/empiricalStats/averaged_fst_{year}.csv"))
        outDict[f"{year}_relatedness"] = _read_matrix(Path(f"../data/empiricalStats/averaged_genRel_{year}.csv"))

        # LD decay (fitted, CLAUDE.md 7.5). The column order is asserted rather than assumed:
        # 7.5e CHECKED that popfile order == specifier order in all three years, but that would
        # break silently if the popfiles were ever regenerated -- and the same class of bug was
        # real in CalcGenRel.py (CLAUDE.md 4).
        s, c, labels = _read_ld_decay(Path(f"../data/empiricalStats/averaged_ldDecay_{year}.csv"))
        spec = _specifier_site_names(year)
        if labels != spec:
            raise ValueError(
                f"{year}: averaged_ldDecay column order != specifier-matrix order.\n"
                f"  ld:        {labels[:3]} ...\n  specifier: {spec[:3]} ...\n"
                f"Every per-subpop matrix in this project is in specifier order (CLAUDE.md 4); "
                f"remap the LD columns before fitting.")
        outDict[f"{year}_ld_sum"], outDict[f"{year}_ld_cnt"] = s, c

    if COMPUTE_FC:
        # The target only exists once ToUseOnBeagles/CalcTemporalFc.py has been RUN on the Beagle
        # machine (CLAUDE.md 7.9.8E). Fitting a statistic with no target is worse than not fitting
        # it, so this raises rather than degrading quietly.
        if not TEMPORAL_FC_OBS.exists():
            raise FileNotFoundError(
                f"COMPUTE_FC is on but the empirical target {TEMPORAL_FC_OBS} does not exist. Run "
                f"ToUseOnBeagles/CalcTemporalFc.py on the Beagle machine, confirm it prints "
                f"fc_common spec {fcc.spec_hash()}, and copy fc_out/averaged_temporalFc.csv here. "
                f"Do NOT hand-write this file.")
        outDict["temporal_fc"] = _read_temporal_fc(TEMPORAL_FC_OBS)
        if fcc.spec_hash() != FC_EMPIRICAL_SPEC:
            raise ValueError(
                f"fc_common spec {fcc.spec_hash()} != {FC_EMPIRICAL_SPEC}, the spec the empirical "
                f"temporal-F_c target was computed under. The two sides would use different sample "
                f"sets or filters and fc_loss would be meaningless. Re-run "
                f"ToUseOnBeagles/CalcTemporalFc.py and set FC_EMPIRICAL_SPEC to the hash it prints.")

    if ldc.spec_hash() != LD_EMPIRICAL_SPEC:
        raise ValueError(
            f"ld_common spec {ldc.spec_hash()} != {LD_EMPIRICAL_SPEC}, the spec the empirical LD "
            f"targets were computed under. The two sides are binned differently and ld_loss would "
            f"be meaningless. Either revert the ld_common change, or re-run "
            f"ToUseOnBeagles/CalculateLD.py (~7.7 h) and update LD_EMPIRICAL_SPEC.")

    return outDict
    

def sample_prior():
    '''
    Sample parameters from the prior distributions.
    '''
    return {
        "m": prior_distributions["m"].rvs(),
        "total_migration": prior_distributions["total_migration"].rvs(),
        "pop": prior_distributions["pop"].rvs(),
        "numClusters": prior_distributions["numClusters"].rvs(),
        # mutation_rate and recombination_rate are both FIXED, not sampled. They stay in the
        # returned dict (and therefore in the CSV) so the output layout is unchanged and every
        # row still records the scale it was run at (CLAUDE.md 10, 7.9.3).
        "mutation_rate": DEFAULT_MUTATION_RATE,
    }


def read_parameters_from_csv(csv_path):
    '''
    Read parameter configurations from a CSV file.
    
    Expected CSV columns: m, pop, numClusters, mutation_rate, recombination_rate
    Each row represents one simulation to run.
    
    :param csv_path: Path to the CSV file with parameters
    :return: List of dictionaries, each containing parameters for one simulation
    '''
    parameters_list = []
    
    try:
        with open(csv_path, mode='r', newline='', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            
            if reader.fieldnames is None:
                raise ValueError(f"CSV file {csv_path} is empty or has no headers")
            
            # Validate that all required columns are present (recombination_rate is optional --
            # fixed at DEFAULT_RECOMBINATION_RATE if absent, 5.4).
            required_cols = {"m", "pop", "numClusters", "mutation_rate"}
            csv_cols = set(reader.fieldnames)
            missing_cols = required_cols - csv_cols
            
            if missing_cols:
                raise ValueError(f"CSV file missing required columns: {missing_cols}. "
                                f"Required columns: {required_cols}")
            
            for row_idx, row in enumerate(reader, start=2):  # start=2 because row 1 is header
                try:
                    parameters = {
                        "m": float(row["m"]),
                        # total_migration is optional; default 0.05 so legacy CSVs without the
                        # column still run.
                        "total_migration": float(row["total_migration"]) if row.get("total_migration") not in (None, "") else 0.05,
                        "pop": int(float(row["pop"])),  # Convert to float first to handle scientific notation
                        "numClusters": int(float(row["numClusters"])),
                        "mutation_rate": float(row["mutation_rate"]),
                        "recombination_rate": float(row["recombination_rate"]) if row.get("recombination_rate") not in (None, "") else DEFAULT_RECOMBINATION_RATE
                    }
                    parameters_list.append(parameters)
                except ValueError as e:
                    print(f"Warning: Row {row_idx} in {csv_path} has invalid values: {e}")
                    continue
        
        if not parameters_list:
            raise ValueError(f"No valid parameter configurations found in {csv_path}")
        
        print(f"Loaded {len(parameters_list)} parameter configuration(s) from {csv_path}")
        return parameters_list
    
    except FileNotFoundError:
        raise FileNotFoundError(f"Input CSV file not found: {csv_path}")


def run_sims_from_csv(input_csv, output_csv="../out/abc_results.csv", simToRun=-1):
    '''
    Run ABC simulations with parameters specified in a CSV file.
    Each row in the input CSV represents one simulation.
    Detailed simulation outputs (diversities and divergences) are saved to detailed_sim_results folder.
    
    :param input_csv: Path to the CSV file with input parameters
    :param output_csv: Path to the output CSV file for results
    :param simToRun: Index of the specific simulation to run (if -1, run all)
    '''
    
    try:
        parameters_list = read_parameters_from_csv(input_csv)
    except Exception as e:
        print(f"Error reading input CSV: {e}")
        return
    
    observed_data = getObservedData()
    
    # Determine if we need to write the header.
    # CHTC's run_code.sh pre-creates this file, so exists() alone would skip the header.
    # Treat a zero-byte file as needing one.
    needs_header = not (Path(output_csv).exists() and Path(output_csv).stat().st_size > 0)
    
    # Create detailed results directory
    output_dir = Path(output_csv).parent
    detailed_results_dir = output_dir / "detailed_sim_results"
    detailed_results_dir.mkdir(parents=True, exist_ok=True)
    print(f"Detailed results will be saved to: {detailed_results_dir}")
    
    # Columns come from CSV_FIELDNAMES -- ONE definition, see the top of this module.
    # No total_loss: the combined standardized distance is built offline by abc_standardize.py.
    fieldnames = CSV_FIELDNAMES
    
    with open(output_csv, mode='a', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        # Write header if the file is new or empty (see needs_header above)
        if needs_header:
            writer.writeheader()
        
        for iteration, parameters in enumerate(parameters_list):
            if simToRun != -1 and iteration != simToRun:
                continue

            print(
                f"Running iteration {iteration + 1}/{len(parameters_list)} "
                f"with m={parameters['m']:.6g}, total_migration={parameters.get('total_migration', 0.05):.4g}, "
                f"pop={int(np.floor(parameters['pop']))}, "
                f"numClusters={parameters['numClusters'] * 33}, "
                f"mutation_rate={parameters['mutation_rate']:.6g}, "
                f"recombination_rate={parameters.get('recombination_rate', DEFAULT_RECOMBINATION_RATE):.6g}..."
            )
            
            try:
                # Run the model
                simulated_data = model(parameters)
                
                # Calculate losses
                losses = calculate_losses(observed_data, simulated_data)
                
                # Copy detailed results for this iteration
                iteration_dir = detailed_results_dir / f"run{iteration + 1}"
                iteration_dir.mkdir(parents=True, exist_ok=True)
                
                # Keep raw features so offline sigma has values, not just losses.
                _copy_raw_features(iteration_dir)
                
                # Built from the shared lists, so a new statistic cannot be silently omitted
                # (CLAUDE.md 10.2). _build_row raises on a missing key rather than letting
                # DictWriter write an empty cell.
                row = _build_row(iteration, parameters, losses)
                
                # Append to CSV
                writer.writerow(row)
                csvfile.flush()  # Ensure data is written immediately
                
                print("  " + _format_losses(losses))
                print(f"  Detailed results saved to: {iteration_dir}")
                
            except Exception as e:
                print(f"  Error in iteration {iteration}: {e}")
                continue
    
    print(f"\n=== CSV-based simulations complete ===")
    print(f"Total iterations: {len(parameters_list)}. Results saved to {output_csv}")
    print(f"Detailed simulation data saved to {detailed_results_dir}")



def run_abc_simulation(num_iterations, output_csv="../out/abc_results.csv"):
    '''
    Run ABC simulations by repeatedly sampling from the prior and computing losses.
    Results are appended to a CSV file.
    
    :param num_iterations: Number of iterations to run
    :param output_csv: Path to the output CSV file
    '''

    observed_data = getObservedData()

    # Determine if we need to write the header.
    # CHTC's run_code.sh pre-creates this file, so exists() alone would skip the header.
    # Treat a zero-byte file as needing one.
    needs_header = not (Path(output_csv).exists() and Path(output_csv).stat().st_size > 0)

    # Raw-feature store for offline standardization; transferred back from CHTC.
    detailed_results_dir = Path(output_csv).parent / "detailed_sim_results"
    detailed_results_dir.mkdir(parents=True, exist_ok=True)

    # Columns come from CSV_FIELDNAMES -- ONE definition, see the top of this module.
    # No total_loss: the combined standardized distance is built offline by abc_standardize.py.
    fieldnames = CSV_FIELDNAMES

    with open(output_csv, mode='a', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        # Write header if the file is new or empty (see needs_header above)
        if needs_header:
            writer.writeheader()

        for iteration in range(num_iterations):
            parameters = sample_prior()
            print(
                f"Running iteration {iteration + 1}/{num_iterations} "
                f"with m={parameters['m']:.6g}, total_migration={parameters.get('total_migration', 0.05):.4g}, "
                f"pop={int(np.floor(parameters['pop']))}, "
                f"numClusters={parameters['numClusters'] * 33}, "
                f"mutation_rate={parameters['mutation_rate']:.6g}, "
                f"recombination_rate={parameters.get('recombination_rate', DEFAULT_RECOMBINATION_RATE):.6g}..."
            )
            
            try:
                # Run the model
                simulated_data = model(parameters)
                
                # Calculate losses
                losses = calculate_losses(observed_data, simulated_data)

                # Keep this trial's raw features for offline standardization.
                iteration_dir = detailed_results_dir / f"run{iteration + 1}"
                iteration_dir.mkdir(parents=True, exist_ok=True)
                _copy_raw_features(iteration_dir)

                # Built from the shared lists, so a new statistic cannot be silently omitted
                # (CLAUDE.md 10.2). _build_row raises on a missing key rather than letting
                # DictWriter write an empty cell.
                row = _build_row(iteration, parameters, losses)

                # Append to CSV
                writer.writerow(row)
                csvfile.flush()  # Ensure data is written immediately

                print("  " + _format_losses(losses))

            except Exception as e:
                print(f"  Error in iteration {iteration}: {e}")
                continue

    print(f"Simulation {iteration + 1} complete. Results saved to {output_csv}")


if __name__ == "__main__":
    # Usage: python ABCAnalysisNoRedis.py <job_id> [num_trials]
    # job_id is a LABEL ONLY -- it names the output file (the submit file remaps it) and appears
    # in the log. It does NOT seed anything, so job ids may repeat or overlap between batches
    # without drawing duplicate parameters.
    # Writes one row per trial to ../out/abc_results.csv plus raw features under
    # ../out/detailed_sim_results/. Afterwards, concatenate the per-job CSVs and run
    # abc_standardize.py. (run_sims_from_csv() is the older CSV-driven path.)
    #
    # RNG: deliberately UNSEEDED (was np.random.seed(job_id) until 2026-09-07). Removed because
    # the reproducibility it offered was never real -- SLiM's forward mating, pyslim.recapitate()
    # and msprime.sim_mutations() are all unseeded, so re-running a job reproduced its PARAMETERS
    # but never its losses. And the parameters are already recorded per row in the output CSV, so
    # nothing recoverable was lost. What it cost was a live footgun: seeds silently collide
    # across batches, and a job id reused from an earlier batch re-draws that batch's exact
    # parameters while every other check looks clean.
    # Independent draws are safe: numpy seeds its global RNG from OS entropy at import, not from
    # the clock, so simultaneously-launched jobs do not collide (verified 2026-09-07).
    # KMeans is unaffected either way -- it uses the fixed Main.KMEANS_SEED.
    if len(sys.argv) < 2:
        print("Usage: python ABCAnalysisNoRedis.py <job_id> [num_trials]")
        sys.exit(1)

    job_id = int(sys.argv[1])
    num_trials = int(sys.argv[2]) if len(sys.argv) > 2 else 100

    # Fixed filename: the submit file remaps it per-process.
    output_csv = "../out/abc_results.csv"
    print(f"Job {job_id}: sampling {num_trials} prior-drawn trials -> {output_csv}")
    run_abc_simulation(num_trials, output_csv=output_csv)

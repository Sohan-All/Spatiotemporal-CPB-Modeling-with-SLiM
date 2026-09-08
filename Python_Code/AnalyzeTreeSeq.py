from pathlib import Path

import tskit
import msprime
import pyslim
import pandas as pd
import math
import numpy as np
import csv

import ld_common as ldc

# LD is ON as of 2026-09-07. Both preconditions are settled: the empirical side has been run so
# the bin scale is measured (CLAUDE.md 7.5.1), and the per-trial cost is 20.6 s -- 0.5% of a
# ~1.2 h trial (7.5.2). ld_loss is a FITTED statistic, so turning this off silently drops it and
# ABCAnalysisNoRedis.model() will raise rather than emit a trial missing it.
# Set COMPUTE_LD=0 in the environment to disable for a diagnostic run that does not need it.
COMPUTE_LD = bool(int(__import__("os").environ.get("COMPUTE_LD", "1")))

# Sites kept for pair enumeration, 1 = all. At the fitted cut (bin_lo >= 562 bp) thin=25 gives
# ~98 bp resolution, which is ample, and takes the LD phase from 5.4 min to 20.6 s (7.5.2).
# Unbiased: ld_loss moves 0.8% between thin=1 and thin=25 (7.5.2), and the empirical side shows
# the same at its own thinning boundary (7.5.1). NOT part of spec_hash -- the two sides may thin
# differently without becoming incomparable.
LD_THIN = 25


def _real_sample_sizes(year):
    '''Diploid individuals actually sequenced per site, in SPECIFIER-MATRIX ROW ORDER.

    Same source and same order as ABCAnalysisNoRedis.get_keep_mask -- popFile{year} counted by
    site name, indexed by specifier row (CLAUDE.md 4). Reimplemented here rather than imported
    because ABCAnalysisNoRedis lazily imports Main, which imports this module.
    Raises rather than silently defaulting if a specifier site is missing from the popfile.
    '''
    year = str(year)
    names = []
    with open(Path(f"../data/Genetic_Data/specifier_matrix_{year}.csv"), encoding="utf-8") as f:
        for line in f:
            if line.strip():
                names.append(line.split(",")[0].strip())

    counts = {}
    with open(Path(f"../data/Genetic_Data/popFile{year}"), encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                counts[parts[1].strip()] = counts.get(parts[1].strip(), 0) + 1

    missing = [s for s in names if s not in counts]
    if missing:
        raise ValueError(f"{year}: specifier sites absent from popFile{year}: {missing}")
    return [counts[s] for s in names]


def _subsample_nodes(ts, pop_idx, time, n_diploid, rng):
    '''Draw n_diploid whole INDIVIDUALS from one deme and return their sample nodes.

    Individuals, not nodes: the empirical unit is a diploid, and drawing loose nodes would mix
    haplotypes from different beetles (same argument as diagnostics/fst_subsample.py).
    Returns every node if the deme holds fewer individuals than requested.
    '''
    nodes = ts.samples(population=pop_idx, time=time)
    by_ind = {}
    for nd in nodes:
        ind = ts.node(nd).individual
        by_ind.setdefault(ind, []).append(nd)
    inds = np.array(sorted(by_ind), dtype=np.int64)
    if len(inds) > n_diploid:
        inds = rng.choice(inds, size=n_diploid, replace=False)
    return np.array(sorted(nd for i in inds for nd in by_ind[i]), dtype=np.int64)


def calculate_ld_decay(ts, genome_indicies, time, year, output_path, rng=None, thin=1):
    '''Per-deme LD decay (mean r^2 vs physical distance), written as bins x demes.

    THE ONE STATISTIC HERE THAT MUST BE SUBSAMPLED. The other four run on whole demes (301-714
    diploids), which 6.6 measured as correct for pi and ~1% on F_st. r^2 is different in kind: its
    small-sample bias is ~1/n_haplotypes (Hill 1981, CLAUDE.md 11), i.e. ~0.06-0.10 at the real
    n = 5-8 diploids -- plausibly larger than the signal. So the simulated demes are cut to each
    site's real n_i and the identical bias appears on both sides (invariant 1). Do NOT repoint the
    other four statistics at these sample sets: 6.6 measured that subsampling makes pi_loss worse.

    Binning, MAF filter and pair enumeration all come from ld_common, which the empirical side
    imports too. Writes ld_{year}.csv with rows = distance bins, columns = demes in
    specifier-matrix row order, mirroring the empirical file's layout (7.5).

    thin: keep every `thin`-th site, applied ONCE to the tree sequence (not per deme).

    COST, measured 2026-09-07 at POPMULT=5000 -- and the ordering is NOT what 7.5c predicted:
                                              2015 (24 demes)   3 years
      per-deme extraction, thin=1  (old)        ~21 min          53 min
      per-deme extraction, thin=25             ~20.5 min        ~52 min   <- thinning alone: 2%
      batched extraction,  thin=1               120 s            5.1 min
      batched extraction,  thin=25              7.3 s            ~18 s
    7.5c framed `thin` as THE cost knob because it timed `decay_one_pop` in isolation. In the real
    call `ts.genotype_matrix()` is 98% of the work (51.3 s vs 1.2 s per deme), and its cost tracks
    SITES AND TREES, not sample count -- so thinning inside the deme loop saved almost nothing and
    the actual fix was to stop paying the extraction once per deme. With batching, thin is a real
    but secondary knob (5.1 min -> 18 s).

    Thinning here is safe and is NOT a spec change:
      * It is UNBIASED for r^2-vs-distance -- which sites you keep does not change the expected
        r^2 at a given separation, only the precision. Measured on the empirical side, where the
        10 kb STAGES boundary drops the pair count 625x and moves mean r^2 by -2.3/-2.5/-2.7%
        (7.5.1), i.e. the decay continuing rather than a step.
      * It is therefore a SIMULATED-SIDE-ONLY knob. It is deliberately not in ld_common and not
        in spec_hash(), so changing it does NOT invalidate the 7.7 h empirical run -- unlike
        touching STAGES, which would (7.5c).
    Simulated site spacing is ~3.9 bp, so thin=T resolves distances down to ~3.9*T bp; at the
    fitted cut of bin_lo >= 562 (7.5.1) thin=25 gives ~98 bp spacing, which is ample.
    '''
    if rng is None:
        # Draws from numpy's global RNG, which the ABC path deliberately leaves UNSEEDED as of
        # 2026-09-07 (see the note in ABCAnalysisNoRedis.__main__). So this varies per trial like
        # SLiM, recapitation and the mutation overlay already did. Diagnostics that need a
        # reproducible subsample pass `rng` explicitly instead.
        rng = np.random.default_rng(np.random.randint(0, 2**31 - 1))

    n_real = _real_sample_sizes(year)
    if len(n_real) != len(genome_indicies):
        raise ValueError(f"{year}: {len(n_real)} specifier rows vs "
                         f"{len(genome_indicies)} demes -- ordering is broken (CLAUDE.md 4)")

    K = len(genome_indicies)
    sum_r2 = np.zeros((ldc.N_BINS, K))
    cnt = np.zeros((ldc.N_BINS, K))

    # THIN THE SITES ONCE, not per deme. genotype_matrix's cost is driven by sites and trees, not
    # by sample count, so thinning inside the deme loop (the obvious place) saves nothing -- see
    # the timing note below.
    if thin > 1:
        keep_sites = np.arange(0, ts.num_sites, thin)
        ts = ts.delete_sites(np.setdiff1d(np.arange(ts.num_sites), keep_sites))

    # ONE genotype_matrix call for the whole year, then slice columns per deme.
    #
    # WHY: measured 2026-09-07, this is 98% of the statistic's cost and the per-deme call was
    # paying it 24/17/20 times over the SAME tree sequence. Extraction is dominated by the site
    # count and the number of trees -- NOT by how many samples are asked for -- so one call over
    # every deme's nodes costs the same as one call over four nodes. Per year at POPMULT=5000:
    #   per-deme extraction, thin=1   53 min for 3 years   <- what this used to do
    #   batched extraction,  thin=1   5.1 min
    #   batched extraction,  thin=25  ~18 s
    # This does NOT contradict 7.5b's "do not call genotype_matrix on the full tree": that is
    # 255k sites x 70,078 samples = 1.8e10. Here it is the ~334 SUBSAMPLED nodes of one year
    # (sum of 2*n_i), i.e. 342 MB at thin=1 and 14 MB at thin=25.
    nodes_by_deme = [_subsample_nodes(ts, idx, time, n_real[k], rng)
                     for k, idx in enumerate(genome_indicies)]
    all_pos = ts.tables.sites.position.astype(np.int64)
    nonempty = [nd for nd in nodes_by_deme if len(nd)]
    if not nonempty:
        raise ValueError(f"{year}: no deme yielded any sample nodes at time={time}")
    H_year = ts.genotype_matrix(samples=np.concatenate(nonempty).astype(np.int32))

    off = 0
    for k, nodes in enumerate(nodes_by_deme):
        H = H_year[:, off:off + len(nodes)]
        off += len(nodes)
        if len(nodes) < 4:                      # need >=2 diploids for any r^2 at all
            continue
        # Per deme, on the sliced columns: a site can be biallelic within one deme and not within
        # another, so this filter must stay per-deme to keep the old semantics exactly.
        biallelic = H.max(axis=1) <= 1          # SLiMMutationModel can stack states at a site
        # Left in its native small dtype; ld_common.standardize() casts. Same reasoning as
        # CalculateLD.haplotypes() -- the cast belongs next to the MAF filter, not here.
        pos_k, H = all_pos[biallelic], H[biallelic]
        # msprime places mutations on integer positions, so several can share one. accumulate()
        # requires ascending UNIQUE positions, and a d = 0 pair is not a distance.
        pos_k, first = np.unique(pos_k, return_index=True)
        ldc.decay_one_pop(pos_k, H[first], sum_r2[:, k], cnt[:, k])

    ldc.write_decay(output_path, sum_r2, cnt, [str(i) for i in range(K)])
    return sum_r2, cnt

def calculate_diversity_and_divergence(ts, genome_indicies, time, output_diversities_path,
                                       output_divergences_path, output_fst_path,
                                       output_relatedness_path):
    '''
    Calculate diversity, divergence, Fst and genetic relatedness for given genome indices at a
    specific time from a tree sequence and output them to csv files.

    Parameters:
    ts: tree sequence object
    genome_indicies: list of population indices to sample from
    time: time point to sample
    output_diversities_path: file path to save diversities (pi, per subpop)
    output_divergences_path: file path to save divergences (d_xy, pairwise) -- diagnostic only
    output_fst_path: file path to save Fst (pairwise, relative differentiation)
    output_relatedness_path: file path to save genetic relatedness (pairwise, per-year centred)
    '''

    #Get all the nodes we are trying to sample from
    pop_samples = []
    for idx in genome_indicies:
        pop_samples.append(ts.samples(population=idx, time=time))

    K = len(pop_samples)
    pairs = [(i, j) for i in range(K) for j in range(K) if i != j]

    #Diversity (pi) per subpop -- all sample sets in a single traversal.
    diversities = np.asarray(ts.diversity(pop_samples), dtype=float)

    #Divergence (d_xy) -- ALL pairs in one traversal via indexes=. Diagonal stays 0.
    #Fst is HUDSON, 1 - Hw/d_xy with Hw = (pi_X + pi_Y)/2. Do NOT switch to ts.Fst: that returns
    #Nei/Slatkin, ~half of Hudson, and the empirical target is pixy's Weir-Cockerham, which
    #matches Hudson (CLAUDE.md 6.7, invariant 9).
    divergences = np.zeros((K, K))
    fsts = np.zeros((K, K))
    if pairs:
        div_flat = ts.divergence(pop_samples, indexes=pairs)
        for (i, j), dv in zip(pairs, div_flat):
            divergences[i, j] = dv
            # dv == 0 only if the pair has no variation at all; 0 beats a NaN in fst_loss.
            fsts[i, j] = 0.0 if dv == 0 else 1.0 - 0.5 * (diversities[i] + diversities[j]) / dv

    #Relatedness is centred across THIS year's subpops -- one call with all sample sets.
    #Never slice a larger matrix to get a smaller one (CLAUDE.md 8#2).
    gr_indexes = [(i, j) for i in range(K) for j in range(K)]
    gr = ts.genetic_relatedness(pop_samples, indexes=gr_indexes)
    relatedness = np.asarray(gr, dtype=float).reshape(K, K)

    #write the data
    with open(output_diversities_path, "w", newline="") as f:
        writer = csv.writer(f)
        for item in diversities:
            writer.writerow([item])

    with open(output_divergences_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(divergences.tolist())

    with open(output_fst_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(fsts.tolist())

    with open(output_relatedness_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(relatedness.tolist())



def analyze_tree_sequence(mutation_rate=None, recombination_rate=None, ancestral_Ne=6700):
    '''
    This function analyzes the tree sequence file generated by the SLiM simulation by
    calculating diversity and divergence statistics before and after recapitation and mutation addition.
    It outputs the results to CSV files.

    mutation_rate, recombination_rate: REQUIRED, no defaults. See the guard below.

    ancestral_Ne: effective size of the panmictic ancestral population used in recapitation.
    Fixed empirical point estimate (6700); exposed here for sensitivity analysis only, NOT
    inferred -- it is confounded with mu via pi = 4*Ne*mu. See CLAUDE.md 5.1.
    '''
    # No defaults: these set the diversity and linkage scale of every output file, and the files
    # record no scale, so a silent fallback is worse than a crash (CLAUDE.md 10.1).
    if mutation_rate is None or recombination_rate is None:
        raise ValueError(
            "analyze_tree_sequence() requires explicit mutation_rate and recombination_rate. "
            "Use ABCAnalysisNoRedis.DEFAULT_MUTATION_RATE (4.646e-7) and "
            "DEFAULT_RECOMBINATION_RATE (2.75e-6). See CLAUDE.md 6.1.1 and 6.3.")

    
    # Load the cluster_data CSV file
    cluster_data = pd.read_csv(Path("../data/cluster_data.csv"))

    assignments_2015 = cluster_data['Genome Assignment 2015']
    assignments_2019 = cluster_data['Genome Assignment 2019']
    assignments_2023 = cluster_data['Genome Assignment 2023']

    total_assignments_2015 = int(max(assignments_2015.dropna()))
    total_assignments_2019 = int(max(assignments_2019.dropna()))
    total_assignments_2023 = int(max(assignments_2023.dropna()))

    genome_indicies_2015 = [-1] * (total_assignments_2015+1)
    genome_indicies_2019 = [-1] * (total_assignments_2019+1)
    genome_indicies_2023 = [-1] * (total_assignments_2023+1)
    
    
    for i in range(len(assignments_2015)):
        if math.isnan(assignments_2015[i]) == False:
            index = int(assignments_2015[i])
            if genome_indicies_2015[index] == -1:
                genome_indicies_2015[index] = i
    for i in range(len(assignments_2019)):
        if math.isnan(assignments_2019[i]) == False:
            index = int(assignments_2019[i])
            if genome_indicies_2019[index] == -1:
                genome_indicies_2019[index] = i
    for i in range(len(assignments_2023)):
        if math.isnan(assignments_2023[i]) == False:
            index = int(assignments_2023[i])
            if genome_indicies_2023[index] == -1:
                genome_indicies_2023[index] = i

    # Load the tree sequence file
    ts = tskit.load(Path("../out/simTreeSeq.trees"))
    

    ts = pyslim.recapitate(ts, recombination_rate=recombination_rate, ancestral_Ne=ancestral_Ne)
    
    
    print("Simplifying tree sequence...")
    samplesToKeep = []
    for idx in genome_indicies_2015:
        samplesToKeep.extend(ts.samples(population=idx, time=16))
    for idx in genome_indicies_2019:
        samplesToKeep.extend(ts.samples(population=idx, time=8))
    for idx in genome_indicies_2023:
        samplesToKeep.extend(ts.samples(population=idx, time=0))
    
    
    # filter_populations=False is REQUIRED. The default renumbers surviving populations, which
    # silently misaligns the ts.samples(population=idx) queries below (CLAUDE.md 2).
    ts = ts.simplify(samples=samplesToKeep, filter_populations=False)
    
    next_id = pyslim.next_slim_mutation_id(ts)
    print("Simulating mutations...")
    ts = msprime.sim_mutations(
            ts,
            rate=mutation_rate,
            model=msprime.SLiMMutationModel(type=0, next_id=next_id),
            keep=True,
    )

    calculate_diversity_and_divergence(
        ts, genome_indicies_2023, time=0,
        output_diversities_path=Path("../data/Output_Data/diversities_2023.csv"),
        output_divergences_path=Path("../data/Output_Data/divergences_2023.csv"),
        output_fst_path=Path("../data/Output_Data/fst_2023.csv"),
        output_relatedness_path=Path("../data/Output_Data/relatedness_2023.csv"))

    calculate_diversity_and_divergence(
        ts, genome_indicies_2019, time=8,
        output_diversities_path=Path("../data/Output_Data/diversities_2019.csv"),
        output_divergences_path=Path("../data/Output_Data/divergences_2019.csv"),
        output_fst_path=Path("../data/Output_Data/fst_2019.csv"),
        output_relatedness_path=Path("../data/Output_Data/relatedness_2019.csv"))

    calculate_diversity_and_divergence(
        ts, genome_indicies_2015, time=16,
        output_diversities_path=Path("../data/Output_Data/diversities_2015.csv"),
        output_divergences_path=Path("../data/Output_Data/divergences_2015.csv"),
        output_fst_path=Path("../data/Output_Data/fst_2015.csv"),
        output_relatedness_path=Path("../data/Output_Data/relatedness_2015.csv"))

    # LD decay (CLAUDE.md 7.5). Kept OFF by default until the empirical side has been run and the
    # bin scale confirmed (6.8) -- and until its per-trial cost is measured, since this runs inside
    # every ABC trial. Enable with COMPUTE_LD=1.
    if COMPUTE_LD:
        print(f"LD decay (ld_common spec {ldc.spec_hash()}) -- "
              f"the empirical side must print the same hash...")
        for year, gidx, t in [("2023", genome_indicies_2023, 0),
                              ("2019", genome_indicies_2019, 8),
                              ("2015", genome_indicies_2015, 16)]:
            s, c = calculate_ld_decay(
                ts, gidx, time=t, year=year,
                output_path=Path(f"../data/Output_Data/ld_{year}.csv"), thin=LD_THIN)
            ldc.report_halfway(s, c, prefix=f"    {year}: ")
    

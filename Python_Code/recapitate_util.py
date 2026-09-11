"""Recapitation that works past msprime's 100-populations-per-event limit (CLAUDE.md 7.9.5).

THE WALL. `pyslim.recapitate(ts, ancestral_Ne=...)` builds a demography in which EVERY SLiM
subpopulation splits from one ancestral population in a SINGLE msprime population_split event.
msprime caps that event at 100 derived populations:

    InputError: Input error in population split: Cannot have more than 100 populations in one
    event. If this is something that you need to do, please open an issue on GitHub

which is why `numClusters` has been capped at 3 (x33 = 99 demes) -- one under the limit. Measured
2026-09-09: 99 demes recapitates fine; 200 demes runs the forward sim happily (112 s) and then
dies here. **The cap is real, it IS recapitation, and it is an msprime API limit rather than a
cost or memory limit.**

THE FIX. Merge in STAGES instead of in one event. Split the demes into groups of <= 99, merge each
group into a throwaway intermediate population at the recapitation time, then merge the
intermediates into the real ancestral population an infinitesimal step later.

WHY THIS IS EXACT, not an approximation. The second merge happens at `np.nextafter` of the first,
i.e. ~1e-13 generations later, and the intermediate populations have size 1.0 (the same value
pyslim assigns the SLiM populations, for the same reason). The probability of any coalescence
inside an intermediate is ~1e-13, so every lineage arrives in the ancestral population at the same
instant it would have under a single split, and the ancestral coalescent that follows is
identical. This is a re-expression of the same demography, not a different model.

Groups nest, so this scales: 99 groups of 99 covers 9,801 demes at two levels, and the loop keeps
adding levels beyond that.

USE. Drop-in for the pyslim call. At <= 99 demes it DELEGATES to pyslim.recapitate unchanged, so
existing behaviour is bit-identical and this cannot perturb any run at the current prior:

    import recapitate_util
    ts = recapitate_util.recapitate(ts, recombination_rate=r, ancestral_Ne=Ne, random_seed=s)

WHAT IT DOES NOT SOLVE. Deme size is `Average Count * POPMULT / numSubpops`, so raising the deme
count SHRINKS every deme -- and SLiM refuses a subpopulation that rounds to zero individuals
(`ERROR (Population::AddSubpopulation): subpopulation p38 empty`). That is a separate, earlier
wall, guarded in Main.main. At 200 demes it needs POPMULT >= ~735; at 400 it needs more than the
prior's floor of 2000. Total simulated N is ~3.33*POPMULT regardless of deme count, so more demes
does NOT mean more individuals -- it means the same individuals cut into finer pieces.
"""

import numpy as np
import msprime
import pyslim

# msprime's hard cap is 100 derived populations per population_split event. Sit one under it, for
# the same reason numClusters=3 (99 demes) has always worked.
MAX_DERIVED_PER_EVENT = 99


def _recap_time(ts):
    """The time all uncoalesced roots sit at. Same check pyslim.recapitate makes, same message."""
    root_times = {ts.node(n).time for t in ts.trees() for n in t.roots}
    if len(root_times) > 1:
        raise ValueError(
            "Not all roots are at the time recapitation expects. Simplifying before recapitating "
            "without keep_input_roots=True is the usual cause (CLAUDE.md 6.4). Observed root "
            f"times: {sorted(root_times)[:5]}")
    return root_times.pop()


def recapitate(ts, ancestral_Ne, **kwargs):
    """pyslim.recapitate, but without the 100-population ceiling.

    Delegates to pyslim below the ceiling so nothing changes for runs at numClusters <= 3.
    `kwargs` are passed to msprime.sim_ancestry exactly as pyslim passes them
    (recombination_rate, random_seed, ...).
    """
    n_pops = ts.num_populations
    if n_pops <= MAX_DERIVED_PER_EVENT:
        return pyslim.recapitate(ts, ancestral_Ne=ancestral_Ne, **kwargs)

    if "demography" in kwargs:
        raise ValueError("pass either ancestral_Ne or demography, not both")

    recap_time = _recap_time(ts)
    demography = msprime.Demography.from_tree_sequence(ts)
    # Size must be > 0 even though every one of these is merged away immediately; pyslim does
    # the same thing for the same reason.
    for pop in demography.populations:
        pop.initial_size = 1.0

    taken = {pop.name for pop in demography.populations}
    ancestral_name = "ancestral"
    while ancestral_name in taken:
        ancestral_name += "_ancestral"

    # The split must be strictly LONGER ago than the roots or it does not apply to them --
    # pyslim's comment, and it is the reason for nextafter rather than equality.
    t = np.nextafter(recap_time, 2 * recap_time)
    level, current = 0, [pop.name for pop in demography.populations]
    while len(current) > MAX_DERIVED_PER_EVENT:
        nxt = []
        for k in range(0, len(current), MAX_DERIVED_PER_EVENT):
            group = current[k:k + MAX_DERIVED_PER_EVENT]
            name = f"_recap_merge_L{level}_{k // MAX_DERIVED_PER_EVENT}"
            while name in taken:
                name += "_"
            taken.add(name)
            demography.add_population(
                name=name, initial_size=1.0,
                description="transient staging population for recapitation; merged away in ~0 time")
            demography.add_population_split(t, derived=group, ancestral=name)
            nxt.append(name)
        current = nxt
        level += 1
        t = np.nextafter(t, 2 * t)

    demography.add_population(
        name=ancestral_name, initial_size=ancestral_Ne,
        description="ancestral population simulated by msprime")
    demography.add_population_split(t, derived=current, ancestral=ancestral_name)

    kwargs["demography"] = demography
    return msprime.sim_ancestry(initial_state=ts, **kwargs)

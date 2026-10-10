"""Samples dropped from the pi / d_xy / F_st / relatedness targets (CLAUDE.md 7.9.16).

Failed genotyping, not biology: observed heterozygosity far below the ~0.27 of every other beetle,
and they read as "kin" of a failed sample in a different field (H15-2023's S212). The other four
H41-2023 beetles, and all seven H41-2015 beetles at the same coordinates, are normal. Dropping the
three keeps H41-2023 in the fit at n=4 (still >= MIN_SUBPOP_N).

Deliberately the three only (Sohan, 2026-10-09). S212 (H15-2023, e=0.631) and S32 (Alsum59-2019,
e=0.456) are milder failures left IN; add them here if that changes. This is NOT
fc_common.EXCLUDED_SAMPLES, which also drops a relative (S27) for reasons specific to F_c, and
whose list is part of the F_c spec hash. The two lists are independent on purpose.

ONE HOME: Python_Code/sample_exclusions.py. COPIED next to ToUseOnBeagles/ because the empirical
side runs on another machine -- keep the two copies identical; both print exclusion_id().
"""

EXCLUDED_SAMPLES = {
    "S221": ("H41-2023", "failed: e=0.954, observed het 0.004"),
    "S222": ("H41-2023", "failed: e=0.816, observed het 0.040"),
    "S223": ("H41-2023", "failed: e=0.558, observed het 0.106"),
}


def exclusion_id():
    """Printable record of the list, so an output log says which sample set produced it."""
    return ",".join(sorted(EXCLUDED_SAMPLES))


def popfile_counts(path):
    """{site: number of diploid samples} from a 'sampleID<TAB>site' popfile, EXCLUDED_SAMPLES
    removed. Both sides take per-site n from here, so an exclusion shrinks n identically."""
    counts = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2 and parts[0].strip() not in EXCLUDED_SAMPLES:
                counts[parts[1].strip()] = counts.get(parts[1].strip(), 0) + 1
    return counts

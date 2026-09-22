"""Is the re-founding toggle's OFF state identical to the production model?

The whole case for a toggle rather than a fork is that OFF must be the CURRENT model, not
merely something close to it -- otherwise every existing result silently changes meaning.
Same SLiM seed, same inputs, so the forward genealogy must match table for table.

Byte comparison is useless here: SLiM records the script path and every -d constant in the
provenance and top-level metadata, so those differ by construction. What must match is the
tables that carry the simulation.
"""
import sys

import numpy as np
import tskit

a = tskit.load(sys.argv[1])
b = tskit.load(sys.argv[2])

print(f"{'':<26}{'production':>14}{'toggle OFF':>14}")
for name in ("num_nodes", "num_edges", "num_individuals", "num_populations",
             "num_sites", "num_mutations", "num_trees", "sequence_length"):
    va, vb = getattr(a, name), getattr(b, name)
    flag = "" if va == vb else "   <-- DIFFERS"
    print(f"  {name:<24}{va:>14}{vb:>14}{flag}")

print("\ntable-by-table, column by column:")
ok = True
checks = [
    ("nodes.time", a.tables.nodes.time, b.tables.nodes.time),
    ("nodes.flags", a.tables.nodes.flags, b.tables.nodes.flags),
    ("nodes.population", a.tables.nodes.population, b.tables.nodes.population),
    ("nodes.individual", a.tables.nodes.individual, b.tables.nodes.individual),
    ("edges.left", a.tables.edges.left, b.tables.edges.left),
    ("edges.right", a.tables.edges.right, b.tables.edges.right),
    ("edges.parent", a.tables.edges.parent, b.tables.edges.parent),
    ("edges.child", a.tables.edges.child, b.tables.edges.child),
    ("individuals.flags", a.tables.individuals.flags, b.tables.individuals.flags),
    ("individuals.location", a.tables.individuals.location, b.tables.individuals.location),
]
for label, xa, xb in checks:
    same = len(xa) == len(xb) and np.array_equal(np.asarray(xa), np.asarray(xb))
    ok &= same
    print(f"  {label:<24}{'IDENTICAL' if same else 'DIFFERS':>14}")

# Individual metadata carries SLiM's pedigree ids and parent ids -- the actual mating record.
ma = [i.metadata for i in a.individuals()]
mb = [i.metadata for i in b.individuals()]
same_meta = ma == mb
ok &= same_meta
print(f"  {'individual metadata':<24}{'IDENTICAL' if same_meta else 'DIFFERS':>14}"
      f"   (SLiM pedigree + parent ids)")

print("\n" + ("PASS - the toggle's OFF state IS the production model."
               if ok else
               "FAIL - OFF is NOT the production model. A toggle is unsafe until this is fixed."))
sys.exit(0 if ok else 1)

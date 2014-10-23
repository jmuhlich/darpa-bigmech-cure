"""Variant of egfr_example where adaptors only bind one EGFR monomer in the
dimer, to make it feasible to visualize the species and reactions."""

from pysb import *
from pysb.core import ComponentSet

import egfr_example

Model(base=egfr_example.model)

rules = ComponentSet()
params = ComponentSet()
for r in L_bind_R, R_R_dimerize, R_Y1068_transphos, R_Y1148_transphos:
    rules.add(r)
    params.add(r.rate_forward)
    if r.rate_reverse:
        params.add(r.rate_reverse)
model.rules -= rules
model.parameters -= params

Initial(EGFR(l=None, r=1, Y1068='p', Y1148='p') %
        EGFR(l=None, r=1, Y1068='u', Y1148='u'),
        EGFR_tot)

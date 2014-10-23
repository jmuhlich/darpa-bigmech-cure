"Print out expanded MAPK cascade rules for egfr_example."

from egfr_example import (kf_bind, kr_bind, kcat_phos, kcat_dephos,
                          Ras, Raf, MEK, ERK, PP2A, MKP)
from pysb.macros import catalyze_state
import pysb.core


# Build handy rate "sets"
klist_bind = [kf_bind, kr_bind]
klist_phos = klist_bind + [kcat_phos]
klist_dephos = klist_bind + [kcat_dephos]

def mapk_single(kinase, pptase, substrate, site):
    """Kinase phos/dephosphorylation."""
    ppt_substrate = substrate()
    if 'k' in ppt_substrate.monomer.sites:
        # Ensure substrates which are themselves kinases don't get
        # dephosphorylated while they are bound to *their* substrate.
        ppt_substrate = ppt_substrate(k=None)
    components = catalyze_state(kinase, 'k',
                                substrate, site, site, 'u', 'p',
                                klist_phos)
    components |= catalyze_state(pptase, 'ppt',
                                 ppt_substrate, site, site, 'p', 'u',
                                 klist_dephos)
    return components

def mapk_double(kinase, pptase, substrate, site1, site2):
    """Distributive + ordered double kinase phos/dephosphorylation."""
    components = mapk_single(kinase, pptase, substrate({site2: 'u'}), site1)
    components |= mapk_single(kinase, pptase, substrate({site1: 'p'}), site2)
    return components


# Disable self-export so we don't get a name conflict error when these rules are
# already in the model.
pysb.core.SelfExporter.do_export = False
mapk_rules = pysb.core.ComponentSet()
mapk_rules |= mapk_single(Ras(S1S2='GTP'), PP2A, Raf, 'X')
mapk_rules |= mapk_double(Raf(X='p'), PP2A, MEK, 'S218', 'S222')
mapk_rules |= mapk_double(MEK(S218='p', S222='p'), MKP, ERK, 'T185', 'Y187')

for r in mapk_rules:
    print r

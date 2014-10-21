"""Simple EGFR -> ERK signaling model.

EGF->Sos elements taken directly from BNG example model "egfr_net.bngl".
Sos->ERK elements adapted from Chen 2008 doi:10.1038/msb.2008.74 .
"""

from pysb import *

Model()


Parameter('EGF_tot', 1.2e6)        # molecule counts
Parameter('EGFR_tot', 1.8e5)       # molecule counts
Parameter('Grb2_tot', 1.0e5)       # molecule counts
Parameter('Shc_tot', 2.7e5)        # molecule counts
Parameter('Sos_tot', 1.3e4)        # molecule counts
Parameter('Grb2_Sos_tot', 4.9e4)   # molecule counts

Parameter("kp1", 1.667e-06)  # ligand-monomer binding (scaled), units: /molecule/s
Parameter("km1", 0.06)       # ligand-monomer dissociation, units: /s
Parameter("kp2", 5.556e-06)  # aggregation of bound monomers (scaled), units: /molecule/s
Parameter("km2", 0.1)        # dissociation of bound monomers, units: /s
Parameter("kp3", 0.5)        # dimer transphosphorylation, units: /s
Parameter("km3", 4.505)      # dimer dephosphorylation, units: /s
Parameter("kp14", 3)         # Shc transphosphorylation, units: /s
Parameter("km14", 0.03)      # Shc dephosphorylation, units: /s
Parameter("km16", 0.005)     # Shc cytosolic dephosphorylation, units: /s
Parameter("kp9", 8.333e-07)  # binding of Grb2 to receptor (scaled), units: /molecule/s
Parameter("km9", 0.05)       # dissociation of Grb2 from receptor, units: /s
Parameter("kp10", 5.556e-06) # binding of Sos to receptor (scaled), units: /molecule/s
Parameter("km10", 0.06)      # dissociation of Sos from receptor, units: /s
Parameter("kp11", 1.25e-06)  # binding of Grb2-Sos to receptor (scaled), units: /molecule/s
Parameter("km11", 0.03)      # diss. of Grb2-Sos from receptor, units: /s
Parameter("kp13", 2.5e-05)   # binding of Shc to receptor (scaled), units: /molecule/s
Parameter("km13", 0.6)       # diss. of Shc from receptor, units: /s
Parameter("kp15", 2.5e-07)   # binding of pShc to receptor (scaled), units: /molecule/s
Parameter("km15", 0.3)       # diss. of pShc from receptor, units: /s
Parameter("kp17", 1.667e-06) # binding of Grb2 to RP-pShc (scaled), units: /molecule/s
Parameter("km17", 0.1)       # diss. of Grb2 from RP-pShc, units: /s
Parameter("kp18", 2.5e-07)   # binding of pShc-Grb2 to receptor (scaled), units: /molecule/s
Parameter("km18", 0.3)       # diss. of pShc-Grb2 from receptor, units: /s
Parameter("kp19", 5.556e-06) # binding of Sos to RP-pShc-Grb2 (scaled), units: /molecule/s
Parameter("km19", 0.0214)    # diss. of Sos from RP-pShc-Grb2, units: /s
Parameter("kp20", 6.667e-08) # binding of pShc-Grb2-Sos to receptor (scaled), units: /molecule/s
Parameter("km20", 0.12)      # diss. of pShc-Grb2-Sos from receptor, units: /s
Parameter("kp24", 5e-06)     # binding of Grb2-Sos to RP-pShc (scaled), units: /molecule/s
Parameter("km24", 0.0429)    # diss. of Grb2-Sos from RP-pShc, units: /s
Parameter("kp21", 1.667e-06) # binding of pShc to Grb2 in cytosol (scaled), units: /molecule/s
Parameter("km21", 0.01)      # diss. of Grb2 and SchP in cytosol, units: /s
Parameter("kp23", 1.167e-05) # binding of pShc to Grb2-Sos in cytosol (scaled), units: /molecule/s
Parameter("km23", 0.1)       # diss. of Grb2-Sos and SchP in cytosol, units: /s
Parameter("kp12", 5.556e-08) # binding of Grb2 to Sos in cytosol (scaled), units: /molecule/s
Parameter("km12", 0.0015)    # diss. of Grb2 and Sos in cytosol, units: /s
Parameter("kp22", 1.667e-05) # binding of pShc-Grb2 to Sos in cytosol (scaled), units: /molecule/s
Parameter("km22", 0.064)     # diss. of pShc-Grb2 and Sos in cytosol, units: /s

# Generic rates for MAPK cascade kinase/phosphatase binding, unbinding and catalysis.
Parameter('kf_bind', 1e-5)
Parameter('kr_bind', 1e-1)
Parameter('kcat_phos', 1e-1)
Parameter('kcat_dephos', 3e-3)


# General monomer site naming conventions:
#   Ynn: tyrosine residue nn, takes state 'u' for unphosphorylated or
#        'p' for phosphorylated
#   Snn: serine, ''
#   Tnn: threonine, ''
#   k: kinase catalytic domain (forms contact with substrate)
#   ppt: phosphatase catalytic domain

# r: receptor binding
Monomer('EGF', ['r'])
# l: ligand binding; r: receptor dimerization
Monomer('EGFR', ['l','r','Y1068','Y1148'], {'Y1068': ['u','p'], 'Y1148': ['u','p']})
# PTB: phosphotyrosine (EGFR) binding; Y*: tyrosines
Monomer('Shc', ['PTB','Y317'], {'Y317': ['u','p']})
# SH2,SH3: binding domains
Monomer('Grb2', ['SH2','SH3'])
# pr: proline-rich (SH3 recognition motif)
Monomer('Sos', ['pr'])
Monomer('Ras', ['k'])
# X: unspecified residue whose phosphorylation we will treat as the regulator of
#    Raf's kinase activity (Ras activation of Raf is actually quite complicated
#    and not perfectly understood)
Monomer('Raf', ['X','k'], {'X': ['u','p']})
Monomer('MEK', ['S218','S222','k'], {'S218': ['u','p'], 'S222': ['u','p']})
Monomer('ERK', ['T185','Y187'], {'T185': ['u','p'], 'Y187': ['u','p']})
Monomer('PP2A', ['ppt'])
Monomer('MKP', ['ppt'])


# Ligand-receptor binding (ligand-monomer)
Rule("L_bind_R", EGFR(l=None,r=None) + EGF(r=None) <> EGFR(l=1,r=None) % EGF(r=1), kp1, km1)

# Receptor-aggregation 
Rule("R_R_dimerize",
     EGFR(l=ANY,r=None) + EGFR(l=ANY,r=None) <>
     EGFR(l=ANY,r=3) % EGFR(l=ANY,r=3), kp2, km2)

# Transphosphorylation of EGFR by RTK
Rule("R_Y1068_transphos", EGFR(r=ANY,Y1068='u') >> EGFR(r=ANY,Y1068='p'), kp3)
Rule("R_Y1148_transphos", EGFR(r=ANY,Y1148='u') >> EGFR(r=ANY,Y1148='p'), kp3)

# Dephosphorylayion
Rule("R_Y1068_dephos", EGFR(Y1068='p') >> EGFR(Y1068='u'), km3)
Rule("R_Y1148_dephos", EGFR(Y1148='p') >> EGFR(Y1148='u'), km3)

# Shc transphosph
Rule("Shc_transphos",
     EGFR(r=ANY,Y1148=('p',1)) % Shc(PTB=1,Y317='u') >>
     EGFR(r=ANY,Y1148=('p',1)) % Shc(PTB=1,Y317='p'), kp14)
Rule("Shc_bound_dephos", Shc(PTB=ANY,Y317='p') >> Shc(PTB=ANY,Y317='u'), km14)

# Y1068 activity
Rule("R_Y1068_bind_Grb2",
     EGFR(Y1068='p') + Grb2(SH2=None,SH3=None) <>
     EGFR(Y1068=('p',1)) % Grb2(SH2=1,SH3=None), kp9, km9)
Rule("R_Y1068_bind_Grb2_Sos",
     EGFR(Y1068='p') + Grb2(SH2=None,SH3=2) % Sos(pr=2) <>
     EGFR(Y1068=('p',1)) % Grb2(SH2=1,SH3=2) % Sos(pr=2), kp11, km11)
Rule("R_Y1068_Grb2_bind_Sos",
     EGFR(Y1068=('p',1)) % Grb2(SH2=1,SH3=None) + Sos(pr=None) <>
     EGFR(Y1068=('p',1)) % Grb2(SH2=1,SH3=2) % Sos(pr=2), kp10, km10)

# Y1148 activity
Rule("R_1148_bind_Shc",
     EGFR(Y1148='p') + Shc(PTB=None,Y317='u') <>
     EGFR(Y1148=('p',1)) % Shc(PTB=1,Y317='u'), kp13, km13)
Rule("R_1148_bind_pShc",
     EGFR(Y1148='p') + Shc(PTB=None,Y317='p') <>
     EGFR(Y1148=('p',1)) % Shc(PTB=1,Y317='p'), kp15, km15)
Rule("R_1148_bind_pShc_Grb2",
     EGFR(Y1148='p') + Shc(PTB=None,Y317=('p',1)) % Grb2(SH2=1,SH3=None) <>
     EGFR(Y1148=('p',2)) % Shc(PTB=2,Y317=('p',1)) % Grb2(SH2=1,SH3=None), kp18, km18)
Rule("R_1148_bind_pShc_Grb2_Sos",
     EGFR(Y1148='p') + Shc(PTB=None,Y317=('p',1)) % Grb2(SH2=1,SH3=3) % Sos(pr=3) <>
     EGFR(Y1148=('p',2)) % Shc(PTB=2,Y317=('p',1)) % Grb2(SH2=1,SH3=3) % Sos(pr=3),
     kp20, km20)
Rule("R_1148_pShc_bind_Grb2",
     EGFR(Y1148=('p',1)) % Shc(PTB=1,Y317='p') + Grb2(SH2=None,SH3=None) <>
     EGFR(Y1148=('p',1)) % Shc(PTB=1,Y317=('p',2)) % Grb2(SH2=2,SH3=None), kp17, km17)
Rule("R_1148_pShc_bind_Grb2_Sos",
     EGFR(Y1148=('p',1)) % Shc(PTB=1,Y317='p') + Grb2(SH2=None,SH3=3) % Sos(pr=3) <> 
     EGFR(Y1148=('p',1)) % Shc(PTB=1,Y317=('p',2)) % Grb2(SH2=2,SH3=3) % Sos(pr=3),
     kp24, km24)
Rule("R_1148_pShc_Grb2_bind_Sos",
     EGFR(Y1148=('p',1)) % Shc(PTB=1,Y317=('p',3)) % Grb2(SH2=3,SH3=None) + Sos(pr=None) <>
     EGFR(Y1148=('p',1)) % Shc(PTB=1,Y317=('p',3)) % Grb2(SH2=3,SH3=4) % Sos(pr=4), kp19, km19)

# Cytosolic 
Rule("pShc_bind_Grb2",
     Shc(PTB=None,Y317='p') + Grb2(SH2=None,SH3=None) <>
     Shc(PTB=None,Y317=('p',1)) % Grb2(SH2=1,SH3=None), kp21, km21)
Rule("pShc_bind_Grb2_Sos",
     Shc(PTB=None,Y317='p') + Grb2(SH2=None,SH3=1) % Sos(pr=1) <>
     Shc(PTB=None,Y317=('p',2)) % Grb2(SH2=2,SH3=1) % Sos(pr=1), kp23, km23)
Rule("Shc_free_dephos", Shc(PTB=None,Y317='p') >> Shc(PTB=None,Y317='u'), km16)
Rule("Grb2_bind_Sos",
     Grb2(SH2=None,SH3=None) + Sos(pr=None) <>
     Grb2(SH2=None,SH3=1) % Sos(pr=1), kp12, km12)
Rule("pShc_Grb2_bind_Sos",
     Shc(PTB=None,Y317=('p',2)) % Grb2(SH2=2,SH3=None) + Sos(pr=None) <> 
     Shc(PTB=None,Y317=('p',2)) % Grb2(SH2=2,SH3=3) % Sos(pr=3), kp22, km22)

# FIXME add Sos activation of Ras

## MAPK section generated by expand_mapk_rules.py .
# (Typically one would just put those macro calls here, but I decided to expand
# the full rule set instead. If people would rather have the macros here we can
# easily change that. -JLM)

# Raf activation
Rule('bind_Ras_Rafu_to_RasRaf', Ras(k=None) + Raf(X='u') <> Ras(k=1) % Raf(X=('u', 1)), kf_bind, kr_bind)
Rule('catalyze_RasRaf_to_Ras_Rafp', Ras(k=1) % Raf(X=('u', 1)) >> Ras(k=None) + Raf(X='p'), kcat_phos)
Rule('bind_PP2A_Rafp_to_PP2ARaf', PP2A(ppt=None) + Raf(X='p', k=None) <> PP2A(ppt=1) % Raf(X=('p', 1), k=None), kf_bind, kr_bind)
Rule('catalyze_PP2ARaf_to_PP2A_Rafu', PP2A(ppt=1) % Raf(X=('p', 1), k=None) >> PP2A(ppt=None) + Raf(X='u', k=None), kcat_dephos)

# MEK activation
Rule('bind_Rafp_MEKuu_to_RafpMEKu', Raf(X='p', k=None) + MEK(S218='u', S222='u') <> Raf(X='p', k=1) % MEK(S218=('u', 1), S222='u'), kf_bind, kr_bind)
Rule('catalyze_RafpMEKu_to_Rafp_MEKpu', Raf(X='p', k=1) % MEK(S218=('u', 1), S222='u') >> Raf(X='p', k=None) + MEK(S218='p', S222='u'), kcat_phos)
Rule('bind_PP2A_MEKpu_to_PP2AMEKu', PP2A(ppt=None) + MEK(S218='p', S222='u', k=None) <> PP2A(ppt=1) % MEK(S218=('p', 1), S222='u', k=None), kf_bind, kr_bind)
Rule('catalyze_PP2AMEKu_to_PP2A_MEKuu', PP2A(ppt=1) % MEK(S218=('p', 1), S222='u', k=None) >> PP2A(ppt=None) + MEK(S218='u', S222='u', k=None), kcat_dephos)
Rule('bind_Rafp_MEKpu_to_RafpMEKp', Raf(X='p', k=None) + MEK(S218='p', S222='u') <> Raf(X='p', k=1) % MEK(S218='p', S222=('u', 1)), kf_bind, kr_bind)
Rule('catalyze_RafpMEKp_to_Rafp_MEKpp', Raf(X='p', k=1) % MEK(S218='p', S222=('u', 1)) >> Raf(X='p', k=None) + MEK(S218='p', S222='p'), kcat_phos)
Rule('bind_PP2A_MEKpp_to_PP2AMEKp', PP2A(ppt=None) + MEK(S218='p', S222='p', k=None) <> PP2A(ppt=1) % MEK(S218='p', S222=('p', 1), k=None), kf_bind, kr_bind)
Rule('catalyze_PP2AMEKp_to_PP2A_MEKpu', PP2A(ppt=1) % MEK(S218='p', S222=('p', 1), k=None) >> PP2A(ppt=None) + MEK(S218='p', S222='u', k=None), kcat_dephos)

# ERK activation
Rule('bind_MEKpp_ERKuu_to_MEKppERKu', MEK(S218='p', S222='p', k=None) + ERK(T185='u', Y187='u') <> MEK(S218='p', S222='p', k=1) % ERK(T185=('u', 1), Y187='u'), kf_bind, kr_bind)
Rule('catalyze_MEKppERKu_to_MEKpp_ERKup', MEK(S218='p', S222='p', k=1) % ERK(T185=('u', 1), Y187='u') >> MEK(S218='p', S222='p', k=None) + ERK(T185='p', Y187='u'), kcat_phos)
Rule('bind_MKP_ERKup_to_MKPERKu', MKP(ppt=None) + ERK(T185='p', Y187='u') <> MKP(ppt=1) % ERK(T185=('p', 1), Y187='u'), kf_bind, kr_bind)
Rule('catalyze_MKPERKu_to_MKP_ERKuu', MKP(ppt=1) % ERK(T185=('p', 1), Y187='u') >> MKP(ppt=None) + ERK(T185='u', Y187='u'), kcat_dephos)
Rule('bind_MEKpp_ERKup_to_MEKppERKp', MEK(S218='p', S222='p', k=None) + ERK(T185='p', Y187='u') <> MEK(S218='p', S222='p', k=1) % ERK(T185='p', Y187=('u', 1)), kf_bind, kr_bind)
Rule('catalyze_MEKppERKp_to_MEKpp_ERKpp', MEK(S218='p', S222='p', k=1) % ERK(T185='p', Y187=('u', 1)) >> MEK(S218='p', S222='p', k=None) + ERK(T185='p', Y187='p'), kcat_phos)
Rule('bind_MKP_ERKpp_to_MKPERKp', MKP(ppt=None) + ERK(T185='p', Y187='p') <> MKP(ppt=1) % ERK(T185='p', Y187=('p', 1)), kf_bind, kr_bind)
Rule('catalyze_MKPERKp_to_MKP_ERKup', MKP(ppt=1) % ERK(T185='p', Y187=('p', 1)) >> MKP(ppt=None) + ERK(T185='p', Y187='u'), kcat_dephos)


Initial(EGF(r=None), EGF_tot)
Initial(Grb2(SH2=None, SH3=None), Grb2_tot)
Initial(Shc(PTB=None, Y317='u'), Shc_tot)
Initial(Sos(pr=None), Sos_tot)
Initial(EGFR(l=None, r=None, Y1068='u', Y1148='u'), EGFR_tot)
Initial(Grb2(SH2=None,SH3=1) % Sos(pr=1), Grb2_Sos_tot)
# FIXME move these parameters to separate lines above
Initial(Ras(k=None), Parameter('Ras_0', 6e4))
Initial(Raf(X='u', k=None), Parameter('Raf_0', 7e4))
Initial(MEK(S218='u', S222='u', k=None), Parameter('MEK_0', 3e6))
Initial(ERK(T185='u', Y187='u'), Parameter('ERK_0', 7e5))
Initial(PP2A(ppt=None), Parameter('PP2A_0', 2e5))
Initial(MKP(ppt=None), Parameter('MKP_0', 1.7e4))


Observable("Dimers", EGFR(r=1) % EGFR(r=1))
Observable("Sos_act",
           Shc(PTB=ANY,Y317=('p',2)) % Grb2(SH2=2,SH3=3) % Sos(pr=3) +
           EGFR(Y1068=('p',1)) % Grb2(SH2=1,SH3=2) % Sos(pr=2))
Observable("RP", EGFR(Y1068=('p',WILD)) + EGFR(Y1148=('p',WILD)))
Observable("Shc_Grb", Shc(Y317=('p',1)) % Grb2(SH2=1))
Observable("Shc_Grb_Sos", Shc(Y317=('p',1)) % Grb2(SH2=1,SH3=2) % Sos(pr=2))
Observable("R_Grb2", EGFR(Y1068=('p',1)) % Grb2(SH2=1))
Observable("R_Shc", EGFR(Y1148=('p',1)) % Shc(PTB=1,Y317='u'))
Observable("R_ShcP", EGFR(Y1148=('p',1)) % Shc(PTB=1,Y317=('p',WILD)))
Observable("ShcP", Shc(Y317=('p',WILD)))
Observable("R_G_S", EGFR(Y1068=('p',1)) % Grb2(SH2=1,SH3=2) % Sos(pr=2))
Observable("R_S_G_S",
           EGFR(Y1148=('p',1)) % Shc(PTB=1,Y317=('p',2)) % Grb2(SH2=2,SH3=3) %
           Sos(pr=3))
Observable("Efgr_total", EGFR)
Observable("Shc_total", Shc)
Observable("Sos_total", Sos)
Observable("Grb2_total", Grb2)
Observable('MEKPP', MEK(S218='p', S222='p'))
Observable('ERKPP', ERK(T185='p', Y187='p'))


if __name__ == '__main__':
    print __doc__, "\n", model
    print "\nNOTE: This model code is designed to be imported and programatically " \
        "manipulated,\nnot executed directly. The above output is merely a " \
        "diagnostic aid."

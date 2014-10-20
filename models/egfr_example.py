"""Simple EGFR -> ERK signaling model.

EGF->Sos elements taken directly from BNG example model "egfr_net.bngl".
Sos->ERK elements adapted from Chen 2008 doi:10.1038/msb.2008.74 .
"""

from pysb import *

Model()


Parameter('egf_tot', 1.2e6)        # molecule counts
Parameter('egfr_tot', 1.8e5)       # molecule counts
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


Monomer('egf', ['r'])
Monomer('egfr', ['l','r','Y1068','Y1148'],
        {'Y1068': ['Y','pY'], 'Y1148': ['Y','pY']})
Monomer('Shc', ['PTB','Y317'], {'Y317': ['Y','pY']})
Monomer('Grb2', ['SH2','SH3'])
Monomer('Sos', ['dom'])


# Ligand-receptor binding (ligand-monomer)
Rule("L_bind_R", egfr(l=None,r=None) + egf(r=None) <> egfr(l=1,r=None) % egf(r=1), kp1, km1)

# Receptor-aggregation 
Rule("R_R_dimerize",
     egfr(l=ANY,r=None) + egfr(l=ANY,r=None) <>
     egfr(l=ANY,r=3) % egfr(l=ANY,r=3), kp2, km2)

# Transphosphorylation of egfr by RTK
Rule("R_Y1068_transphos", egfr(r=ANY,Y1068="Y") >> egfr(r=ANY,Y1068="pY"), kp3)
Rule("R_Y1148_transphos", egfr(r=ANY,Y1148="Y") >> egfr(r=ANY,Y1148="pY"), kp3)

# Dephosphorylayion
Rule("R_Y1068_dephos", egfr(Y1068="pY") >> egfr(Y1068="Y"), km3)
Rule("R_Y1148_dephos", egfr(Y1148="pY") >> egfr(Y1148="Y"), km3)

# Shc transphosph
Rule("Shc_transphos",
     egfr(r=ANY,Y1148=("pY",1)) % Shc(PTB=1,Y317="Y") >>
     egfr(r=ANY,Y1148=("pY",1)) % Shc(PTB=1,Y317="pY"), kp14)
Rule("Shc_bound_dephos", Shc(PTB=ANY,Y317="pY") >> Shc(PTB=ANY,Y317="Y"), km14)

# Y1068 activity
Rule("R_Y1068_bind_Grb2",
     egfr(Y1068="pY") + Grb2(SH2=None,SH3=None) <>
     egfr(Y1068=("pY",1)) % Grb2(SH2=1,SH3=None), kp9, km9)
# FIXME why not explicitly mention Sos rather than use ANY?
Rule("R_Y1068_bind_Grb2_Sos",
     egfr(Y1068="pY") + Grb2(SH2=None,SH3=ANY) <>
     egfr(Y1068=("pY",1)) % Grb2(SH2=1,SH3=ANY), kp11, km11)
Rule("R_Y1068_Grb2_bind_Sos",
     egfr(Y1068=("pY",1)) % Grb2(SH2=1,SH3=None) + Sos(dom=None) <>
     egfr(Y1068=("pY",1)) % Grb2(SH2=1,SH3=2) % Sos(dom=2), kp10, km10)

# Y1148 activity
Rule("R_1148_bind_Shc",
     egfr(Y1148="pY") + Shc(PTB=None,Y317="Y") <>
     egfr(Y1148=("pY",1)) % Shc(PTB=1,Y317="Y"), kp13, km13)
Rule("R_1148_bind_pShc",
     egfr(Y1148="pY") + Shc(PTB=None,Y317="pY") <>
     egfr(Y1148=("pY",1)) % Shc(PTB=1,Y317="pY"), kp15, km15)
Rule("R_1148_bind_pShc_Grb2",
     egfr(Y1148="pY") + Shc(PTB=None,Y317=("pY",1)) % Grb2(SH2=1,SH3=None) <>
     egfr(Y1148=("pY",2)) % Shc(PTB=2,Y317=("pY",1)) % Grb2(SH2=1,SH3=None), kp18, km18)
Rule("R_1148_bind_pShc_Grb2_Sos",
     egfr(Y1148="pY") + Shc(PTB=None,Y317=("pY",1)) % Grb2(SH2=1,SH3=3) % Sos(dom=3) <>
     egfr(Y1148=("pY",2)) % Shc(PTB=2,Y317=("pY",1)) % Grb2(SH2=1,SH3=3) % Sos(dom=3),
     kp20, km20)
Rule("R_1148_pShc_bind_Grb2",
     egfr(Y1148=("pY",1)) % Shc(PTB=1,Y317="pY") + Grb2(SH2=None,SH3=None) <>
     egfr(Y1148=("pY",1)) % Shc(PTB=1,Y317=("pY",2)) % Grb2(SH2=2,SH3=None), kp17, km17)
Rule("R_1148_pShc_bind_Grb2_Sos",
     egfr(Y1148=("pY",1)) % Shc(PTB=1,Y317="pY") + Grb2(SH2=None,SH3=3) % Sos(dom=3) <> 
     egfr(Y1148=("pY",1)) % Shc(PTB=1,Y317=("pY",2)) % Grb2(SH2=2,SH3=3) % Sos(dom=3),
     kp24, km24)
# FIXME why not explicitly mention egfr-pY1148 rather than use ANY?
Rule("R_1148_pShc_Grb2_bind_Sos",
     Shc(PTB=ANY,Y317=("pY",2)) % Grb2(SH2=2,SH3=None) + Sos(dom=None) <>
     Shc(PTB=ANY,Y317=("pY",2)) % Grb2(SH2=2,SH3=3) % Sos(dom=3), kp19, km19)

# Cytosolic 
Rule("pShc_bind_Grb2",
     Shc(PTB=None,Y317="pY") + Grb2(SH2=None,SH3=None) <>
     Shc(PTB=None,Y317=("pY",1)) % Grb2(SH2=1,SH3=None), kp21, km21)
# FIXME why not explicitly mention Sos rather than use ANY?
Rule("pShc_bind_Grb2_Sos",
     Shc(PTB=None,Y317="pY") + Grb2(SH2=None,SH3=ANY) <>
     Shc(PTB=None,Y317=("pY",1)) % Grb2(SH2=1,SH3=ANY), kp23, km23)
Rule("Shc_free_dephos", Shc(PTB=None,Y317="pY") >> Shc(PTB=None,Y317="Y"), km16)
Rule("Grb2_bind_Sos",
     Grb2(SH2=None,SH3=None) + Sos(dom=None) <>
     Grb2(SH2=None,SH3=1) % Sos(dom=1), kp12, km12)
Rule("pShc_Grb2_bind_Sos",
     Shc(PTB=None,Y317=("pY",2)) % Grb2(SH2=2,SH3=None) + Sos(dom=None) <> 
     Shc(PTB=None,Y317=("pY",2)) % Grb2(SH2=2,SH3=3) % Sos(dom=3), kp22, km22)


Initial(egf(r=None), egf_tot)
Initial(Grb2(SH2=None, SH3=None), Grb2_tot)
Initial(Shc(PTB=None, Y317='Y'), Shc_tot)
Initial(Sos(dom=None), Sos_tot)
Initial(egfr(l=None, r=None, Y1068='Y', Y1148='Y'), egfr_tot)
Initial(Grb2(SH2=None,SH3=1) % Sos(dom=1), Grb2_Sos_tot)


# Observables
"""
    Molecules     Dimers       egfr(r!1).egfr(r!1)
    Molecules     Sos_act      Shc(PTB!+,Y317~pY!2).Grb2(SH2!2,SH3!3).Sos(dom!3), egfr(Y1068~pY!1).Grb2(SH2!1,SH3!2).Sos(dom!2)
    Molecules     RP           egfr(Y1068~pY!?), egfr(Y1148~pY!?)
    Molecules     Shc_Grb      Shc(Y317~pY!1).Grb2(SH2!1)
    Molecules     Shc_Grb_Sos  Shc(Y317~pY!1).Grb2(SH2!1,SH3!2).Sos(dom!2)
    Molecules     R_Grb2       egfr(Y1068~pY!1).Grb2(SH2!1)
    Molecules     R_Shc        egfr(Y1148~pY!1).Shc(PTB!1,Y317~Y)
    Molecules     R_ShcP       egfr(Y1148~pY!1).Shc(PTB!1,Y317~pY!?)
    Molecules     ShcP         Shc(Y317~pY!?)
    Molecules     R_G_S        egfr(Y1068~pY!1).Grb2(SH2!1,SH3!2).Sos(dom!2)
    # Strong differences are seen for R_G_S in comparison with path model
    Molecules     R_S_G_S      egfr(Y1148~pY!1).Shc(PTB!1,Y317~pY!2).Grb2(SH2!2,SH3!3).Sos(dom!3)

    Molecules     Efgr_total  egfr
    Molecules     Shc_total   Shc
    Molecules     Sos_total   Sos
    Molecules     Grb2_total  Grb2
"""


if __name__ == '__main__':
    print __doc__, "\n", model
    print "\nNOTE: This model code is designed to be imported and programatically " \
        "manipulated,\nnot executed directly. The above output is merely a " \
        "diagnostic aid."
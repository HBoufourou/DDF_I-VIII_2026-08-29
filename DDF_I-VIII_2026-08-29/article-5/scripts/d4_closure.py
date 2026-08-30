import numpy as np

# ============ ENTREES (toutes issues du corpus) ============
hbar_c_eV_um = 0.1973          # eV.um
ell_um   = 5.81                # Article II : longueur du puits lineaire
eps_eV   = 24e-3               # Article II/IV : epsilon = F*ell = 24 meV = mu (ombre)
MS_GeV   = 6000.0              # Article IV : M_S ~ 3.5-9.8 TeV, mediane
MPl_GeV  = 1.22e19
eta_B    = 6.1e-10             # Planck
m_p      = 0.938
MX_GeV   = 5.0                 # Article VII : fenetre benchmark 4-6 GeV

# ============ D3 : mu_w = (d_y Phi)/M  =  F / M ============
F_eV2 = eps_eV * (hbar_c_eV_um/ell_um)     # pente du puits en eV^2
print(f"D3  pente du puits F = eps/ell = {F_eV2:.3e} eV^2  (= {eps_eV*1e3:.0f} meV x {hbar_c_eV_um/ell_um*1e3:.1f} meV)")
mu_w_meV  = F_eV2/eps_eV                    # normalisation a l'echelle meV (le secteur lui-meme)
mu_w_TeV  = F_eV2/(MS_GeV*1e9)              # normalisation a l'echelle du mur
print(f"    mu_w [kappa ~ 1/mu(meV)] = F/eps = 1/ell = {mu_w_meV:.3e} eV")
print(f"    mu_w [kappa ~ 1/M_S]     = {mu_w_TeV:.3e} eV")

# ============ D2 : gel du transfert winding<->B ============
# (a) operateur dim-6 (courant.courant)/M^2 : Gamma ~ T^5/M^4 ; gel Gamma=H
Tf_a = (MS_GeV**4/ MPl_GeV)**(1/3)
print(f"\nD2a gel par Gamma=H (op. 1/M_S^2)      : T_f = (M_S^4/M_Pl)^(1/3) = {Tf_a*1e3:.0f} MeV")
# (b) VERROU NUCLEONIQUE : le transfert DeltaB=1 doit passer par un etat plus lourd
#     que le proton (sinon p -> X_leger + l en ~ns). Le gel est alors de Boltzmann :
#     T_f ~ M_X / ln(...) ~ M_X/25  (gel standard d'un processus a seuil)
Tf_b = MX_GeV/25
print(f"D2b gel de Boltzmann sous le seuil M_X : T_f ~ M_X/25 = {Tf_b*1e3:.0f} MeV")

# verification du danger nucleon si PAS de seuil :
Gam_p = m_p**5/ (MS_GeV**4)                 # GeV
tau_s = 6.58e-25/Gam_p
print(f"    [sans seuil] duree de vie du proton via op 1/M_S^2 : {tau_s:.1e} s  << 1e34 ans  -> INTERDIT")

# ============ D1+D4 : fermeture ============
# eta_B ~ b_w * chi * mu_w / n_gamma |_(T_f)  ;  chi ~ T^2/6, n_gamma = 0.244 T^3
#   =>  b_w * mu_w ~ eta_B * (0.244*6) * T_f  ~  1.5 * eta_B * T_f
for name,Tf in [("D2a (Gamma=H)",Tf_a), ("D2b (seuil X)",Tf_b)]:
    need = 1.5*eta_B*Tf*1e9    # eV
    print(f"\nD4  [{name}]  b_w*mu_w requis = 1.5*eta_B*T_f = {need:.2e} eV")
    for lab,mw in [("kappa ~ 1/mu (meV)",mu_w_meV), ("kappa ~ 1/M_S (TeV)",mu_w_TeV)]:
        b_req = need/mw
        verdict = "ORDRE UN  -> fenetre NON VIDE" if 1e-3<b_req<1e3 else f"deficit 10^{np.log10(b_req):.0f} -> MORT"
        print(f"      {lab:22s}: b_w requis = {b_req:.2e}   {verdict}")

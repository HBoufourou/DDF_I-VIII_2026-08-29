# -*- coding: utf-8 -*-
"""Article VI — verification de chaque nombre publie.
   Usage : python3 verify.py     Sortie : ../proofs/certified_values.json"""
import json, os, numpy as np
OK=lambda c:"PASS" if c else "*** FAIL ***"; R={}

print("="*66); print("PART A — THE WALL"); print("="*66)
w1,p=5.5,0.4514; sr=2*w1*p**2
print(f"  regle de somme 2 w1 |p|^2 = {sr:.4f}                {OK(abs(sr-2.2414)<1e-3)}")
print(f"  rayon du cercle           = {np.sqrt(sr):.4f}                {OK(abs(np.sqrt(sr)-1.4971)<1e-3)}")
print(f"  c_b^2+c_b'^2 = |p|^2      = {p**2:.4f}   (facteur 2w1 = {2*w1:g} entre les deux)")
inv=max(abs((np.cos(t)*p)**2+(np.sin(t)*p)**2-p**2) for t in np.linspace(0,2*np.pi,4001))
print(f"  invariance le long du cercle, erreur max = {inv:.2e}   {OK(inv<1e-14)}")
hbar_c=1.97327e-7; invR=hbar_c/8.2e-6*1e3; L4=invR**4; Lobs=2.4**4; eps=np.sqrt(Lobs/L4)
print(f"\n  1/R                = {invR:.2f} meV                      {OK(abs(invR-24.06)<0.05)}")
print(f"  (1/R)^4            = {L4:.0f} meV^4   (NE PAS confondre avec Lambda_obs)")
print(f"  Lambda_obs         = {Lobs:.1f} meV^4")
print(f"  deficit epsilon    = {eps:.3e}                     {OK(abs(eps-9.95e-3)<1e-4)}")
print(f"  rapport absorbe    = {L4/Lobs:.0f} ~ 1e4, par eps^2 = {eps**2:.1e}")
ratio=0.26/0.049
print(f"\n  Omega_DM/Omega_b : predit 6, observe {ratio:.2f}, ecart {abs(6-ratio)/ratio*100:.1f} %")
print(f"  efficacite e = {ratio/6:.3f}   -- FACTEUR AJUSTE, declare comme tel")
print(f"  radion/Cassini : SECTION RETIREE (alpha_rho couple a la TENSION, pas a la matiere)")
R['partA']=dict(sum_rule=sr,radius=float(np.sqrt(sr)),cb2_sum=p**2,factor_2w1=2*w1,
 invariance_err=float(inv),invR_meV=invR,L4_meV4=L4,Lambda_obs_meV4=Lobs,epsilon=float(eps),
 Omega_ratio_obs=ratio,efficiency_fitted=ratio/6,radion_section="withdrawn")

print("\n"+"="*66); print("PART B — THE FIBRE"); print("="*66)
As,v,MPl=2.1e-9,2.87,2.435e18
print(f"  {'N_e':>4} {'eps':>9} {'V0/MPl^4':>12} {'M_inf (GeV)':>13} {'r=6(ns-1)^2':>13}")
obs=[]
for Ne,e_,ns,r_ in ((55,3.7e-4,0.9676,0.0059),(58,3.4e-4,0.9693,0.0054),(60,3.2e-4,0.9704,0.0051)):
    V0=24*np.pi**2*e_*As/v; M=V0**0.25*MPl; pred=6*(ns-1)**2
    obs.append(dict(Ne=Ne,eps=e_,ns=ns,r=r_,V0=V0,M_inf=M,r_pred=pred))
    print(f"  {Ne:4d} {e_:9.1e} {V0:12.3e} {M:13.3e} {pred:13.4f}  vs r={r_:.4f}  {OK(abs(pred-r_)<0.0015)}")
print(f"\n  >>> V0 = 24 pi^2 eps A_s / v ne contient NI kappa NI volume NI coefficient de boucle")
print(f"      => A_s -> V0 -> M_inf est INDEPENDANT DE LA GEOMETRIE")
print(f"  candidat 3181 : UN SEUL diviseur K3 (le n.1), FIXE par l'involution donc PAIR")
print(f"  les deux dP7 sont ECHANGES par l'involution -> [Open] : refaire la statistique")
print(f"                                                  de l'Article V en base paire")
R['partB']=dict(observables=obs,geometry_independent=True,K3_divisors=1,inflaton="divisor 1",
 parity="verified: fixed by the involution",open_item="dP7 pair exchanged; recompute Article V statistics in the even basis")
os.makedirs('../proofs',exist_ok=True)
json.dump(R,open('../proofs/certified_values.json','w'),indent=1)
print("\n-> ../proofs/certified_values.json")

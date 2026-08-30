"""c_b — LA DERIVATION FERMEE. Le terme de flux F0^2 est AUSSI de rang un ;
deux rang-un dans un espace 3D laissent un probleme 2x2 a UN SEUL parametre. 07/08/2026."""
import numpy as np, sympy as sp
sq=sp.sqrt
a1=1/(4*sq(7)); b1=-sq(7)/4; a2=sq(sp.Rational(5,28)); b2=-1/sq(35)
print("="*74); print("1. LE VECTEUR DU TERME DE FLUX F0^2 (meme comptage que la tension)")
print("="*74)
print("  IIA massive, frame d'Einstein 10D : terme (1/2) e^{5 phi/2} F0^2 , F0 = 0-forme")
print("  => poids dilaton 5/2 (contre 5/4 pour la tension D8 : (p-3)/4)")
print("  terme de VOLUME (il remplit les 10 dims) => il porte AUSSI e^{beta1 rho} du segment")
wf_phi=sp.Rational(5,2); wf_chi=4*a2+5*b2; wf_rho=9*a1+b1
w_phi=sp.Rational(5,4); w_chi=4*a2+5*b2; w_rho=9*a1
v=[sp.simplify(x*sq(2)) for x in (wf_phi,wf_chi,wf_rho)]
w=[sp.simplify(x*sq(2)) for x in (w_phi,w_chi,w_rho)]
print(f"  w (tension) = ({float(w[0]):.4f}, {float(w[1]):.4f}, {float(w[2]):.4f}) , |w|^2 = {sp.simplify(sum(x**2 for x in w))}")
print(f"  v (flux)    = ({float(v[0]):.4f}, {float(v[1]):.4f}, {float(v[2]):.4f}) , |v|^2 = {sp.simplify(sum(x**2 for x in v))}")
dot=sp.simplify(sum(a*b for a,b in zip(w,v)))
print(f"\n  >>> DEUX THEOREMES DE NORME : |w|^2 = 6 (deja connu) et |v|^2 = {sp.simplify(sum(x**2 for x in v))} (NOUVEAU)")
print(f"  >>> et le produit scalaire est ENTIER : w.v = {dot}")
cos=sp.simplify(dot/(sq(6)*sq(14)))
print(f"  >>> cos(angle) = {dot}/sqrt(84) = {sp.nsimplify(cos)} = {float(cos):.6f}  ->  angle = {float(sp.deg(sp.acos(cos))):.2f} deg")
print("\n"+"="*74); print("2. LA STRUCTURE : DEUX RANG-UN => UN PROBLEME 2x2 A UN PARAMETRE")
print("="*74)
print("  M^2 = kappa w w^T + delta v v^T   (les deux sources, les deux de rang un)")
print("  => la direction orthogonale a w ET a v est de masse NULLE : ce n'est pas Phi")
print("     (Phi doit tirer sa masse du puits). Phi vit donc dans le plan span{w, v},")
print("     et c'est l'etat le PLUS LEGER des deux massifs.")
print("  => c_b ne depend plus que de UN nombre : r = delta/kappa.")
W=np.array([float(x) for x in w]); V=np.array([float(x) for x in v])
print(f"\n  {'r = delta/kappa':>15} {'c_b':>9} {'fenetre [0,30;0,61]'}")
sols=[]
for r in (0.01,0.1,0.25,0.5,1.0,2.0,4.0,8.0,16.0,50.0,200.0):
    M=np.outer(W,W)+r*np.outer(V,V)
    ev,U=np.linalg.eigh(M); n=U[:,1]              # le plus leger des deux massifs
    cb=abs(n@W); sols.append((r,cb))
    print(f"  {r:15.2f} {cb:9.4f} {'OUI' if 0.30<=cb<=0.61 else '-'}")
print("\n  -> recherche des r qui donnent exactement les bords de la fenetre :")
from scipy.optimize import brentq
def cb_of(r):
    M=np.outer(W,W)+r*np.outer(V,V); ev,U=np.linalg.eigh(M); return abs(U[:,1]@W)
for target in (0.30,0.61):
    try:
        r=brentq(lambda x: cb_of(np.exp(x))-target,-8,8); print(f"     c_b = {target:.2f}  <=>  r = {np.exp(r):.4f}")
    except Exception as e: print(f"     c_b = {target}: pas de solution dans la plage")

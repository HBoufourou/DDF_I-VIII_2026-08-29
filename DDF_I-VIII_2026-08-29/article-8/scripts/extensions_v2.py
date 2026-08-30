"""EXTENSIONS v2 (05/08) : (1) k=4 exponentielles (limitation levee) ;
(2) premier vide explicite racetrack a 2 modules (Etapes 3-4, machinerie complete)."""
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import minimize_scalar, minimize
what=np.array([5*np.sqrt(2)/4,np.sqrt(10/7)]); what/=np.linalg.norm(what)
LO,HI=73.4,81.9
print("=== (1) k = 4 EXPONENTIELLES — methode : racine sur l'angle du ker 2D ===")
rng=np.random.default_rng(20260805)
S=np.array([1,-1,2,-2,np.sqrt(2),-np.sqrt(2),1/np.sqrt(2),-1/np.sqrt(2),np.sqrt(7)/2,-np.sqrt(7)/2,2/np.sqrt(7),-2/np.sqrt(7)])
th=[]
for _ in range(2500):
    A=rng.uniform(0.5,2.0,4)*rng.choice([1,-1],4); E=rng.choice(S,(4,2))
    ns=null_space(E.T)
    if ns.shape[1]!=2: continue
    M=np.column_stack([E,-np.ones(4)])
    def res(a):
        u=np.cos(a)*ns[:,0]+np.sin(a)*ns[:,1]; r=u/A
        if np.any(r<=0): return 1e3
        sol,rr,rk,_=np.linalg.lstsq(M,np.log(r),rcond=None)
        return rr[0] if len(rr) else 1e3
    best=(1e3,None)
    for a0 in np.linspace(0,2*np.pi,4,endpoint=False):
        m=minimize_scalar(res,bracket=(a0,a0+0.3)) if False else minimize(lambda v:res(v[0]),[a0],method='Nelder-Mead',options={'xatol':1e-8,'fatol':1e-12,'maxiter':120})
        if m.fun<best[0]: best=(m.fun,m.x[0])
    if best[0]>1e-12: continue
    a=best[1]; u=np.cos(a)*ns[:,0]+np.sin(a)*ns[:,1]
    H=(E.T*u)@E; ev,evec=np.linalg.eigh(H)
    if ev[0]<=0 or ev[0]/ev[1]>0.1: continue
    th.append(np.degrees(np.arccos(min(1,abs(evec[:,0]@what)))))
th=np.array(th); iso=(HI-LO)/90
if len(th)>100:
    f=np.mean((th>=LO)&(th<=HI))
    print(f"  minima valides k=4 : {len(th)} ; fraction fenetre = {f:.1%} (isotrope {iso:.1%}, ratio {f/iso:.2f})")
    print(f"  vs k=3 mesure cordes : 16-19 %  =>  verdict S2 {'INCHANGE' if 0.5>f>0.03 else 'A REVOIR'} a k=4")
else: print(f"  minima valides : {len(th)} — statistique insuffisante, limitation MAINTENUE")
print()
print("=== (2) PREMIER VIDE EXPLICITE : racetrack a 2 modules (K = -ln(S+Sb) - 3ln(T+Tb)) ===")
A_,a_,B_,b_,W0=1.0,2*np.pi/10,1.2,2*np.pi/12,-1e-3
def W(s,t): return W0+A_*np.exp(-a_*t)+B_*np.exp(-b_*s)
def Ws(s,t): return -b_*B_*np.exp(-b_*s)
def Wt(s,t): return -a_*A_*np.exp(-a_*t)
def V(x):
    s,t=x
    if s<=0.05 or t<=0.05: return 1e3
    K=-np.log(2*s)-3*np.log(2*t)
    DS=Ws(s,t)+(-1/(2*s))*W(s,t); DT=Wt(s,t)+(-3/(2*t))*W(s,t)
    return np.exp(K)*((2*s)**2*DS**2/1+ (2*t)**2/3*DT**2 - 3*W(s,t)**2)
best=None
for s0 in (5,10,20):
    for t0 in (5,10,20):
        m=minimize(V,[s0,t0],method='Nelder-Mead',options={'xatol':1e-12,'fatol':1e-18,'maxiter':4000})
        if best is None or m.fun<best.fun: best=m
s0,t0=best.x
print(f"  minimum : Re S = {s0:.3f}, Re T = {t0:.3f} ; V_min = {best.fun:.3e} (<0 : AdS, uplift requis — standard, declare)")
# Hessienne en champs CANONIQUES : sigma = ln(2s)/sqrt2 ; tau = sqrt(3/2) ln(2t)
h=1e-5
def Vc(y):
    sg,ta=y; return V([np.exp(np.sqrt(2)*sg)/2, np.exp(ta*np.sqrt(2/3))/2])
y0=np.array([np.log(2*s0)/np.sqrt(2), np.sqrt(3/2)*np.log(2*t0)])
H=np.zeros((2,2))
for i in range(2):
    for j in range(2):
        e1=np.zeros(2);e1[i]=h; e2=np.zeros(2);e2[j]=h
        H[i,j]=(Vc(y0+e1+e2)-Vc(y0+e1-e2)-Vc(y0-e1+e2)+Vc(y0-e1-e2))/(4*h*h)
ev,evec=np.linalg.eigh(H)
print(f"  Hessienne canonique (sigma, tau) : valeurs propres = {ev[0]:.3e}, {ev[1]:.3e} ; hierarchie {ev[0]/ev[1]:.2e}")
nl=evec[:,0]
print(f"  ETAT LEGER : n = {nl[0]:+.3f} sigma + {nl[1]:+.3f} tau   (sigma ~ dilaton 4D, tau ~ volume)")
import math
for lab,d in (("pure sigma (dilaton)",[1,0]),("pure tau (volume)",[0,1]),("diagonale (1,1)/sqrt2",[1/np.sqrt(2),1/np.sqrt(2)])):
    print(f"  angle vs direction de tension candidate [{lab:24s}] : {math.degrees(math.acos(min(1,abs(np.dot(nl,d))))):.1f} deg")
print("  [le theta DEFINITIF attend le DICTIONNAIRE (S,T)<->(phi,chi) — la marche critique")
print("   declaree de l'Etape 3 ; l'etat propre est livre pour que tout dictionnaire s'applique]")

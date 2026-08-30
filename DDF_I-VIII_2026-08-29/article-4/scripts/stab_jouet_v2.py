"""E3a-3 v2 — ROBUSTESSE du verdict S2 : le rapport fenetre/isotrope est-il stable
quand on change (i) la mesure de tirage des exposants, (ii) le nombre de termes k,
(iii) la coupure de hierarchie ? Methode exacte (ker E^T + systeme lineaire). Seed fixe."""
import numpy as np
from scipy.linalg import null_space
what=np.array([5*np.sqrt(2)/4,np.sqrt(10/7)]); what/=np.linalg.norm(what)
LO,HI=73.4,81.9; ISO=(HI-LO)/90
def run(measure,k,cut,N,seed):
    rng=np.random.default_rng(seed); th=[]
    stringset=np.array([1,-1,2,-2,np.sqrt(2),-np.sqrt(2),1/np.sqrt(2),-1/np.sqrt(2),np.sqrt(7)/2,-np.sqrt(7)/2,2/np.sqrt(7),-2/np.sqrt(7)])
    for _ in range(N):
        A=rng.uniform(0.5,2.0,k)*rng.choice([1,-1],k)
        if measure=="unif": E=rng.uniform(-2,2,(k,2))
        elif measure=="gauss": E=rng.normal(0,1,(k,2))
        else: E=rng.choice(stringset,(k,2))          # exposants "de cordes" (rationnels sur sqrt2, sqrt7)
        ns=null_space(E.T)
        if ns.shape[1]!=k-2: continue
        # k=3 : ker dim1 -> 2 signes ; k=4 : ker dim2 -> echantillonner 6 directions
        dirs=[ns[:,0],-ns[:,0]] if k==3 else [np.cos(a)*ns[:,0]+np.sin(a)*ns[:,1] for a in np.linspace(0,2*np.pi,6,endpoint=False)]
        for u_dir in dirs:
            r=u_dir/A
            if np.any(r<=0): continue
            M=np.column_stack([E,-np.ones(k)])
            sol,res,rk,_=np.linalg.lstsq(M,np.log(r),rcond=None)
            if rk<3 or (len(res) and res[0]>1e-16): continue    # k=4 : systeme surdetermine, exiger residu nul
            H=(E.T*u_dir)@E
            ev,evec=np.linalg.eigh(H)
            if ev[0]<=0 or ev[0]/ev[1]>cut: continue
            th.append(np.degrees(np.arccos(min(1,abs(evec[:,0]@what))))); break
    th=np.array(th); n=len(th)
    if n<200: return n,None,None
    f=np.mean((th>=LO)&(th<=HI)); return n,f,np.sqrt(f*(1-f)/n)
print(f"{'mesure':8s} {'k':>2s} {'coupure':>8s} {'minima':>7s} {'fenetre':>8s} {'ratio/iso':>9s}")
for meas in ("unif","gauss","cordes"):
    for k,cut,N in ((3,0.1,150000),(3,0.02,150000),(4,0.1,60000)):
        n,f,se=run(meas,k,cut,N,20260805)
        print(f"{meas:8s} {k:2d} {cut:8.2f} {n:7d} " + (f"{f:7.1%} {f/ISO:9.2f}" if f is not None else "  (stat. insuffisante)"))
print(f"\nattente isotrope : {ISO:.1%} ; verdict S2 ROBUSTE si tous les ratios restent ~0.8-1.5")

"""AUDIT DE MINIMUM v2 — corrige les 3 defauts de la v1 sur les potentiels de supergravite.
Diagnostic du cas E1 : min_eig ~ 1e-10 alors que V ~ 1e-7 => on lit du BRUIT, pas la courbure.
06/08/2026."""
import numpy as np

def audit_minimum_v2(V, x, echelles=None, verbose=True):
    """x : point teste. echelles : ordre de grandeur naturel de chaque variable
    (par defaut |x_i|). Travaille en variables LOG => pas relatifs, tous les moduli
    sur le meme pied. Test SANS DIMENSION : gradient logarithmique et valeurs propres
    comparees a l'echelle naturelle |V|, jamais a zero absolu."""
    x=np.asarray(x,float); n=len(x)
    sc=np.abs(x) if echelles is None else np.asarray(echelles,float)
    sc=np.where(sc>0,sc,1.0)
    V0=V(x); scaleV=abs(V0) if V0!=0 else 1.0
    # pas optimaux : eps^(1/3) pour le gradient, eps^(1/4) pour la Hessienne (relatifs)
    eps=np.finfo(float).eps
    hg=eps**(1/3)*sc; hh=eps**(1/4)*sc
    g=np.zeros(n)
    for i in range(n):
        e=np.zeros(n); e[i]=hg[i]
        g[i]=(V(x+e)-V(x-e))/(2*hg[i])
    glog=g*sc/scaleV                      # gradient LOGARITHMIQUE, sans dimension
    H=np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            e1=np.zeros(n); e1[i]=hh[i]; e2=np.zeros(n); e2[j]=hh[j]
            H[i,j]=(V(x+e1+e2)-V(x+e1-e2)-V(x-e1+e2)+V(x-e1-e2))/(4*hh[i]*hh[j])
    Hlog=H*np.outer(sc,sc)/scaleV         # Hessienne sans dimension
    Hlog=(Hlog+Hlog.T)/2
    ev=np.linalg.eigvalsh(Hlog)
    # plancher de bruit estime : erreur relative sur V amplifiee par la difference seconde
    bruit=eps**0.5*max(1.0,np.max(np.abs(Hlog)))
    d=dict(grad_log=np.linalg.norm(glog), eigs_log=ev, bruit=bruit,
           stationnaire=np.linalg.norm(glog)<1e-4,
           minimum=bool(np.linalg.norm(glog)<1e-4 and ev.min()>bruit),
           direction_plate=bool(abs(ev.min())<bruit),
           hierarchie=float(ev.min()/ev.max()) if ev.max()>0 else np.nan)
    if verbose:
        print(f"  |grad log| = {d['grad_log']:.2e}   (stationnaire si < 1e-4)")
        print(f"  valeurs propres log = {np.array2string(ev,precision=3)}")
        print(f"  plancher de bruit   = {bruit:.2e}")
        print(f"  -> stationnaire : {d['stationnaire']} | minimum net : {d['minimum']}"
              f" | direction plate (physique ou bruit) : {d['direction_plate']}")
    return d

if __name__=="__main__":
    print("="*70); print("DEMO : pourquoi la v1 echouait — un potentiel de type LVS"); print("="*70)
    # jouet reproduisant la pathologie : V ~ 1e-7, variables d'echelles 1 et 1e3, direction plate
    def Vtoy(z):
        s,tb=z
        return 1e-7*((s-4.0)**2 + 1e-6*(np.log(tb/1000.))**2)
    x=[4.0,1000.]
    print("\n--- v1 (pas fixe h=1e-6, seuils absolus) ---")
    h=1e-6
    g=np.array([(Vtoy(np.array(x)+h*np.eye(2)[i])-Vtoy(np.array(x)-h*np.eye(2)[i]))/(2*h) for i in range(2)])
    H=np.zeros((2,2))
    for i in range(2):
        for j in range(2):
            e1=h*np.eye(2)[i]; e2=h*np.eye(2)[j]
            H[i,j]=(Vtoy(np.array(x)+e1+e2)-Vtoy(np.array(x)+e1-e2)-Vtoy(np.array(x)-e1+e2)+Vtoy(np.array(x)-e1-e2))/(4*h*h)
    ev=np.linalg.eigvalsh(H)
    print(f"  grad_norm = {np.linalg.norm(g):.3e} ; min_eig = {ev.min():.3e}  -> 'est_minimum = False'")
    print("  (exactement le motif de son audit E1 : eigenvalue au niveau du bruit)")
    print("\n--- v2 (pas relatifs, tests sans dimension) ---")
    audit_minimum_v2(Vtoy,x,echelles=[1.0,1000.])
    print("\n  => le meme point : STATIONNAIRE, avec une direction quasi plate IDENTIFIEE")
    print("     comme telle au lieu d'etre lue comme 'pas un minimum'.")

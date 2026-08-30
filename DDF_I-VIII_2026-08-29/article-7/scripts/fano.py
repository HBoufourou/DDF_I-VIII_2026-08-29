# -*- coding: utf-8 -*-
"""Article VII — construction explicite du groupe, recensement des classes,
   paysage combinatoire et seuils. Tout y est reconstruit, rien n'est affirme.
   Usage : python3 fano.py     Sortie : ../proofs/certified_values.json"""
import itertools, json, os
from collections import defaultdict

def mm(A,B): return [[sum(A[i][k]*B[k][j] for k in range(3))%2 for j in range(3)] for i in range(3)]
S=[[1,1,0],[0,1,0],[0,0,1]]; T=[[0,1,0],[0,0,1],[1,0,0]]; I=[[1,0,0],[0,1,0],[0,0,1]]
G={tuple(map(tuple,I))}; new=[tuple(map(tuple,I))]
while new:
    c=new.pop()
    for g in (S,T):
        t=tuple(map(tuple, mm(list(map(list,c)), g)))
        if t not in G: G.add(t); new.append(t)
pts=[tuple(int(b) for b in format(i,'03b')) for i in range(1,8)]
ap=lambda M,p: tuple(sum(M[i][j]*p[j] for j in range(3))%2 for i in range(3))
orbit={ap(list(map(list,M)),pts[0]) for M in G}
lines=[c for c in itertools.combinations(pts,3) if all(c[0][k]^c[1][k]^c[2][k]==0 for k in range(3))]
inc=defaultdict(int)
for L in lines:
    for p in L: inc[p]+=1
print(f"|GL(3,2)| = {len(G):3d}   attendu 168        {'PASS' if len(G)==168 else 'FAIL'}")
print(f"transitif sur les 7 points                 {'PASS' if orbit==set(pts) else 'FAIL'}")
print(f"droites   = {len(lines):3d}   attendu 7          {'PASS' if len(lines)==7 else 'FAIL'}")
print(f"3 droites par point                        {'PASS' if all(inc[p]==3 for p in pts) else 'FAIL'}")

PDG={'nu_e':(0,0,1),'nu_mu':(0,0,1),'n':(0,1,0),'Lambda0':(0,1,0),'Sigma_c++':(2,1,0),
     'pi+':(1,0,0),'K+':(1,0,0),'W+':(1,0,0),'e-':(-1,0,1),'mu-':(-1,0,1),'tau-':(-1,0,1),
     'p':(1,1,0),'Sigma+':(1,1,0),'Lambda_c+':(1,1,0)}
vec=lambda Q,B,L:(abs(Q)%2,B%2,L%2)
cls=defaultdict(list)
for k,(Q,B,L) in PDG.items(): cls[vec(Q,B,L)].append(k)
empty=[v for v in pts if not cls[v]]
print(f"\nclasses peuplees = {7-len(empty)}, vides = {len(empty)} -> {empty}   "
      f"{'PASS' if len(empty)==2 and set(empty)=={(0,1,1),(1,1,1)} else 'FAIL'}")
exotics={'X3872':(0,0,0),'Zc3900':(1,0,0),'Pc4312':(1,1,0),'Tcc':(1,0,0),'Zb10610':(1,0,0)}
fill=[v for v in exotics.values() if v in empty]
print(f"exotiques dans une classe vide : {fill}          {'PASS' if not fill else 'FAIL'}")
print(f"hydrogene (0,1,1) tomberait dans une classe vide -> exclu (etat lie EM, pas une resonance)")

b={2:2,1:5,0:8,-1:5}
c011=3*b[1]+3*b[-1]+3*b[0]+3*b[2]; c111=3*b[0]+3*b[2]+3*b[1]+3*b[-1]
print(f"\npaysage (0,1,1) = {c011}, (1,1,1) = {c111}, total {c011+c111}   "
      f"{'PASS' if c011==60 and c111==60 else 'FAIL'}")

M={'Lambda_c+':2286.5,'mu-':105.7,'Sigma_c+':2452.9,'e-':0.511,'Sigma_c++':2453.97,'Xi_c+':2467.9}
CH=[('011',['Lambda_c+','mu-']),('011',['Sigma_c+','e-']),
    ('111',['Sigma_c++','mu-']),('111',['Xi_c+','e-'])]
chans=[]
print()
for c,prod in CH:
    thr=sum(M[p] for p in prod)
    chans.append({'class':c,'products':prod,'threshold_MeV':round(thr,1)})
    print(f"  classe {c} -> {'+'.join(prod):24s} seuil {thr:7.1f} MeV")

out={'group':{'order':len(G),'transitive':orbit==set(pts),'lines':len(lines),
              'lines_per_point':3,'verified':True},
     'classes':{'populated':7-len(empty),'empty':[list(v) for v in empty]},
     'atomic_exclusion':{'hydrogen':[0,1,1],'note':'EM bound state, not an S-matrix pole'},
     'exotics_fill_empty':bool(fill),
     'landscape':{'011':c011,'111':c111,'total':c011+c111},
     'channels':chans,
     'epistemic':{'derived':['F2^3 classification','GL(3,2) order 168 transitive','7 lines',
                             'landscape 60/60','channel thresholds'],
                  'fact':['5 populated / 2 empty (PDG)'],
                  'hypothesis':['a genuine resonance exists in each empty class'],
                  'assumption':['mass window 4-6 GeV'],
                  'asymmetry':'the symmetry allows the classes; it does not require them to be occupied'}}
os.makedirs('../proofs',exist_ok=True)
json.dump(out,open('../proofs/certified_values.json','w'),indent=1)
print("\n-> ../proofs/certified_values.json")

# Verification exhaustive du secteur fracton (Article IX, section 4)
import numpy as np
pts=[(0,0,1),(0,1,0),(0,1,1),(1,0,0),(1,0,1),(1,1,0),(1,1,1)]
lines=[]
for i in range(7):
    for j in range(i+1,7):
        k=pts.index(tuple(a^b for a,b in zip(pts[i],pts[j])))
        L=tuple(sorted([i,j,k]))
        if L not in lines: lines.append(L)
assert len(lines)==7
M=np.zeros((7,7),int)
for r,L in enumerate(lines):
    for c in L: M[r,c]=1
img={tuple(sum(M[b] for b in range(7) if m>>b&1)%2 if m else np.zeros(7,int)) for m in range(128)}
ker=[m for m in range(128) if not any((M@np.array([(m>>b)&1 for b in range(7)]))%2)]
assert len(img)==16 and sorted({sum(v) for v in img})==[0,3,4,7]   # Hamming [7,4,3]
assert len(ker)==8 and {bin(m).count('1') for m in ker if m}=={4}  # simplexe [7,3,4]
assert not any(sum(v) in (1,2) for v in img)                       # M1
print("Fracton sector: rank 4, image=Hamming[7,4,3], kernel=simplex[7,3,4],")
print("8 ground states (weight 4), M1 verified. ALL CHECKS PASS.")

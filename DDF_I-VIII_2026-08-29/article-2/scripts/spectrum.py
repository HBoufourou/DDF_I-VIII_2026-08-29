
# parity correction (v9): a single field carries ONE branch
from scipy.special import ai_zeros as _z
_a,_b,_,_=_z(3)
print("even-branch ratios :"," : ".join(f"{-x/-_b[0]:.2f}" for x in _b))
print("odd-branch ratios  :"," : ".join(f"{-x/-_a[0]:.2f}" for x in _a))
print("interleaved (two-field realisation) : 1 : 2.29 : 3.19 : 4.01")

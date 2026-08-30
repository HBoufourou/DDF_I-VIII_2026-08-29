
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os, csv, json

HERE = os.path.dirname(os.path.abspath(__file__))
ART  = os.path.dirname(HERE)
FIG, DAT, PRF = (os.path.join(ART, s) for s in ("figures","data","proofs"))
for d in (FIG, DAT, PRF): os.makedirs(d, exist_ok=True)

eta_B=6.1e-10; mu_meV=34e-3; mu_TeV=1.4e-16; tau_obs=1e34*3.156e7
F_eV2 = 24e-3*(0.197327/5.81*1e3)*1e-3
def b_w(Tf,mu): return 1.5*eta_B*Tf/mu
bw_thr, bw_GH, bw_TeV = b_w(5e9/25,mu_meV), b_w(47e6,mu_meV), b_w(5e9/25,mu_TeV)
print(f"F = {F_eV2:.2e} eV2 (attendu 8.15e-4) ; b_w = {bw_thr:.1f} / {bw_GH:.1f} / {bw_TeV:.1e}")
assert abs(bw_thr-5.4)<0.1 and abs(bw_GH-1.3)<0.1 and 1e14<bw_TeV<1e16

MX=np.linspace(2,8,200); bw=1.5*eta_B*(MX*1e9/25)/mu_meV
fig,ax=plt.subplots(figsize=(9,5.5))
ax.semilogy(MX,bw,'b-',lw=2,label='meV branch: $b_w=1.5\\,\\eta_B T_f/\\mu_w$, $T_f=M_X/25$')
ax.axhspan(1,10,color='g',alpha=.15); ax.text(2.1,7,'order-one window',color='g')
ax.axhline(bw_TeV,color='r',ls='--',lw=2,label=f'TeV branch: $b_w\\simeq{bw_TeV:.1e}$ (15 orders short)')
ax.axvspan(4,6,color='grey',alpha=.2); ax.text(4.1,3e2,'benchmark 4-6 GeV',rotation=90,va='bottom')
ax.plot([5],[bw_thr],'ko',ms=8); ax.annotate(f'threshold: $b_w={bw_thr:.1f}$',(5,bw_thr),textcoords='offset points',xytext=(10,5))
ax.plot([5],[bw_GH],'ks',ms=8); ax.annotate(f'$\\Gamma=H$: $b_w={bw_GH:.1f}$',(5,bw_GH),textcoords='offset points',xytext=(10,-14))
ax.set_yscale('log'); ax.set_ylim(0.5,1e16)
ax.set_xlabel('$M_X$ [GeV]'); ax.set_ylabel('$b_w$ required')
ax.set_title('The closure: meV branch closes at order one; TeV branch fails by fifteen orders')
ax.legend(loc='lower right'); ax.grid(alpha=.3,which='both')
fig.tight_layout(); fig.savefig(os.path.join(FIG,'fig1_closure.png'),dpi=300); plt.close(fig)

MS=np.logspace(2,17,200); tau=1.2e-9*(MS/1e3)**4; MS_x=1e3*(tau_obs/1.2e-9)**0.25
fig,ax=plt.subplots(figsize=(9,5.5))
ax.loglog(MS,tau,'b-',lw=2,label='$\\tau_p\\sim1.2\\times10^{-9}\\,{\\rm s}\\,(M_S/{\\rm TeV})^4$')
ax.axhline(tau_obs,color='g',ls='--',lw=2,label='observed $\\tau_p>10^{34}$ yr')
ax.plot([1e3],[1.2e-9],'ro',ms=9)
ax.annotate('wall scale: $10^{-9}$ s\n(51 orders short)',(1e3,1.2e-9),textcoords='offset points',xytext=(12,8))
ax.plot([MS_x],[tau_obs],'g^',ms=10)
ax.annotate(f'would need $M_S\\simeq{MS_x:.0e}$ GeV\n$\\Rightarrow$ lock $m_X>m_p$ mandatory',(MS_x,tau_obs),textcoords='offset points',xytext=(-230,30))
ax.set_xlabel('$M_S$ [GeV]'); ax.set_ylabel('$\\tau_p$ [s]')
ax.set_title('The nucleon lock: a TeV-scale contact transfer is dead on arrival')
ax.legend(loc='lower right'); ax.grid(alpha=.3,which='both')
fig.tight_layout(); fig.savefig(os.path.join(FIG,'fig2_lock.png'),dpi=300); plt.close(fig)

fig,ax=plt.subplots(figsize=(9,4.8),subplot_kw={'projection':'mollweide'})
lon,lat=np.radians(264-180),np.radians(48)
ax.scatter([lon],[lat],marker='*',s=260,c='r',edgecolors='k',zorder=5)
ax.annotate('residual = CMB dipole\n$(l,b)=(264^\\circ,48^\\circ)$, $A=\\beta=1.23\\times10^{-3}$',(lon,lat),textcoords='offset points',xytext=(-40,-40))
ax.set_xticklabels([]); ax.grid(alpha=.3)
ax.set_title('Prediction I: direction fixed a priori (no sky scan, no trials factor)')
fig.tight_layout(); fig.savefig(os.path.join(FIG,'fig3_dipole.png'),dpi=300); plt.close(fig)

w=csv.writer(open(os.path.join(DAT,'d8_closure.csv'),'w',newline=''))
w.writerows([['route','Tf_MeV','bw_mu_required_eV','mu_w_eV','bw_required'],
 ['threshold_at_MX',200,1.8e-1,34e-3,5.4],['Gamma_eq_H',47,4.3e-2,34e-3,1.3],
 ['TeV_branch',200,1.8e-1,1.4e-16,1.3e15]])
w=csv.writer(open(os.path.join(DAT,'d8_lock.csv'),'w',newline=''))
w.writerows([['MS_GeV','tau_p_s']]+[[f'{m:.1e}',f'{1.2e-9*(m/1e3)**4:.1e}'] for m in [1e3,1e6,1e9,1e12,MS_x,1e16]])
w=csv.writer(open(os.path.join(DAT,'d8_dipole.csv'),'w',newline=''))
w.writerows([['l_deg','b_deg','A'],[264,48,1.23e-3]])
json.dump({'F_eV2':8.15e-4,'mu_w_meV_eV':34e-3,'mu_w_TeV_eV':1.4e-16,'eta_B':eta_B,
 'bw_threshold':5.4,'bw_GammaH':1.3,'bw_TeV':1.3e15,'tau_p_TeV_s':1.2e-9,
 'MS_thresholdless_GeV':float(f'{MS_x:.2e}'),'dipole_l':264,'dipole_b':48,'A':1.23e-3},
 open(os.path.join(PRF,'d8_certified.json'),'w'),indent=2)
print("OK : figures + data + proofs régénérés dans", ART)

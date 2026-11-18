"""Block-average analysis of Kremer-Grest SLLOD runs (LJ units).
sigma = -P ; N1 = sxx-syy ; N2 = syy-szz ; eta = sxy/rate ; Psi = N/rate^2."""
import numpy as np, glob, json, re
res=[]
for f in sorted(glob.glob('stress_r*.dat')):
    tag=f[7:-4]; rate={'r002':0.02,'r005':0.05,'r01':0.1,'r02':0.2}[tag]
    d=np.loadtxt(f,comments='#'); 
    n=len(d); d=d[n//3:]                     # discard first third (transient)
    pxx,pyy,pzz,pxy=d[:,1],d[:,2],d[:,3],d[:,4]
    q={'eta':-pxy/rate,'N1':pyy-pxx,'N2':pzz-pyy}
    out={'rate':rate,'n_samples_used':len(d)}
    B=5; m=len(d)//B
    for k,v in q.items():
        bl=np.array([v[i*m:(i+1)*m].mean() for i in range(B)])
        out[k]=bl.mean(); out[k+'_se']=bl.std(ddof=1)/np.sqrt(B)
    out['Psi1']=out['N1']/rate**2; out['Psi2']=out['N2']/rate**2
    out['Psi1_se']=out['N1_se']/rate**2; out['Psi2_se']=out['N2_se']/rate**2
    out['N2_over_N1']=out['N2']/out['N1']
    res.append(out)
json.dump(res,open('kg_results.json','w'),indent=1)
print("rate  eta(se)  N1(se)  N2(se)  Psi2(se)  N2/N1")
for r in res: print(f"{r['rate']:.3f} {r['eta']:.3f}({r['eta_se']:.3f}) {r['N1']:.3f}({r['N1_se']:.3f}) {r['N2']:.3f}({r['N2_se']:.3f}) {r['Psi2']:.2f}({r['Psi2_se']:.2f}) {r['N2_over_N1']:.3f}")

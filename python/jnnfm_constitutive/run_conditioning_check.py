"""Conditioning and wall-curvature checks on the approach to the regularity boundary.

Usage:  python run_conditioning_check.py -0.3 -0.8 -1.2 -1.5 -1.6 -1.64 -1.66
        (first argument is K, the rest are the shrinking rates lambda; S = 1.5)

Three diagnostics are reported along the negative branch, and a table is written
to data/conditioning_check.csv.

1. Bordered Newton-Jacobian: smallest singular value and condition number.
   The wall row of the momentum equation is replaced by the boundary condition
   F(0) = lambda, so the singular coefficient 1 - 2KF never enters this matrix
   at the wall. The test therefore shows only that no progressive deterioration
   of the bordered Jacobian is observed; it does not probe the singular
   coefficient directly. For K = -0.3, S = 1.5 the condition number stays between
   about 1.9e7 and 1.7e7 while 1 - 2K*lambda falls from 0.52 to 0.004.

2. Coefficient-level diagnostic. In the closed momentum system the matrix that
   multiplies (F'', G'') is (1 - 2KF) times the 2x2 identity, so its smallest
   singular value is min_eta |1 - 2KF|. This is the object that becomes singular;
   it is reported together with the location where the minimum is attained (the
   wall, eta = 0, on the negative branch).

3. Wall curvature F''(0) at N = 160, 200 and 240. The curvature is resolved to
   better than 0.2% for 1 - 2K*lambda >= 0.04 (F''(0) = -40.4 at lambda = -1.6);
   closer to the boundary it keeps steepening but is no longer converged in N
   (values between about -74 and -77 at lambda = -1.66, 1 - 2K*lambda = 0.004).
"""
import numpy as np, sys
from cheb_base import base_state_from_bvp
def jac(K,lam,S,N=100,einf=30.0):
    r=base_state_from_bvp(K,lam,S,N=N,einf=einf)
    eta,F,G,H,De,D2=r; n=N+1; I=np.eye(n)
    Fp,Gp=De@F,De@G; Fpp,Gpp=D2@F,D2@G; Dl=1-2*K*F
    J11=np.diag(Dl)@D2+np.diag(-2*K*Fpp-2*F)-np.diag(H)@De+2*K*np.diag(Fp)@De
    J12=np.diag(2*G)-2*K*np.diag(Gp)@De; J13=-np.diag(Fp)
    J21=np.diag(-2*K*Gpp-2*G)+2*K*np.diag(Gp)@De
    J22=np.diag(Dl)@D2-np.diag(2*F)-np.diag(H)@De+2*K*np.diag(Fp)@De; J23=-np.diag(Gp)
    J=np.block([[J11,J12,J13],[J21,J22,J23],[2*I,np.zeros((n,n)),De]])
    w,f=0,N
    for row in (w,n+w,2*n+w,f,n+f): J[row,:]=0; J[row,row]=1.0
    s=np.linalg.svd(J,compute_uv=False)
    return s[-1], s[0]/s[-1], 1-2*K*lam

def wall_diag(K,lam,S,N):
    eta,F,G,H,De,D2=base_state_from_bvp(K,lam,S,N=N,einf=30.0)
    Dl=1-2*K*F
    i=int(np.argmin(np.abs(Dl)))
    return np.abs(Dl[i]), eta[i], (D2@F)[0]

if __name__=="__main__":
    import csv, os
    K=float(sys.argv[1]) if len(sys.argv)>1 else -0.3; S=1.5; Ns=(160,200,240)
    rows=[]
    print("%7s %9s %11s %11s %9s %11s   F''(0) at N=160/200/240"%("lam","1-2K*lam","sigma_min(J)","kappa(J)","eta_min","min|1-2KF|"))
    for lam in [float(x) for x in sys.argv[2:]]:
        smin,kap,fac=jac(K,lam,S)
        d=[wall_diag(K,lam,S,N) for N in Ns]
        c=[x[2] for x in d]
        spread=(max(c)-min(c))/abs(c[1])
        print("%7.4f %9.4f %11.3e %11.3e %9.3f %11.4f   %.3f / %.3f / %.3f  (spread %.2f%%)"%(lam,fac,smin,kap,d[1][1],d[1][0],*c,100*spread))
        rows.append([K,S,lam,fac,smin,kap,d[1][0],d[1][1],*c,spread])
    out=os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..","data","conditioning_check.csv")
    with open(out,"w",newline="") as fh:
        w=csv.writer(fh)
        w.writerow(["K","S","lambda","one_minus_2K_lambda","sigma_min_bordered_jacobian","cond_bordered_jacobian",
                    "min_abs_1_minus_2KF","eta_at_min","Fpp0_N160","Fpp0_N200","Fpp0_N240","rel_spread_Fpp0"])
        w.writerows(rows)
    print("wrote",os.path.normpath(out))

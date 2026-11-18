import numpy as np
from rr_solver import solve, continue_lambda

def chebD(N):
    x=np.cos(np.pi*np.arange(N+1)/N)
    c=np.hstack([2.,np.ones(N-1),2.])*(-1)**np.arange(N+1)
    X=np.tile(x,(N+1,1)).T; dX=X-X.T
    D=np.outer(c,1./c)/(dX+np.eye(N+1)); D-=np.diag(D.sum(axis=1))
    return D,x

def gamma1_momentum(K,lam,S,N=110,einf=25.0,sol=None):
    """Smallest momentum eigenvalue via Chebyshev collocation, H eliminated algebraically."""
    if sol is None:
        # continue from lambda=0 to target
        s=None; l=0.0; step=0.05 if lam>0 else -0.05
        while abs(l-lam)>1e-9:
            nxt = lam if abs(lam-l)<abs(step) else l+step
            t=solve(s,K,nxt,S,einf=einf)
            if t.status==0: s,l=t,nxt
            else:
                step/=2
                if abs(step)<1e-7: return None
        sol=s
    D,xi=chebD(N); eta=einf*(1-xi)/2
    De=D*(-2/einf); D2=De@De
    y=sol.sol(eta); F0,F0p,G0,G0p,H0=y
    Dl=1-2*K*F0
    F0pp=(F0**2-G0**2+H0*F0p-K*(F0p**2-G0p**2))/Dl
    G0pp=(2*F0*G0+H0*G0p-2*K*F0p*G0p)/Dl
    # H-operator: H' = -2F, H(0)=0  ->  H = P F
    D2h=De.copy(); D2h[0,:]=0; D2h[0,0]=1
    Z=np.eye(N+1); Z[0,0]=0
    P=-2.0*np.linalg.solve(D2h,Z)
    I=np.eye(N+1); O=np.zeros((N+1,N+1))
    A_FF=(np.diag(Dl)@D2+2*K*np.diag(F0pp)-2*K*np.diag(F0p)@De-np.diag(H0)@De+2*np.diag(F0)-np.diag(F0p)@P)
    A_FG=2*K*np.diag(G0p)@De-2*np.diag(G0)
    A_GF=-4*K*np.diag(G0pp)+2*K*np.diag(G0p)@De-2*np.diag(G0)-np.diag(G0p)@P
    A_GG=(np.diag(Dl)@D2+2*K*np.diag(F0pp)+2*K*np.diag(F0p)@De-np.diag(H0)@De-2*np.diag(F0))
    A=np.block([[A_FF,A_FG],[A_GF,A_GG]]); Bm=np.block([[I,O],[O,I]])
    Af=-A.copy(); Bf=Bm.copy(); m=N+1
    for r,cix in [(0,0),(m,m),(m+N,m+N)]:
        Af[r,:]=0; Af[r,cix]=1; Bf[r,:]=0
    Af[N,:]=0; Af[N,0:m]=De[N,:]; Bf[N,:]=0
    gs=np.linspace(0.0,1.2,900)
    sv=np.array([np.linalg.svd(Af-g*Bf,compute_uv=False)[-1] for g in gs])
    # local minima
    idx=[i for i in range(1,len(gs)-1) if sv[i]<sv[i-1] and sv[i]<sv[i+1]]
    if not idx: return None,sol
    return gs[idx[0]],sol

"""Illustrative application: rotating-disk contactor, no-radial-pressure-gradient closure, lambda=0, S=0. Parameters are assumed, not measured."""
import numpy as np
from scipy.integrate import solve_bvp
def solve(K,Z=40.0,guess=None):
    def rhs(e,y):
        F,Fp,G,Gp,H=y; D=1-2*K*F
        return np.vstack([Fp,(F**2-G**2+H*Fp-K*(Fp**2-Gp**2))/D,Gp,(2*F*G+H*Gp-2*K*Fp*Gp)/D,-2*F])
    bc=lambda a,b: np.array([a[0],a[2]-1,a[4],b[0],b[2]])
    if guess is None:
        x=np.linspace(0,Z,800);d=np.exp(-x);y=np.vstack([0.2*x*d,0.2*(1-x)*d,d,-d,-0.88*(1-d)])
    else:x,y=guess.x,guess.y
    return solve_bvp(rhs,bc,x,y,tol=1e-10,max_nodes=200000)
s0=solve(0.0); print("K=0:",s0.status,"F'(0)=%.6f G'(0)=%.6f H(inf)=%.6f"%(s0.y[1,0],s0.y[3,0],s0.y[4,-1]))
res={}
for K in (-0.0314,-0.1,-0.3142):
    s=s0
    for k in np.linspace(0,K,12)[1:]: s=solve(k,guess=s); assert s.status==0,(k,s.status)
    res[K]=s
    print("K=%.4f: status %d F'(0)=%.6f G'(0)=%.6f H(inf)=%.6f min(1-2KF)=%.4f max F=%.4f"%(K,s.status,s.y[1,0],s.y[3,0],s.y[4,-1],(1-2*K*s.y[0]).min(),s.y[0].max()))
# dimensional example
rho,eta0,nu=1000.0,0.01,1e-5; R=0.02; n=300; Om=2*np.pi*n/60; Psi2=-1e-5
K=Psi2*Om/eta0
print("Omega=%.4f rad/s K=Psi2*Om/eta0=%.5f r*sqrt(Om/nu) at rim=%.3f Re_R=%.1f"%(Om,K,R*np.sqrt(Om/nu),Om*R**2/nu))
s=res[-0.0314]
for name,sol,Kv in (("Newtonian",s0,0.0),("RR K=%.4f"%K,solve(K,guess=s),K)):
    Gp=sol.y[3,0]; Hinf=sol.y[4,-1]
    M=np.pi/2*eta0*Om*np.sqrt(Om/nu)*abs(Gp)*R**4
    Fp=sol.y[1,0]; gam=R*Om*np.sqrt(Om/nu)*np.hypot(Fp,Gp)  # resultant wall shear rate
    tau=eta0*gam
    Wi2=abs(Kv)*R*np.sqrt(Om/nu)*np.hypot(Fp,Gp)
    Q=np.pi*R**2*np.sqrt(nu*Om)*abs(Hinf)
    print(f"{name}: -G'(0)={abs(Gp):.5f} torque/face M={M:.6e} N m  rim shear rate={gam:.1f} 1/s  rim wall shear stress={tau:.3f} Pa  rim Wi2={Wi2:.3f}  entrainment Q={Q:.6e} m3/s  N2(rim)={Psi2*gam**2 if Kv else 0:.3f} Pa")

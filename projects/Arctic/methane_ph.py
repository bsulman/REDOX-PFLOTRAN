from pylab import *

phs=arange(3.0,8.5,0.25) 

def pHdep(phs,monod_k,inhib_k,norm=False):
    H_concs=10**-phs
    out=H_concs/(H_concs+10**(-monod_k))*10**(-inhib_k)/(H_concs+10**-inhib_k)
    # Function looks symmetrical around two k values based on some optimizations
    if norm:
        out=out/pHdep(0.5*(monod_k+inhib_k),monod_k,inhib_k,norm=False)
    return out

# Data from Kotsyurbenko et al 2007 at 15 C (they also have other temperature points)
kotsy_ph=array([6.0,4.8,3.8])
kotsy_acet=array([22,17,1.7])
kotsy_acet_err=array([7,6,0.9])
kotsy_total=array([33,26,14])
kotsy_total_err=array([9,7,5])
kotsy_acetfrac=array([0.68,0.65,0.12])
kotsy_acetfrac_err=array([0.26,0.29,0.07])

def lsqerr(x):
    acet=pHdep(kotsy_ph,x[0],x[1],norm=True)
    hydro=pHdep(kotsy_ph,x[2],x[3],norm=True)*x[4]
    frac=acet/(acet+hydro)

    # return ((frac-kotsy_acetfrac)**2).sum()
    # return ((acet+hydro-kotsy_total/kotsy_total.max())**2).sum()
    return (concatenate([acet+hydro-kotsy_total/kotsy_total.max(),frac-kotsy_acetfrac])**2).sum()

from scipy.optimize import fmin
x=fmin(lsqerr,[5,5,5.0,5.0,0.5])
print(x)

f,axs=subplots(nrows=2,num='pH effects',clear=True)
acetaclastic=pHdep(phs,x[0],x[1],norm=True)
hydrogenotrophic=pHdep(phs,x[2],x[3],norm=True)*x[4]
axs[0].plot(phs,acetaclastic,label='Acetaclastic')
axs[0].plot(phs,hydrogenotrophic,label='Hydrogenotrophic')
axs[0].plot(phs,hydrogenotrophic+acetaclastic,label='Total')
axs[0].errorbar(kotsy_ph,kotsy_acet/kotsy_total.max(),yerr=kotsy_acet_err/kotsy_total.max(),c='C0',ls='None',marker='o')
axs[0].errorbar(kotsy_ph,kotsy_total/kotsy_total.max(),yerr=kotsy_total_err/kotsy_total.max(),c='C2',ls='None',marker='o')
axs[0].errorbar(kotsy_ph,(kotsy_total-kotsy_acet)/kotsy_total.max(),yerr=(kotsy_total_err**2+kotsy_acet_err**2)**0.5/kotsy_total.max(),c='C3',ls='None',marker='o')
axs[0].legend()

axs[1].plot(phs,acetaclastic/(acetaclastic+hydrogenotrophic))
axs[1].errorbar(kotsy_ph,kotsy_acetfrac,yerr=kotsy_acetfrac_err,c='C0',ls='None',marker='o')
axs[0].set(ylabel='Relative methanogenesis rate',title='Relative methanogenesis rate',xlabel='pH')
axs[1].set(ylabel='Acetaclastic fraction',title='Acetaclastic fraction of total methanogenesis',xlabel='pH')

axs[0].set_title('(a)',loc='left')
axs[1].set_title('(b)',loc='left')


show()
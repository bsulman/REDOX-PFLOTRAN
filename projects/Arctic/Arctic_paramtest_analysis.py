import plot_Arctic_Fe as pa
import glob
import xarray
from numpy import *
from matplotlib import pyplot

sims=['organic_trough','organic_nottrough','mineral_trough','mineral_nottrough']
stats=[ 'RMSE'  ,'nRMSE' ,'MAE' ,'R' ,'R2']
vars=['pH','Fe_II','CO2flux','CH4flux','log(DOM)','log(Acetate)']

def get_stats(filenames):
    groups=pa.get_groups(filenames[0])

    ac_scales=unique([float(g.split('_')[3]) for g in groups])
    frates=unique([float(g.split('_')[-1]) for g in groups])
    allstats=xarray.Dataset(coords={'acetate_scale':ac_scales,'ferm_rate':frates,'sim':sims,'stat':stats})
    for var in vars:
        allstats[var]=xarray.DataArray(nan,allstats.coords)

    for group in groups:
        print(group)
        mod_obs=pa.mod_obs_comp(filenames,group[group.find('_acetate'):])
        ac_scale=float(group.split('_')[3])
        frate=float(group.split('_')[-1])
        for sim in mod_obs:
            statsdata=pa.getstats(mod_obs[sim])
            for stat in statsdata:
                for var in statsdata[stat].index:
                    allstats[var].loc[{'acetate_scale':ac_scale,'ferm_rate':frate,'sim':sim,'stat':stat}]=statsdata[stat][var]

def plot_stats(allstats,stat,figlabel=None):
    f,a=pyplot.subplots(nrows=len(sims),ncols=len(vars),num=figlabel,clear=True,figsize=(15,10),constrained_layout=False)
    x=allstats['ferm_rate']
    y=allstats['acetate_scale']
    for varnum in range(len(vars)):
        for simnum in range(len(sims)):
            sim=sims[simnum]
            var=vars[varnum]
            h=a[simnum,varnum].pcolormesh(x,y,allstats[var].sel(sim=sim,stat=stat),shading='auto')
            cb=pyplot.colorbar(h,ax=a[simnum,varnum])
            cb.set_label(stat)
            a[simnum,varnum].set(yscale='log',xlabel='Fermentation rate',ylabel='Acetate scale',title=var,ylim=(y[0]*0.9,y[-1]*1.1))
            if var=='pH':
                a[simnum,varnum].set_title('%s: pH'%sim)
            
    return f

if __name__=='__main__':
    allstats=xarray.open_dataset('Arctic_Fe_output/paramtest_stats_2021-05-06.nc')
    for stat in stats:
        print(stat)
        f=plot_stats(allstats,stat)
        f.tight_layout()
        f.savefig('Arctic_Fe_output/paramtest_%s.png'%stat)
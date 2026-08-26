from plot_ELM_alquimia_result import plot_vars
import xarray
import matplotlib.pyplot as plt
import numpy as np

yr=1999
data=xarray.open_dataset('/lustre/or-scratch/cades-ccsi/b0u/Alaska_alquimia_7cell_AK-BEO_ICB20TRCNRDCTCBC/run/Alaska_alquimia_7cell_AK-BEO_ICB20TRCNRDCTCBC.elm.h0.%d-01-01-00000.nc'%yr)

landcover_types=[
  'Trough',
  'LCP center',
  'LCP transition',
  'HCP center',
  'HCP transition',
  'Rim',
  'Average',
]

c=plt.get_cmap('RdBu_r')(np.linspace(0,1,6))

vars=['VWC','O2','Fe2','FeOxide','CH4','Sulfide','acetate','temperature','pH']

f,a=plt.subplots(ncols=len(vars),num='Profiles',clear=True,figsize=(9.2,4.4),sharey=True)
for num in range(6):
    plot_vars(data.isel(lndgrid=num).sel(time=slice(f'{yr:d}-06-01',f'{yr:d}-10-01')),vars=vars,figsize=(6,8.5),maxdepth=2.0,a_profile=a,
              vmin={'FeOxide':450,'Fe2':-0.01,'CH4':-0.01,'pH':2.0},vmax={'Fe2':1.75,'CH4':14.0,'pH':12.0},profile_color=c[num],quantiles=[])
    
for num in range(1,len(vars)):
    a[num].set_ylabel('')
for num in range(len(vars)):
    a[num].xaxis.label.set_fontsize('large')

a[0].legend(labels=landcover_types[:6],fontsize='large')
a[0].yaxis.label.set_fontsize('large')

for n in range(len(landcover_types)):
    plot_vars(data.isel(lndgrid=n),plotname=landcover_types[n],
              vars=['VWC','O2','frozen','temperature','acetate','CH4'],vmax={'Fe2':1.0,'CH4':10.0},
              figsize=(6,8.5),maxdepth=1.0,snapshots=[24*365//2],mean_profile=False)

f,a=plt.subplots(num='Carbon time series',clear=True,nrows=3)
# (data['TOTSOMC']/1000).plot(ax=a[0],label='Total SOM C')
for n in range(len(landcover_types)):
    (data.isel(lndgrid=n)['TOTVEGC']/1000).plot(ax=a[0],label=landcover_types[n])
    a[0].set(title='C pools',xlabel='Time (year)',ylabel='C stock (kg C m$^{-2}$')
    a[0].legend()
    a[0].set_ylim(bottom=0.0)

    # (data['HR']*3600*24).plot(ax=a[1],label='HR')
    # Should add a smoothed curve
    (data.isel(lndgrid=n)['HR']*3600*24).resample(time='7D').mean().plot(ax=a[1],ls='-',c='C%d'%n)
    (data.isel(lndgrid=n)['GPP']*3600*24).resample(time='7D').mean().plot(ax=a[1],ls='--',c='C%d'%n)
    (data.isel(lndgrid=n)['NEE']*3600*24).resample(time='7D').mean().plot(ax=a[1],ls=':',c='C%d'%n)
    a[1].set(title='C fluxes',xlabel='Time (year)',ylabel='C flux (g C m$^{-2}$ day$^{-1}$)')
    # a[1].legend()
    a[1].axhline(0.0,c='k',ls=':',lw=0.5)
    a[1].legend(labels=['HR','GPP','NEE'])

    (data.isel(lndgrid=n)['CH4FLUX_ALQUIMIA']*3600*24).resample(time='7D').mean().plot(ax=a[2],label=landcover_types[n])
    a[2].set(title='CH$_4$ flux',xlabel='Time (year)',ylabel='C flux (g C m$^{-2}$ day$^{-1}$)')
    a[2].axhline(0.0,c='k',ls=':',lw=0.5)

surfdata_multicell=xarray.open_dataset('/nfs/data/ccsi/proj-shared/b0u/NGEE_ELM/BEO_surfdata_multicell.nc')
tide_data_multicell=xarray.open_dataset('/nfs/data/ccsi/proj-shared/b0u/NGEE_ELM/BEO_hydro_BC_multicell.nc')
f,a=plt.subplots(num='Water heights',clear=True,nrows=1)
a.fill_between(np.arange(len(landcover_types)),np.zeros(len(landcover_types))-.1,surfdata_multicell['ht_above_stream'],ls='-',color='brown',label='Soil surface',step='mid')
a.axhspan(-0.1,tide_data_multicell['tide_height'].mean(),color='b',alpha=0.5,label='Water level')
plt.xticks(ticks=np.arange(len(landcover_types)),labels=landcover_types)
a.set_ylabel('Height (m)')
a.set(xlim=(0,len(landcover_types)-1.5),ylim=(-0.1,0.23),title='Polygon landform levels')

plt.show()
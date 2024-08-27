import plot_ELM_alquimia_result as pltELM
import xarray
import matplotlib.pyplot as plt
import sys,os
import numpy as np
import pandas

yr=1999

# SETx=xarray.open_dataset(f'/home/b0u/ELM_outputs/marsh_alquimia_saltmarshrest_US-SETX_ICB20TRCNRDCTCBC/run/marsh_alquimia_saltmarshrest_US-SETX_ICB20TRCNRDCTCBC.elm.h0.{yr:d}-01-01-00000.nc')
# SETx_lowFe=xarray.open_dataset(f'/home/b0u/ELM_outputs/marsh_alquimia_saltmarshrest_lowFe_US-SETX_ICB20TRCNRDCTCBC/run/marsh_alquimia_saltmarshrest_lowFe_US-SETX_ICB20TRCNRDCTCBC.elm.h0.{yr:d}-01-01-00000.nc')

# SETx=xarray.open_dataset(f'/lustre/or-scratch/cades-ccsi/b0u/marsh_alquimia_saltmarshrest_lowerdrain_US-SETX_ICB20TRCNRDCTCBC/run/marsh_alquimia_saltmarshrest_lowerdrain_US-SETX_ICB20TRCNRDCTCBC.elm.h0.{yr:d}-01-01-00000.nc')

SETx=xarray.open_dataset(f'/lustre/or-scratch/cades-ccsi/b0u/marsh_alquimia_saltmarshrest_2023_US-SETX_ICB20TRCNRDCTCBC/run/marsh_alquimia_saltmarshrest_2023_US-SETX_ICB20TRCNRDCTCBC.elm.h0.{yr:d}-01-01-00000.nc')
SETx_lowFe=xarray.open_dataset(f'/lustre/or-scratch/cades-ccsi/b0u/marsh_alquimia_saltmarshrest_lowFe_2023_US-SETX_ICB20TRCNRDCTCBC/run/marsh_alquimia_saltmarshrest_lowFe_2023_US-SETX_ICB20TRCNRDCTCBC.elm.h0.{yr:d}-01-01-00000.nc')
SETx_floodresponse=xarray.open_dataset(f'/lustre/or-scratch/cades-ccsi/b0u/marsh_alquimia_saltmarshrest_2023_floodresponse_US-SETX_ICB20TRCNRDCTCBC/run/marsh_alquimia_saltmarshrest_2023_floodresponse_US-SETX_ICB20TRCNRDCTCBC.elm.h0.{yr:d}-01-01-00000.nc')

landcover=['Fresh marsh',
              'Flooded forest/swamp',
              'Woody estuarine',
              'Salt marsh',
              'Degraded salt marsh',
              'Restored salt marsh',
              'Freshwater shrub',
              'Farmed wetland',
            ]

marsh_healthy=SETx.isel(lndgrid=landcover.index('Salt marsh'))
marsh_degraded=SETx.isel(lndgrid=landcover.index('Degraded salt marsh'))
marsh_degraded_lowFe=SETx_lowFe.isel(lndgrid=landcover.index('Degraded salt marsh'))
marsh_restored=SETx.isel(lndgrid=landcover.index('Restored salt marsh'))

Hoch_obs=pandas.read_excel('/home/b0u/SETx/DOE_Marsh_Master_Dataset_all.xlsx',sheet_name='JDMWMA sites',skiprows=2,skipfooter=1)

# vars=['CH4','sulfide','Fe2','H2OSFC','CH4flux']
vars=['temperature','salinity','pH','Fe2','Sulfide','Sulfate','CH4']
start=f'{yr:d}-10-01'
end=f'{yr:d}-11-01'
vmax={'CH4':2,'sulfide':25,'Fe2':0.25,'sulfate':40.0,'salinity':32,'pH':7.0}
vmin={'pH':5.0,'salinity':-1}

f,a=plt.subplots(nrows=1,ncols=len(vars),num='Marsh restoration comp',clear=True,figsize=(16,5))
pltELM.plot_vars(marsh_healthy.sel(time=slice(start,end)),vars,a_profile=a,maxdepth=0.3,profile_color='green',vmax=vmax,vmin=vmin,mean_profile=True,quantiles=[0.2,0.8])
pltELM.plot_vars(marsh_degraded.sel(time=slice(start,end)),vars,a_profile=a,maxdepth=0.3,profile_color='purple',vmax=vmax,vmin=vmin,mean_profile=True,quantiles=[0.2,0.8])
pltELM.plot_vars(marsh_restored.sel(time=slice(start,end)),vars,a_profile=a,maxdepth=0.3,profile_color='orange',vmax=vmax,vmin=vmin,mean_profile=True,quantiles=[0.2,0.8])
pltELM.plot_vars(marsh_degraded_lowFe.sel(time=slice(start,end)),vars,a_profile=a,maxdepth=0.3,profile_color='magenta',vmax=vmax,vmin=vmin,mean_profile=True,quantiles=[0.2,0.8])

Hoch_means=Hoch_obs.groupby(Hoch_obs['site-station'].str[0]).mean()
Hoch_errs=Hoch_obs.groupby(Hoch_obs['site-station'].str[0]).sem()

a[vars.index('temperature')].errorbar(Hoch_means['temperature']['H'],0.15,xerr=Hoch_errs['temperature']['H'],c='green',marker='o')
a[vars.index('temperature')].errorbar(Hoch_means['temperature']['D'],0.17,xerr=Hoch_errs['temperature']['D'],c='purple',marker='o')
a[vars.index('temperature')].errorbar(Hoch_means['temperature']['R'],0.19,xerr=Hoch_errs['temperature']['R'],c='orange',marker='o')

a[vars.index('salinity')].errorbar(Hoch_means['salinity']['H'],0.15,xerr=Hoch_errs['salinity']['H'],c='green',marker='o')
a[vars.index('salinity')].errorbar(Hoch_means['salinity']['D'],0.17,xerr=Hoch_errs['salinity']['D'],c='purple',marker='o')
a[vars.index('salinity')].errorbar(Hoch_means['salinity']['R'],0.19,xerr=Hoch_errs['salinity']['R'],c='orange',marker='o')

a[vars.index('pH')].errorbar(Hoch_means['pH']['H'],0.15,xerr=Hoch_errs['pH']['H'],c='green',marker='o')
a[vars.index('pH')].errorbar(Hoch_means['pH']['D'],0.17,xerr=Hoch_errs['pH']['D'],c='purple',marker='o')
a[vars.index('pH')].errorbar(Hoch_means['pH']['R'],0.19,xerr=Hoch_errs['pH']['R'],c='orange',marker='o')

a[vars.index('Fe2')].errorbar(Hoch_means['ferrous']['H']*1e-3,0.15,xerr=Hoch_errs['ferrous']['H']*1e-3,c='green',marker='o')
a[vars.index('Fe2')].errorbar(Hoch_means['ferrous']['D']*1e-3,0.17,xerr=Hoch_errs['ferrous']['D']*1e-3,c='purple',marker='o')
a[vars.index('Fe2')].errorbar(Hoch_means['ferrous']['R']*1e-3,0.19,xerr=Hoch_errs['ferrous']['R']*1e-3,c='orange',marker='o')

a[vars.index('Sulfide')].errorbar(Hoch_means['sulfide']['H']*1e-3,0.15,xerr=Hoch_errs['sulfide']['H']*1e-3,c='green',marker='o')
a[vars.index('Sulfide')].errorbar(Hoch_means['sulfide']['D']*1e-3,0.17,xerr=Hoch_errs['sulfide']['D']*1e-3,c='purple',marker='o')
a[vars.index('Sulfide')].errorbar(Hoch_means['sulfide']['R']*1e-3,0.19,xerr=Hoch_errs['sulfide']['R']*1e-3,c='orange',marker='o')

a[vars.index('Sulfate')].errorbar(Hoch_means['sulfate']['H'],0.15,xerr=Hoch_errs['sulfate']['H'],c='green',marker='o')
a[vars.index('Sulfate')].errorbar(Hoch_means['sulfate']['D'],0.17,xerr=Hoch_errs['sulfate']['D'],c='purple',marker='o')
a[vars.index('Sulfate')].errorbar(Hoch_means['sulfate']['R'],0.19,xerr=Hoch_errs['sulfate']['R'],c='orange',marker='o')

# How do I calculate porosity or total water volume from Matt's data? He has bulk density and water weight % but not volume %
# For now, assuming H2O weight % approximates volume % but should check with Matt
# CH4 inventory is in umol/m2. Need to convert to mmol/L H2O
a[vars.index('CH4')].cla()
pltELM.plot_var(marsh_healthy['CH4_vr'].T.squeeze()[:10,:]/12.011,profileax=a[vars.index('CH4')],vmax=vmax.get('CH4',5.0),vmin=vmin.get('CH4',0.0),label='CH$_4$ concentration\n(mmol C/L H$_2$O)',
                     profile_color='green',mean_profile=True,quantiles=[0.2,0.8],
                     maxdepth=0.3,title='Soil CH$_4$',axlabel='CH$_4$ (mmol C/m$^3$)')
pltELM.plot_var(marsh_degraded['CH4_vr'].T.squeeze()[:10,:]/12.011,profileax=a[vars.index('CH4')],vmax=vmax.get('CH4',5.0),vmin=vmin.get('CH4',0.0),label='CH$_4$ concentration\n(mmol C/L H$_2$O)',
                     profile_color='purple',mean_profile=True,quantiles=[0.2,0.8],
                     maxdepth=0.3,title='Soil CH$_4$',axlabel='CH$_4$ (mmol C/m$^3$)')
pltELM.plot_var(marsh_degraded_lowFe['CH4_vr'].T.squeeze()[:10,:]/12.011,profileax=a[vars.index('CH4')],vmax=vmax.get('CH4',5.0),vmin=vmin.get('CH4',0.0),label='CH$_4$ concentration\n(mmol C/L H$_2$O)',
                     profile_color='magenta',mean_profile=True,quantiles=[0.2,0.8],
                     maxdepth=0.3,title='Soil CH$_4$',axlabel='CH$_4$ (mmol C/m$^3$)')
pltELM.plot_var(marsh_restored['CH4_vr'].T.squeeze()[:10,:]/12.011,profileax=a[vars.index('CH4')],vmax=vmax.get('CH4',5.0),vmin=vmin.get('CH4',0.0),label='CH$_4$ concentration\n(mmol C/L H$_2$O)',
                     profile_color='orange',mean_profile=True,quantiles=[0.2,0.8],
                     maxdepth=0.3,title='Soil CH$_4$',axlabel='CH$_4$ (mmol C/m$^3$)')
a[vars.index('CH4')].errorbar(Hoch_means['CH4 Inventory']['H']/0.25*1e-3,0.15,xerr=Hoch_errs['CH4 Inventory']['H']/0.25*1e-3,c='green',marker='o')
a[vars.index('CH4')].errorbar(Hoch_means['CH4 Inventory']['D']/0.25*1e-3,0.17,xerr=Hoch_errs['CH4 Inventory']['D']/0.25*1e-3,c='purple',marker='o')
a[vars.index('CH4')].errorbar(Hoch_means['CH4 Inventory']['R']/0.25*1e-3,0.19,xerr=Hoch_errs['CH4 Inventory']['R']/0.25*1e-3,c='orange',marker='o')

a[1].legend(['Healthy marsh','Degraded marsh','Restored marsh'],loc='lower right')

f,a=plt.subplots(nrows=2,num='Timeseries',clear=True)
pltELM.plot_vars(marsh_healthy.sel(time=slice(start,end)).resample(time='1D').mean(),['H2OSFC','CH4flux'],a_contour=a,profile_color='green',vmax=vmax,vmin=vmin,mean_profile=True,quantiles=[0.1,0.9],WT_thresh=0.85)
pltELM.plot_vars(marsh_degraded.sel(time=slice(start,end)).resample(time='1D').mean(),['H2OSFC','CH4flux'],a_contour=a,profile_color='purple',vmax=vmax,vmin=vmin,mean_profile=True,quantiles=[0.1,0.9],WT_thresh=0.85)
pltELM.plot_vars(marsh_restored.sel(time=slice(start,end)).resample(time='1D').mean(),['H2OSFC','CH4flux'],a_contour=a,profile_color='orange',vmax=vmax,vmin=vmin,mean_profile=True,quantiles=[0.1,0.9],WT_thresh=0.85)
pltELM.plot_vars(marsh_degraded_lowFe.sel(time=slice(start,end)).resample(time='1D').mean(),['H2OSFC','CH4flux'],a_contour=a,profile_color='magenta',vmax=vmax,vmin=vmin,mean_profile=True,quantiles=[0.1,0.9],WT_thresh=0.85)

a[0].legend(labels=['Healthy marsh','Degraded marsh','Restored marsh'],handles=a[0].lines[::2],loc='lower right')
a[1].errorbar(marsh_restored['time'].sel(time='%d-10-15 12:00'%yr).item(),Hoch_means['CH4 flux']['H']*12e-6,yerr=Hoch_errs['CH4 flux']['H']*12e-6,color='green',marker='o')
a[1].errorbar(marsh_restored['time'].sel(time='%d-10-18 12:00'%yr).item(),Hoch_means['CH4 flux']['D']*12e-6,yerr=Hoch_errs['CH4 flux']['D']*12e-6,color='purple',marker='o')
a[1].errorbar(marsh_restored['time'].sel(time='%d-10-20 12:00'%yr).item(),Hoch_means['CH4 flux']['R']*12e-6,yerr=Hoch_errs['CH4 flux']['R']*12e-6,color='orange',marker='o')



f,a=plt.subplots(nrows=2,num='Forcing',clear=True)
pltELM.plot_vars(marsh_healthy,['H2OSFC_TIDE'],a_contour=a,profile_color='green',vmax=vmax,vmin=vmin,mean_profile=True,quantiles=[0.1,0.9],WT_thresh=0.85)
pltELM.plot_vars(marsh_degraded,['H2OSFC_TIDE'],a_contour=a,profile_color='purple',vmax=vmax,vmin=vmin,mean_profile=True,quantiles=[0.1,0.9],WT_thresh=0.85)
pltELM.plot_vars(marsh_restored,['H2OSFC_TIDE'],a_contour=a,profile_color='orange',vmax=vmax,vmin=vmin,mean_profile=True,quantiles=[0.1,0.9],WT_thresh=0.85)
a[1].plot(marsh_healthy['time'],marsh_healthy['SALINITY'],c='green')
a[1].plot(marsh_degraded['time'],marsh_healthy['SALINITY'],c='purple')
a[1].plot(marsh_restored['time'],marsh_healthy['SALINITY'],c='orange')
a[1].set(xlabel='Time',ylabel='Salinity (psu)',title='Salinity')
a[1].xaxis.set_major_formatter(pltELM.format_nc_time('%b-%d'))
# a[0].legend(labels=['Healthy marsh','Degraded marsh','Restored marsh'],handles=a[0].lines[::2],loc='lower right')
a[0].set_xlabel('Time')

vars=['H2OSFC','VWC','DOC','CH4','FeOxide','Sulfide','Sulfate','salinity']
# vars=['H2OSFC','VWC','O2','DIC','pH']
pltELM.plot_vars(marsh_healthy,plotname='Healthy marsh',
              vars=vars,vmax=vmax,
              figsize=(6,8.5),maxdepth=1.0,snapshots=(np.linspace(0.25,0.75,3)*24*365).astype(int),mean_profile=False)

pltELM.plot_vars(marsh_degraded,plotname='Degraded marsh',
              vars=vars,vmax=vmax,
              figsize=(6,8.5),maxdepth=1.0,snapshots=[24*365//2],mean_profile=False)

pltELM.plot_vars(marsh_degraded_lowFe,plotname='Degraded marsh (low Fe)',
              vars=vars,vmax=vmax,
              figsize=(6,8.5),maxdepth=1.0,snapshots=[24*365//2],mean_profile=False)

pltELM.plot_vars(marsh_restored,plotname='Restored marsh',
              vars=vars,vmax=vmax,
              figsize=(6,8.5),maxdepth=1.0,snapshots=[24*365//2],mean_profile=False)


# (data['TOTSOMC']/1000).plot(ax=a[0],label='Total SOM C')
def plot_fluxes(a,data,name,color):
    (data['TOTVEGC']/1000).plot(ax=a[0],label=name,c=color)
    # (data['TOTSOMC']/1000).plot(ax=a[0],c=color,ls='--')
    a[0].set(title='C pools',xlabel='Time (year)',ylabel='C stock (kg C m$^{-2}$')
    a[0].legend()
    a[0].set_ylim(bottom=0.0)

    # (data['HR']*3600*24).plot(ax=a[1],label='HR')
    # Should add a smoothed curve
    (data['HR']*3600*24).resample(time='7D').mean().plot(ax=a[1],ls='-',c=color)
    (data['GPP']*3600*24).resample(time='7D').mean().plot(ax=a[1],ls='--',c=color)
    (data['NEE']*3600*24).resample(time='7D').mean().plot(ax=a[1],ls=':',c=color)
    a[1].set(title='C fluxes',xlabel='Time (year)',ylabel='C flux (g C m$^{-2}$ day$^{-1}$)')
    # a[1].legend()
    a[1].axhline(0.0,c='k',ls=':',lw=0.5)
    a[1].legend(labels=['HR','GPP','NEE'])

    (data['CH4FLUX_ALQUIMIA']*3600*24).resample(time='7D').mean().plot(ax=a[2],label=name,c=color)
    a[2].set(title='CH$_4$ flux',xlabel='Time (year)',ylabel='C flux (g C m$^{-2}$ day$^{-1}$)')
    a[2].axhline(0.0,c='k',ls=':',lw=0.5)

f,a=plt.subplots(num='Carbon time series flood response',clear=True,nrows=3)
plot_fluxes(a,marsh_healthy,'Healthy marsh','green')
plot_fluxes(a,marsh_degraded,'Degraded marsh','purple')
plot_fluxes(a,marsh_degraded_lowFe,'Degraded marsh (low Fe)','magenta')
plot_fluxes(a,marsh_restored,'Restored marsh','orange')

a[2].errorbar(marsh_restored['time'].sel(time='%d-10-15 12:00'%yr).item(),Hoch_means['CH4 flux']['H']*12e-6,yerr=Hoch_errs['CH4 flux']['H']*12e-6,color='green',marker='o')
a[2].errorbar(marsh_restored['time'].sel(time='%d-10-18 12:00'%yr).item(),Hoch_means['CH4 flux']['D']*12e-6,yerr=Hoch_errs['CH4 flux']['D']*12e-6,color='purple',marker='o')
a[2].errorbar(marsh_restored['time'].sel(time='%d-10-20 12:00'%yr).item(),Hoch_means['CH4 flux']['R']*12e-6,yerr=Hoch_errs['CH4 flux']['R']*12e-6,color='orange',marker='o')

(SETx_floodresponse.isel(lndgrid=landcover.index('Salt marsh'))['GPP']*3600*24).resample(time='7D').mean().plot(ax=a[1],ls='--',c='green')
(SETx_floodresponse.isel(lndgrid=landcover.index('Degraded salt marsh'))['GPP']*3600*24).resample(time='7D').mean().plot(ax=a[1],ls='--',c='purple')
(SETx_floodresponse.isel(lndgrid=landcover.index('Restored salt marsh'))['GPP']*3600*24).resample(time='7D').mean().plot(ax=a[1],ls='--',c='orange')

(SETx_floodresponse.isel(lndgrid=landcover.index('Salt marsh'))['TOTVEGC']/1000).resample(time='7D').mean().plot(ax=a[0],ls='--',c='green')
(SETx_floodresponse.isel(lndgrid=landcover.index('Degraded salt marsh'))['TOTVEGC']/1000).resample(time='7D').mean().plot(ax=a[0],ls='--',c='purple')
(SETx_floodresponse.isel(lndgrid=landcover.index('Restored salt marsh'))['TOTVEGC']/1000).resample(time='7D').mean().plot(ax=a[0],ls='--',c='orange')

plt.show()
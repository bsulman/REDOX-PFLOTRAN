import plot_ELM_alquimia_result as pltELM
import xarray
import matplotlib.pyplot as plt
import sys,os
import numpy as np

if len(sys.argv)>1:
    yrs=np.arange(int(sys.argv[1]),int(sys.argv[-1])+1)
else:
    yrs=[1881]

datadir=os.path.expanduser('~/ELM_outputs/')

saline_plantsal=xarray.open_mfdataset([datadir+'alquimia_salinity_US-PHM_ICB20TRCNRDCTCBC/run/alquimia_salinity_US-PHM_ICB20TRCNRDCTCBC.elm.h0.%s-01-01-00000.nc'%yr for yr in yrs]).squeeze()
saline=xarray.open_mfdataset([datadir+'alquimia_salinity_noplantsal_US-PHM_ICB20TRCNRDCTCBC/run/alquimia_salinity_noplantsal_US-PHM_ICB20TRCNRDCTCBC.elm.h0.%s-01-01-00000.nc'%yr for yr in yrs]).squeeze()

fresh=xarray.open_mfdataset([datadir+'alquimia_nosalinity_US-PHM_ICB20TRCNRDCTCBC/run/alquimia_nosalinity_US-PHM_ICB20TRCNRDCTCBC.elm.h0.%s-01-01-00000.nc'%yr for yr in yrs]).squeeze()

import pandas,matplotlib.dates
WT_shad=pandas.read_csv('/home/b0u/PIE_porewater_WT_data/knb-lter-pie.569.1/MAR-RO-Wtable-Shad-2019.csv',
                        parse_dates=[['Date','Time']],na_values='NA')
marsh_height_shad={'S101':1.109,'S102':1.125,'S103':1.095,'S104':1.049}

WT_nelson=pandas.read_csv('/home/b0u/PIE_porewater_WT_data/knb-lter-pie.568.1/MAR-RO-Wtable-Nel-2019.csv',
                            parse_dates=[['Date','Time']],na_values='NA')
marsh_height_nelson={'N201':1.368,'N202':1.360,'N203':1.474,'N204':1.697,'N205':1.812}

porewater1=pandas.read_csv('/home/b0u/PIE_porewater_WT_data/knb-lter-pie.71.6/MAR-PR-Porewater.csv',
                        parse_dates=['Date'],na_values='NA')
porewater2=pandas.read_csv('/home/b0u/PIE_porewater_WT_data/knb-lter-pie.34.20/PIE_MP_porewater_2022.csv',
                        parse_dates=[['YEAR','MONTH','DAY']],na_values='NA')
fluxdata=pandas.read_excel('/home/b0u/PIE_porewater_WT_data/knb-lter-pie.503.2 (Flux data)/MAR-SH-EddyFluxTower-2015_v2.xlsx',sheet_name=1,index_col=0)

vars2=['H2OSFC','salinity','temperature']
vars=['Sulfate','Sulfide','DOC','DIC','CH4','pH']
vmax={'sulfate':25.0,'CH4':1.5,'sulfide':1.0,'salinity':30,'DIC':10.0,'pH':8.8}
vmin={'salinity':0.0,'pH':6.6,'sulfate':0.0}
f,a=plt.subplots(ncols=3,nrows=len(vars),num='Salinity contour comp',clear=True,figsize=(9,9.5))
f2,a2=plt.subplots(ncols=2,nrows=len(vars2),num='Salinity contour comp 2',clear=True,figsize=(9,8))
start=f'{yrs[0]:d}-04-30'
end=f'{yrs[-1]:d}-11-01'
pltELM.plot_vars(saline.sel(time=slice(start,end)),vars,a_contour=a[:,0],a_profile=a[:,2],profile_color='orange',vmax=vmax,vmin=vmin,do_snapshots=False,mean_profile=True,quantiles=[0.1,0.9])
pltELM.plot_vars(fresh.sel(time=slice(start,end)),vars,a_contour=a[:,1],a_profile=a[:,2],profile_color='blue',vmax=vmax,vmin=vmin,do_snapshots=False,mean_profile=True,quantiles=[0.1,0.9])
pltELM.plot_vars(saline.sel(time=slice(start,end)),vars2,a_contour=a2[:,0],a_profile=a2[:,1],profile_color='orange',vmax=vmax,vmin=vmin,do_snapshots=False,mean_profile=True,quantiles=[0.1,0.9])
# pltELM.plot_vars(fresh.isel(time=slice(start,end)),'fresh',vars2,a=a2[:,[1,2]],profile_color='blue',vmax=vmax,do_snapshots=False,mean_profile=True,quantiles=[0.1,0.9])

Shad4=porewater1.set_index('Location').loc['Shad-4m']
for date in Shad4['Date'].unique():
    sample=Shad4.set_index(['Date']).loc[date]
    a2[vars2.index('salinity'),-1].plot(sample['Sal'].astype(float),sample['Depth']/100,'o',ms=0.75,c='k')
    a[vars.index('DOC'),2].plot(sample['DOC'].astype(float)/1000,sample['Depth']/100,'o',ms=0.75,c='k')
    a[vars.index('Sulfide'),2].plot(sample['H2S'].astype(float),sample['Depth']/100,'o',ms=0.75,c='k') # Per Anne Giblin, sulfide concentrations in the LTER dataset are actually in mM, not uM
Shad10=porewater1.set_index('Location').loc['Shad-10m']
for date in Shad10['Date'].unique():
    sample=Shad10.set_index(['Date']).loc[date]
    a2[vars2.index('salinity'),-1].plot(sample['Sal'].astype(float),sample['Depth']/100,'s',ms=0.75,c='k')
    a[vars.index('DOC'),2].plot(sample['DOC'].astype(float)/1000,sample['Depth']/100,'s',ms=0.75,c='k')
    a[vars.index('Sulfide'),2].plot(sample['H2S'].astype(float),sample['Depth']/100,'s',ms=0.75,c='k')

a[1,2].legend(labels=['Saline','Fresh','Measured'])
a2[0,1].legend(labels=['Saline','Fresh','Measured'])
a2[0,0].set_xlim(a2[1,0].get_xlim())

pltELM.letter_label(a)
pltELM.letter_label(a2)

for z in [2,10,20,40]:
    mid=fluxdata[f'TS{z:d}'].mean()
    left=fluxdata[f'TS{z:d}'].quantile(0.1)
    right=fluxdata[f'TS{z:d}'].quantile(0.9)
    a2[vars2.index('temperature'),-1].errorbar(mid,z/100,xerr=np.array([[mid-left],[right-mid]]),c='k',marker='o')

a2[vars2.index('temperature'),-1].set_xlim(right=28)

f,a=plt.subplots(ncols=2,nrows=3,num='Iron',clear=True,figsize=(6.4,7))
vmax['FeS']=saline['soil_FeS'].T.squeeze()[:8,:].max().compute()
vmax['Fe2']=(fresh['soil_Fe2'].T.squeeze()[:9,:]/(saline['watsat'].T.squeeze()[:9,:].to_masked_array())).compute().quantile(0.95)
# vmax['FeS']=480
vmin['FeS']=0
# vmax['FeOxide']=550
vmin['FeOxide']=0.0
pltELM.plot_vars(saline.sel(time=slice(start,end)),['Fe2'],a_profile=a[0,1],profile_color='orange',vmax=vmax,vmin=vmin,do_snapshots=False,mean_profile=True,quantiles=[0.1,0.9])
pltELM.plot_vars(fresh.sel(time=slice(start,end)),['Fe2'],a_contour=a[0,0],a_profile=a[0,1],profile_color='blue',vmax=vmax,vmin=vmin,do_snapshots=False,mean_profile=True,quantiles=[0.1,0.9])
a[0,1].legend(['Saline','Fresh'])
pltELM.plot_vars(saline.sel(time=slice(start,end)),['FeOxide','FeS'],a_profile=a[[1,2],1],profile_color='orange',vmax=vmax,vmin=vmin,do_snapshots=False,mean_profile=True,quantiles=[0.1,0.9])
pltELM.plot_vars(fresh.sel(time=slice(start,end)),['FeOxide','FeS'],a_profile=a[[1,2],1],profile_color='blue',vmax=vmax,vmin=vmin,do_snapshots=False,mean_profile=True,quantiles=[0.1,0.9])
# (fresh['soil_FeOxide'].T.squeeze()[:10,:].diff(dim='time')*1e3).mean(dim='time').plot(ax=a[1,0],vmin=-3e-2,vmax=3e-2,cmap='RdBu_r',cbar_kwargs={'label':'Fe Oxide change (mmol m$^{-3}$ day$^{-1}$)'})
a[1,0].plot((fresh['soil_FeOxide'].T.squeeze()[:10,:].diff(dim='time')*1e3*24).mean(dim='time'),fresh['levdcmp'][:10],c='blue')
a[1,0].plot((saline['soil_FeOxide'].T.squeeze()[:10,:].diff(dim='time')*1e3*24).mean(dim='time'),fresh['levdcmp'][:10],c='orange')
a[1,0].set(title='Fe oxide change',ylim=(1.5,0),xlabel='Fe Oxide change (mmol Fe m$^{-3}$ day$^{-1}$)',ylabel='Soil depth (m)',xlim=(-2.1e-2,3e-2))
a[1,0].axvline(0,c='k',lw=0.5,ls=':')

a[2,0].plot((fresh['soil_FeS'].T.squeeze()[:10,:].diff(dim='time')*1e3*24).mean(dim='time'),fresh['levdcmp'][:10],c='blue')
a[2,0].plot((saline['soil_FeS'].T.squeeze()[:10,:].diff(dim='time')*1e3*24).mean(dim='time'),fresh['levdcmp'][:10],c='orange')
a[2,0].set(title='Fe sulfide change',ylim=(1.5,0),xlabel='Fe Sulfide change (mmol Fe m$^{-3}$ day$^{-1}$)',ylabel='Soil depth (m)')
a[2,0].axvline(0,c='k',lw=0.5,ls=':')

pltELM.letter_label(a)

# start_zoom='%d-04-15'%yrs[-1] #int(8760/12*5)
# end_zoom='%d-08-31'%yrs[-1] #int(8760/12*7)
# f,a=plt.subplots(ncols=3,nrows=len(vars),num='Salinity contour comp zoomed',clear=True,figsize=(8,9))
# pltELM.plot_vars(saline.sel(time=slice(start_zoom,end_zoom)),'saline',vars,a=a[:,[0,2]],profile_color='orange',vmax=vmax,do_snapshots=False,mean_profile=True,quantiles=[0.1,0.9])
# pltELM.plot_vars(fresh.sel(time=slice(start_zoom,end_zoom)),'fresh',vars,a=a[:,[1,2]],profile_color='blue',vmax=vmax,do_snapshots=False,mean_profile=True,quantiles=[0.1,0.9])

# a[-1,2].legend(labels=['Saline','Fresh'])

CH4flux_obs=pandas.read_excel('/home/b0u/PIE_porewater_WT_data/PLM_CH4_gap-filled.xlsx',sheet_name=0,index_col=0)
import datetime
CH4flux_resampled=(CH4flux_obs['gap-filled'].resample('7D').mean()*1e-3)

f,a=plt.subplots(num='GHG fluxes',clear=True,nrows=4,figsize=(5.6,8.5))
# (saline['CH4FLUX_ALQUIMIA']/12.011*1e6).plot(ax=a,c='orange',label='Saline')
# (fresh['CH4FLUX_ALQUIMIA']/12.011*1e6).plot(ax=a,c='blue',label='Fresh')
(saline['CH4FLUX_ALQUIMIA']/12.011*1e6).resample(time='7D').mean().plot(ax=a[0],c='orange',label='Saline')
(fresh['CH4FLUX_ALQUIMIA']/12.011*1e6).resample(time='7D').mean().plot(ax=a[0],c='blue',label='Fresh')
(saline_plantsal['CH4FLUX_ALQUIMIA']/12.011*1e6).resample(time='7D').mean().plot(ax=a[0],c='red',label='Saline + plant response')

(saline['CH4FLUX_ALQUIMIA']/12.011*1e6).resample(time='7D').mean().plot(ax=a[1],c='orange',label='Saline')
(fresh['CH4FLUX_ALQUIMIA']/12.011*1e6).resample(time='7D').mean().plot(ax=a[1],c='blue',label='Fresh')
(saline_plantsal['CH4FLUX_ALQUIMIA']/12.011*1e6).resample(time='7D').mean().plot(ax=a[1],c='red',label='Saline + plant response')


a[0].plot(fresh['time'][0].item()+(CH4flux_resampled.index-datetime.datetime(2019,1,1)).to_pytimedelta(),
        CH4flux_resampled,c='red',label='Saline obs',linestyle='--')

a[0].set(title='Methane flux',xlabel='Time',ylabel='Methane flux ($\mu$mol m$^{-2}$ s$^{-1}$)')
a[0].legend()

a[1].plot(fresh['time'][0].item()+(CH4flux_resampled.index-datetime.datetime(2019,1,1)).to_pytimedelta(),
        CH4flux_resampled,c='red',label='Saline obs',linestyle='--')

a[1].set(title='Methane flux',xlabel='Time',ylabel='Methane flux ($\mu$mol m$^{-2}$ s$^{-1}$)',ylim=(-0.001,0.005))
a[1].legend()
a[1].axhline(0.0,c='k',lw=0.5,ls=':')

(saline['HR']/12.011*1e6).resample(time='7D').mean().plot(ax=a[2],c='orange',label='Saline')
(fresh['HR']/12.011*1e6).resample(time='7D').mean().plot(ax=a[2],c='blue',label='Fresh')
(saline_plantsal['HR']/12.011*1e6).resample(time='7D').mean().plot(ax=a[2],c='red',label='Saline + plant response')
a[2].set(title='Soil CO$_2$ flux',xlabel='Time',ylabel='CO$_2$ flux ($\mu$mol m$^{-2}$ s$^{-1}$)')
a[2].legend()

(saline['DIC_RUNOFF']/12.011*1e6).resample(time='7D').mean().plot(ax=a[3],c='orange',ls='-')
(fresh['DIC_RUNOFF']/12.011*1e6).resample(time='7D').mean().plot(ax=a[3],c='blue',ls='-')
(saline_plantsal['DIC_RUNOFF']/12.011*1e6).resample(time='7D').mean().plot(ax=a[3],c='red',ls='-')
a[3].set(title='Soil lateral DIC flux',xlabel='Time',ylabel='DIC flux ($\mu$mol C m$^{-2}$ s$^{-1}$)')
# a[3].legend()
a[3].set_ylim(bottom=-0.01)

pltELM.letter_label(a)

f,a=plt.subplots(num='SOC',clear=True)
pltELM.plot_vars(fresh,['soilC'],'SOC',a_profile=[a],profile_color='blue',do_snapshots=False,mean_profile=True,quantiles=[0.1,0.9],maxdepth=1.0)
pltELM.plot_vars(saline,['soilC'],'SOC',a_profile=[a],profile_color='orange',do_snapshots=False,mean_profile=True,quantiles=[0.1,0.9],maxdepth=1.0)
pltELM.plot_vars(saline_plantsal,['soilC'],'SOC',a_profile=[a],profile_color='red',do_snapshots=False,mean_profile=True,quantiles=[0.1,0.9],maxdepth=1.0)
# pltELM.letter_label(a)
soilprofile=pandas.read_csv('/home/b0u/PIE_porewater_WT_data/Spivak_soil_profiles/dataset-827298_bulk-soil-properties__v1.tsv',delimiter='\t',na_values=['nd'])
profile_gr=soilprofile[(soilprofile['Location']=='MARSH')].groupby(['Depth_min'])
a.errorbar(profile_gr['Carbon_Density'].mean()*1e3/1e6,profile_gr['Depth_min'].mean()/100,xerr=profile_gr['Carbon_Density'].sem()*1e3/1e6)
a.legend(labels=['Fresh','Saline','Saline + reduced GPP','Obs'],handles=a.lines)

start_zoom='%d-07-15 01:00'%yrs[-1] #int(8760/12*5)
end_zoom='%d-07-16 03:00'%yrs[-1] #int(8760/12*7)
f,a=plt.subplots(num='Tides and fluxes',clear=True,nrows=2,sharex=True,figsize=(7.88, 6.92))
(fresh['H2OSOI'].T.squeeze()/fresh['watsat'].T.squeeze()[:10,:]).sel(time=slice(start_zoom,end_zoom)).plot(ax=a[0],cmap='Blues',cbar_kwargs={'label':'Volumetric water\n(fraction of saturation)'})
(-fresh['H2OSFC']/1000).where(fresh['H2OSFC']>0).sel(time=slice(start_zoom,end_zoom)).plot(ax=a[0])
# a[0].set(title='Water level',xlabel='',ylabel='Water level (mm)')
a[0].set(title='Soil water',ylim=(0.7,-.75),xlabel='',ylabel='Soil depth (m)')

# (saline['GPP']/12.011*1e6).sel(time=slice(start_zoom,end_zoom)).plot(ax=a[2],c='orange',label='Saline')
# (saline_plantsal['GPP']/12.011*1e6).sel(time=slice(start_zoom,end_zoom)).plot(ax=a[2],c='red',label='Saline (plant response)')
# (fresh['GPP']/12.011*1e6).sel(time=slice(start_zoom,end_zoom)).plot(ax=a[2],c='blue',label='Fresh')
(saline['HR']/12.011*1e6).sel(time=slice(start_zoom,end_zoom)).plot(ax=a[1],c='orange',label='Soil CO$_2$ flux (Saline)',ls='--')
(saline_plantsal['HR']/12.011*1e6).sel(time=slice(start_zoom,end_zoom)).plot(ax=a[1],c='red',label='Soil CO$_2$ flux (Saline + lower GPP)',ls='--')
(fresh['HR']/12.011*1e6).sel(time=slice(start_zoom,end_zoom)).plot(ax=a[1],c='blue',label='Soil CO$_2$ flux (Fresh)',ls='--')

(saline['CH4FLUX_ALQUIMIA']/12.011*1e6).sel(time=slice(start_zoom,end_zoom)).plot(ax=a[1],c='orange',label='Methane flux (Saline)')
(fresh['CH4FLUX_ALQUIMIA']/12.011*1e6).sel(time=slice(start_zoom,end_zoom)).plot(ax=a[1],c='blue',label='Methane flux (Fresh)')
(saline_plantsal['CH4FLUX_ALQUIMIA']/12.011*1e6).sel(time=slice(start_zoom,end_zoom)).plot(ax=a[1],c='red',label='Methane flux (Saline + lower GPP)')

(saline['DIC_RUNOFF']/12.011*1e6).sel(time=slice(start_zoom,end_zoom)).plot(ax=a[1],c='orange',label='Lateral DIC flux (Saline)',ls=':')
(saline_plantsal['DIC_RUNOFF']/12.011*1e6).sel(time=slice(start_zoom,end_zoom)).plot(ax=a[1],c='red',label='Lateral DIC flux (Saline + lower GPP)',ls=':')
(fresh['DIC_RUNOFF']/12.011*1e6).sel(time=slice(start_zoom,end_zoom)).plot(ax=a[1],c='blue',label='Lateral DIC flux (Fresh)',ls=':')

a[1].set(title='Soil carbon fluxes',xlabel='Time',ylabel='Carbon flux ($\mu$ mol C m$^{-2}$ s$^{-1}$)')
a[1].legend(ncol=1,loc='upper left')
a[1].xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b-%d'))

pltELM.letter_label(a)

import pandas,matplotlib.dates
WT_shad=pandas.read_csv('/home/b0u/PIE_porewater_WT_data/knb-lter-pie.569.1/MAR-RO-Wtable-Shad-2019.csv',
                        parse_dates=[['Date','Time']],na_values='NA')
marsh_height_shad={'S101':1.109,'S102':1.125,'S103':1.095,'S104':1.049}

WT_nelson=pandas.read_csv('/home/b0u/PIE_porewater_WT_data/knb-lter-pie.568.1/MAR-RO-Wtable-Nel-2019.csv',
                            parse_dates=[['Date','Time']],na_values='NA')
marsh_height_nelson={'N201':1.368,'N202':1.360,'N203':1.474,'N204':1.697,'N205':1.812}

porewater1=pandas.read_csv('/home/b0u/PIE_porewater_WT_data/knb-lter-pie.71.6/MAR-PR-Porewater.csv',
                        parse_dates=['Date'],na_values='NA')
porewater2=pandas.read_csv('/home/b0u/PIE_porewater_WT_data/knb-lter-pie.34.20/PIE_MP_porewater_2022.csv',
                        parse_dates=[['YEAR','MONTH','DAY']],na_values='NA')

marshcol=saline
f,a=plt.subplots(num='PIE data',clear=True,nrows=3,ncols=4,figsize=(10,4))
ax=a[0,0]
for wellnum in marsh_height_shad.keys():
    h=ax.plot(WT_shad['Date_Time'],WT_shad['Logger '+wellnum+'  Dc (m)']-marsh_height_shad[wellnum],label=wellnum)[0]
ax.axhline(marsh_height_shad[wellnum]-marsh_height_shad[wellnum],ls=':',lw=1.0,c='k')

# ax.plot(WT_shad['Date_Time'],WT_shad['Logger Sref (tide gauge)   Dc (m)'],'k--',label='Sref')
ax.legend()
ax.set_xlim((matplotlib.dates.datestr2num('2019-07-01'),matplotlib.dates.datestr2num('2019-07-04')))
ax.set(title='Shad water table',xlabel='Water table (m)',ylabel='Date',ylim=(-1,1))
ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b-%d'))

ax=a[1,0]
for wellnum in marsh_height_nelson.keys():
    h=ax.plot(WT_nelson['Date_Time'],WT_nelson['Logger '+wellnum+'  Dc (m)']-marsh_height_nelson[wellnum],label=wellnum)[0]
ax.axhline(marsh_height_nelson[wellnum]-marsh_height_nelson[wellnum],ls=':',lw=1.0,c='k')

# ax.plot(WT_nelson['Date_Time'],WT_nelson['Logger Nref (tide gauge)   Dc (m)'],'k--',label='Sref')
ax.legend()
ax.set_xlim((matplotlib.dates.datestr2num('2019-07-01'),matplotlib.dates.datestr2num('2019-07-04')))
ax.set(title='Nelson water table',xlabel='Water table (m)',ylabel='Date',ylim=(-1,1))
ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b-%d'))

yr=marshcol.time[0].item().year
timeslice=slice(f'{yr}-07-01',f'{yr}-07-03')
H2OSFC=marshcol['H2OSFC'].sel(time=timeslice)/1000
WTD=marshcol['ZWT'].sel(time=timeslice)
a[2,0].plot(H2OSFC.time,H2OSFC.where(H2OSFC>0,-WTD))
a[2,0].plot(H2OSFC.time,marshcol['H2OSFC_TIDE'].sel(time=timeslice)/1000,ls='--')
a[2,0].axhline(0,ls=':',c='k')
a[2,0].xaxis.set_major_formatter(pltELM.format_nc_time('%b-%d'))
a[2,0].set(title='Model water table',xlabel='Water table (m)',ylabel='Date',ylim=(-1,1))

Nelson4=porewater1.set_index('Location').loc['Nelson-4m']
for date in Nelson4['Date'].unique():
    sample=Nelson4.set_index(['Date']).loc[date]
    a[1,1].plot(sample['Sal'].astype(float),sample['Depth'],'o')
    a[1,2].plot(sample['DOC'].astype(float),sample['Depth'],'o')
    a[1,3].plot(sample['H2S'].astype(float)*1000,sample['Depth'],'o') # Per Anne Giblin, units are wrong in data file for H2S and it should be in mM
Nelson10=porewater1.set_index('Location').loc['Nelson-10m']
for date in Nelson10['Date'].unique():
    sample=Nelson10.set_index(['Date']).loc[date]
    a[1,1].plot(sample['Sal'].astype(float),sample['Depth'],'s')
    a[1,2].plot(sample['DOC'].astype(float),sample['Depth'],'s')
    a[1,3].plot(sample['H2S'].astype(float)*1000,sample['Depth'],'s')

a[1,1].set(title='Nelson salinity',xlabel='Salinity (ppt)',ylabel='Depth (cm)',ylim=(100,0))
a[1,2].set(title='Nelson DOC',xlabel='DOC concentration ($\mu$M)',ylabel='Depth (cm)',ylim=(100,0))
a[1,3].set(title='Nelson sulfide',xlabel='Sulfide concentration ($\mu$M)',ylabel='Depth (cm)',ylim=(100,0))

Shad4=porewater1.set_index('Location').loc['Shad-4m']
for date in Shad4['Date'].unique():
    sample=Shad4.set_index(['Date']).loc[date]
    a[0,1].plot(sample['Sal'].astype(float),sample['Depth'],'o')
    a[0,2].plot(sample['DOC'].astype(float),sample['Depth'],'o')
    a[0,3].plot(sample['H2S'].astype(float)*1000,sample['Depth'],'o')
Shad10=porewater1.set_index('Location').loc['Shad-10m']
for date in Shad10['Date'].unique():
    sample=Shad10.set_index(['Date']).loc[date]
    a[0,1].plot(sample['Sal'].astype(float),sample['Depth'],'s')
    a[0,2].plot(sample['DOC'].astype(float),sample['Depth'],'s')
    a[0,3].plot(sample['H2S'].astype(float)*1000,sample['Depth'],'s')

a[0,1].set(title='Shad salinity',xlabel='Salinity (ppt)',ylabel='Depth (cm)',ylim=(100,0))
a[0,2].set(title='Shad DOC',xlabel='DOC concentration ($\mu$M)',ylabel='Depth (cm)',ylim=(100,0))
a[0,3].set(title='Shad sulfide',xlabel='Sulfide concentration ($\mu$M)',ylabel='Depth (cm)',ylim=(100,0))

porosity=marshcol['watsat'].sel(time=slice(start,end)).T.squeeze()[:10,:].to_masked_array()
mod_sulfide_monthly=marshcol['soil_sulfide'].resample(time='1M').mean()[:,:10]*1e3/porosity[:10,0]
mod_DOC_monthly=marshcol['DOC_vr'].resample(time='1M').mean()[:,:10]*1e3/porosity[:10,0]/12.011
mod_salinity_monthly=marshcol['soil_salinity'].resample(time='1M').mean()[:,:10]
for m in range(5,10):
    a[2,3].plot(mod_sulfide_monthly.isel(time=m),mod_sulfide_monthly.levdcmp*100)
    a[2,2].plot(mod_DOC_monthly.isel(time=m),mod_sulfide_monthly.levdcmp*100)
    a[2,1].plot(mod_salinity_monthly.isel(time=m),mod_salinity_monthly.levdcmp*100)

a[2,1].set(title='Model salinity',xlabel='Salinity (ppt)',ylabel='Depth (cm)',ylim=(100,0))
a[2,3].set(title='Model sulfide',xlabel='Sulfide concentration ($\mu$M)',ylabel='Depth (cm)',ylim=(100,0))
a[2,2].set(title='Model DOC',xlabel='DOC concentration ($\mu$M)',ylabel='Depth (cm)',ylim=(100,0))

pltELM.letter_label(a)

f,a=plt.subplots(nrows=4,ncols=1,clear=True,num='Tide demo',figsize=(5,9),sharex=True)
pltELM.plot_vars(saline.isel(time=slice(int(8760/12*6.1),int(8760/12*7.2))),a_contour=a,
                 vars=['H2OSFC','VWC','salinity','oxygen'],maxdepth=0.45,vmax={'salinity':35})
pltELM.letter_label(a)

def save_all_figs(dirname,format='png',**kwargs):
    for fname in plt.get_figlabels():
        fname_fixed=fname.replace('/','-')
        print(fname_fixed)
        plt.figure(fname_fixed).savefig('{dirname:s}/{fname:s}.{format}'.format(dirname=dirname,format=format,fname=fname_fixed),**kwargs)



plt.show()
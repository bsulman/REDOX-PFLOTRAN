import xarray
import matplotlib.pyplot as plt
import matplotlib.colors
import sys
import numpy as np
from plot_ELM_alquimia_result import format_nc_time

d=xarray.open_dataset(sys.argv[1])

marshcol=d.isel(lndgrid=0)

def plot_tides(a,start=0,end=8760,order=['H2OSFC','vertflow','SWC','oxygen','sulfate','sulfide','salinity','DOC','DIC','CH4'],salinity_vmax=80,sulfate_vmax=None,methane_vmax=0.5,zmax=1.5):
    for num,var in enumerate(order):
        VWC=(marshcol['H2OSOI']/marshcol['watsat']).T.squeeze()[:10,start:end]
        porosity=marshcol['watsat'].T.squeeze()[:10,start:end].to_masked_array()
        if var=='SWC':
            VWC.plot(ax=a[num],vmin=0,vmax=1.0,cmap='Blues',cbar_kwargs={'label':'Volumetric water\n(fraction of saturation)'})
            a[num].set(ylim=(zmax,0),title='Soil water content',ylabel='Depth (m)',xlabel='Time')
        elif var=='H2OSFC':
            marshcol['H2OSFC'][start:end].plot(ax=a[num])
            a[num].set(title='Surface water depth',ylabel='Surface water depth\n(mm)',xlabel='Time')
        elif var=='latflow':
            marshcol['QFLX_LAT_AQU'][start:end].plot(ax=a[num])
            a[num].set(title='Lateral water flow into column',xlabel='Time',ylabel='Water flow rate\n(mm s$^{-1}$)')
            a[num].axhline(0,lw=0.5,ls=':',c='k')
        elif var=='drain':
            marshcol['QDRAI'][start:end].plot(ax=a[num])
            a[num].set(title='Drainage flow',xlabel='Time',ylabel='Water flow rate\n(mm s$^{-1}$)')
            a[num].axhline(0,lw=0.5,ls=':',c='k')
        elif var=='vertflow':
            marshcol['QFLX_ADV'].T[:10,start:end].plot(ax=a[num],vmax=0.005)
            a[num].set(ylim=(zmax,0),title='Vertical flux',ylabel='Depth (m)',xlabel='Time')
        elif var=='drain_vr':
            (marshcol['QDRAI_VR'].T[:10,start:end]/3600).plot(ax=a[num],cbar_kwargs={'label':'Drainage flux\n(mm s$^{-1}$)'})
            a[num].set(ylim=(zmax,0),title='Subsurface drainage by layer',ylabel='Depth (m)',xlabel='Time')
        elif var=='salinity':
            (marshcol['soil_salinity'].T[:10,start:end]).plot(ax=a[num],cbar_kwargs={'label':'Salinity (ppt)'},vmin=0,vmax=salinity_vmax)
            a[num].set(ylim=(zmax,0),title='Soil salinity',ylabel='Depth (m)',xlabel='Time')
        elif var=='oxygen':
            # 32 converting from mol/m3 to mg/L = 16 g/mol * 2 (O2) *1000mg/g * 1/1000 m3/L
            (marshcol['soil_O2'].T[:10,start:end]/porosity[:10,:]).plot(ax=a[num],cbar_kwargs={'label':'O$_2$ concentration\n(mmol L$^{-1}$)'})
            a[num].set(ylim=(0.3,0),title='Soil O$_2$ concentration',ylabel='Depth (m)',xlabel='Time')
        elif var=='DOC':
            # Convert to umol/L 1e6 umol/mol * 1e-3 m3/L
            (marshcol['DOC_vr'].T[:10,start:end]/porosity[:10,:]/12.011).plot(ax=a[num],cbar_kwargs={'label':'DOC concentration\n(mmol L$^{-1}$)'})
            a[num].set(ylim=(zmax,0),title='Soil DOC concentration',ylabel='Depth (m)',xlabel='Time')
        elif var=='DIC':
            (marshcol['DIC_vr'].T[:10,start:end]/porosity[:10,:]/12.011).plot(ax=a[num],cbar_kwargs={'label':'DIC concentration\n(mmol L$^{-1}$)'})
            a[num].set(ylim=(zmax,0),title='Soil DIC concentration',ylabel='Depth (m)',xlabel='Time')
        elif var=='CH4':
            (marshcol['CH4_vr'].T[:10,start:end]/porosity[:10,:]*1/12.011).plot(ax=a[num],cbar_kwargs={'label':'CH4 concentration\n(mmol L$^{-1}$)'},vmax=methane_vmax)
            a[num].set(ylim=(zmax,0),title='Soil CH4 concentration',ylabel='Depth (m)',xlabel='Time')
        elif var=='sulfate':
            (marshcol['soil_sulfate'].T[:10,start:end]/porosity[:10,:]).plot(ax=a[num],cbar_kwargs={'label':'Sulfate concentration\n(mmol L$^{-1}$)'},vmax=sulfate_vmax)
            a[num].set(ylim=(zmax,0),title='Soil sulfate concentration',ylabel='Depth (m)',xlabel='Time')
        elif var=='sulfide':
            (marshcol['soil_sulfide'].T[:10,start:end]/porosity[:10,:]).plot(ax=a[num],cbar_kwargs={'label':'Sulfide concentration\n(mmol L$^{-1}$)'})
            a[num].set(ylim=(zmax,0),title='Soil sulfide concentration',ylabel='Depth (m)',xlabel='Time')
        elif var=='sulfide_ratio':
            (marshcol['soil_sulfide'].T[:10,start:end]/marshcol['soil_sulfate'].T[:10,start:end]).plot(ax=a[num],cbar_kwargs={'label':'Sulfide:sulfate ratio'})
            a[num].set(ylim=(zmax,0),title='Soil sulfide:sulfate ratio',ylabel='Depth (m)',xlabel='Time')
        elif var=='pH':
            (marshcol['soil_pH'].T[:10,start:end]).plot(ax=a[num],cbar_kwargs={'label':'pH'})
            a[num].set(ylim=(zmax,0),title='Soil pH',ylabel='Depth (m)',xlabel='Time')
        elif var=='NEE':
            (marshcol['NEE'][start:end]*1e6/12).plot(ax=a[num])
            a[num].set(title='Net ecosystem CO$_2$ exchange',xlabel='Time',ylabel='NEE (umol C m$^{-2}$ s$^{-1}$)')
            a[num].axhline(0.0,color='k',ls=':',lw=0.5)
        elif var=='CH4flux' or var=='CH4FLUX_ALQUIMIA':
            (marshcol['CH4FLUX_ALQUIMIA'][start:end]*1e6/12).plot(ax=a[num])
            a[num].set(title='Surface CH$_4$ emission',xlabel='Time',ylabel='CH$_4$ emission\n(umol C m$^{-2}$ s$^{-1}$)')

        else:
            # raise ValueError('Plot type %s not implemented'%var)
            marshcol[var].isel(time=slice(start,end)).plot(ax=a[num])

f,a=plt.subplots(nrows=5,ncols=2,clear=True,num='Tidal result: Year',figsize=(8,9))
plot_tides(a.ravel())

f,a=plt.subplots(nrows=5,ncols=2,clear=True,num='Tidal result: Months',figsize=(8,9))
plot_tides(a.ravel(),start=int(8760/12*6),end=int(8760/12*7.5))

f,a=plt.subplots(nrows=4,ncols=1,clear=True,num='Tide demo',figsize=(5,9),sharex=True)
plot_tides(a,start=int(8760/12*6.1),end=int(8760/12*7.2),order=['H2OSFC','SWC','salinity','oxygen'],zmax=0.45,sulfate_vmax=200)
a[-1].xaxis.set_major_formatter(format_nc_time('%b-%d'))

# forcing=xarray.open_dataset('/home/b0u/tidesPLM_2018_gapfill.nc',decode_times=False)
start=int(8760/12*5)
end=int(8760/12*6.5)
start=0
end=-1
f,a=plt.subplots(num='Salinity forcing',clear=True)
# a.plot(d['time'][forcing['time'].values][start:end],forcing['tide_salinity'].squeeze()[start:end])
marshcol['SALINITY'][start:end].plot(ax=a)
a.set(title='Channel salinity forcing',xlabel='Time',ylabel='Salinity (ppt)')

f,a=plt.subplots(nrows=5,ncols=2,clear=True,num='Tidal result: Days',figsize=(8,9))
start=195
end=197
if len(d.lndgrid)==2:
    (d.isel(lndgrid=1)['H2OSFC'][int(8760/365*start):int(8760/365*end)]-800).plot(ax=a[0,0],ls='--')
else:
    (d.isel(lndgrid=0)['H2OSFC_TIDE'][int(8760/365*start):int(8760/365*end)]).plot(ax=a[0,0],ls='--')
plot_tides(a.ravel(),start=int(8760/365*start),end=int(8760/365*end))
a[0,0].legend(labels=['Tidal channel','Marsh'])



z=marshcol['levdcmp'][:10]
dz=np.array([float(n) for n in '1.7512817916255204E-002   2.7578969259676251E-002   4.5470033242413201E-002   7.4967410986208557E-002  0.12360036510228053       0.20378255101043175       0.33598062644843263       0.55393840536868488       0.91329003158906108        1.5057607013992766        2.4825796969813321        4.0930819526214002        6.7483512780057175        11.126150294204420        13.851152141963599'.split()])[:10]

start=int(8760/365*195)
end=int(8760/365*(end+1))
maxtide=marshcol['H2OSFC'][start:end].argmax()
mintide=marshcol['H2OSFC'][start:end].argmin()
sal=marshcol['soil_salinity'][start:end,:10]

# a.plot(sal[maxtide],z,label='High tide')
# a.plot(sal[maxtide-12],z,label='Low tide')
# a.plot(sal[sal.mean(dim='levdcmp').argmax()],z,label='High salinity')
# a.plot(sal[sal.mean(dim='levdcmp').argmin()],z,label='Low salinity')


f,a=plt.subplots(nrows=5,num='Salinity time series',clear=True,figsize=(4,8))

cm=plt.get_cmap('Blues')
norm=matplotlib.colors.Normalize(0.5,1.0)
for n in range(12):
    x=sal.mean(dim='levdcmp')[:-12].argmin()
    a[0].plot(sal[x+n],z,c=cm(norm((marshcol['H2OSOI']/marshcol['watsat']).T.squeeze()[:7,start+n].mean())))
    a[1].plot(marshcol['soil_O2'][start:end,:10][x+n]*32,z,c=cm(norm((marshcol['H2OSOI']/marshcol['watsat']).T.squeeze()[:7,start+n].mean())),label='Moisture = %1.1f'%((marshcol['H2OSOI']/marshcol['watsat']).T.squeeze()[:7,start+n].mean()))
    a[3].plot((marshcol['H2OSOI']/marshcol['watsat']).T.squeeze()[:10,start+n],z,c=cm(norm((marshcol['H2OSOI']/marshcol['watsat']).T.squeeze()[:7,start+n].mean())))
    a[2].plot(marshcol['DOC_vr'][start:end,:10][x+n]/12,z,c=cm(norm((marshcol['H2OSOI']/marshcol['watsat']).T.squeeze()[:7,start+n].mean())),label='Moisture = %1.1f'%((marshcol['H2OSOI']/marshcol['watsat']).T.squeeze()[:7,start+n].mean()))
    a[4].plot(marshcol['DIC_vr'][start:end,:10][x+n]/12,z,c=cm(norm((marshcol['H2OSOI']/marshcol['watsat']).T.squeeze()[:7,start+n].mean())),label='Moisture = %1.1f'%((marshcol['H2OSOI']/marshcol['watsat']).T.squeeze()[:7,start+n].mean()))
    a[4].plot(marshcol['CH4_vr'][start:end,:10][x+n]/12,z,c=cm(norm((marshcol['H2OSOI']/marshcol['watsat']).T.squeeze()[:7,start+n].mean())),label='Moisture = %1.1f'%((marshcol['H2OSOI']/marshcol['watsat']).T.squeeze()[:7,start+n].mean()),ls='--')

# a[0].legend()
zmax=0.5
a[0].set(ylim=(zmax,0),title='Salinity profiles',xlabel='Salinity (ppt)',ylabel='Depth (m)')
a[1].set(ylim=(zmax,0),title='Oxygen profiles',xlabel='Oxygen (mg L$^{-1}$)',ylabel='Depth (m)')
a[1].legend()
a[3].set(ylim=(zmax,0),title='Moisture profiles',xlabel='Moisture (fraction of saturation)',ylabel='Depth (m)')
a[2].set(ylim=(zmax,0),title='DOC profiles',xlabel='DOC (mmol C L$^{-1}$)',ylabel='Depth (m)')
a[4].set(ylim=(zmax,0),title='DIC profiles',xlabel='DIC (mmol C L$^{-1}$)',ylabel='Depth (m)')



f,a=plt.subplots(num='Soil C profile',clear=True)
a.plot((marshcol['LITR1C_vr']+marshcol['LITR2C_vr']+marshcol['LITR3C_vr']+marshcol['SOIL1C_vr']+marshcol['SOIL2C_vr']+marshcol['SOIL3C_vr']+marshcol['SOIL4C_vr']).mean(dim='time')[:10],z,'k-',label='Total')
a.plot((marshcol['LITR1C_vr']+marshcol['LITR2C_vr']+marshcol['LITR3C_vr']+marshcol['SOIL1C_vr']+marshcol['SOIL2C_vr']+marshcol['SOIL3C_vr']+marshcol['SOIL4C_vr']).max(dim='time')[:10],z,'k--')
a.plot((marshcol['LITR1C_vr']+marshcol['LITR2C_vr']+marshcol['LITR3C_vr']+marshcol['SOIL1C_vr']+marshcol['SOIL2C_vr']+marshcol['SOIL3C_vr']+marshcol['SOIL4C_vr']).min(dim='time')[:10],z,'k--')

for var in ['LITR1C','LITR2C','LITR3C','SOIL1C','SOIL2C','SOIL3C','SOIL4C']:
    l=a.plot(marshcol[var+'_vr'].mean(dim='time')[:10],z,label=var)[0]
    a.plot(marshcol[var+'_vr'].max(dim='time')[:10],z,ls='--',color=l.get_color())
    a.plot(marshcol[var+'_vr'].min(dim='time')[:10],z,ls='--',color=l.get_color())

a.set(title='Soil organic matter pools',xlabel='C concentration (g C m$^{-3}$)',ylabel='Depth (m)',ylim=(zmax,0))
a.legend()

f,a=plt.subplots(num='Demo plot',clear=True,nrows=4)
plot_tides(a.ravel(),start=int(8760/12*6),end=int(8760/12*7.5),order=['NEE','CH4FLUX_ALQUIMIA','H2OSFC','salinity'])

f,a=plt.subplots(num='Hydro plot',clear=True,nrows=6,sharex=True)
# (d.isel(lndgrid=1)['H2OSFC'][int(8760/12*6):int(8760/12*7.0)]-800).plot(ax=a[0],ls='--',c='navy')
plot_tides(a.ravel(),start=int(8760/12*6),end=int(8760/12*6.1),order=['H2OSFC','latflow','drain','drain_vr','SWC','salinity'],zmax=1.25,salinity_vmax=40)
# a[0].legend(labels=['Tidal channel','Marsh'])
a[3].xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b-%d'))

# Plum Island water table
import pandas
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

f,a=plt.subplots(num='PIE data',clear=True,nrows=3,ncols=4,figsize=(10,4))
ax=a[0,0]
for wellnum in marsh_height_shad.keys():
    h=ax.plot(WT_shad['Date_Time'],WT_shad['Logger '+wellnum+'  Dc (m)']-marsh_height_shad[wellnum],label=wellnum)[0]
ax.axhline(marsh_height_shad[wellnum]-marsh_height_shad[wellnum],ls=':',lw=1.0,c='k')

# ax.plot(WT_shad['Date_Time'],WT_shad['Logger Sref (tide gauge)   Dc (m)'],'k--',label='Sref')
ax.legend()
ax.set_xlim((matplotlib.dates.datestr2num('2019-07-01'),matplotlib.dates.datestr2num('2019-07-04')))
ax.set(title='Shad water table',xlabel='Water table (m)',ylabel='Date')
ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b-%d'))

ax=a[1,0]
for wellnum in marsh_height_nelson.keys():
    h=ax.plot(WT_nelson['Date_Time'],WT_nelson['Logger '+wellnum+'  Dc (m)']-marsh_height_nelson[wellnum],label=wellnum)[0]
ax.axhline(marsh_height_nelson[wellnum]-marsh_height_nelson[wellnum],ls=':',lw=1.0,c='k')

# ax.plot(WT_nelson['Date_Time'],WT_nelson['Logger Nref (tide gauge)   Dc (m)'],'k--',label='Sref')
ax.legend()
ax.set_xlim((matplotlib.dates.datestr2num('2019-07-01'),matplotlib.dates.datestr2num('2019-07-04')))
ax.set(title='Nelson water table',xlabel='Water table (m)',ylabel='Date')
ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b-%d'))

yr=marshcol.time[0].item().year
timeslice=slice(f'{yr}-07-01',f'{yr}-07-03')
H2OSFC=marshcol['H2OSFC'].sel(time=timeslice)/1000
WTD=marshcol['ZWT'].sel(time=timeslice)
a[2,0].plot(H2OSFC.time,H2OSFC.where(H2OSFC>0,-WTD))
a[2,0].axhline(0,ls=':',c='k')
a[2,0].xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b-%d'))
a[2,0].set(title='Model water table',xlabel='Water table (m)',ylabel='Date')

Nelson4=porewater1.set_index('Location').loc['Nelson-4m']
for date in Nelson4['Date'].unique():
    sample=Nelson4.set_index(['Date']).loc[date]
    a[1,1].plot(sample['Sal'].astype(float),sample['Depth'])
    a[1,2].plot(sample['DOC'].astype(float),sample['Depth'])
    a[1,3].plot(sample['H2S'].astype(float),sample['Depth'])
Nelson10=porewater1.set_index('Location').loc['Nelson-10m']
for date in Nelson10['Date'].unique():
    sample=Nelson10.set_index(['Date']).loc[date]
    a[1,1].plot(sample['Sal'].astype(float),sample['Depth'],ls='--')
    a[1,2].plot(sample['DOC'].astype(float),sample['Depth'],ls='--')
    a[1,3].plot(sample['H2S'].astype(float),sample['Depth'],ls='--')

a[1,1].set(title='Nelson salinity',xlabel='Salinity (ppt)',ylabel='Depth (cm)',ylim=(100,0))
a[1,2].set(title='Nelson DOC',xlabel='DOC concentration ($\mu$M)',ylabel='Depth (cm)',ylim=(100,0))
a[1,3].set(title='Nelson sulfide',xlabel='Sulfide concentration ($\mu$M)',ylabel='Depth (cm)',ylim=(100,0))

Shad4=porewater1.set_index('Location').loc['Shad-4m']
for date in Shad4['Date'].unique():
    sample=Shad4.set_index(['Date']).loc[date]
    a[0,1].plot(sample['Sal'].astype(float),sample['Depth'])
    a[0,2].plot(sample['DOC'].astype(float),sample['Depth'])
    a[0,3].plot(sample['H2S'].astype(float),sample['Depth'])
Shad10=porewater1.set_index('Location').loc['Shad-10m']
for date in Shad10['Date'].unique():
    sample=Shad10.set_index(['Date']).loc[date]
    a[0,1].plot(sample['Sal'].astype(float),sample['Depth'],ls='--')
    a[0,2].plot(sample['DOC'].astype(float),sample['Depth'],ls='--')
    a[0,3].plot(sample['H2S'].astype(float),sample['Depth'],ls='--')

a[0,1].set(title='Shad salinity',xlabel='Salinity (ppt)',ylabel='Depth (cm)',ylim=(100,0))
a[0,2].set(title='Shad DOC',xlabel='DOC concentration ($\mu$M)',ylabel='Depth (cm)',ylim=(100,0))
a[0,3].set(title='Shad sulfide',xlabel='Sulfide concentration ($\mu$M)',ylabel='Depth (cm)',ylim=(100,0))

porosity=marshcol['watsat'].T.squeeze()[:10,start:end].to_masked_array()
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

plt.show()
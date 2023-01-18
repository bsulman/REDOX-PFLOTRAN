import xarray
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
import sys
import numpy as np


data=xarray.open_mfdataset(sys.argv[1:])


VWC=data['H2OSOI'].T.squeeze()[:10,:]
O2=data['soil_O2'].T.squeeze()[:10,:]
DOC=data['DOC_vr'].T.squeeze()[:10,:]/12
DIC=data['DIC_vr'].T.squeeze()[:10,:]/12
CH4=data['CH4_vr'].T.squeeze()[:10,:]/12
porosity=data['watsat'].T.squeeze()[:10,:]
z=data['levdcmp'][:10]
t=data['time']
dz=np.array([float(n) for n in '1.7512817916255204E-002   2.7578969259676251E-002   4.5470033242413201E-002   7.4967410986208557E-002  0.12360036510228053       0.20378255101043175       0.33598062644843263       0.55393840536868488       0.91329003158906108        1.5057607013992766        2.4825796969813321        4.0930819526214002        6.7483512780057175        11.126150294204420        13.851152141963599'.split()])[:10]

maxdepth=1.5

inds=VWC.isel(levgrnd=6).load().argsort()
snapshots=[inds[len(inds)//10].item(),inds[len(inds)//10*9].item(),inds[len(inds)//5].item()]
snapshot_styles=['--',':','-']

f,a=plt.subplots(num='Water, oxygen, and carbon',clear=True,nrows=6,ncols=2,gridspec_kw={'width_ratios':[1,0.5]},figsize=(6,8))

VWC.plot(ax=a[0,0],vmax=1,cbar_kwargs={'label':'Volumetric water\n(fraction of saturation)'})
watervol=(porosity).to_masked_array()
(O2/watervol).plot(ax=a[1,0],vmax=(O2/watervol).isel(levdcmp=0).compute().quantile(0.95),cbar_kwargs={'label':'Oxygen concentration\n(mmol/L H$_2$O)'})
# alqO2['soil_O2'].T.sel(time=slice('0096-01-01',None,None)).plot(ax=a[1],vmax=.2,cbar_kwargs={'label':'Oxygen concentration\n(mol/m$^3$)'})
(DOC/watervol/1000).plot(ax=a[2,0],cbar_kwargs={'label':'DOC concentration\n(mol C/L H$_2$O)'})
(DIC/watervol/1000).plot(ax=a[3,0],vmax=(DIC/watervol/1000).compute().quantile(0.95),cbar_kwargs={'label':'DIC concentration\n(mol C/L H$_2$O)'})
(CH4/watervol/1000).plot(ax=a[4,0],vmax=(CH4/watervol/1000).compute().quantile(0.95),cbar_kwargs={'label':'CH4 concentration\n(mol C/L H$_2$O)'})
soilC=(data['SOIL1C_vr']+data['SOIL2C_vr']+data['SOIL3C_vr']+data['SOIL4C_vr']+data['LITR1C_vr']+data['LITR2C_vr']+data['LITR3C_vr']).T.squeeze()[:10,:]
soilC.plot(ax=a[5,0],vmax=soilC.compute().quantile(0.95),cbar_kwargs={'label':'SOC concentration\n(g C/m$^3$)'})

# a[0,1].plot(VWC.mean(dim='time'),z)
# a[1,1].plot((O2/VWC.to_masked_array()).mean(dim='time'),z)
# a[2,1].plot((DOC/VWC.to_masked_array()/1000).mean(dim='time'),z)
# a[3,1].plot((DIC/VWC.to_masked_array()/1000).mean(dim='time'),z)

for num in range(len(snapshots)):
    for ax in a[:,0]:
        ax.axvline(ax.xaxis.convert_units(t[snapshots[num]].item()),ls=snapshot_styles[num],c='gray',lw=2.0)
    a[0,1].plot(VWC.isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')
    a[1,1].plot((O2/watervol).isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')
    a[2,1].plot((DOC/watervol/1000).isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')
    a[3,1].plot((DIC/watervol/1000).isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')
    a[4,1].plot((CH4/watervol/1000).isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')
    a[5,1].plot((soilC).isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')


a[0,0].set(title='Soil water content',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
a[1,0].set(title='Soil O$_2$',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
a[2,0].set(title='Soil DOC',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
a[3,0].set(title='Soil DIC',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
a[4,0].set(title='Soil CH4',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
a[5,0].set(title='Soil organic C',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')

a[0,1].set(title='Soil water content',ylim=(maxdepth,0),xlabel='VWC (frac of sat)',ylabel='Soil depth (m)')
a[1,1].set(title='Soil O$_2$',ylim=(maxdepth,0),xlabel='O$_2$ (mmol/L)',ylabel='Soil depth (m)')
a[2,1].set(title='Soil DOC',ylim=(maxdepth,0),xlabel='DOC (mol C/L)',ylabel='Soil depth (m)')
a[3,1].set(title='Soil DIC',ylim=(maxdepth,0),xlabel='DIC (mol C/L)',ylabel='Soil depth (m)')
a[4,1].set(title='Soil CH4',ylim=(maxdepth,0),xlabel='CH4 (mol C/L)',ylabel='Soil depth (m)')
a[5,1].set(title='Soil organic C',ylim=(maxdepth,0),xlabel='SOC (g C/m$^3$)',ylabel='Soil depth (m)')

Fe2=data['soil_Fe2'].T.squeeze()[:10,:]
FeOxide=data['soil_FeOxide'].T.squeeze()[:10,:]
pH=data['soil_pH'].T.squeeze()[:10,:]

f,a=plt.subplots(num='Fe',clear=True,nrows=3,ncols=2,gridspec_kw={'width_ratios':[1,0.5]},figsize=(6,8))
Fe2.plot(ax=a[0,0],vmax=Fe2.load()[:7,:].quantile(0.95),cbar_kwargs={'label':'Fe(II) concentration (mol Fe/m$^3$)'})
# a[0,1].plot(Fe2.mean(dim='time'),z)
a[0,0].set(title='Fe(II) concentration',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
a[0,1].set(title='Fe(II) concentration',ylim=(maxdepth,0),xlabel='Fe(II) (mol Fe/m$^3$)',ylabel='Soil depth (m)')

(FeOxide).plot(ax=a[1,0],cbar_kwargs={'label':'Fe oxide concentration (mol Fe/m$^3$)'})
# a[1,1].plot((FeOxide).mean(dim='time'),z)
a[1,0].set(title='Fe oxide concentration',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
a[1,1].set(title='Fe oxide concentration',ylim=(maxdepth,0),xlabel='Fe oxide (mol Fe/m$^3$)',ylabel='Soil depth (m)')

pH.plot(ax=a[2,0],vmin=2,cbar_kwargs={'label':'pH'})
# a[2,1].plot(pH.mean(dim='time'),z)
a[2,0].set(title='Soil pH',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
a[2,1].set(title='Soil pH',ylim=(maxdepth,0),xlabel='Soil pH',ylabel='Soil depth (m)')

for num in range(len(snapshots)):
    for ax in a[:,0]:
        ax.axvline(ax.xaxis.convert_units(t[snapshots[num]].item()),ls=snapshot_styles[num],c='gray',lw=2.0)
    a[0,1].plot(Fe2.isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')
    a[1,1].plot(FeOxide.isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')
    a[2,1].plot(pH.isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')

sulfate=data['soil_sulfate'].T.squeeze()[:10,:]
sulfide=data['soil_sulfide'].T.squeeze()[:10,:]
f,a=plt.subplots(num='Sulfur',clear=True,nrows=2,ncols=2,gridspec_kw={'width_ratios':[1,0.5]},figsize=(6,8))
sulfate.plot(ax=a[0,0],vmax=sulfate.load()[:7,:].quantile(0.95),cbar_kwargs={'label':'Sulfate concentration (mol/m$^3$)'})
# a[0,1].plot(Fe2.mean(dim='time'),z)
a[0,0].set(title='Sulfate concentration',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
a[0,1].set(title='Sulfate concentration',ylim=(maxdepth,0),xlabel='Sulfate (mol/m$^3$)',ylabel='Soil depth (m)')

(sulfide).plot(ax=a[1,0],cbar_kwargs={'label':'Sulfide concentration (mol/m$^3$)'})
# a[1,1].plot((FeOxide).mean(dim='time'),z)
a[1,0].set(title='Sulfide concentration',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
a[1,1].set(title='Sulfide concentration',ylim=(maxdepth,0),xlabel='Fe oxide (mol Fe/m$^3$)',ylabel='Soil depth (m)')


for num in range(len(snapshots)):
    for ax in a[:,0]:
        ax.axvline(ax.xaxis.convert_units(t[snapshots[num]].item()),ls=snapshot_styles[num],c='gray',lw=2.0)
    a[0,1].plot(sulfate.isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')
    a[1,1].plot(sulfide.isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')


if 'SIC_vr' in data:
    calcite=data['SIC_vr'].T.squeeze()[:10,:]
    # f,a=plt.subplots(num='Calcite',clear=True,nrows=1,ncols=2,gridspec_kw={'width_ratios':[1,0.5]},figsize=(6,3))
    f,a=plt.subplot_mosaic(
        '''
        AB
        C.
        ''',
        gridspec_kw={'width_ratios':[1,0.5],'height_ratios':[1,0.5]},num='Soil inorganic C',clear=True,
    )
    (calcite).plot(ax=a['A'],cbar_kwargs={'label':'Soil inorganic C concentration (mol C/m$^3$)'})
    a['B'].plot((calcite).mean(dim='time'),z)
    a['A'].set(title='Soil inorganic C concentration',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
    a['B'].set(title='Soil inorganic C concentration',ylim=(maxdepth,0),xlabel='Soil inorganic C (mol C/m$^3$)',ylabel='Soil depth (m)')
    (calcite*dz[:,None]).sum(dim='levdcmp').plot(ax=a['C'])
    a['C'].set(title='Column total soil inorganic C',xlabel='Time (year)',ylabel='Total soil inorganic C (g C/m$^2$)',xlim=a['A'].get_xlim())


f,a=plt.subplots(num='Carbon time series',clear=True,nrows=2)
data['TOTSOMC'].plot(ax=a[0],label='Total SOM C')
data['TOTLITC'].plot(ax=a[0],label='Total litter C')
data['TOTVEGC'].plot(ax=a[0],label='Total vegetation C')
a[0].set(title='C pools',xlabel='Time (year)',ylabel='C stock (g C m$^{-2}$')
a[0].legend()

data['HR'].plot(ax=a[1],label='HR')
# Should add a smoothed curve
data['HR'].resample(time='7D').mean().plot(ax=a[1])
data['GPP'].plot(ax=a[1],label='GPP')
a[1].set(title='C fluxes (+ is to atmosphere)',xlabel='Time (year)',ylabel='C flux (g C m$^{-2}$ s$^{-1}$)')

f,a=plt.subplots(num='Time step',clear=True)
data['chem_dt'].plot(ax=a)
a.axhline(data['chem_dt'].mean(),c='k',ls='--')

plt.show()
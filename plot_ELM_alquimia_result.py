import xarray
import matplotlib.pyplot as plt
import matplotlib.dates
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



inds=VWC.isel(levgrnd=6).load().argsort()
snapshots=[inds[len(inds)//10].item(),inds[len(inds)//10*9].item(),inds[len(inds)//5].item()]
snapshot_styles=['--',':','-']

watervol=(porosity).to_masked_array()

def plot_vars(plotname,vars,figsize=(4,4),maxdepth=1.5,vmax={}):
    f,a=plt.subplots(num=plotname,clear=True,nrows=len(vars),ncols=2,gridspec_kw={'width_ratios':[1,0.5]},figsize=figsize,squeeze=False)
    for varnum,var in enumerate(vars):
        if var=='VWC':
            (VWC/porosity).plot(ax=a[varnum,0],vmax=1,cbar_kwargs={'label':'Volumetric water\n(fraction of saturation)'},cmap='Blues')
            for num in range(len(snapshots)):
                a[varnum,0].axvline(a[varnum,0].xaxis.convert_units(t[snapshots[num]].item()),ls=snapshot_styles[num],c='gray',lw=2.0)
                a[varnum,1].plot(VWC.isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')
            a[varnum,0].set(title='Soil water content',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
            a[varnum,1].set(title='Soil water content',ylim=(maxdepth,0),xlabel='VWC (frac of sat)',ylabel='Soil depth (m)')
        if var=='frozen':
            frozenfrac= (data['SOILICE']/(data['SOILLIQ']+data['SOILICE'])).T.squeeze()[:10,:]
            frozenfrac.plot(ax=a[varnum,0],vmax=1,cbar_kwargs={'label':'Frozen water fraction'},cmap='Blues_r')
            for num in range(len(snapshots)):
                a[varnum,0].axvline(a[varnum,0].xaxis.convert_units(t[snapshots[num]].item()),ls=snapshot_styles[num],c='gray',lw=2.0)
                a[varnum,1].plot(frozenfrac.isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')
            a[varnum,0].set(title='Frozen water fraction',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
            a[varnum,1].set(title='Frozen water fraction',ylim=(maxdepth,0),xlabel='Frozen water fraction',ylabel='Soil depth (m)')
        elif var=='O2':
            if 'O2' in vmax:
                vmax_O2=vmax['O2']
            else:
                vmax_O2=(O2/watervol).isel(levdcmp=0).compute().quantile(0.95)
            (O2/watervol).plot(ax=a[varnum,0],vmax=vmax_O2,cbar_kwargs={'label':'Oxygen concentration\n(mmol/L H$_2$O)'})
            for num in range(len(snapshots)):
                a[varnum,0].axvline(a[varnum,0].xaxis.convert_units(t[snapshots[num]].item()),ls=snapshot_styles[num],c='gray',lw=2.0)
                a[varnum,1].plot((O2/watervol).isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')
            a[varnum,0].set(title='Soil O$_2$',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
            a[varnum,1].set(title='Soil O$_2$',ylim=(maxdepth,0),xlabel='O$_2$ (mmol/L)',ylabel='Soil depth (m)',xlim=(0,vmax_O2))
        elif var=='DOC':
            if 'DOC' in vmax:
                vmax_DOC=vmax['DOC']
            else:
                vmax_DOC=(DOC/watervol).max().compute()
            (DOC/watervol).plot(ax=a[varnum,0],vmax=vmax_DOC,cbar_kwargs={'label':'DOC concentration\n(mmol C/L H$_2$O)'})
            for num in range(len(snapshots)):
                a[varnum,0].axvline(a[varnum,0].xaxis.convert_units(t[snapshots[num]].item()),ls=snapshot_styles[num],c='gray',lw=2.0)
                a[varnum,1].plot((DOC/watervol).isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')
            a[varnum,0].set(title='Soil DOC',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
            a[varnum,1].set(title='Soil DOC',ylim=(maxdepth,0),xlabel='DOC (mmol C/L)',ylabel='Soil depth (m)',xlim=(0,vmax_DOC))
        elif var=='DIC':
            if 'DIC' in vmax:
                vmax_DIC=vmax['DIC']
            else:
                vmax_DIC=(DOC/watervol).max().compute()
            (DIC/watervol).plot(ax=a[varnum,0],cbar_kwargs={'label':'DIC concentration\n(mmol C/L H$_2$O)'},vmax=vmax_DIC)
            for num in range(len(snapshots)):
                a[varnum,0].axvline(a[varnum,0].xaxis.convert_units(t[snapshots[num]].item()),ls=snapshot_styles[num],c='gray',lw=2.0)
                a[varnum,1].plot((DOC/watervol).isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')
            a[varnum,0].set(title='Soil DIC',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
            a[varnum,1].set(title='Soil DIC',ylim=(maxdepth,0),xlabel='DIC (mmol C/L)',ylabel='Soil depth (m)',xlim=(0,vmax_DIC))       
        elif var=='CH4':
            if 'CH4' in vmax:
                vmax_CH4=vmax['CH4']
            else:
                vmax_CH4=(CH4/watervol*1000).max().compute()
            (CH4/watervol*1000).plot(ax=a[varnum,0],cbar_kwargs={'label':'CH4 concentration\n($\mu$mol C/L H$_2$O)'},vmax=vmax_CH4)
            for num in range(len(snapshots)):
                a[varnum,0].axvline(a[varnum,0].xaxis.convert_units(t[snapshots[num]].item()),ls=snapshot_styles[num],c='gray',lw=2.0)
                a[varnum,1].plot((CH4/watervol*1000).isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')
            a[varnum,0].set(title='Soil CH4',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
            a[varnum,1].set(title='Soil CH4',ylim=(maxdepth,0),xlabel='CH$_4$ ($\mu$mol C/L)',ylabel='Soil depth (m)',xlim=(0,vmax_CH4))        
        elif var=='soilC':
            soilC=(data['SOIL1C_vr']+data['SOIL2C_vr']+data['SOIL3C_vr']+data['SOIL4C_vr']+data['LITR1C_vr']+data['LITR2C_vr']+data['LITR3C_vr']).T.squeeze()[:10,:]/1000
            soilC.plot(ax=a[varnum,0],vmax=soilC.compute().quantile(0.95),cbar_kwargs={'label':'SOC concentration\n(kg C/m$^3$)'})
            for num in range(len(snapshots)):
                a[varnum,0].axvline(a[varnum,1].xaxis.convert_units(t[snapshots[num]].item()),ls=snapshot_styles[num],c='gray',lw=2.0)
                a[varnum,1].plot((soilC).isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')
            a[varnum,0].set(title='Soil organic C',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
            a[varnum,1].set(title='Soil organic C',ylim=(maxdepth,0),xlabel='SOC (g C/m$^3$)',ylabel='Soil depth (m)')
        elif var=='Fe2':
            Fe2=data['soil_Fe2'].T.squeeze()[:10,:]/watervol*1000
            if 'Fe2' in vmax:
                vmax_Fe2=vmax['Fe2']
            else:
                vmax_Fe2=Fe2.load()[:7,:].quantile(0.95)
            Fe2.plot(ax=a[varnum,0],vmax=vmax_Fe2,cbar_kwargs={'label':'Fe(II) concentration\n($\mu$mol Fe/L H$_2$O)'})
            a[varnum,0].set(title='Fe(II) concentration',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
            a[varnum,1].set(title='Fe(II) concentration',ylim=(maxdepth,0),xlabel='Fe(II) ($\mu$mol Fe/L H$_2$O)',ylabel='Soil depth (m)',xlim=(0,vmax_Fe2))
            for num in range(len(snapshots)):
                a[varnum,0].axvline(a[varnum,0].xaxis.convert_units(t[snapshots[num]].item()),ls=snapshot_styles[num],c='gray',lw=2.0)
                a[varnum,1].plot(Fe2.isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')
        elif var=='FeOxide':
            FeOxide=data['soil_FeOxide'].T.squeeze()[:10,:]
            if 'FeOxide' in vmax:
                vmax_FeOxide=vmax['FeOxide']
            else:
                vmax_FeOxide=FeOxide.max().compute()
            (FeOxide).plot(ax=a[varnum,0],cbar_kwargs={'label':'Fe oxide concentration\n(mol Fe/m$^3$)'},vmax=vmax_FeOxide)
            a[varnum,0].set(title='Fe oxide concentration',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
            a[varnum,1].set(title='Fe oxide concentration',ylim=(maxdepth,0),xlabel='Fe oxide (mol Fe/m$^3$)',ylabel='Soil depth (m)',xlim=(FeOxide.min().compute(),vmax_FeOxide))
            for num in range(len(snapshots)):
                a[varnum,0].axvline(a[varnum,0].xaxis.convert_units(t[snapshots[num]].item()),ls=snapshot_styles[num],c='gray',lw=2.0)
                a[varnum,1].plot(FeOxide.isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')
        elif var=='pH':
            pH=data['soil_pH'].T.squeeze()[:10,:]
            if 'pH' in vmax:
                vmax_pH=vmax['pH']
            else:
                vmax_pH=pH.max().compute()
            pH.plot(ax=a[varnum,0],vmin=3.5,vmax=vmax_pH,cbar_kwargs={'label':'pH'})
            a[varnum,0].set(title='Soil pH',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
            a[varnum,1].set(title='Soil pH',ylim=(maxdepth,0),xlabel='Soil pH',ylabel='Soil depth (m)',xlim=(3.5,vmax_pH))
            for num in range(len(snapshots)):
                a[varnum,1].plot(pH.isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')
                a[varnum,0].axvline(a[varnum,0].xaxis.convert_units(t[snapshots[num]].item()),ls=snapshot_styles[num],c='gray',lw=2.0)
        elif var=='Sulfate':
            sulfate=data['soil_sulfate'].T.squeeze()[:10,:]/watervol*1000
            if 'sulfate' in vmax:
                vmax_sulfate=vmax['sulfate']
            else:
                vmax_sulfate=sulfate.load()[:7,:].quantile(0.95)
            sulfate.plot(ax=a[varnum,0],vmax=vmax_sulfate,cbar_kwargs={'label':'Sulfate concentration\n($\mu$mol/L H$_2$O)'})
            a[varnum,0].set(title='Sulfate concentration',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
            a[varnum,1].set(title='Sulfate concentration',ylim=(maxdepth,0),xlabel='Sulfate ($\mu$mol/L H$_2$O)',ylabel='Soil depth (m)',xlim=(0,vmax_sulfate))
            for num in range(len(snapshots)):
                a[varnum,0].axvline(a[varnum,0].xaxis.convert_units(t[snapshots[num]].item()),ls=snapshot_styles[num],c='gray',lw=2.0)
                a[varnum,1].plot(sulfate.isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')
        elif var=='Sulfide':
            sulfide=data['soil_sulfide'].T.squeeze()[:10,:]/watervol*1000
            if 'sulfide' in vmax:
                vmax_sulfide=vmax['sulfide']
            else:
                vmax_sulfide=sulfide.load()[:7,:].quantile(0.95)            
            (sulfide).plot(ax=a[varnum,0],cbar_kwargs={'label':'Sulfide concentration\n(mol/m$^3$)'},vmax=vmax_sulfide)
            a[varnum,0].set(title='Sulfide concentration',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
            a[varnum,1].set(title='Sulfide concentration',ylim=(maxdepth,0),xlabel='Sulfide ($\mu$mol/L H$_2$O)',ylabel='Soil depth (m)',xlim=(0,vmax_sulfide))
            for num in range(len(snapshots)):
                a[varnum,1].plot(sulfide.isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')
                a[varnum,0].axvline(a[varnum,0].xaxis.convert_units(t[snapshots[num]].item()),ls=snapshot_styles[num],c='gray',lw=2.0)
        elif var=='vertflow':
            vertflow=data['QFLX_ADV'].squeeze().T[:10,:]*3600
            vmax_vertflow=vmax.get('vertflow',0.005*3600)
            vertflow.plot(vmax=vmax_vertflow,ax=a[varnum,0],cbar_kwargs={'label':'Vertical flow rate (mm/hour)'})
            a[varnum,0].set(title='Vertical water flow',ylim=(maxdepth,0),xlabel='Time (year)',ylabel='Soil depth (m)')
            a[varnum,1].set(title='Vertical water flow',ylim=(maxdepth,0),xlabel='Water flow (mm/hour)',ylabel='Soil depth (m)',xlim=(0,vmax_vertflow))
            for num in range(len(snapshots)):
                a[varnum,1].plot(vertflow.isel(time=snapshots[num]),z,ls=snapshot_styles[num],c='gray')
                a[varnum,0].axvline(a[varnum,0].xaxis.convert_units(t[snapshots[num]].item()),ls=snapshot_styles[num],c='gray',lw=2.0)

        a[varnum,0].xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b-%d'))

plot_vars('Water and oxygen',['VWC','O2','frozen','vertflow'],figsize=(6,8.5),maxdepth=2.2,vmax={'vertflow':1e-1})
plot_vars('Carbon',['soilC','DOC','DIC','CH4'],figsize=(6,8.5),vmax={'DOC':10.0,'DIC':10.0},maxdepth=2.2)
plot_vars('Redox',['Fe2','FeOxide','Sulfate','Sulfide','pH'],figsize=(6,8.5),maxdepth=2.2,vmax={'pH':9.0})



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


f,a=plt.subplots(num='Carbon time series',clear=True,nrows=3)
# (data['TOTSOMC']/1000).plot(ax=a[0],label='Total SOM C')
(data['TOTLITC']/1000).plot(ax=a[0],label='Total litter C')
(data['TOTVEGC']/1000).plot(ax=a[0],label='Total vegetation C')
a[0].set(title='C pools',xlabel='Time (year)',ylabel='C stock (kg C m$^{-2}$')
a[0].legend()

# (data['HR']*3600*24).plot(ax=a[1],label='HR')
# Should add a smoothed curve
(data['HR']*3600*24).resample(time='7D').mean().plot(ax=a[1],label='HR')
(data['GPP']*3600*24).resample(time='7D').mean().plot(ax=a[1],label='GPP')
(data['NEE']*3600*24).resample(time='7D').mean().plot(ax=a[1],label='NEE')
a[1].set(title='C fluxes (+ is to atmosphere)',xlabel='Time (year)',ylabel='C flux (g C m$^{-2}$ day$^{-1}$)')
a[1].legend()

(data['CH4FLUX_ALQUIMIA']*3600*24).plot(ax=a[2])
a[2].set(title='CH$_4$ flux',xlabel='Time (year)',ylabel='C flux (g C m$^{-2}$ day$^{-1}$)')

f,a=plt.subplots(num='Time step',clear=True)
data['chem_dt'].plot(ax=a)
a.axhline(data['chem_dt'].mean(),c='k',ls='--')

plt.show()
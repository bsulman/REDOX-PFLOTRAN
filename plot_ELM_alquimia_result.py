import xarray
import matplotlib.pyplot as plt
import matplotlib.dates
import sys
import numpy as np
import nc_time_axis,cftime


dz=np.array([float(n) for n in '1.7512817916255204E-002   2.7578969259676251E-002   4.5470033242413201E-002   7.4967410986208557E-002  0.12360036510228053       0.20378255101043175       0.33598062644843263       0.55393840536868488       0.91329003158906108        1.5057607013992766        2.4825796969813321        4.0930819526214002        6.7483512780057175        11.126150294204420        13.851152141963599'.split()])[:10]
snapshot_styles=['--',':','-']



def plot_var(vardata,contourax=None,profileax=None,vmax=None,vmin=None,label=None,cmap=None,
            snapshots=[],mean_profile=False,profile_color='gray',quantiles=[],
            title=None,axlabel=None,maxdepth=None):
    if contourax is not None:
        vardata.plot(ax=contourax,vmax=vmax,vmin=vmin,cbar_kwargs={'label':label},cmap=cmap)
        contourax.set(title=title,ylim=(maxdepth,0),xlabel='Time',ylabel='Soil depth (m)')
    if profileax is not None:
        if 'levdcmp' in vardata.dims:
            z=vardata['levdcmp']
        else:
            z=vardata['levgrnd']
        for num in range(len(snapshots)):
            if contourax is not None:
                contourax.axvline(contourax.xaxis.convert_units(vardata['time'][snapshots[num]].item()),ls=snapshot_styles[num],c='gray',lw=2.0)
            profileax.plot(vardata.isel(time=snapshots[num]),z,ls=snapshot_styles[num],c=profile_color)
        if mean_profile:
            profileax.plot(vardata.mean(dim='time'),z,ls='-',c=profile_color,lw=1.5,label=label)
        if len(quantiles) == 2:
            profileax.fill_betweenx(z,vardata.load().quantile(quantiles[0],dim='time'),vardata.load().quantile(quantiles[1],dim='time'),fc=profile_color,alpha=0.3)
        profileax.set(title=title,ylim=(maxdepth,0),xlabel=axlabel,ylabel='Soil depth (m)')
        profileax.set_xlim(left=vmin,right=vmax)


def plot_vars(data,vars,plotname=None,figsize=(4,4),maxdepth=1.5,vmax={},vmin={},a_contour=None,a_profile=None,profile_color='gray',mean_profile=True,snapshots=[],quantiles=[0.1,0.9],WT_thresh=1.0,time_format='%b-%d'):
    if a_contour is None and a_profile is None:
        f,a=plt.subplots(num=plotname,clear=True,nrows=len(vars),ncols=2,gridspec_kw={'width_ratios':[1,0.5]},figsize=figsize,squeeze=False)
        a_contour=a[:,0]
        a_profile=a[:,1]
    inds=data['H2OSOI'].T.squeeze()[:10,:].isel(levgrnd=6).load().argsort()
    porosity=data['watsat'].T.squeeze()[:10,:]
    watervol=(porosity).to_masked_array()
    contourax=None
    profileax=None
    for varnum,var in enumerate(vars):
        if a_contour is not None:
            contourax = np.atleast_1d(a_contour)[varnum]
        if a_profile is not None:
            profileax = np.atleast_1d(a_profile)[varnum]
        if var=='VWC' or var=='SWC':
            plot_var((data['H2OSOI'].T.squeeze()[:10,:]/porosity),contourax,profileax,vmax=1.0,vmin=0.0,
                    label='Volumetric water\n(fraction of saturation)',cmap='Blues',
                    snapshots=snapshots,profile_color=profile_color,mean_profile=mean_profile,quantiles=quantiles,
                    title='Soil water content',axlabel='VWC (frac of sat)',maxdepth=maxdepth)
        if var=='frozen':
            if 'SOILICE' in data:
                frozenfrac= (data['SOILICE']/(data['SOILLIQ']+data['SOILICE'])).T.squeeze()[:10,:]
                plot_var(frozenfrac,contourax,profileax,vmax=1.0,vmin=0.0,label='Frozen water fraction',cmap='Blues_r',
                         snapshots=snapshots,profile_color=profile_color,mean_profile=mean_profile,quantiles=quantiles,
                         title='Frozen water fraction',axlabel='Frozen water fraction',maxdepth=maxdepth)
        elif var=='O2' or var=='oxygen':
            if 'O2' in vmax:
                vmax_O2=vmax['O2']
            else:
                vmax_O2=(data['soil_O2'].T.squeeze()[:10,:]/watervol).isel(levdcmp=0).compute().quantile(0.95)
            plot_var(data['soil_O2'].T.squeeze()[:10,:]/watervol,contourax,profileax,vmax=vmax_O2,
                     label='Oxygen concentration\n(mmol/L H$_2$O)',snapshots=snapshots,profile_color=profile_color,
                     mean_profile=mean_profile,quantiles=quantiles,
                     maxdepth=maxdepth,title='Soil O$_2$',axlabel='O$_2$ (mmol/L)')
        elif var=='DOC':
            if 'DOC' in vmax:
                vmax_DOC=vmax['DOC']
            else:
                vmax_DOC=(data['DOC_vr'].T.squeeze()[:10,:]/watervol/12.011).max().compute()
            plot_var(data['DOC_vr'].T.squeeze()[:10,:]/watervol/12.011,contourax,profileax,vmax=vmax_DOC,
                     label='DOC concentration\n(mmol/L H$_2$O)',snapshots=snapshots,profile_color=profile_color,mean_profile=mean_profile,
                     quantiles=quantiles,maxdepth=maxdepth,title='Soil DOC',axlabel='DOC (mmol C/L)')
        elif var=='DIC':
            if 'DIC' in vmax:
                vmax_DIC=vmax['DIC']
            else:
                vmax_DIC=(data['DIC_vr'].T.squeeze()[:10,:]/watervol/12.011).max().compute()
            plot_var(data['DIC_vr'].T.squeeze()[:10,:]/watervol/12.011,contourax,profileax,vmax=vmax_DIC,label='DIC concentration\n(mmol/L H$_2$O)',
                     snapshots=snapshots,profile_color=profile_color,mean_profile=mean_profile,quantiles=quantiles,
                     maxdepth=maxdepth,title='Soil DIC',axlabel='DIC (mmol C/L)')
        elif var=='acetate':
            if 'acetate' in vmax:
                vmax_acetate=vmax['acetate']
            else:
                vmax_acetate=(data['soil_acetate'].T.squeeze()[:10,:]/watervol).max().compute()
            plot_var(data['soil_acetate'].T.squeeze()[:10,:]/watervol,contourax,profileax,vmax=vmax_acetate,label='Acetate concentration\n(mmol/L H$_2$O)',
                     snapshots=snapshots,profile_color=profile_color,mean_profile=mean_profile,quantiles=quantiles,
                     maxdepth=maxdepth,title='Soil Acetate',axlabel='Acetate (mmol/L)')
        elif var=='CH4':
            if 'CH4' in vmax:
                vmax_CH4=vmax['CH4']
            else:
                vmax_CH4=(data['CH4_vr'].T.squeeze()[:10,:]/watervol/12.011).max().compute()
            plot_var(data['CH4_vr'].T.squeeze()[:10,:]/watervol/12.011,contourax,profileax,vmax=vmax_CH4,vmin=vmin.get('CH4',0.0),label='CH$_4$ concentration\n(mmol C/L H$_2$O)',
                     snapshots=snapshots,profile_color=profile_color,mean_profile=mean_profile,quantiles=quantiles,
                     maxdepth=maxdepth,title='Soil CH$_4$',axlabel='CH$_4$ (mmol C/L)')      
        elif var=='soilC':
            soilC=(data['SOIL1C_vr']+data['SOIL2C_vr']+data['SOIL3C_vr']+data['SOIL4C_vr']+data['LITR1C_vr']+data['LITR2C_vr']+data['LITR3C_vr']).T.squeeze()[:10,:]/100**3
            plot_var(soilC,contourax,profileax,vmax=None,label='SOC concentration\n(kg C/m$^3$)',snapshots=snapshots,profile_color=profile_color,
                     mean_profile=mean_profile,quantiles=quantiles,
                     maxdepth=maxdepth,title='Soil organic C',axlabel='SOC (g C/cm$^3$)')
        elif var=='Fe2':
            Fe2=data['soil_Fe2'].T.squeeze()[:10,:]/watervol
            if 'Fe2' in vmax:
                vmax_Fe2=vmax['Fe2']
            else:
                vmax_Fe2=Fe2.load()[:7,:].quantile(0.95)
            plot_var(Fe2,contourax,profileax,vmax=vmax_Fe2,label='Fe(II) concentration\n(mmol Fe/L H$_2$O)',vmin=vmin.get('Fe2',0.0),
                     snapshots=snapshots,profile_color=profile_color,mean_profile=mean_profile,quantiles=quantiles,
                     maxdepth=maxdepth,title='Fe(II) concentration',axlabel='Fe(II) (mmol Fe/L H$_2$O)')
        elif var=='FeOxide':
            FeOxide=data['soil_FeOxide'].T.squeeze()[:10,:]
            if 'FeOxide' in vmax:
                vmax_FeOxide=vmax['FeOxide']
            else:
                vmax_FeOxide=FeOxide.max().compute()
            plot_var(FeOxide,contourax,profileax,vmax=vmax_FeOxide,vmin=vmin.get('FeOxide',FeOxide.min().compute()),label='Fe oxide concentration\n(mol Fe/m$^3$)',
                     snapshots=snapshots,profile_color=profile_color,mean_profile=mean_profile,quantiles=quantiles,
                     maxdepth=maxdepth,title='Fe oxide concentration',axlabel='Fe oxide (mol Fe/m$^3$)')
        elif var=='FeS':
            plot_var(data['soil_FeS'].T.squeeze()[:10,:],contourax,profileax,vmax=vmax.get('FeS',data['soil_FeS'].T.squeeze()[:10,:].max().compute()),
                     vmin=vmin.get('FeS',data['soil_FeS'].T.squeeze()[:10,:].min().compute()),
                     label='Fe sulfide concentration\n(mol Fe/m$^3$)',
                     snapshots=snapshots,profile_color=profile_color,mean_profile=mean_profile,quantiles=quantiles,
                     maxdepth=maxdepth,title='Fe sulfide concentration',axlabel='Fe sulfide (mol Fe/m$^3$)')
        elif var=='pH':
            pH=data['soil_pH'].T.squeeze()[:10,:]
            if 'pH' in vmax:
                vmax_pH=vmax['pH']
            else:
                vmax_pH=pH.max().compute()
            plot_var(pH,contourax,profileax,vmin=vmin.get('pH',pH.min().compute()),vmax=vmax_pH,label='pH',snapshots=snapshots,
                     profile_color=profile_color,mean_profile=mean_profile,quantiles=quantiles,
                     maxdepth=maxdepth,title='Soil pH',axlabel='Soil pH')
        elif var=='Sulfate' or var=='sulfate':
            sulfate=data['soil_sulfate'].T.squeeze()[:10,:]/watervol
            if 'sulfate' in vmax:
                vmax_sulfate=vmax['sulfate']
            else:
                vmax_sulfate=sulfate.load()[:7,:].quantile(0.95)
            plot_var(sulfate,contourax,profileax,vmax=vmax_sulfate,vmin=vmin.get('sulfate',0.0),label='Sulfate concentration\n(mmol/L H$_2$O)',
                     snapshots=snapshots,profile_color=profile_color,mean_profile=mean_profile,quantiles=quantiles,
                     maxdepth=maxdepth,title='Sulfate concentration',axlabel='Sulfate (mmol/L H$_2$O)')
        elif var=='Sulfide' or var=='sulfide':
            sulfide=data['soil_sulfide'].T.squeeze()[:10,:]/watervol
            if 'sulfide' in vmax:
                vmax_sulfide=vmax['sulfide']
            else:
                vmax_sulfide=sulfide.load()[:7,:].quantile(0.95)            
            plot_var(sulfide,contourax,profileax,vmax=vmax_sulfide,label='Sulfide concentration\n(mmol/L H$_2$O)',
                     snapshots=snapshots,profile_color=profile_color,mean_profile=mean_profile,quantiles=quantiles,
                     maxdepth=maxdepth,title='Sulfide concentration',axlabel='Sulfide (mmol/L H$_2$O)',)
        elif var=='vertflow':
            vertflow=data['QFLX_ADV'].squeeze().T[:10,:]*3600
            vmax_vertflow=vmax.get('vertflow',0.005*3600)
            plot_var(vertflow,contourax,profileax,vmax=vmax_vertflow,label='Vertical flow rate (mm/hour)',
                     snapshots=snapshots,profile_color=profile_color,mean_profile=mean_profile,quantiles=quantiles,
                     maxdepth=maxdepth,title='Vertical water flow',axlabel='Water flow (mm/hour)')
        elif var=='salinity':
            plot_var(data['soil_salinity'].squeeze().T[:10,:],contourax,profileax,vmax=vmax.get('salinity',50),vmin=vmin.get('salinity',None),label='Salinity (ppt)',
                     snapshots=snapshots,profile_color=profile_color,mean_profile=mean_profile,quantiles=quantiles,
                     maxdepth=maxdepth,title='Soil salinity',axlabel='Salinity (ppt)')
        elif var=='temperature':
            plot_var((data['TSOI'].squeeze().T[:10,:]-273.15),contourax,profileax,label='Temperature (C)',
                     snapshots=snapshots,profile_color=profile_color,mean_profile=mean_profile,quantiles=quantiles,
                     maxdepth=maxdepth,title='Soil temperature',axlabel='Temperature (C)')
        if contourax is not None:
            if var=='NPP' or var=='GPP' or var=='NEE':
                (data[var].squeeze()/12e-6).plot(ax=contourax,c=profile_color)
                contourax.set(title=var,ylabel=var+'\n($\mu$ mol m$^{-2}$ s$^{-1}$)',xlabel='Time (year)')
                if profileax is not None:
                    profileax.set_visible(False)
            elif var=='TOTVEGC':
                (data[var].squeeze()).plot(ax=contourax,c=profile_color)
                contourax.set(title='Total vegetation biomass C',ylabel='Biomass C\n(g C m$^{-2}$ s$^{-1}$)',xlabel='Time (year)')
                if profileax is not None:
                    profileax.set_visible(False)
            elif var=='H2OSFC':
                H2OSFC=data['H2OSFC'].squeeze()*1e-3
                # ZWT=data['ZWT'].squeeze()
                VWC=(data['H2OSOI']/data['watsat']).T.squeeze()[:10,:]
                lev1=(VWC>=WT_thresh).argmax(dim='levgrnd').compute()
                lev2=lev1.where(lev1==0,lev1-1).compute()
                ZWT=(data['levgrnd'][lev2]+(data['levgrnd'][lev1]-data['levgrnd'][lev2])/(VWC[lev1]-VWC[lev2])*(WT_thresh-VWC[lev2])).where(lev1>0,0.0)

                (H2OSFC.where(H2OSFC>0,-ZWT)).plot(ax=contourax,c=profile_color)
                contourax.set(title='Water level',ylabel='Water level (m)',xlabel='Time (year)')
                contourax.axhline(0.0,c='k',lw=0.5,ls=':')
                if profileax is not None:
                    profileax.set_visible(False)
            elif var=='H2OSFC_TIDE':
                if 'H2OSFC_TIDE' in data:
                    (data[var].squeeze()*1e-3).plot(ax=contourax)
                    contourax.set(title='Tide water level',ylabel='Water level (m)',xlabel='Time (year)')
                    if profileax is not None:
                        profileax.set_visible(False)
                    contourax.axhline(0.0,c='k',lw=0.5,ls=':')
            elif var=='latflow':
                if 'QFLX_LAT_AQU' in data:
                    (data['QFLX_LAT_AQU'].squeeze()).plot(ax=contourax)
                    contourax.set(title='Lateral water flow into column',xlabel='Time',ylabel='Water flow rate\n(mm s$^{-1}$)')
                    if profileax is not None:
                        profileax.set_visible(False)
                    contourax.axhline(0.0,c='k',lw=0.5,ls=':')
            elif var=='drain':
                (data['QDRAI'].squeeze()).plot(ax=contourax)
                contourax.set(title='Drainage flow',xlabel='Time',ylabel='Water flow rate\n(mm s$^{-1}$)')
                if profileax is not None:
                    profileax.set_visible(False)
                contourax.axhline(0.0,c='k',lw=0.5,ls=':')
            elif var=='CH4flux':
                (data['CH4FLUX_ALQUIMIA'].squeeze()/12.011*1e6).plot(ax=contourax,c=profile_color)
                if profileax is not None:
                    profileax.set_visible(False)
                contourax.set(title='Methane flux',xlabel='Time',ylabel='Methane flux ($\mu$mol m$^{-2}$ s$^{-1}$)')
        
        if((data['time'][-1].item()-data['time'][0].item()).days<365 and contourax is not None):
            contourax.xaxis.set_major_formatter(format_nc_time(time_format))

# Fix an issue with dates from model netCDF output not formatting correctly
class format_nc_time(matplotlib.dates.ticker.Formatter):
    def __init__(self,fmt,calendar='noleap'):
        self.calendar = calendar
        self.fmt=fmt
    def __call__(self, x, pos=0):
        return cftime.num2date(x, nc_time_axis._TIME_UNITS,calendar=self.calendar).strftime(self.fmt)

def letter_label(a):
    from string import ascii_lowercase
    num=0
    for ax in a.ravel():
        if ax.get_visible():
            label='('+ascii_lowercase[num]+')'
            ax.set_title(label,loc='left')
            num = num + 1

if __name__ == '__main__':
    if not sys.argv[-1].endswith('.nc'):
        gridnum=int(sys.argv[-1])
        files=sys.argv[1:-1]
    else:
        gridnum=0
        files=sys.argv[1:]

    data=xarray.open_mfdataset(files).isel(lndgrid=gridnum)

    plot_vars(data,plotname='Water and oxygen',vars=['VWC','O2','frozen','temperature'],figsize=(6,8.5),maxdepth=2.2,vmax={'vertflow':1e-1})
    plot_vars(data,plotname='Carbon',vars=['soilC','DOC','DIC','CH4','acetate'],figsize=(6,8.5),vmax={'DOC':10.0,'DIC':10.0,'CH4':2.5},maxdepth=2.2)
    plot_vars(data,plotname='Redox',vars=['Fe2','FeOxide','Sulfate','Sulfide','pH'],figsize=(6,8.5),maxdepth=2.2,vmax={'pH':9.0})
    if 'QFLX_LAT_AQU' in data:
        plot_vars(data,plotname='Hydro',vars=['H2OSFC','salinity','drain','latflow'])
        if 'H2OSFC_TIDE' in data:
            (data['H2OSFC_TIDE']/1000).plot(ax=plt.figure('Hydro').axes[0],c='k',ls='--')

    f,a=plt.subplots(num='Carbon time series',clear=True,nrows=3)
    # (data['TOTSOMC']/1000).plot(ax=a[0],label='Total SOM C')
    (data['TOTLITC']/1000).plot(ax=a[0],label='Total litter C')
    (data['TOTVEGC']/1000).plot(ax=a[0],label='Total vegetation C')
    a[0].set(title='C pools',xlabel='Time (year)',ylabel='C stock (kg C m$^{-2}$')
    a[0].legend()
    a[0].set_ylim(bottom=0.0)

    # (data['HR']*3600*24).plot(ax=a[1],label='HR')
    # Should add a smoothed curve
    (data['HR']*3600*24).resample(time='7D').mean().plot(ax=a[1],label='HR')
    (data['GPP']*3600*24).resample(time='7D').mean().plot(ax=a[1],label='GPP')
    (data['NEE']*3600*24).resample(time='7D').mean().plot(ax=a[1],label='NEE')
    a[1].set(title='C fluxes',xlabel='Time (year)',ylabel='C flux (g C m$^{-2}$ day$^{-1}$)')
    a[1].legend()
    a[1].axhline(0.0,c='k',ls=':',lw=0.5)

    (data['CH4FLUX_ALQUIMIA']*3600*24).plot(ax=a[2])
    a[2].set(title='CH$_4$ flux',xlabel='Time (year)',ylabel='C flux (g C m$^{-2}$ day$^{-1}$)')
    a[2].axhline(0.0,c='k',ls=':',lw=0.5)

    f,a=plt.subplots(num='Time step',clear=True)
    data['chem_dt'].plot(ax=a)
    a.axhline(data['chem_dt'].mean(),c='k',ls='--')

    plt.show()
from pylab import *
import Arctic_Fe_CH4 as arctic
import decomp_network
import xarray
import sys
import pandas
from string import ascii_lowercase


def plot_result(result,SOM_ax=None,pH_ax=None,Fe_ax=None,CO2flux_ax=None,CH4flux_ax=None,porewater_ax=None,cation_ax=None,
                SOM_args={},pH_args={},Fe_args={},CO2flux_args={},CH4flux_args={},porewater_args={},cation_args={},ls='-',
                do_legend=False,gdrywt=False,BD=None,SOC_pct=None,cellulose_SOC_frac=1.0,obs=None,obs_marker=None,color='C0'):
    t=result['time']
    if obs is not None:
        obs_mod=pandas.DataFrame(index=obs.Incubation_Time)
        obs_mod_pts=array([abs(t-obst).argmin() for obst in obs.Incubation_Time])
        obs_mod_pts[obs_mod_pts==0]=1
    else:
        obs_mod=None
    if gdrywt:
        if SOC_pct is None and BD is None:
            raise TypeError('SOC_pct or BD must be a number if gdrywt is True')
        SOC_mol_m3=result['Total Sorbed cellulose'].isel(time=0)/cellulose_SOC_frac # mol SOC/m3
        SOC_gC_cm3=SOC_mol_m3*12/100**3
        SOC_gC_gdwt=SOC_pct/100
        cm3_to_dwt=SOC_gC_gdwt/SOC_gC_cm3 # Conversion from /cm3 to /gdwt. This is 1/ bulk density in g/cm3
        if BD is not None:
            cm3_to_dwt=1.0/BD
        # print('cm3_to_dwt = %1.3f'%cm3_to_dwt)
    if SOM_ax is not None:
        if gdrywt:
            # Plot in %, i.e. gC/gdrywt *100.
            l=SOM_ax.plot(t,result['Total Sorbed cellulose']/SOC_mol_m3*SOC_pct+SOC_pct*(1-cellulose_SOC_frac),label='SOM',ls=ls,color=color,**SOM_args)[0]
            SOM_ax.set_ylabel('Concentration\n(SOC %)')
        else:
            l=SOM_ax.plot(t,result['Total Sorbed cellulose']*1e-3,label='SOM',ls=ls,**SOM_args)[0]
            SOM_ax.set_ylabel('Concentration\n(mmol C/cm$^{-3}$)')

        SOM_ax.set_title('SOM remaining')
        SOM_ax.set_xlabel('Time (days)')

    if obs is not None:
        obs_mod['pH_obs']=obs['pH'].set_axis(obs.Incubation_Time)
        obs_mod['pH_mod']=-log10(result['Free H+'])[obs_mod_pts]
    if pH_ax is not None:
        pH_ax.plot(t,-log10(result['Free H+']),c=color,ls=ls)

        pH_ax.set_title('pH')
        pH_ax.set_ylabel('pH')
        pH_ax.set_xlabel('Time (days)')
        
        if obs is not None:
            if obs_marker is not None:
                markers=['o','s','^','+','<','x']
                obs_cats=unique(obs[obs_marker])
                print(obs_cats)
                for x in range(len(obs_cats)):
                    xx=obs[obs_marker]==obs_cats[x]
                    pH_ax.plot(obs.Incubation_Time[xx],obs['pH'][xx],marker=markers[x],linestyle='None',label='Measured',color='k')
            else:
                pH_ax.plot(obs.Incubation_Time,obs['pH'],marker='o',linestyle='None',label='Measured',color='k',**pH_args)
    
    if obs is not None:
        obs_mod['Fe_II_obs']=obs['Fe_II'].set_axis(obs.Incubation_Time)
        obs_mod['Fe_II_mod']=result['Total Fe++'][obs_mod_pts]*result['Porosity']*1e3*cm3_to_dwt
    if Fe_ax is not None:
        molar_volume=34.3600 # From database. cm3/mol
        molar_weight = 106.8690
        # Add sorbed Fe to this
        if gdrywt:
            # Assume we can use SOC % to convert from volume to dry weight
            l=Fe_ax.plot(t,result['Fe(OH)3 VF']/molar_volume*1e6*cm3_to_dwt   ,label='Iron oxide',c=color,ls='-',**Fe_args)[0]
            Fe_ax.set_ylabel('Concentration\n($\mu$mol g dwt$^{-1}$)')
            l=Fe_ax.plot(t,result['Total Fe+++']*result['Porosity']*1e3*cm3_to_dwt   ,label='Fe$^{3\!\!+}$',ls='--',c=color,**Fe_args)[0]
            
            l=Fe_ax.plot(t,result['Total Fe++']*result['Porosity']*1e3*cm3_to_dwt ,ls=':' ,c=color ,label='Fe$^{2\!\!+}$',**Fe_args)[0]
            Fe_ax.plot(t,(result['Total Sorbed Fe++']+result['Total Sorbed Fe+++'])*1e6/100**3*cm3_to_dwt,ls='dashdot',c=color,**Fe_args,label='Sorbed Fe')
            if obs is not None:
                if obs_marker is not None:
                    markers=['o','s','^','+','<','x']
                    obs_cats=unique(obs[obs_marker])
                    print(obs_cats)
                    for x in range(len(obs_cats)):
                        xx=obs[obs_marker]==obs_cats[x]
                        Fe_ax.plot(t,obs.Incubation_Time[xx],obs['Fe_II'][xx],marker=markers[x],linestyle='None',label='Measured Fe$^{2\!\!+}$',color='k')
                else:
                    Fe_ax.plot(obs.Incubation_Time,obs['Fe_II'],marker='o',linestyle='None',label='Measured Fe$^{2\!\!+}$',color='k')
        else:
            l=Fe_ax.plot(t,result['Fe(OH)3 VF']/molar_volume*1e6   ,label='Iron oxide (solid)')[0]
            Fe_ax.set_ylabel('Concentration\n($\mu$mol cm$^{-3}$)')
        
            # M/L to umol/cm3: 1e6/1e3=1e3
            l=Fe_ax.plot(t,result['Total Fe+++']*result['Porosity']*1e3   ,label='Fe+++',ls='-')[0]
            
            l=Fe_ax.plot(t,result['Total Fe++']*result['Porosity']*1e3 ,ls='--'  ,label='Fe++')[0]
        
        Fe_ax.set_title('Fe species')
        
        Fe_ax.set_xlabel('Time (days)')
        if do_legend:
            Fe_ax.legend(fontsize='small')
    
    if obs is not None:
        CO2mean=obs[['CO2_1','CO2_2','CO2_3']].astype(float).diff().mean(axis=1)/obs['Incubation_Time'].diff()
        obs_mod['CO2flux_obs']=CO2mean.set_axis(obs.Incubation_Time)
        obs_mod['CO2flux_mod']=(diff(result['Total Tracer']*result['Porosity'])/diff(result['time'])*1e3*cm3_to_dwt)[obs_mod_pts]
        obs_mod['CO2flux_obs_err']=(obs[['CO2_1','CO2_2','CO2_3']].astype(float).diff().sem(axis=1) /obs['Incubation_Time'].diff()).set_axis(obs['Incubation_Time'])
    if CO2flux_ax is not None:
        # gasflux_ax.set_yscale('log')
        if gdrywt:
            
            l=CO2flux_ax.plot(t[:-1],diff(result['Total Tracer']*result['Porosity'])/diff(result['time'])*1e3*cm3_to_dwt,label='CO2',color=color,ls=ls,**CO2flux_args)[0]
            CO2flux_ax.set_ylabel('CO$_2$ flux rate\n($\mu$mol g dwt$^{-1}$ day$^{-1}$)')
            
            if obs is not None:
                CO2mean=obs[['CO2_1','CO2_2','CO2_3']].astype(float).diff().mean(axis=1)/obs['Incubation_Time'].diff()
                CO2std=obs[['CO2_1','CO2_2','CO2_3']].astype(float).diff().std(axis=1) /obs['Incubation_Time'].diff()
                

                if obs_marker is not None:
                    markers=['o','s','^','+','<','x']
                    obs_cats=unique(obs[obs_marker])
                    print(obs_cats)
                    for x in range(len(obs_cats)):
                        xx=obs[obs_marker]==obs_cats[x]
                        CO2flux_ax.errorbar(obs.Incubation_Time[xx],CO2mean[xx],yerr=CO2std[xx],marker=markers[x],linestyle='None',label='Measured CO$_2$',color='C0')
                else:
                    CO2flux_ax.errorbar(obs.Incubation_Time,CO2mean,yerr=CO2std,marker='o',linestyle='None',label='Measured CO$_2$',color='k')
                
        else:
            
            l=CO2flux_ax.plot(result['time'][:-1],diff(result['Total Tracer']*result['Porosity'])/diff(result['time'])*1e3,label='CO2',c=color,**CO2flux_args)[0]
            CO2flux_ax.set_ylabel('CO$_2$ flux rate\n($\mu$mol cm$^{-3}$ day$^{-1}$)')

        CO2flux_ax.set_title('CO$_2$ flux')
        
        CO2flux_ax.set_xlabel('Time (days)')
        # if do_legend:
        #     CO2flux_ax.legend(fontsize='small')

    if obs is not None:
        CH4mean=obs[['CH4_1','CH4_2','CH4_3']].astype(float).diff().mean(axis=1)/obs['Incubation_Time'].diff()
        obs_mod['CH4flux_obs']=CH4mean.set_axis(obs.Incubation_Time)
        obs_mod['CH4flux_obs_err']=(obs[['CH4_1','CH4_2','CH4_3']].astype(float).diff().sem(axis=1)/obs['Incubation_Time'].diff()).set_axis(obs['Incubation_Time'])
        obs_mod['CH4flux_mod']=(diff(result['Total CH4(aq)']*result['Porosity'])/diff(result['time'])*1e3*cm3_to_dwt)[obs_mod_pts]
    if CH4flux_ax is not None:
        if gdrywt:
            l=CH4flux_ax.plot(t[:-1],diff(result['Total CH4(aq)']*result['Porosity'])/diff(result['time'])*1e3*cm3_to_dwt,label='CH4',c=color,ls=ls,**CH4flux_args)[0]
            CH4flux_ax.set_ylabel('CH$_4$ Flux rate\n($\mu$mol g dwt$^{-1}$ day$^{-1}$)')
            
            if obs is not None:
                CH4mean=obs[['CH4_1','CH4_2','CH4_3']].astype(float).diff().mean(axis=1)/obs['Incubation_Time'].diff() 
                CH4std=obs[['CH4_1','CH4_2','CH4_3']].astype(float).diff().sem(axis=1) /obs['Incubation_Time'].diff() 
                
                if obs_marker is not None:
                    markers=['o','s','^','+','<','x']
                    obs_cats=unique(obs[obs_marker])
                    print(obs_cats)
                    for x in range(len(obs_cats)):
                        xx=obs[obs_marker]==obs_cats[x]
                        CH4flux_ax.errorbar(obs.Incubation_Time[xx],CH4mean[xx],yerr=CH4std[xx],marker=markers[x],linestyle='None',label='Measured CH$_4$',color='k')
                else:
                    CH4flux_ax.errorbar(obs.Incubation_Time,CH4mean,yerr=CH4std,marker='o',linestyle='None',label='Measured CH$_4$',color='k')
                
                
        else:
            l=CH4flux_ax.plot(t[:-1],diff(result['Total CH4(aq)']*result['Porosity'])/diff(result['time'])*1e3,label='CH4',**CH4flux_args,c=color)[0]
            
            CH4flux_ax.set_ylabel('CH$_4$ flux rate\n($\mu$mol cm$^{-3}$ day$^{-1}$)')

        CH4flux_ax.set_title('CH$_4$ flux')
        
        CH4flux_ax.set_xlabel('Time (days)')
        # if do_legend:
        #     CH4flux_ax.legend(fontsize='small')
        
    if obs is not None:
        obs_mod['DOM_obs']=obs['WEOC'].set_axis(obs.Incubation_Time)*1e-6/BD*1000/result['Porosity'][0].item()
        obs_mod['DOM_mod']=result['Total DOM1'][obs_mod_pts]
        obs_mod['Acetate_obs']=obs['TOAC'].set_axis(obs.Incubation_Time)*1e-6/BD*1000/result['Porosity'][0].item()
        obs_mod['Acetate_mod']=result['Total Acetate-'][obs_mod_pts]
    if porewater_ax is not None:
        porewater_ax.set_yscale('log')
        porewater_ax.plot(t,result['Total DOM1'],ls='-',label='DOM',c=color)
        porewater_ax.plot(t,result['Total Acetate-'],ls='--',label='Acetate',c=color)
        # porewater_ax.plot(t,result['Total O2(aq)'],'--',label='O2',c='C4')
        # porewater_ax.plot(t,result['Total Fe+++'],'--',label='Fe+++',c='C1')
        # porewater_ax.plot(t,result['Free Fe+++'],':',label='Fe+++',c='C1')
        # porewater_ax.plot(t,result['Total Fe++'],':',label='Fe++',c='C2')
        
        porewater_ax.set_title('Porewater concentrations')
        porewater_ax.set_ylabel('Concentration (M)')
        porewater_ax.set_xlabel('Time (days)')
        
        
        if obs is not None:
            if obs_marker is not None:
                markers=['o','s','^','+','<','x']
                obs_cats=unique(obs[obs_marker])
                print(obs_cats)
                for x in range(len(obs_cats)):
                    xx=obs[obs_marker]==obs_cats[x]
                    porewater_ax.plot(obs.Incubation_Time[xx],obs['TOAC'][xx]*1e-6/BD*1000/result['Porosity'][0].item(),marker=markers[x],linestyle='None',label='Measured TOAC',color='C3',**porewater_args)
                    porewater_ax.plot(obs.Incubation_Time[xx],obs['WEOC'][xx]*1e-6/BD*1000/result['Porosity'][0].item(),marker=markers[x],linestyle='None',label='Measured WEOC',color='C0',**porewater_args)
            else:
                porewater_ax.plot(obs.Incubation_Time,obs['TOAC']*1e-6/BD*1000/result['Porosity'][0].item(),marker='o',linestyle='None',label='Measured TOAC',color='k')
                porewater_ax.plot(obs.Incubation_Time,obs['WEOC']*1e-6/BD*1000/result['Porosity'][0].item(),marker='s',linestyle='None',label='Measured WEOC',color='k')


        
        if do_legend:
            porewater_ax.legend(fontsize='small')
            
    if cation_ax is not None:
        for n,cation in enumerate(decomp_network.pools_list_to_dict(arctic.make_reactions())['Cation exchange']['cations'].keys() ):
            if cation == 'H+':
                continue
            cation_ax.plot(t,result['Total Sorbed '+cation]/(BD*1e-3*100**3)*1000,c='C'+str(n),label=cation,**cation_args)
        cation_ax.plot(t,result['CEC H+']/(BD*1e-3*100**3)*1000,label='H+',c='C'+str(n+1),ls=ls,**cation_args)
        
        cation_ax.set_ylabel('Exch conc\n(mmol/kg)')
        
        if do_legend:
            cation_ax.legend(fontsize='small')

    return obs_mod
        

def plot_timeseries(data_list,group_suffix='',color='C0',fig=None,cation_fig=None):
    # if isinstance(fname,str):
    #     fnames=[fname]
    # else:
    #     fnames=fname

    result_organic_trough=None
    result_organic_nottrough=None
    result_mineral_trough=None
    result_mineral_nottrough=None
    result_highO2_organic=None

    # for f in fnames:
    #     if 'organic_trough'+group_suffix in get_groups(f):
    #         result_organic_trough=xarray.open_dataset(f,group='organic_trough'+group_suffix)
    #     if 'organic_nottrough'+group_suffix in get_groups(f):
    #         result_organic_nottrough=xarray.open_dataset(f,group='organic_nottrough'+group_suffix)
    #     if 'mineral_trough'+group_suffix in get_groups(f):
    #         result_mineral_trough=xarray.open_dataset(f,group='mineral_trough'+group_suffix)
    #     if 'mineral_nottrough'+group_suffix in get_groups(f):
    #         result_mineral_nottrough=xarray.open_dataset(f,group='mineral_nottrough'+group_suffix)
    #     if 'highO2_organic'+group_suffix in get_groups(f):
    #         result_highO2_organic=xarray.open_dataset(f,group='highO2_organic'+group_suffix)


    result_organic_trough=data_list['organic_trough'+group_suffix]
    result_organic_nottrough=data_list['organic_nottrough'+group_suffix]
    result_mineral_trough=data_list['mineral_trough'+group_suffix]
    result_mineral_nottrough=data_list['mineral_nottrough'+group_suffix]
    result_highO2_organic=data_list['highO2_organic'+group_suffix]

    if result_organic_trough is None:
        raise ValueError('result_organic_trough'+group_suffix+' not found')
    if result_organic_nottrough is None:
        raise ValueError('result_organic_nottrough'+group_suffix+' not found')
    if result_mineral_trough is None:
        raise ValueError('result_mineral_trough'+group_suffix+' not found')
    if result_mineral_nottrough is None:
        raise ValueError('result_mineral_nottrough'+group_suffix+' not found')
    if result_highO2_organic is None:
        raise ValueError('result_highO2_organic'+group_suffix+' not found')

    import os
    # fig,axes=subplots(3,1,num='Organic horizon Anoxic',figsize=(6,8.4),clear=True)
    if cation_fig is None:
        cation_fig,cation_axes=subplots(ncols=5,num='Cations',clear=True,figsize=(13,4))
    else:
        cation_axes=cation_fig.axes
    if fig is None:
        fig,axs=subplots(nrows=5,ncols=5,num='Time series',clear=True,figsize=(16,8),sharex=False)
        newplot=True
    else:
        axs=np.reshape(fig.axes,(5,5))
        newplot=False
    from string import ascii_lowercase
    for x in range(5):
        for y in range(5):
            axs[x,y].set_title('('+ascii_lowercase[x*5+y]+')',loc='left')

    axes=axs[:,1]
    obs=arctic.get_layer(9,'Organic')
    obs_mod=plot_result(result_organic_trough,Fe_ax=axes[2],CO2flux_ax=axes[0],CH4flux_ax=axes[1],pH_ax=axes[3],porewater_ax=axes[4],cation_ax=cation_axes[1],do_legend=False,BD=arctic.BD_organic_trough,
                gdrywt=True,SOC_pct=arctic.SOC_layermean['Organic'],cellulose_SOC_frac=arctic.cellulosefrac,#BD=BD_layerest2_trough['Organic',True],
                # obs=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Anoxic')&(Barrow_synthesis_data['Soil_layer'].str.capitalize()=='Organic')
                        # &(Barrow_synthesis_data['Incubation_Temperature']>4)&(Barrow_synthesis_data['Microtopography'].str.lower()=='trough')])
                        obs=obs,color=color)
    # axes[0].set_ylim(bottom=-0.01)
    obs_mod['Fe_II_obs_err']=pandas.Series([5.71,3.4,16.0,7.6],index=[0,17,30,57])/np.sqrt(2)
    obs_mod['DOM_obs_err']=pandas.Series([309-218,376-262],index=[30,57])/np.sqrt(3)*1e-6/arctic.BD_organic_trough*1000/result_organic_trough['Porosity'][0].item()
    obs_mod['TOAC_obs_err']=pandas.Series([50,10],index=[30,57])/np.sqrt(3)*1e-6/arctic.BD_organic_trough*1000/result_organic_trough['Porosity'][0].item()

    axes[2].errorbar(obs_mod.index,obs_mod['Fe_II_obs'],yerr=obs_mod['Fe_II_obs_err'],marker='None',c='k')
    axes[4].errorbar(obs_mod.index,obs_mod['DOM_obs'],yerr=obs_mod['DOM_obs_err'],marker='None',c='k')
    axes[4].errorbar(obs_mod.index,obs_mod['Acetate_obs'],yerr=obs_mod['TOAC_obs_err'],marker='None',c='k')

    # axes[0].set_ylim(bottom=0.9e-11,top=15)
    axes[0].set_yscale('linear')
    

    # xmax=axes[0].get_xlim()[1]
    # xmax=max(arctic.simlength,90)
    xmax=85
    if newplot:
        axes[0].text(0.5,1.3,'%s Organic horizon (Anoxic)'%(obs['Microtopography'].iloc[0]),fontsize='large',transform=axes[0].transAxes,ha='center')
        for ax in axes:
            ax.axvspan(arctic.initfrac*arctic.simlength,xmax,color='b',alpha=0.1)
            ax.set_xlim(left=-5,right=xmax)
        
    axes=axs[:,2]
    obs=arctic.get_layer(5,'Organic')
    obs_mod=plot_result(result_organic_nottrough,Fe_ax=axes[2],CO2flux_ax=axes[0],CH4flux_ax=axes[1],pH_ax=axes[3],porewater_ax=axes[4],cation_ax=cation_axes[2],do_legend=False,
                gdrywt=True,SOC_pct=arctic.SOC_layermean['Organic'],cellulose_SOC_frac=arctic.cellulosefrac,BD=arctic.BD_organic_nottrough,
                # obs=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Anoxic')&(Barrow_synthesis_data['Soil_layer'].str.capitalize()=='Organic')
                        # &(Barrow_synthesis_data['Incubation_Temperature']>4)&(Barrow_synthesis_data['Microtopography'].str.lower()!='trough')])
                        obs=obs,color=color)
    # axes[0].set_ylim(bottom=-0.01)
    obs_mod['Fe_II_obs_err']=pandas.Series([1.04,3.9,5.7,4.2],index=[0,15,29,55])/np.sqrt(2)
    obs_mod['DOM_obs_err']=pandas.Series([147-84],index=[29])/np.sqrt(3)*1e-6/arctic.BD_organic_nottrough*1000/result_organic_nottrough['Porosity'][0].item()
    obs_mod['TOAC_obs_err']=pandas.Series([5],index=[29])/np.sqrt(3)*1e-6/arctic.BD_organic_nottrough*1000/result_organic_nottrough['Porosity'][0].item()

    axes[2].errorbar(obs_mod.index,obs_mod['Fe_II_obs'],yerr=obs_mod['Fe_II_obs_err'],marker='None',c='k')
    axes[4].errorbar(obs_mod.index,obs_mod['DOM_obs'],yerr=obs_mod['DOM_obs_err'],marker='None',c='k')
    axes[4].errorbar(obs_mod.index,obs_mod['Acetate_obs'],yerr=obs_mod['TOAC_obs_err'],marker='None',c='k')

    # axes[0].set_ylim(bottom=0.9e-11,top=15)
    axes[0].set_yscale('linear')


    # xmax=axes[0].get_xlim()[1]
    # xmax=90
    if newplot:
        axes[0].text(0.5,1.3,'%s Organic horizon (Anoxic)'%(obs['Microtopography'].iloc[0]),fontsize='large',transform=axes[0].transAxes,ha='center')
        for ax in axes:
            ax.axvspan(arctic.initfrac*arctic.simlength,xmax,color='b',alpha=0.1)
            ax.set_xlim(left=-5,right=xmax)


    # fig,axes=subplots(3,1,num='Mineral horizon Anoxic',figsize=(6,8.4),clear=True)
    axes=axs[:,3]
    obs=arctic.get_layer(9,'Mineral')
    obs_mod=plot_result(result_mineral_trough,Fe_ax=axes[2],CO2flux_ax=axes[0],CH4flux_ax=axes[1],pH_ax=axes[3],porewater_ax=axes[4],cation_ax=cation_axes[3],do_legend=False,
        gdrywt=True,SOC_pct=arctic.SOC_layermean['Mineral'],cellulose_SOC_frac=arctic.cellulosefrac,BD=arctic.BD_mineral_trough,
        # obs=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Anoxic')&(Barrow_synthesis_data['Soil_layer'].str.capitalize()=='Mineral')
            # &(Barrow_synthesis_data['Incubation_Temperature']>4)&(Barrow_synthesis_data['Microtopography'].str.lower()=='trough')])
            obs=obs,color=color)
    # axes[0].set_ylim(bottom=-0.01)

    obs_mod['Fe_II_obs_err']=pandas.Series([4.5,2.2,4.7,1.5],index=[0,17,30,57])/np.sqrt(2)
    obs_mod['DOM_obs_err']=pandas.Series([239-156,213-89],index=[30,57])/np.sqrt(3)*1e-6/arctic.BD_mineral_trough*1000/result_mineral_trough['Porosity'][0].item()
    obs_mod['TOAC_obs_err']=pandas.Series([5,2],index=[30,57])/np.sqrt(3)*1e-6/arctic.BD_mineral_trough*1000/result_mineral_trough['Porosity'][0].item()


    axes[2].errorbar(obs_mod.index,obs_mod['Fe_II_obs'],yerr=obs_mod['Fe_II_obs_err'],marker='None',c='k')
    axes[4].errorbar(obs_mod.index,obs_mod['DOM_obs'],yerr=obs_mod['DOM_obs_err'],marker='None',c='k')
    axes[4].errorbar(obs_mod.index,obs_mod['Acetate_obs'],yerr=obs_mod['TOAC_obs_err'],marker='None',c='k')

    # axes[0].set_ylim(bottom=0.9e-11,top=3.0)
    axes[0].set_yscale('linear')

    # xmax=axes[0].get_xlim()[1]
    # xmax=90
    if newplot:
        for ax in axes:
            ax.axvspan(arctic.initfrac*arctic.simlength,xmax,color='b',alpha=0.1)
            ax.set_xlim(left=-5,right=xmax)
        axes[0].text(0.5,1.3,'%s Mineral horizon (Anoxic)'%(obs['Microtopography'].iloc[0]),fontsize='large',transform=axes[0].transAxes,ha='center')


    axes=axs[:,4]
    obs=arctic.get_layer(5,'Mineral')
    obs_mod=plot_result(result_mineral_nottrough,Fe_ax=axes[2],CO2flux_ax=axes[0],CH4flux_ax=axes[1],pH_ax=axes[3],porewater_ax=axes[4],cation_ax=cation_axes[4],do_legend=False,
        gdrywt=True,SOC_pct=arctic.SOC_layermean['Mineral'],cellulose_SOC_frac=arctic.cellulosefrac,BD=arctic.BD_mineral_nottrough,
        # obs=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Anoxic')&(Barrow_synthesis_data['Soil_layer'].str.capitalize()=='Mineral')&
            # (Barrow_synthesis_data['Incubation_Temperature']>4)&(Barrow_synthesis_data['Microtopography'].str.lower()!='trough')])
            obs=obs,color=color)
    # axes[0].set_ylim(bottom=-0.01)

    obs_mod['Fe_II_obs_err']=pandas.Series([0.38,3.6,4.5,2.3],index=[0,15,29,55])/np.sqrt(2)
    obs_mod['DOM_obs_err']=pandas.Series([22,28+16],index=[29,55])/np.sqrt(3)*1e-6/arctic.BD_mineral_nottrough*1000/result_mineral_nottrough['Porosity'][0].item()
    obs_mod['TOAC_obs_err']=pandas.Series([8,2],index=[29,55])/np.sqrt(3)*1e-6/arctic.BD_mineral_nottrough*1000/result_mineral_nottrough['Porosity'][0].item()

    axes[2].errorbar(obs_mod.index,obs_mod['Fe_II_obs'],yerr=obs_mod['Fe_II_obs_err'],marker='None',c='k')
    axes[4].errorbar(obs_mod.index,obs_mod['DOM_obs'],yerr=obs_mod['DOM_obs_err'],marker='None',c='k')
    axes[4].errorbar(obs_mod.index,obs_mod['Acetate_obs'],yerr=obs_mod['TOAC_obs_err'],marker='None',c='k')

    # axes[0].set_ylim(bottom=0.9e-11,top=3.0)
    axes[0].set_yscale('linear')

    # xmax=axes[0].get_xlim()[1]
    # xmax=90
    if newplot:
        for ax in axes:
            ax.axvspan(arctic.initfrac*arctic.simlength,xmax,color='b',alpha=0.1)
            ax.set_xlim(left=-5,right=xmax)
        axes[0].text(0.5,1.3,'%s Mineral horizon (Anoxic)'%(obs['Microtopography'].iloc[0]),fontsize='large',transform=axes[0].transAxes,ha='center')


    # Oxic
    # fig,axes=subplots(3,1,num='Organic horizon Oxic',figsize=(6,8.4),clear=True)
    axes=axs[:,0]
    obs=arctic.get_layer(3,'Organic')
    plot_result(result_highO2_organic,Fe_ax=axes[2],CO2flux_ax=axes[0],CH4flux_ax=axes[1],pH_ax=axes[3],porewater_ax=axes[4],cation_ax=cation_axes[0],do_legend=False,
            gdrywt=True,SOC_pct=arctic.SOC_layermean['Organic'],cellulose_SOC_frac=arctic.cellulosefrac,BD=arctic.BD_atmoO2_organic,
            # obs=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Oxic')&(Barrow_synthesis_data['Soil_layer'].str.capitalize()=='Organic')&(Barrow_synthesis_data['Incubation_Temperature']>4)])
            obs=obs,color=color)
    # axes[0].set_ylim(bottom=-0.01,top=30)
    axes[1].set_ylim(bottom=-0.01,top=0.2)

    axes[0].set_yscale('linear')
    # axes[0].set_ylim(-5e-2,155e-2) # This leaves out Core NGADG0003 which has 10x higher flux for some reason
    # axes[1].set_ylim(bottom=0.9e-11)

    # axes[0].set_xlim(right=71)
    for ax in fig.axes:
        ax.set_xlim(left=-5,right=xmax)
    for ax in axs[3,:]:
        ax.set_ylim(3.8,6.1)
    for ax in axs[4,:]:
        ax.set_ylim(bottom=1e-5)

    if newplot:
        leg=axs[2,4].legend(loc='upper right',ncol=2,edgecolor='k')
        cation_axes[0].legend()
        axes[4].legend(loc='upper right')
        axes[0].text(0.5,1.3,'%s Organic horizon (Oxic)'%(obs['Microtopography'].iloc[0]),fontsize='large',transform=axes[0].transAxes,ha='center')

    return 


def get_groups(filename):
    import netCDF4
    with netCDF4.Dataset(filename) as d:
        groups=list(d.groups.keys())
    return groups


def mod_obs_comp(data_list,group_suffix=''):
    # if isinstance(fname,str):
    #     fnames=[fname]
    # else:
    #     fnames=fname

    result_organic_trough=None
    result_organic_nottrough=None
    result_mineral_trough=None
    result_mineral_nottrough=None
    result_highO2_organic=None

    # for f in fnames:
    #     if 'organic_trough'+group_suffix in get_groups(f):
    #         result_organic_trough=xarray.open_dataset(f,group='organic_trough'+group_suffix)
    #     if 'organic_nottrough'+group_suffix in get_groups(f):
    #         result_organic_nottrough=xarray.open_dataset(f,group='organic_nottrough'+group_suffix)
    #     if 'mineral_trough'+group_suffix in get_groups(f):
    #         result_mineral_trough=xarray.open_dataset(f,group='mineral_trough'+group_suffix)
    #     if 'mineral_nottrough'+group_suffix in get_groups(f):
    #         result_mineral_nottrough=xarray.open_dataset(f,group='mineral_nottrough'+group_suffix)
    #     if 'highO2_organic'+group_suffix in get_groups(f):
    #         result_highO2_organic=xarray.open_dataset(f,group='highO2_organic'+group_suffix)

    result_organic_trough=data_list['organic_trough'+group_suffix]
    result_organic_nottrough=data_list['organic_nottrough'+group_suffix]
    result_mineral_trough=data_list['mineral_trough'+group_suffix]
    result_mineral_nottrough=data_list['mineral_nottrough'+group_suffix]
    result_highO2_organic=data_list['highO2_organic'+group_suffix]

    obs_mods={}
    if result_organic_trough is not None:
            obs_mods['organic_trough']=plot_result(result_organic_trough,BD=arctic.BD_organic_trough,
                gdrywt=True,SOC_pct=arctic.SOC_layermean['Organic'],cellulose_SOC_frac=arctic.cellulosefrac,obs=arctic.get_layer(9,'Organic'))
    
    if result_organic_nottrough is not None:
            obs_mods['organic_nottrough']=plot_result(result_organic_nottrough,gdrywt=True,SOC_pct=arctic.SOC_layermean['Organic'],cellulose_SOC_frac=arctic.cellulosefrac,BD=arctic.BD_organic_nottrough,
                        obs=arctic.get_layer(5,'Organic'))
    if result_mineral_trough is not None:
            obs_mods['mineral_trough']=plot_result(result_mineral_trough,gdrywt=True,SOC_pct=arctic.SOC_layermean['Mineral'],cellulose_SOC_frac=arctic.cellulosefrac,
                                    BD=arctic.BD_mineral_trough,obs=arctic.get_layer(9,'Mineral'))
    if result_mineral_nottrough is not None:
            obs_mods['mineral_nottrough']=plot_result(result_mineral_nottrough,gdrywt=True,SOC_pct=arctic.SOC_layermean['Mineral'],cellulose_SOC_frac=arctic.cellulosefrac,
                         BD=arctic.BD_mineral_nottrough,obs=arctic.get_layer(5,'Mineral'))
    if result_highO2_organic is not None:
            obs_mods['highO2_organic']=plot_result(result_highO2_organic,gdrywt=True,SOC_pct=arctic.SOC_layermean['Organic'],cellulose_SOC_frac=arctic.cellulosefrac,
                        BD=arctic.BD_atmoO2_organic,obs=arctic.get_layer(3,'Organic'))

    obs_mods['organic_trough']['Fe_II_obs_err']=pandas.Series([5.71,3.4,16.0,7.6],index=[0,17,30,57])/np.sqrt(2) # Chowdhury et al 2015 Table 1 and 4
    obs_mods['mineral_trough']['Fe_II_obs_err']=pandas.Series([4.5,2.2,4.7,1.5],index=[0,17,30,57])/np.sqrt(2)
    obs_mods['organic_nottrough']['Fe_II_obs_err']=pandas.Series([1.04,3.9,5.7,4.2],index=[0,15,29,55])/np.sqrt(2)
    obs_mods['mineral_nottrough']['Fe_II_obs_err']=pandas.Series([0.38,3.6,4.5,2.3],index=[0,15,29,55])/np.sqrt(2)
    obs_mods['highO2_organic']['Fe_II_obs_err']=np.nan

    # Herndon et al. 2015
    obs_mods['organic_trough']['DOM_obs_err']=pandas.Series([309-218,376-262],index=[30,57])/np.sqrt(3)*1e-6/arctic.BD_organic_trough*1000/result_organic_trough['Porosity'][0].item()
    obs_mods['mineral_trough']['DOM_obs_err']=pandas.Series([239-156,213-89],index=[30,57])/np.sqrt(3)*1e-6/arctic.BD_mineral_trough*1000/result_mineral_trough['Porosity'][0].item()
    obs_mods['organic_nottrough']['DOM_obs_err']=pandas.Series([147-84],index=[29])/np.sqrt(3)*1e-6/arctic.BD_organic_nottrough*1000/result_organic_nottrough['Porosity'][0].item()
    obs_mods['mineral_nottrough']['DOM_obs_err']=pandas.Series([22,28+16],index=[29,55])/np.sqrt(3)*1e-6/arctic.BD_mineral_nottrough*1000/result_mineral_nottrough['Porosity'][0].item()
    obs_mods['highO2_organic']['DOM_obs_err']=np.nan

    obs_mods['organic_trough']['TOAC_obs_err']=pandas.Series([50,10],index=[30,57])/np.sqrt(3)*1e-6/arctic.BD_organic_trough*1000/result_organic_trough['Porosity'][0].item()
    obs_mods['mineral_trough']['TOAC_obs_err']=pandas.Series([5,2],index=[30,57])/np.sqrt(3)*1e-6/arctic.BD_mineral_trough*1000/result_mineral_trough['Porosity'][0].item()
    obs_mods['organic_nottrough']['TOAC_obs_err']=pandas.Series([5],index=[29])/np.sqrt(3)*1e-6/arctic.BD_organic_nottrough*1000/result_organic_nottrough['Porosity'][0].item()
    obs_mods['mineral_nottrough']['TOAC_obs_err']=pandas.Series([8,2],index=[29,55])/np.sqrt(3)*1e-6/arctic.BD_mineral_nottrough*1000/result_mineral_nottrough['Porosity'][0].item()
    obs_mods['highO2_organic']['TOAC_obs_err']=np.nan

    return obs_mods
 
def getstats(obs_mod,logged=['DOM','Acetate']):
    specs=['pH','Fe_II','CO2flux','CH4flux','DOM','Acetate']
    out=pandas.DataFrame(index=specs,columns=['RMSE','nRMSE','MAE','R','R2'])
    for spec in specs:
        mod=obs_mod[spec+'_mod']
        obs=obs_mod[spec+'_obs']
        if spec in logged:
            mod=log10(mod)
            obs=log10(obs)
        out['RMSE'][spec] = sqrt(((obs-mod)**2).mean())
        out['nRMSE'][spec] = out['RMSE'][spec]/obs.mean()
        out['MAE'][spec]  = abs(obs-mod).mean()
        if spec in logged:
            out['R'][spec]    = log10(obs_mod[[spec+'_obs',spec+'_mod']]).corr().iloc[0,1]
        else:
            out['R'][spec] = obs_mod[[spec+'_obs',spec+'_mod']].corr().iloc[0,1]
        out['R2'][spec]   = out['R'][spec]**2

    return out.rename({spec:'log(%s)'%spec for spec in logged})

def save_all_figs(dirname,format='png',**kwargs):
    for fname in get_figlabels():
        fname_fixed=fname.replace('/','-')
        print(fname_fixed)
        figure(fname_fixed).savefig('{dirname:s}/{fname:s}.{format}'.format(dirname=dirname,format=format,fname=fname_fixed),**kwargs)


if __name__ == '__main__':
    pos={'Fe(OH)3': (212.0, 348.0),
        'Fe+++': (212.0, 306.0),
        'cellulose': (1303.0, 594.0),
        'DOM1': (1230.0, 450.0),
        'O2(aq)': (909.0, 306.0),
        'HCO3-': (1009.0, 162.0),
        'Fe++': (216.0, 162.0),
        'CH4(aq)': (726.0, 18.0),
        'Acetate-': (689.0, 306.0),
        'H2(aq)': (1275.0, 306.0),
        'Mg++': (106.0, 18.0),
        'Ca++': (190.0, 18.0),
        'Na+': (266.0, 18.0),
        'K+': (27.0, 18.0),
        'Aerobic decomposition': (1377.0, 450.0),
        'Hydrolysis': (1265.0, 522.0),
        'fermentation': (1220.0, 378.0),
        'DOM aerobic respiration': (1211.0, 234.0),
        'Fe(II) microbial oxidation': (339.0, 90.0),
        'Cation exchange': (148.0, 90.0),
        'Hydrogenotrophic methanogenesis': (1009.0, 90.0),
        'Acetate aerobic respiration': (807.0, 234.0),
        'Hydrogen oxidation': (1012.0, 234.0),
        'Fe(III) reduction': (333.0, 234.0),
        'Acetaclastic methanogenesis': (571.0, 234.0)}
    networkfig2=figure('Reaction network (with reactions)',clear=True,figsize=(9.4,6.6))
    drawn=decomp_network.draw_network_with_reactions(decomp_network.decomp_network(pools=arctic.pools,reactions=arctic.make_reactions()),node_size=1200,arrowstyle='->',arrowsize=8.5,#edge_color='gray',
                omit=['NH4+','Rock(s)','gas','surf_complex','secondary','H+','implicit','Fe(II) abiotic oxidation','Fe(OH)2','Calcite'],
                namechanges={'cellulose':'Cellulose','DOM1':'DOM','O2(aq)':'O$_2$','CH4(aq)':'CH$_4$','HCO3-':'CO$_2$',
                            'Fe(OH)2':'Fe(OH)$_2$','Fe(OH)3':'Fe(OH)$_3$','Fe++':r'Fe$^\mathrm{+\!\!+}$','Fe+++':r'Fe$^\mathrm{+\!\!+\!\!\!+}$',
                            'Acetate-':'Acetate','fermentation':'Fermentation','Fe(II) microbial oxidation':r'Fe$^\mathrm{+\!\!+}$'+'\noxidation',
                            'DOM aerobic respiration':'DOM aerobic\nrespiration','Hydrogen oxidation':'Hydrogen\noxidation','Aerobic decomposition':'Aerobic\ndecomposition',
                            'Acetaclastic methanogenesis':'Acetoclastic\nmethanogenesis','Acetate aerobic respiration':'Acetate\naerobic respiration',
                            'Hydrogenotrophic methanogenesis':'Hydrogenotrophic\nmethanogenesis','Gas':'Dissolved gas','Primary aqueous':'Dissolved ion',
                            'Unknown':'Cation exchange','Fe(III) reduction':r'Fe$^\mathrm{+\!\!+\!\!\!+}$'+'\nreduction','H2(aq)':'H$_2$','SOMdecomp Reaction':'SOM Decomposition'
                            },markers={'mineral':'8'},pos=pos)

    colors={'Anaerobic':'C1','Periodic':'C2','Low Fe':'C3','Aerobic':'C0'}

    if len(sys.argv)<2:
        # fname='Arctic_Fe_output/results_2021-03-31.nc'
        import glob
        fname=sorted(glob.glob('Arctic_Fe_output/output_2022-06-30_0?.nc'))
    else:
        fname=sys.argv[1:]

    data_list={}
    data_list_noinhib={}
    data_list_noFe={}

    import netCDF4
    for filename in fname:
        with netCDF4.Dataset(filename) as d:
            groups=list(d.groups.keys())

        for g in groups:
            print(g,flush=True)
            
            d=xarray.open_dataset(filename,group=g,decode_times=False)
            if 'nperiods' in g:
                gs=g.split('_')
                d['ndryperiods']=int(gs[gs.index('nperiods')+1])
                d['pH']=float(gs[gs.index('pH')+1])
                if 'Fescale' in gs:
                    d['Fescale']=float(gs[gs.index('Fescale')+1])
                else:
                    d['Fescale']=1.0
            
                newdims=['pH','ndryperiods','Fescale']
                d=d.expand_dims(newdims).set_coords(newdims)
            
            if 'noFeCH4Inhibition' in g:
                data_list_noinhib[g]=d
            elif 'noFe' in g:
                data_list_noFe[g]=d
            else:
                data_list[g]=d

    results_periodic=xarray.combine_by_coords(data_list[xx] for xx in data_list if 'pH' in xx and xx.startswith('organic'))
    results_periodic_mineral=xarray.combine_by_coords(data_list[xx] for xx in data_list if 'pH' in xx and xx.startswith('mineral'))
    results_periodic_noinhib=xarray.combine_by_coords(data_list_noinhib[xx] for xx in data_list_noinhib if 'pH' in xx and xx.startswith('organic'))
    results_periodic_noinhib_mineral=xarray.combine_by_coords(data_list_noinhib[xx] for xx in data_list_noinhib if 'pH' in xx and xx.startswith('mineral'))
    results_periodic_noFe=xarray.combine_by_coords(data_list_noFe[xx] for xx in data_list_noFe if 'pH' in xx and xx.startswith('organic'))
    results_periodic_noFe_mineral=xarray.combine_by_coords(data_list_noFe[xx] for xx in data_list_noFe if 'pH' in xx and xx.startswith('mineral'))


    # plot_timeseries(data_list,'_BD_Bockheim',ls=':')
    # # plot_timeseries(fname,'_BD_porosity')
    # plot_timeseries(data_list_noinhib,'_BD_Bockheim_noFeCH4Inhibition',ls='-',fig=figure('Time series'),cation_fig=figure('Cations'))
    # plot_timeseries(data_list_noFe,'_BD_Bockheim_noFe',ls='--',fig=figure('Time series'),cation_fig=figure('Cations'))

    cm=get_cmap('plasma')
    norm=matplotlib.colors.Normalize(0,len(results_periodic_mineral['Fescale']))
    plot_timeseries(data_list_noinhib,'_BD_Bockheim_noFeCH4Inhibition_Fescale_2.0',color=cm(norm(3)))
    plot_timeseries(data_list_noinhib,'_BD_Bockheim_noFeCH4Inhibition',fig=figure('Time series'),cation_fig=figure('Cations'),color=cm(norm(2)))
    plot_timeseries(data_list_noinhib,'_BD_Bockheim_noFeCH4Inhibition_Fescale_0.5',fig=figure('Time series'),cation_fig=figure('Cations'),color=cm(norm(1)))

    obs_mods=mod_obs_comp(data_list,'_BD_Bockheim')
    obs_mods_noinhib=mod_obs_comp(data_list_noinhib,'_BD_Bockheim_noFeCH4Inhibition')
    obs_mods_noFe=mod_obs_comp(data_list_noFe,'_BD_Bockheim_noFe')



    simnames={
       'organic_trough':'Trough organic',
       'organic_nottrough':'Rim organic',
       'mineral_trough':'Trough mineral',
       'mineral_nottrough':'Rim mineral',
       'highO2_organic':'Center organic (oxic)',
    }
    markers=['o','s','^','+','<','x']
    norm=matplotlib.colors.Normalize(0,max([obs_mods[sim].index.max() for sim in obs_mods]))
    cmap=plt.get_cmap()

    Fescale_unit=30.0 #nmol/L/s

    f,axs=subplots(4,2,clear=True,num='1-1 plots',figsize=(7,10))
    # f2,axs2=subplots(4,2,clear=True,num='1-1 plots no inhib',figsize=(4,12))
    Feax=axs[0,0]
    pHax=axs[0,1]
    CH4noinhib=axs[1,0]
    CH4inhib=axs[1,1]
    CH4noFe=axs[2,0]
    CO2ax=axs[2,1]
    DOMax=axs[3,0]
    orgacidax=axs[3,1]
    for num,sim in enumerate(obs_mods.keys()):
        CO2ax.scatter(obs_mods[sim]['CO2flux_obs'],obs_mods[sim]['CO2flux_mod'],marker=markers[num],s=10,c=obs_mods[sim].index,label=simnames[sim],norm=norm)
        CH4inhib.scatter(obs_mods[sim]['CH4flux_obs'],obs_mods[sim]['CH4flux_mod'],marker=markers[num],s=10,c=obs_mods[sim].index,label=simnames[sim],norm=norm)
        for nn in range(len(obs_mods[sim]['CO2flux_obs'])):
            CO2ax.errorbar(obs_mods[sim]['CO2flux_obs'].iloc[nn],obs_mods[sim]['CO2flux_mod'].iloc[nn],xerr=obs_mods[sim]['CO2flux_obs_err'].iloc[nn],marker='None',c=cmap(norm(obs_mods[sim].index[nn])))
            CH4inhib.errorbar(obs_mods[sim]['CH4flux_obs'].iloc[nn],obs_mods[sim]['CH4flux_mod'].iloc[nn],xerr=obs_mods[sim]['CH4flux_obs_err'].iloc[nn],marker='None',c=cmap(norm(obs_mods[sim].index[nn])))
            Feax.errorbar(obs_mods[sim]['Fe_II_obs'].iloc[nn],obs_mods[sim]['Fe_II_mod'].iloc[nn],xerr=obs_mods[sim]['Fe_II_obs_err'].iloc[nn],marker='None',c=cmap(norm(obs_mods[sim].index[nn])))
            DOMax.errorbar(obs_mods[sim]['DOM_obs'].iloc[nn],obs_mods[sim]['DOM_mod'].iloc[nn],xerr=obs_mods[sim]['DOM_obs_err'].iloc[nn],marker='None',c=cmap(norm(obs_mods[sim].index[nn])))
            orgacidax.errorbar(obs_mods[sim]['Acetate_obs'].iloc[nn],obs_mods[sim]['Acetate_mod'].iloc[nn],xerr=obs_mods[sim]['TOAC_obs_err'].iloc[nn],marker='None',c=cmap(norm(obs_mods[sim].index[nn])))
        DOMax.scatter(obs_mods[sim]['DOM_obs'],obs_mods[sim]['DOM_mod'],        marker=markers[num],s=10,c=obs_mods[sim].index,label=simnames[sim],norm=norm)
        orgacidax.scatter(obs_mods[sim]['Acetate_obs'],obs_mods[sim]['Acetate_mod'],marker=markers[num],s=10,c=obs_mods[sim].index,label=simnames[sim],norm=norm)
        Feax.scatter(obs_mods[sim]['Fe_II_obs'],obs_mods[sim]['Fe_II_mod'],    marker=markers[num],s=10,c=obs_mods[sim].index,label=simnames[sim],norm=norm)
        pHax.scatter(obs_mods[sim]['pH_obs'],obs_mods[sim]['pH_mod'],          marker=markers[num],s=10,c=obs_mods[sim].index,label=simnames[sim],norm=norm)

        CH4noinhib.scatter(obs_mods_noinhib[sim]['CH4flux_obs'],obs_mods_noinhib[sim]['CH4flux_mod'],marker=markers[num],s=10,c=obs_mods_noinhib[sim].index,label=simnames[sim],norm=norm)
        for nn in range(len(obs_mods[sim]['CO2flux_obs'])):
            CH4noinhib.errorbar(obs_mods_noinhib[sim]['CH4flux_obs'].iloc[nn],obs_mods_noinhib[sim]['CH4flux_mod'].iloc[nn],xerr=obs_mods_noinhib[sim]['CH4flux_obs_err'].iloc[nn],marker='None',c=cmap(norm(obs_mods[sim].index[nn])))
        CH4noFe.scatter(obs_mods_noFe[sim]['CH4flux_obs'],obs_mods_noFe[sim]['CH4flux_mod'],marker=markers[num],s=10,c=obs_mods_noFe[sim].index,label=simnames[sim],norm=norm)
        for nn in range(len(obs_mods[sim]['CO2flux_obs'])):
            CH4noFe.errorbar(obs_mods_noFe[sim]['CH4flux_obs'].iloc[nn],obs_mods_noFe[sim]['CH4flux_mod'].iloc[nn],xerr=obs_mods_noFe[sim]['CH4flux_obs_err'].iloc[nn],marker='None',c=cmap(norm(obs_mods[sim].index[nn])))


    CO2ax.set(title='CO$_2$ flux',xlabel='Obs CO$_2$ flux rate\n($\mu$mol g dwt$^{-1}$ day$^{-1}$)',ylabel='Mod CO$_2$ flux rate\n($\mu$mol g dwt$^{-1}$ day$^{-1}$)')
    CH4inhib.set(title='CH$_4$ flux (Fe(III) inhib)',xlabel='Obs CH$_4$ flux rate\n($\mu$mol g dwt$^{-1}$ day$^{-1}$)',ylabel='Mod CH$_4$ flux rate\n($\mu$mol g dwt$^{-1}$ day$^{-1}$)')
    CH4noinhib.set(title='CH$_4$ flux (no Fe(III) inhib)',xlabel='Obs CH$_4$ flux rate\n($\mu$mol g dwt$^{-1}$ day$^{-1}$)',ylabel='Mod CH$_4$ flux rate\n($\mu$mol g dwt$^{-1}$ day$^{-1}$)')
    CH4noFe.set(title='CH$_4$ flux (no Fe)',xlabel='Obs CH$_4$ flux rate\n($\mu$mol g dwt$^{-1}$ day$^{-1}$)',ylabel='Mod CH$_4$ flux rate\n($\mu$mol g dwt$^{-1}$ day$^{-1}$)')
    DOMax.set(title='DOM concentration',xlabel='Obs DOM (mol L$^{-1}$)',ylabel='Mod DOM (mol L$^{-1}$)')
    orgacidax.set(title='Organic acid concentration',xlabel='Obs organic acids (mol L$^{-1}$)',ylabel='Mod acetate (mol L$^{-1}$)')
    Feax.set(title='Fe(II) concentration',xlabel='Obs Fe(II) ($\mu$mol g dwt$^{-1}$)',ylabel='Mod Fe(II) ($\mu$mol g dwt$^{-1}$)')
    pHax.set(title='pH',xlabel='Obs pH',ylabel='Mod pH')
    axs[0,0].legend(loc='lower right')
    cb=colorbar(ax=axs[:],mappable=axs[0,1].collections[0])
    cb.set_label('Incubation day')
    for num,ax in enumerate(axs.ravel()):
        ax.axline((mean(ax.get_xlim()),mean(ax.get_xlim())),slope=1.0,c='k',lw=0.4,ls='--')
        ax.set_title('('+ascii_lowercase[num]+')',loc='left')

    DOMax.set(yscale='log',xscale='log')
    orgacidax.set(yscale='log',xscale='log')
    fitstats=getstats(pandas.concat(obs_mods))
    CO2ax.text(0.05,0.95,f"R$^2$ = {fitstats['R2']['CO2flux']:1.2f}",transform=CO2ax.transAxes,ha='left',va='top')
    CH4inhib.text(0.05,0.95,f"R$^2$ = {fitstats['R2']['CH4flux']:1.2f}",transform=CH4inhib.transAxes,ha='left',va='top')
    CH4noinhib.text(0.05,0.95,f"R$^2$ = {getstats(pandas.concat(obs_mods_noinhib))['R2']['CH4flux']:1.2f}",transform=CH4noinhib.transAxes,ha='left',va='top')
    CH4noFe.text(0.05,0.95,f"R$^2$ = {getstats(pandas.concat(obs_mods_noFe))['R2']['CH4flux']:1.2f}",transform=CH4noFe.transAxes,ha='left',va='top')
    Feax.text(0.05,0.95,f"R$^2$ = {fitstats['R2']['Fe_II']:1.2f}",transform=Feax.transAxes,ha='left',va='top')
    pHax.text(0.05,0.95,f"R$^2$ = {fitstats['R2']['pH']:1.2f}",transform=pHax.transAxes,ha='left',va='top')
    DOMax.text(0.05,0.95,f"R$^2$ = {fitstats['R2']['log(DOM)']:1.2f}",transform=DOMax.transAxes,ha='left',va='top')
    orgacidax.text(0.05,0.95,f"R$^2$ = {fitstats['R2']['log(Acetate)']:1.2f}",transform=orgacidax.transAxes,ha='left',va='top')

    # axs[3,1].set_visible(False)

    # results_periodic=xarray.open_dataset(fname,group='periodic_sims')

    # fig,axes=subplots(4,1,num='Periodic inundation',figsize=(6,8.4),clear=True)
    # # CH4ax=axes[0].twinx()
    # plot_result(results_periodic.squeeze().isel(ndryperiods=1,pH=2),Fe_ax=axes[2],CO2flux_ax=axes[0],CH4flux_ax=axes[1],pH_ax=axes[3],#porewater_ax=axes[3],
    #         gdrywt=True,SOC_pct=arctic.SOC_layermean['Organic'],cellulose_SOC_frac=arctic.cellulosefrac,BD=arctic.BD_layerest2_trough['Organic',False])
    # plot_result(results_periodic_noinhib.squeeze().isel(ndryperiods=1,pH=2),Fe_ax=axes[2],CO2flux_ax=axes[0],CH4flux_ax=axes[1],pH_ax=axes[3],#porewater_ax=axes[3],
    #         gdrywt=True,SOC_pct=arctic.SOC_layermean['Organic'],cellulose_SOC_frac=arctic.cellulosefrac,BD=arctic.BD_layerest2_trough['Organic',False])
    # # axes[0].set_ylim(bottom=-0.01)
    # axes[2].legend(labels=['Fe(OH)$_3$','Fe$^{3\!\!+}$','Fe$^{2\!\!+}$','Sorbed Fe'],loc=(0.83,0.2),fontsize='small')

    # axes[0].set_ylim(bottom=0.9e-11)
    # # axes[0].legend()
    # # CH4ax.legend()
    # # axes[3].legend(loc='lower left')
    # t=fig.axes[0].set_title('(a)',loc='left')  
    # t=fig.axes[1].set_title('(b)',loc='left')  
    # t=fig.axes[2].set_title('(c)',loc='left')  
    # t=fig.axes[3].set_title('(d)',loc='left')  

    cmap_O2=matplotlib.colors.LinearSegmentedColormap.from_list('cmap_O2',['blue','white'])
    # xmax=axes[0].get_xlim()[1]
    # for ax in axes:
    #     ax.pcolormesh(results_periodic['time'],ax.get_ylim(),results_periodic.isel(ndryperiods=1,pH=2)['Total O2(aq)'].to_masked_array()[:-1,None].T,cmap=cmap_O2,zorder=-1,shading='flat',alpha=0.15)
    #     ax.set_xlim(right=xmax,left=0)




    def cumulative_flux(data,var,pH,horizon,Fescale=1):
        return array([data.sel(ndryperiods=n,pH=pH,Fescale=Fescale)[var].isel(time=arctic.simlength*24-1) for n in data['ndryperiods'].values[1:]])*porosity/arctic.BD_layerest2_trough[horizon,False]

    cm=get_cmap('plasma')
    # norm=matplotlib.colors.Normalize(min(results_periodic['pH']),max(results_periodic['pH']))
    norm=matplotlib.colors.Normalize(0,len(results_periodic_mineral['Fescale']))
    porosity=results_periodic['Porosity'].isel(time=0).mean().item()
    pH=5.0
    fig,axs=subplots(ncols=2,nrows=2,num='Periodic by dry periods',clear=True,squeeze=False,figsize=(8,4))
    x=(150)/results_periodic_mineral['ndryperiods'].values[1:]
    for n,Fescale in enumerate(results_periodic['Fescale'].values):
        c=cm(norm(n))
        # axs[0,0].plot(results_periodic_noinhib['ndryperiods'].values[:],array([results_periodic_noinhib.sel(ndryperiods=n,pH=pH,Fescale=1)['Total Tracer'].isel(time=arctic.simlength*24-1) for n in results_periodic_noinhib['ndryperiods'].values[:]])*porosity/arctic.BD_layerest2_trough['Organic',False],'o-',label=Fescale,c=c)
        axs[0,0].plot(x,cumulative_flux(results_periodic_noinhib_mineral,'Total Tracer',pH,'Mineral',Fescale),'o-',c=c,label='%d nM s$^{-1}$'%int(Fescale*Fescale_unit))
        axs[0,1].plot(x,cumulative_flux(results_periodic_noinhib_mineral,'Total Tracer',pH,'Mineral',Fescale)-cumulative_flux(results_periodic_noFe_mineral,'Total Tracer',pH,'Mineral',Fescale),'o-',c=c,label='%d nM s$^{-1}$'%int(Fescale*Fescale_unit))
        axs[0,1].plot(x,cumulative_flux(results_periodic_mineral,'Total Tracer',pH,'Mineral',Fescale=Fescale)-cumulative_flux(results_periodic_noFe_mineral,'Total Tracer',pH,'Mineral',Fescale),'^:',c=c,label='%d nM s$^{-1}$'%int(Fescale*Fescale_unit))
        
        # axs[0,0].plot(results_periodic_noFe['ndryperiods'].values[:],array([results_periodic_noFe.sel(ndryperiods=n,pH=pH,Fescale=1)['Total Tracer'].isel(time=arctic.simlength*24-1) for n in results_periodic_noFe['ndryperiods'].values[:]])*porosity/arctic.BD_layerest2_trough['Organic',False],'+--',label=Fescale,c=c)
        # axs[0,0].plot(results_periodic_noFe_mineral['ndryperiods'].values[:],cumulative_flux(results_periodic_noFe_mineral,'Total Tracer',pH,'Mineral'),'+--',c=c)
        
        # axs[1,0].plot(results_periodic['ndryperiods'].values[:],array([results_periodic.sel(ndryperiods=n,pH=pH,Fescale=1)['Total CH4(aq)'].isel(time=arctic.simlength*24-1) for n in results_periodic['ndryperiods'].values[:]])*porosity/arctic.BD_layerest2_trough['Organic',False],'^:',label=Fescale,c=c)
        axs[1,0].plot(x,cumulative_flux(results_periodic_mineral,'Total CH4(aq)',pH,'Mineral',Fescale),'^:',c=c)
        axs[1,1].plot(x,cumulative_flux(results_periodic_mineral,'Total CH4(aq)',pH,'Mineral',Fescale)-cumulative_flux(results_periodic_noFe_mineral,'Total CH4(aq)',pH,'Mineral',Fescale),'^:',c=c)
        
        # axs[1,0].plot(results_periodic_noinhib['ndryperiods'].values[:],array([results_periodic_noinhib.sel(ndryperiods=n,pH=pH,Fescale=1)['Total CH4(aq)'].isel(time=arctic.simlength*24-1) for n in results_periodic_noinhib['ndryperiods'].values[:]])*porosity/arctic.BD_layerest2_trough['Organic',False],'o-',label=Fescale,c=c)
        axs[1,0].plot(x,cumulative_flux(results_periodic_noinhib_mineral,'Total CH4(aq)',pH,'Mineral',Fescale),'o-',c=c)
        axs[1,1].plot(x,cumulative_flux(results_periodic_noinhib_mineral,'Total CH4(aq)',pH,'Mineral',Fescale)-cumulative_flux(results_periodic_noFe_mineral,'Total CH4(aq)',pH,'Mineral',Fescale),'o-',c=c)

        # axs[1,0].plot(results_periodic_noFe['ndryperiods'].values[:],array([results_periodic_noFe.sel(ndryperiods=n,pH=pH,Fescale=1)['Total CH4(aq)'].isel(time=arctic.simlength*24-1) for n in results_periodic_noFe['ndryperiods'].values[:]])*porosity/arctic.BD_layerest2_trough['Organic',False],'+--',label=pH,c=c)
        # axs[1,0].plot(results_periodic_noFe_mineral['ndryperiods'].values[:],cumulative_flux(results_periodic_noFe_mineral,'Total CH4(aq)',pH,'Mineral'),'+--',c=cm(norm(Fescale)))


    axs[0,0].set_xlabel('Length of oxic-anoxic cycles (days)')
    axs[0,1].set_xlabel('Length of oxic-anoxic cycles (days)')
    axs[1,0].set_xlabel('Length of oxic-anoxic cycles (days)')
    axs[1,1].set_xlabel('Length of oxic-anoxic cycles (days)')

    axs[0,0].set_ylabel('Cumulative CO$_2$ flux\n(mmol g dwt$^{-1}$)')
    axs[0,1].set_ylabel('Cumulative CO$_2$ flux diff.\n(mmol g dwt$^{-1}$)')
    axs[1,0].set_ylabel('Cumulative CH$_4$ flux\n(mmol g dwt$^{-1}$)')
    axs[1,1].set_ylabel('Cumulative CH$_4$ flux diff.\n(mmol g dwt$^{-1}$)')

    axs[0,0].set_title('Cumulative CO$_2$ flux')
    axs[0,1].set_title('Fe effect on cumulative CO$_2$ flux')
    axs[1,0].set_title('Cumulative CH$_4$ flux')
    axs[1,1].set_title('Fe effect on cumulative CH$_4$ flux')

    
    # axs[0,0].set_ylim(bottom=0)
    # axs[0,1].set_ylim(bottom=0)
    axs[1,0].set_ylim(bottom=0)
    # axs[1,1].set_ylim(bottom=0)
    
    axs[0,0].legend(title='Fermentation rate')
    # axs[0,0].set_xticks(results_periodic['ndryperiods'].values[:])
    # axs[0,1].set_xticks(results_periodic['ndryperiods'].values[:])
    # axs[1,0].set_xticks(results_periodic['ndryperiods'].values[:])
    # axs[1,1].set_xticks(results_periodic['ndryperiods'].values[:])

    axs[1,0].legend(handles=[axs[1,1].lines[0],axs[1,1].lines[1]],labels=['Fe(III) inhibition','No inhibition'])

    axs[0,0].set_title('(a)',loc='left')
    axs[0,1].set_title('(b)',loc='left')
    axs[1,0].set_title('(c)',loc='left')
    axs[1,1].set_title('(d)',loc='left')

    axs[1,1].axhline(0,c='k',ls=':',lw=0.5)


    ndry_plotted=[4]
    # fig,axs=subplots(ncols=len(ndry_plotted),nrows=6,num='Periodic comparison organic',clear=True,figsize=(10,8.5),squeeze=False)
    fig=figure(num='Periodic comparison organic',clear=True,figsize=(10,8.5),constrained_layout=False)
    gs=matplotlib.gridspec.GridSpec(nrows=6,ncols=1,figure=fig,left=0.1,right=0.9,bottom=0.1,top=0.9,hspace=0.1)
    gs2=matplotlib.gridspec.GridSpecFromSubplotSpec(2,1,gs[1,0],hspace=0.1)

    axs=atleast_2d([fig.add_subplot(gs[0,0]),fig.add_subplot(gs2[1,0]),fig.add_subplot(gs[2,0]),
        fig.add_subplot(gs[3,0]),fig.add_subplot(gs[4,0]),fig.add_subplot(gs[5,0]),fig.add_subplot(gs2[0,0])]).T

    Fescales=[1,3]
    pH=5.0
    for nn,ndry in enumerate(ndry_plotted):
        
        # for nnn,Fescale in enumerate(results_periodic['Fescale'][1:]):
        for nnn in Fescales:
            Fescale=results_periodic['Fescale'][nnn]
            ls=':'
            c=cm(matplotlib.colors.Normalize(0,len(results_periodic_mineral['Fescale']))(nnn))
            d=results_periodic.sel(ndryperiods=ndry,pH=pH,Fescale=Fescale)
            plot_result(d,CH4flux_ax=axs[0,nn],gdrywt=True,CO2flux_ax=axs[1,nn],ls=ls,color=c,
                    BD=arctic.BD_layerest2_trough['Organic',False],SOC_pct=arctic.SOC_layermean['Organic'],CH4flux_args={},CO2flux_args={})
            plot_result(d,gdrywt=True,CO2flux_ax=axs[6,nn],ls=ls,color=c,
                    BD=arctic.BD_layerest2_trough['Organic',False],SOC_pct=arctic.SOC_layermean['Organic'],CH4flux_args={},CO2flux_args={})
            axs[2,nn].plot(results_periodic['time'][:-1],(d['Total Fe++'].diff(dim='time')*24*1e3*porosity+d['Total Sorbed Fe++'].diff(dim='time')*24*1e6/100**3)/arctic.BD_layerest2_trough['Organic',False],ls,c=c,label=str(ndry))

            axs[3,nn].plot(results_periodic['time'],d['Fe(OH)3 VF']/arctic.molar_volume_FeOH3*1e3/arctic.BD_layerest2_trough['Organic',False],ls,c=c)
            axs[4,nn].plot(results_periodic['time'],-log10(d['Free H+']),ls,c=c)
            axs[5,nn].plot(results_periodic['time'],d['Total Acetate-']*1e3,c=c,ls=ls)
            # axs[1,nn].plot(results_periodic['time'][:-1],diff(results_periodic['Total Tracer'].sel(ndryperiods=ndry,pH=pH)*results_periodic['Porosity'].sel(ndryperiods=ndry,pH=pH,time=0))/diff(results_periodic['time'])*1e3/arctic.BD_layerest2_trough['Organic',False],ls,c='b')
            CO2flux=diff(d['Total Tracer']*d['Porosity'].sel(time=0))/diff(results_periodic['time'])*1e3/arctic.BD_layerest2_trough['Organic',False]
            CH4flux=diff(d['Total CH4(aq)']*d['Porosity'].sel(time=0))/diff(results_periodic['time'])*1e3/arctic.BD_layerest2_trough['Organic',False]
            # axs[1,nn].plot(results_periodic['time'][:-1],CO2flux/CH4flux,ls=ls,c=c)

            ls='-'
            d=results_periodic_noinhib.sel(ndryperiods=ndry,pH=pH,Fescale=Fescale)
            plot_result(d,CH4flux_ax=axs[0,nn],gdrywt=True,CO2flux_ax=axs[1,nn],ls=ls,color=c,
                    BD=arctic.BD_layerest2_trough['Organic',False],SOC_pct=arctic.SOC_layermean['Organic'],CH4flux_args={},CO2flux_args={})
            plot_result(d,gdrywt=True,CO2flux_ax=axs[6,nn],ls=ls,color=c,
                    BD=arctic.BD_layerest2_trough['Organic',False],SOC_pct=arctic.SOC_layermean['Organic'],CH4flux_args={},CO2flux_args={})
            axs[2,nn].plot(d['time'][:-1],(d['Total Fe++'].diff(dim='time')*24*1e3*porosity+d['Total Sorbed Fe++'].diff(dim='time')*24*1e6/100**3)/arctic.BD_layerest2_trough['Organic',False],ls,c=c,label=str(ndry))

            axs[3,nn].plot(results_periodic_noinhib['time'],d['Fe(OH)3 VF']/arctic.molar_volume_FeOH3*1e3/arctic.BD_layerest2_trough['Organic',False],ls,c=c)
            axs[4,nn].plot(results_periodic_noinhib['time'],-log10(d['Free H+']),ls,c=c)
            axs[5,nn].plot(results_periodic_noinhib['time'],d['Total Acetate-']*1e3,c=c,ls=ls)
            # axs[1,nn].plot(results_periodic_noinhib['time'][:-1],diff(results_periodic_noinhib['Total Tracer'].sel(ndryperiods=ndry,pH=pH)*results_periodic_noinhib['Porosity'].sel(ndryperiods=ndry,pH=pH,time=0))/diff(results_periodic_noinhib['time'])*1e3/arctic.BD_layerest2_trough['Organic',False],ls,c='b')
            CO2flux=diff(d['Total Tracer']*d['Porosity'].sel(time=0))/diff(results_periodic_noinhib['time'])*1e3/arctic.BD_layerest2_trough['Organic',False]
            CH4flux=diff(d['Total CH4(aq)']*d['Porosity'].sel(time=0))/diff(results_periodic_noinhib['time'])*1e3/arctic.BD_layerest2_trough['Organic',False]
            # axs[1,nn].plot(results_periodic_noinhib['time'][:-1],CO2flux/CH4flux,ls=ls,c=c)

            ls='--'
            d=results_periodic_noFe.sel(ndryperiods=ndry,pH=pH,Fescale=Fescale)
            plot_result(d,CH4flux_ax=axs[0,nn],gdrywt=True,CO2flux_ax=axs[1,nn],ls=ls,color=c,
                    BD=arctic.BD_layerest2_trough['Organic',False],SOC_pct=arctic.SOC_layermean['Organic'],CH4flux_args={},CO2flux_args={})
            plot_result(d,gdrywt=True,CO2flux_ax=axs[6,nn],ls=ls,color=c,
                    BD=arctic.BD_layerest2_trough['Organic',False],SOC_pct=arctic.SOC_layermean['Organic'],CH4flux_args={},CO2flux_args={})
            axs[2,nn].plot(results_periodic_noFe['time'][:-1],(d['Total Fe++'].diff(dim='time')*24*1e3*porosity+d['Total Sorbed Fe++'].diff(dim='time')*24*1e6/100**3)/arctic.BD_layerest2_trough['Organic',False],ls,c=c,label=str(ndry))

            axs[3,nn].plot(results_periodic_noFe['time'],d['Fe(OH)3 VF']/arctic.molar_volume_FeOH3*1e3/arctic.BD_layerest2_trough['Organic',False],ls,c=c)
            axs[4,nn].plot(results_periodic_noFe['time'],-log10(d['Free H+']),ls,c=c)
            axs[5,nn].plot(results_periodic_noFe['time'],d['Total Acetate-']*1e3,c=c,ls=ls)
            # axs[1,nn].plot(results_periodic_noFe['time'][:-1],diff(results_periodic_noFe['Total Tracer'].sel(ndryperiods=ndry,pH=pH)*results_periodic_noFe['Porosity'].sel(ndryperiods=ndry,pH=pH,time=0))/diff(results_periodic_noFe['time'])*1e3/arctic.BD_layerest2_trough['Organic',False],ls,c='b')
            CO2flux=diff(d['Total Tracer']*d['Porosity'].sel(time=0))/diff(results_periodic_noFe['time'])*1e3/arctic.BD_layerest2_trough['Organic',False]
            CH4flux=diff(d['Total CH4(aq)']*d['Porosity'].sel(time=0))/diff(results_periodic_noFe['time'])*1e3/arctic.BD_layerest2_trough['Organic',False]
            # axs[1,nn].plot(results_periodic_noFe['time'][:-1],CO2flux/CH4flux,ls=ls,c=c)



        axs[0,nn].set_xlim(0,arctic.simlength)
        # axs[1,nn].set_xlim(0,arctic.simlength)
        # axs[1,nn].set_title('CO$_2$:CH$_4$ flux ratio')
        # axs[1,nn].set_ylabel('CO$_2$:CH$_4$ ratio')
        # axs[1,nn].set_ylim(-2,300)
        # axs[1,nn].set_xlabel('Time (days)')
        axs[0,nn].set_title('CH$_4$ flux rate')

        axs[1,nn].set_ylim(0,18)
        axs[1,nn].set_ylabel(axs[1,nn].get_ylabel(),y=1.0)
        axs[6,nn].set_ylim(40,150)
        axs[6,nn].set_ylabel('')
        axs[6,nn].get_xaxis().set_visible(False)
        axs[6,nn].spines['bottom'].set_linestyle('--')
        axs[1,nn].spines['top'].set_linestyle('--')

        axs[2,nn].set_ylim(bottom=0)
        axs[2,nn].set_xlabel('Time (days)')
        axs[2,nn].set_ylabel('Fe(II) production rate\n($\mu$mol g dwt$^{-1}$ day$^{-1}$)')
        axs[2,nn].set_title('Fe(II) production rate')

        axs[3,nn].set_xlabel('Time (days)')
        axs[3,nn].set_ylabel('Fe oxide \n(mmol g dwt$^{-1}$)')
        axs[3,nn].set_title('Reactive Fe oxide minerals')

        axs[4,nn].set_xlabel('Time (days)')
        axs[4,nn].set_ylabel('pH')
        axs[4,nn].set_title('pH')

        axs[5,nn].set_ylim(bottom=1e-5)
        axs[5,nn].set_xlabel('Time (days)')
        axs[5,nn].set_ylabel('Concentration (mM)')
        axs[5,nn].set_title('Acetate concentration')

        for x in range(axs.shape[0]):
            # for num in range(ndry):
                # axs[x,nn].axvspan(arctic.simlength*arctic.oxicfrac/ndry+num*arctic.simlength/ndry,(num+1)*arctic.simlength/ndry,color='b',alpha=0.1,zorder=-1)
            axs[x,nn].pcolormesh(results_periodic['time'],axs[x,nn].get_ylim(),d['Total O2(aq)'].to_masked_array()[:-1,None].T,cmap=cmap_O2,zorder=-1,shading='flat',alpha=0.15)
            axs[x,nn].set_xlim(0,arctic.simlength)
            axs[x,nn].set_title('('+ascii_lowercase[x*len(ndry_plotted)+nn]+')',loc='left')
        axs[6,0].set_title('(b)',loc='left')
            
    l1=axs[2,0].legend(labels=[int(x*Fescale_unit) for x in results_periodic['Fescale'][Fescales].values],
                       handles=(axs[0,0].lines[1],axs[0,0].lines[4]),
                       title='Ferm. rate (nM s$^{-1}$)',ncol=3,loc=(0.1,1.05))
    l2=axs[2,0].legend(labels=['Fe(III) inhibition','No inhibition','No Fe'],title='Fe(III) inhibition of CH$_4$',
                        handles=[axs[0,0].lines[0],axs[0,0].lines[1],axs[0,0].lines[2]],ncol=3,loc=(0.565,1.05))
    axs[2,0].add_artist(l1)
    gs.update(hspace=0.85,bottom=0.05,top=0.96,left=0.08,right=0.98)



    # fig,axs=subplots(ncols=len(ndry_plotted),nrows=6,num='Periodic comparison mineral',clear=True,figsize=(10,8.5),squeeze=False)
    fig=figure(num='Periodic comparison mineral',clear=True,figsize=(10,8.5),constrained_layout=False)
    gs=matplotlib.gridspec.GridSpec(nrows=6,ncols=1,figure=fig,left=0.1,right=0.9,bottom=0.1,top=0.9,hspace=0.1)
    gs2=matplotlib.gridspec.GridSpecFromSubplotSpec(2,1,gs[1,0],hspace=0.1)

    axs=atleast_2d([fig.add_subplot(gs[0,0]),fig.add_subplot(gs2[1,0]),fig.add_subplot(gs[2,0]),
        fig.add_subplot(gs[3,0]),fig.add_subplot(gs[4,0]),fig.add_subplot(gs[5,0]),fig.add_subplot(gs2[0,0])]).T

    from string import ascii_lowercase

    pH=5.0
    for nn,ndry in enumerate(ndry_plotted):
        # for nnn,Fescale in enumerate(results_periodic_mineral['Fescale'][1:]):
        for nnn in Fescales:
            Fescale=results_periodic_mineral['Fescale'][nnn]
            ls=':'
            c=cm(matplotlib.colors.Normalize(0,len(results_periodic_mineral['Fescale']))(nnn))
            d=results_periodic_mineral.sel(ndryperiods=ndry,pH=pH,Fescale=Fescale)
            plot_result(d,CH4flux_ax=axs[0,nn],gdrywt=True,CO2flux_ax=axs[1,nn],ls=ls,color=c,
                    BD=arctic.BD_layerest2_trough['Mineral',False],SOC_pct=arctic.SOC_layermean['Mineral'],CH4flux_args={},CO2flux_args={})
            plot_result(d,gdrywt=True,CO2flux_ax=axs[6,nn],ls=ls,color=c,
                    BD=arctic.BD_layerest2_trough['Mineral',False],SOC_pct=arctic.SOC_layermean['Mineral'],CH4flux_args={},CO2flux_args={})
            axs[2,nn].plot(d['time'][:-1],(d['Total Fe++'].diff(dim='time')*24*1e3*porosity+d['Total Sorbed Fe++'].diff(dim='time')*24*1e6/100**3)/arctic.BD_layerest2_trough['Mineral',False],ls,c=c,label=str(ndry))

            axs[3,nn].plot(d['time'],d['Fe(OH)3 VF']/arctic.molar_volume_FeOH3*1e3/arctic.BD_layerest2_trough['Mineral',False],ls,c=c)
            axs[4,nn].plot(d['time'],-log10(d['Free H+']),ls,c=c)
            axs[5,nn].plot(d['time'],d['Total Acetate-']*1e3,c=c,ls=ls)
            # axs[1,nn].plot(results_periodic_mineral['time'][:-1],diff(results_periodic_mineral['Total Tracer'].sel(ndryperiods=ndry,pH=pH)*results_periodic_mineral['Porosity'].sel(ndryperiods=ndry,pH=pH,time=0))/diff(results_periodic_mineral['time'])*1e3/arctic.BD_layerest2_trough['Mineral',False],ls,c='b')
            CO2flux=diff(d['Total Tracer']*d['Porosity'])/diff(d['time'])*1e3/arctic.BD_layerest2_trough['Mineral',False]
            CH4flux=diff(d['Total CH4(aq)']*d['Porosity'])/diff(d['time'])*1e3/arctic.BD_layerest2_trough['Mineral',False]
            # axs[1,nn].plot(results_periodic_mineral['time'][:-1],CO2flux/CH4flux,ls=ls,c=c)

            ls='-'
            d=results_periodic_noinhib_mineral.sel(ndryperiods=ndry,pH=pH,Fescale=Fescale)
            plot_result(d,CH4flux_ax=axs[0,nn],gdrywt=True,CO2flux_ax=axs[1,nn],ls=ls,color=c,
                    BD=arctic.BD_layerest2_trough['Mineral',False],SOC_pct=arctic.SOC_layermean['Mineral'],CH4flux_args={},CO2flux_args={})
            plot_result(d,gdrywt=True,CO2flux_ax=axs[6,nn],ls=ls,color=c,
                    BD=arctic.BD_layerest2_trough['Mineral',False],SOC_pct=arctic.SOC_layermean['Mineral'],CH4flux_args={},CO2flux_args={})
            axs[2,nn].plot(d['time'][:-1],(d['Total Fe++'].diff(dim='time')*24*1e3*porosity+d['Total Sorbed Fe++'].diff(dim='time')*24*1e6/100**3)/arctic.BD_layerest2_trough['Mineral',False],ls,c=c,label=str(ndry))

            axs[3,nn].plot(d['time'],d['Fe(OH)3 VF']/arctic.molar_volume_FeOH3*1e3/arctic.BD_layerest2_trough['Mineral',False],ls,c=c)
            axs[4,nn].plot(d['time'],-log10(d['Free H+']),ls,c=c)
            axs[5,nn].plot(d['time'],d['Total Acetate-']*1e3,c=c,ls=ls)
            # axs[1,nn].plot(results_periodic_mineral['time'][:-1],diff(results_periodic_mineral['Total Tracer'].sel(ndryperiods=ndry,pH=pH)*results_periodic_mineral['Porosity'].sel(ndryperiods=ndry,pH=pH,time=0))/diff(results_periodic_mineral['time'])*1e3/arctic.BD_layerest2_trough['Mineral',False],ls,c='b')
            CO2flux=diff(d['Total Tracer']*d['Porosity'])/diff(d['time'])*1e3/arctic.BD_layerest2_trough['Mineral',False]
            CH4flux=diff(d['Total CH4(aq)']*d['Porosity'])/diff(d['time'])*1e3/arctic.BD_layerest2_trough['Mineral',False]
            # axs[1,nn].plot(results_periodic_mineral['time'][:-1],CO2flux/CH4flux,ls=ls,c=c)

            ls='--'
            d=results_periodic_noFe_mineral.sel(ndryperiods=ndry,pH=pH,Fescale=Fescale)
            plot_result(d,CH4flux_ax=axs[0,nn],gdrywt=True,CO2flux_ax=axs[1,nn],ls=ls,color=c,
                    BD=arctic.BD_layerest2_trough['Mineral',False],SOC_pct=arctic.SOC_layermean['Mineral'],CH4flux_args={},CO2flux_args={})
            plot_result(d,gdrywt=True,CO2flux_ax=axs[6,nn],ls=ls,color=c,
                    BD=arctic.BD_layerest2_trough['Mineral',False],SOC_pct=arctic.SOC_layermean['Mineral'],CH4flux_args={},CO2flux_args={})        
            axs[2,nn].plot(d['time'][:-1],(d['Total Fe++'].diff(dim='time')*24*1e3*porosity+d['Total Sorbed Fe++'].diff(dim='time')*24*1e6/100**3)/arctic.BD_layerest2_trough['Mineral',False],ls,c=c,label=str(ndry))

            axs[3,nn].plot(d['time'],d['Fe(OH)3 VF']/arctic.molar_volume_FeOH3*1e3/arctic.BD_layerest2_trough['Mineral',False],ls,c=c)
            axs[4,nn].plot(d['time'],-log10(d['Free H+']),ls,c=c)
            axs[5,nn].plot(d['time'],d['Total Acetate-']*1e3,c=c,ls=ls)
            # axs[1,nn].plot(results_periodic_mineral['time'][:-1],diff(results_periodic_mineral['Total Tracer'].sel(ndryperiods=ndry,pH=pH)*results_periodic_mineral['Porosity'].sel(ndryperiods=ndry,pH=pH,time=0))/diff(results_periodic_mineral['time'])*1e3/arctic.BD_layerest2_trough['Mineral',False],ls,c='b')
            CO2flux=diff(d['Total Tracer']*d['Porosity'])/diff(d['time'])*1e3/arctic.BD_layerest2_trough['Mineral',False]
            CH4flux=diff(d['Total CH4(aq)']*d['Porosity'])/diff(d['time'])*1e3/arctic.BD_layerest2_trough['Mineral',False]
            # axs[1,nn].plot(results_periodic_mineral['time'][:-1],CO2flux/CH4flux,ls=ls,c=c)


        axs[0,nn].set_xlim(0,arctic.simlength)
        # axs[1,nn].set_xlim(0,arctic.simlength)
        # axs[1,nn].set_title('CO$_2$:CH$_4$ flux ratio')
        # axs[1,nn].set_ylabel('CO$_2$:CH$_4$ ratio')
        # axs[1,nn].set_ylim(-2,300)
        # axs[1,nn].set_xlabel('Time (days)')
        axs[0,nn].set_title('CH$_4$ flux rate')

        axs[1,nn].set_ylim(0,3.5)
        axs[1,nn].set_ylabel(axs[1,nn].get_ylabel(),y=1.0)
        axs[6,nn].set_ylim(12,38)
        axs[6,nn].set_ylabel('')
        axs[6,nn].get_xaxis().set_visible(False)
        axs[6,nn].spines['bottom'].set_linestyle('--')
        axs[1,nn].spines['top'].set_linestyle('--')

        axs[2,nn].set_ylim(bottom=0)
        axs[2,nn].set_xlabel('Time (days)')
        axs[2,nn].set_ylabel('Fe(II) prod. rate\n($\mu$mol g dwt$^{-1}$ day$^{-1}$)')
        axs[2,nn].set_title('Fe(II) production rate')

        axs[3,nn].set_xlabel('Time (days)')
        axs[3,nn].set_ylabel('Fe oxide \n(mmol g dwt$^{-1}$)')
        axs[3,nn].set_title('Reactive Fe oxide minerals')

        axs[4,nn].set_xlabel('Time (days)')
        axs[4,nn].set_ylabel('pH')
        axs[4,nn].set_title('pH')

        axs[5,nn].set_ylim(bottom=1e-5)
        axs[5,nn].set_xlabel('Time (days)')
        axs[5,nn].set_ylabel('Concentration (mM)')
        axs[5,nn].set_title('Acetate concentration')

        for x in range(axs.shape[0]):
            # for num in range(ndry):
                # axs[x,nn].axvspan(arctic.simlength*arctic.oxicfrac/ndry+num*arctic.simlength/ndry,(num+1)*arctic.simlength/ndry,color='b',alpha=0.1,zorder=-1)
            axs[x,nn].pcolormesh(results_periodic['time'],axs[x,nn].get_ylim(),d['Total O2(aq)'].to_masked_array()[:-1,None].T,cmap=cmap_O2,zorder=-1,shading='flat',alpha=0.15)
            axs[x,nn].set_xlim(0,arctic.simlength)
            axs[x,nn].set_title('('+ascii_lowercase[x*len(ndry_plotted)+nn]+')',loc='left')
        axs[6,0].set_title('(b)',loc='left')
            
    l1=axs[2,0].legend(labels=[int(x*Fescale_unit) for x in results_periodic['Fescale'][Fescales].values],
                       handles=(axs[0,0].lines[1],axs[0,0].lines[4]),
                       title='Ferm. rate (nM s$^{-1}$)',ncol=3,loc=(0.1,1.05))
    l2=axs[2,0].legend(labels=['Fe(III) inhibition','No inhibition','No Fe'],title='Fe(III) inhibition of CH$_4$',
                        handles=[axs[0,0].lines[0],axs[0,0].lines[1],axs[0,0].lines[2]],ncol=3,loc=(0.7,1.05))
    axs[2,0].add_artist(l1)
    gs.update(hspace=0.85,bottom=0.05,top=0.96,left=0.08,right=0.98)


    f,a=subplots(num='DOM concentrations',nrows=3,clear=True)
    results_periodic_noinhib_mineral['Total Acetate-'].sel(ndryperiods=3,pH=5.0,Fescale=0.5).plot(ls='-',label='Acetate',color='C0',ax=a[0])
    results_periodic_noinhib_mineral['Total Acetate-'].sel(ndryperiods=3,pH=5.0,Fescale=1.0).plot(ls='-',color='C0',lw=0.5,ax=a[0])
    results_periodic_noFe_mineral['Total Acetate-'].sel(ndryperiods=3,pH=5.0,Fescale=0.5).plot(ls='--',label='Acetate noFe',color='C0',ax=a[0])
    results_periodic_noFe_mineral['Total Acetate-'].sel(ndryperiods=3,pH=5.0,Fescale=1.0).plot(ls='--',color='C0',lw=0.5,ax=a[0])

    results_periodic_noinhib_mineral['Total DOM1'].sel(ndryperiods=3,pH=5.0,Fescale=0.5).plot(ls='-',label='DOM1',color='C1',ax=a[0])
    results_periodic_noinhib_mineral['Total DOM1'].sel(ndryperiods=3,pH=5.0,Fescale=1.0).plot(ls='-',color='C1',lw=0.5,ax=a[0])
    results_periodic_noFe_mineral['Total DOM1'].sel(ndryperiods=3,pH=5.0,Fescale=0.5).plot(ls='--',label='DOM1 noFe',color='C1',ax=a[0])
    results_periodic_noFe_mineral['Total DOM1'].sel(ndryperiods=3,pH=5.0,Fescale=1.0).plot(ls='--',color='C1',lw=0.5,ax=a[0])

    results_periodic_noinhib_mineral['Total Sorbed cellulose'].sel(ndryperiods=3,pH=5.0,Fescale=0.5).plot(ls='-',label='DOM1',color='C1',ax=a[1])
    results_periodic_noinhib_mineral['Total Sorbed cellulose'].sel(ndryperiods=3,pH=5.0,Fescale=1.0).plot(ls='-',color='C1',lw=0.5,ax=a[1])
    results_periodic_noFe_mineral['Total Sorbed cellulose'].sel(ndryperiods=3,pH=5.0,Fescale=0.5).plot(ls='--',label='DOM1 noFe',color='C1',ax=a[1])
    results_periodic_noFe_mineral['Total Sorbed cellulose'].sel(ndryperiods=3,pH=5.0,Fescale=1.0).plot(ls='--',color='C1',lw=0.5,ax=a[1])

    results_periodic_noinhib_mineral['Total Tracer'].sel(ndryperiods=3,pH=5.0,Fescale=0.5).plot(ls='-',label='DOM1',color='C1',ax=a[2])
    results_periodic_noinhib_mineral['Total Tracer'].sel(ndryperiods=3,pH=5.0,Fescale=1.0).plot(ls='-',color='C1',lw=0.5,ax=a[2])
    results_periodic_noFe_mineral['Total Tracer'].sel(ndryperiods=3,pH=5.0,Fescale=0.5).plot(ls='--',label='DOM1 noFe',color='C1',ax=a[2])
    results_periodic_noFe_mineral['Total Tracer'].sel(ndryperiods=3,pH=5.0,Fescale=1.0).plot(ls='--',color='C1',lw=0.5,ax=a[2])


    a[0].legend()

    # Comparing fermentation rates with organic acid accumulation from Herndon et al. 2015
    # Convert from umol C/g SOC/30 days to nmol C/L/s
    print(' Fermentation rate values from Herndon et al. 2015, Figure 5') 
    print('organic center',200e-6/(30*24*3600)*arctic.BD_organic_nottrough/arctic.SOC_organic_nottrough*1000*arctic.porosity_organic_nottrough*1e9)
    print('mineral center',20e-6/(30*24*3600)*arctic.BD_mineral_nottrough/arctic.SOC_mineral_nottrough*1000*arctic.porosity_mineral_nottrough*1e9)
    print('organic trough',75e-6/(30*24*3600)*arctic.BD_organic_trough/arctic.SOC_organic_trough*1000*arctic.porosity_organic_trough*1e9)
    print('mineral trough',10e-6/(30*24*3600)*arctic.BD_mineral_trough/arctic.SOC_mineral_trough*1000*arctic.porosity_mineral_trough*1e9)
    
    # Values from "Acetate production" in Table 4
    print(' Fermentation rate values from Herndon et al. 2015, Table 4 "Acetate production"') 
    print('organic center',310e-6/(30*24*3600)*arctic.BD_organic_nottrough/arctic.SOC_organic_nottrough*1000*arctic.porosity_organic_nottrough*1e9)
    print('mineral center',60e-6/(30*24*3600)*arctic.BD_mineral_nottrough/arctic.SOC_mineral_nottrough*1000*arctic.porosity_mineral_nottrough*1e9)
    print('organic trough',66e-6/(30*24*3600)*arctic.BD_organic_trough/arctic.SOC_organic_trough*1000*arctic.porosity_organic_trough*1e9)
    print('mineral trough',40e-6/(30*24*3600)*arctic.BD_mineral_trough/arctic.SOC_mineral_trough*1000*arctic.porosity_mineral_trough*1e9)
    
    save_all_figs('Arctic_Fe_output/figs/2022-09-08',format='jpg',dpi=300)

    show()
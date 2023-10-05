import matplotlib
import matplotlib.pyplot as plt
# matplotlib.use('Agg')
import xarray
import sys
import pandas
import numpy as np

filenames=sys.argv[1:-1]

if sys.argv[-1].endswith('.nc'):
    filenames.append(sys.argv[-1])
    figdir='Mn_output/Figures'
else:
    figdir=sys.argv[-1]


# result_files = [
# 'Mn_ph35_saved_2020-08-31.nc',
# 'Mn_ph40_saved_2020-08-31.nc',
#   'Mn_ph45_saved_2020-08-31.nc',
#   'Mn_ph50_saved_2020-08-31.nc',
#   'Mn_ph55_saved_2020-08-31.nc',
# ]

leafC_mass_conv=1.0 # For simulations after I started applying the C to dry mass conversion in the sims. Otherwise it should be 0.4
molar_volume_birnessite = 251.1700 # cm3/mol Divide by 7 below because birnessite is defined in database as 7 mols of Mn
Mn_molarmass=54.94        #g/mol


def getval_from_str(s,subs,valtype=float):
    return valtype(s[s.index(subs)+len(subs):].partition('_')[0])

def getgroups(filename):
    import netCDF4
    with netCDF4.Dataset(filename) as d:
        groups=list(d.groups.keys())
    return groups

def load_data(filenames,chunks={},readgroups=False):
    data_list=[]
    incubation_list=[]

    allfilenames=[]
    groupnames=[]

    from Manganese_profile import setup_sims,make_gname,molar_mass
    expected_sims=setup_sims()
    expected_sims['simnum']=expected_sims.index

    for filenum,filename in enumerate(sorted(filenames)):
        print('Reading file %s'%filename)

        if readgroups:
            groups=getgroups(filename)
        else:
            groups=[]
            for n in range(filenum,len(expected_sims),len(filenames)):
                name=make_gname(expected_sims['pH'][n],
                        expected_sims['Ndep'][n]*1000/molar_mass['N']/100**2/(365*24*3600),
                        expected_sims['warming'][n],expected_sims['anox_freq'][n],expected_sims['anox_len'][n],
                        expected_sims['birn_rate'][n])
                groups.append(name )
                # if (expected_sims['Ndep'][n]==0) and (expected_sims['warming'][n]==0) and (expected_sims['anox_len'][n]==0.5):
                #     groups.append(name+'_incubation')

        for g in groups:
            print(g,flush=True)

            d=xarray.open_dataset(filename,group=g,chunks=chunks,decode_times=False)
            if 'soil_pH' not in d:
                d['soil_pH']=getval_from_str(g,'pH')
                d['Ndep']=getval_from_str(g,'Ndep',int)
                d['warming']=getval_from_str(g,'warming')
                d['redox_cycles']=getval_from_str(g,'anox_freq',int)
                d['anox_lenscales']=getval_from_str(g,'anox_len')
                newdims=['soil_pH','Ndep','warming','redox_cycles']
                d=d.expand_dims(newdims).set_coords(newdims)
            if 'birnrate' not in d:
                d['birnrate']=getval_from_str(g,'birnrate')
                d=d.expand_dims('birnrate').set_coords('birnrate')
            if 'incubation' in g:
                incubation_list.append(d)
            else:
                data_list.append(d)
                allfilenames.append(filename)
                groupnames.append(g)

    pH_sims=np.array([d['soil_pH'].item() for d in data_list])
    Ndep_sims=np.array([d['Ndep'].item() for d in data_list])
    warming_sims=np.array([d['warming'].item() for d in data_list])
    anox_freq_sims=np.array([d['redox_cycles'].item() for d in data_list])
    anox_len_sims=np.array([d['anox_lenscales'].item() for d in data_list])
    birnrate_sims=np.array([d['birnrate'].item() for d in data_list])
    allsims=pandas.DataFrame({'Ndep':Ndep_sims,'warming':warming_sims,'anox_freq':anox_freq_sims,'anox_len':anox_len_sims,'pH':pH_sims,'birnrate':birnrate_sims})

    comp=allsims.merge(expected_sims,how='outer',indicator=True)
    missing=comp[comp['_merge']!='both']
    if len(missing)>0:
        print('Warning: Not all expected sims were found:')
        print(missing)
        print('Generating all nan data for missing simulations.')
        for num in missing.index:
            d=(data_list[0]*np.nan).assign_coords(
                soil_pH=np.atleast_1d(missing.loc[num]['pH']),
                Ndep=np.atleast_1d(missing.loc[num]['Ndep']),
                warming=np.atleast_1d(missing.loc[num]['warming']),
                redox_cycles=np.atleast_1d(missing.loc[num]['anox_freq']),
                anox_lenscales=np.atleast_1d(missing.loc[num]['anox_len']))
            data_list.append(d)

        pH_sims=np.array([d['soil_pH'].item() for d in data_list])
        Ndep_sims=np.array([d['Ndep'].item() for d in data_list])
        warming_sims=np.array([d['warming'].item() for d in data_list])
        anox_freq_sims=np.array([d['redox_cycles'].item() for d in data_list])
        anox_len_sims=np.array([d['anox_lenscales'].item() for d in data_list])
        birnrate_sims=np.array([d['birnrate'].item() for d in data_list])

    x=(anox_len_sims==1.0)&(anox_freq_sims>0)
    # if len(pH_sims[x])%len(np.unique(pH_sims)) != 0:
    #     raise ValueError("Simulations don't appear to divide evenly by pH and birnessite rate!")
    if not x.any():
        raise ValueError('No simulations meet criteria for alldata')
    alldata=xarray.combine_by_coords(data_list[xx] for xx in np.nonzero(x)[0])
    x=(Ndep_sims==0)&(warming_sims==0)&(birnrate_sims.round(5)==4.0)&(pH_sims==5.0)
    if not x.any():
        raise ValueError('No simulations meet criteria for redoxdata')
    redoxdata=xarray.combine_by_coords(data_list[xx] for xx in np.nonzero(x)[0])

    print('Finished combining data',flush=True)
    return alldata,redoxdata,incubation_list

if len(filenames)>0:
    alldata,redoxdata,incubation_list=load_data(filenames,chunks={},readgroups=False) # Use chunks={'time':3650} for long datasets
controldata=alldata.isel(Ndep=0,warming=0)

def getdata(soil_pH,Ndep,warming,redox_cycles=50,anox_lenscale=0.5):
    # return data_list[flatnonzero((pH_sims==soil_pH)&(Ndep_sims==Ndep)&(warming_sims==warming)&(anox_freq_sims==redox_cycles)&(anox_len_sims==anox_lenscale))[0]].squeeze()
    return alldata.sel(soil_pH=soil_pH,Ndep=Ndep,warming=warming,anox_lenscales=anox_lenscale,redox_cycles=redox_cycles)

def save_one_fig(f,dirname=figdir,format='png',**kwargs):
    fname_fixed=f.get_label().replace('/','-')
    savename='{dirname:s}/{fname:s}.{format}'.format(dirname=dirname,format=format,fname=fname_fixed)
    print(savename)
    f.savefig(savename,**kwargs)

def letter_label(ax,label=None,x=0.03,y=1.03,title=True,**kwargs):
    from string import ascii_lowercase
    if isinstance(label,int):
        label='('+ascii_lowercase[label]+')'
    if title:
        return ax.set_title(label,loc='left')
    else:
        return ax.text(x,y,label,transform=ax.transAxes,**kwargs)

dt=(controldata['time'][1]-controldata['time'][0]).item()
oneyr=int(365/dt)

anox_baseline_i=0
anox_baseline=redoxdata['anox_lenscales'][anox_baseline_i].load().item()
birnrate_baseline_i=2

# Multiplying by surface area of 100 m2/m3 and converting to units of mmol/m3/year
birnrate_labels=[str(num) for num in (1e-12*2.0**np.arange(-2,8,2)*100*1e3*3600*24*365).round(1)]

f,axs=plt.subplots(num='Figure S1 Saturated time fraction',clear=True)
O2=redoxdata['Total O2(aq)'].sel(soil_pH=5.0).squeeze()
satfrac=(O2<O2.max()*0.9).mean(dim='time').T
h=axs.pcolormesh(satfrac['anox_lenscales'],-satfrac['depth'],satfrac*100,cmap='RdBu',vmin=0,vmax=100,shading='auto')
cb=f.colorbar(h,ax=axs)
cb.set_label('Saturated fraction of time (%)')
cb.set_ticks([0,25,50,75,100])
axs.set(title='Water-saturated fraction of time',xlabel='Drainage time scale (days)',ylabel='Depth (cm)')

save_one_fig(f)

total_litter_all=((alldata['Total Sorbed Cellulose']+alldata['Total Sorbed Lignin']).coarsen(time=oneyr,boundary='trim').mean()*12e-3*(alldata.z_bottom-alldata.z_top)/100).sum(dim='depth')

results_byph=controldata.isel(birnrate=birnrate_baseline_i)
results_bybirnrate=controldata.squeeze()
pHs=controldata['soil_pH'].to_masked_array()


total_MAOM=((results_bybirnrate['Total DOM3'].coarsen(time=oneyr,boundary='trim').mean()*results_bybirnrate.saturation*results_bybirnrate.Porosity.mean())*12e-3/1000*100**3*(results_bybirnrate.z_bottom-results_bybirnrate.z_top)/100).sum(dim='depth')
total_SOC=(((results_bybirnrate['Total DOM3']+results_bybirnrate['Total DOM1']+results_bybirnrate['Total DOM2']).\
            coarsen(time=oneyr,boundary='trim').mean()*results_bybirnrate.saturation*results_bybirnrate.Porosity.mean())\
                *12e-3/1000*100**3*(results_bybirnrate.z_bottom-results_bybirnrate.z_top)/100).sum(dim='depth')
total_CO2=((results_bybirnrate['Total Tracer'].coarsen(time=oneyr,boundary='trim').mean()*results_bybirnrate.saturation*results_bybirnrate.Porosity.mean())*12e-3/1000*100**3*(results_bybirnrate.z_bottom-results_bybirnrate.z_top)/100).sum(dim='depth')

# kg C/m2
total_litter=((results_bybirnrate['Total Sorbed Cellulose']+results_bybirnrate['Total Sorbed Lignin']).load().coarsen(time=oneyr,boundary='trim').mean()*12e-3*(results_bybirnrate.z_bottom-results_bybirnrate.z_top)/100).sum(dim='depth').load()

#z_bottom and z_top are in cm. This calculation results in mol/m2 (first converts mol/L to mol/cm3, then multiplies by depth in cm to get to mol/m2 in each layer)
totalMn=((results_bybirnrate['Total Mn++']+results_bybirnrate['Total Mn+++']+results_bybirnrate['Total chelated_Mn+++'])/1000*results_bybirnrate.Porosity.mean()*results_bybirnrate.saturation + \
            (results_bybirnrate['Birnessite2 VF']*7/molar_volume_birnessite) + # mol/cm3 \
            results_bybirnrate['Total Sorbed Mn++']/100**3).load().coarsen(time=oneyr,boundary='trim').mean()*100**2*(results_bybirnrate.z_bottom-results_bybirnrate.z_top)
birnessite=results_bybirnrate['Birnessite2 VF'].load().coarsen(time=oneyr,boundary='trim').mean()*7/molar_volume_birnessite*100**2*(results_bybirnrate.z_bottom-results_bybirnrate.z_top)

bioavail_Mn=(totalMn-birnessite).sum(dim='depth')

Figure3,axs=plt.subplots(nrows=4,num='Figure 3',clear=True,figsize=(5,8))
markers=['o','x','+','^']
yrs=[39,19,9,29]
n_birnrate=len(total_litter.birnrate)-1
cmap=plt.get_cmap('Blues')
norm=matplotlib.colors.LogNorm(vmin=total_litter.birnrate.min(),vmax=total_litter.birnrate.max())
for n,t in enumerate(yrs[:1]):
    for birnrate in range(n_birnrate):
        l0=axs[0].plot(bioavail_Mn.isel(time=t,birnrate=birnrate)*1e3,total_litter.isel(time=t,birnrate=birnrate),markers[n],c='k') #=cmap(norm(total_litter.birnrate[birnrate])))#='%1.1f'%(1-t/39))
        l1=axs[3].plot(bioavail_Mn.isel(time=t,birnrate=birnrate)*1e3,total_MAOM.isel(time=t,birnrate=birnrate),markers[n],c='k') #=cmap(norm(total_litter.birnrate[birnrate])))#,c='%1.1f'%(1-t/39))
        l1=axs[1].plot(bioavail_Mn.isel(time=t,birnrate=birnrate)*1e3,(total_CO2.isel(time=t,birnrate=birnrate)-total_CO2.isel(time=t-1,birnrate=birnrate))*1000,markers[n],c='k') #=cmap(norm(total_litter.birnrate[birnrate])))#,c='%1.1f'%(1-t/39))
        l2=axs[2].plot(bioavail_Mn.isel(time=t,birnrate=birnrate)*1e3,(total_SOC+total_litter).isel(time=t,birnrate=birnrate),markers[n],c='k') #=cmap(norm(total_litter.birnrate[birnrate])))#,c='%1.1f'%(1-t/39))
        if birnrate==n_birnrate-1:
            l0[0].set_label(f'Year {t+1}')

print('litter stock max = %1.2f, min = %1.2f'%(total_litter.isel(time=t,birnrate=birnrate_baseline_i).max(),total_litter.isel(time=t,birnrate=birnrate_baseline_i).min()))

axs[0].set(title='Particulate Organic C',xlabel='Bioavailable soil Mn (mmol m$^{-2}$)',ylabel='C stock (kg C m$^{-2}$)')
axs[3].set(title='Mineral-associated Organic C',xlabel='Bioavailable soil Mn (mmol m$^{-2}$)',ylabel='C stock (kg C m$^{-2}$)')
axs[1].set(title='Annual CO$_2$ production',xlabel='Bioavailable soil Mn (mmol m$^{-2}$)',ylabel='C flux (g C m$^{-2}$ year$^{-1}$)')
axs[2].set(title='Total soil and forest floor C',xlabel='Bioavailable soil Mn (mmol m$^{-2}$)',ylabel='C stock (kg C m$^{-2}$)')
# axs[0].legend()

# cb=plt.colorbar(ax=axs,mappable=plt.cm.ScalarMappable(cmap=cmap,norm=norm))
# cb.set_label('Relative birnessite dissolution rate')

for num in range(len(axs)):
    letter_label(axs[num],num)

save_one_fig(Figure3)

Figure3_expanded,axs=plt.subplots(nrows=3,ncols=4,num='Figure S5 C stocks in all scenarios',clear=True,figsize=(12,10),sharey='row')
d=alldata.isel(warming=0,Ndep=0)
# MAOM_control=((d['Total DOM3'].coarsen(time=oneyr,boundary='trim').mean()*d.saturation*d.Porosity.mean())*12e-3/1000*100**3*(d.z_bottom-d.z_top)/100).sum(dim='depth')
# totC_control=(d['Total Sorbed Cellulose'].coarsen(time=oneyr,boundary='trim').mean()*12e-3*(d.z_bottom-d.z_top)/100).sum(dim='depth')+\
#             (d['Total Sorbed Lignin'].coarsen(time=oneyr,boundary='trim').mean()*12e-3*(d.z_bottom-d.z_top)/100).sum(dim='depth') +\
#                 (((d['Total DOM3']+d['Total DOM1']+d['Total DOM2']).\
#             coarsen(time=oneyr,boundary='trim').mean()*d.saturation*d.Porosity.mean())\
#                 *12e-3/1000*100**3*(d.z_bottom-d.z_top)/100).sum(dim='depth')
for num in range(4):
    w=[1,2,0,0][num]
    n=[0,0,1,2][num]
    d=alldata.isel(warming=w,Ndep=n)
    total_cellulose=(d['Total Sorbed Cellulose'].load().coarsen(time=oneyr,boundary='trim').mean()*12e-3*(d.z_bottom-d.z_top)/100).sum(dim='depth').squeeze()
    total_lignin=(d['Total Sorbed Lignin'].load().coarsen(time=oneyr,boundary='trim').mean()*12e-3*(d.z_bottom-d.z_top)/100).sum(dim='depth').squeeze()

    total_MAOM=((d['Total DOM3'].load().coarsen(time=oneyr,boundary='trim').mean()*d.saturation*d.Porosity.mean())*12e-3/1000*100**3*(d.z_bottom-d.z_top)/100).sum(dim='depth').squeeze()
    total_SOC=(((d['Total DOM3']+d['Total DOM1']+d['Total DOM2']).load().\
            coarsen(time=oneyr,boundary='trim').mean()*d.saturation*d.Porosity.mean())\
                *12e-3/1000*100**3*(d.z_bottom-d.z_top)/100).sum(dim='depth')
    MAOM_control=total_MAOM*0.0
    totC_control=total_SOC*0.0

    totalMn=((d['Total Mn++']+d['Total Mn+++']+d['Total chelated_Mn+++'])/1000*d.Porosity.mean()*d.saturation + \
            (d['Birnessite2 VF']*7/molar_volume_birnessite) + \
            d['Total Sorbed Mn++']/100**3).squeeze().load().coarsen(time=oneyr,boundary='trim').mean()*100**2*(d.z_bottom-d.z_top)
    birnessite=d['Birnessite2 VF'].squeeze().load().coarsen(time=oneyr,boundary='trim').mean()*7/molar_volume_birnessite*100**2*(d.z_bottom-d.z_top)

    bioavail_Mn=(totalMn-birnessite).sum(dim='depth').squeeze()


    for t in [9,19,29,39]:
        # print(num,t)
        axs[0,num].plot(bioavail_Mn.isel(time=t).to_masked_array().ravel()*1e3,(total_cellulose+total_lignin).isel(time=t).to_masked_array().ravel(),'o',label='Year %d'%(t+1),c='%1.1f'%(1-t/39))
        axs[1,num].plot(bioavail_Mn.isel(time=t).to_masked_array().ravel()*1e3,total_MAOM.isel(time=t).to_masked_array().ravel(),'o',label='Year %d'%(t+1),c='%1.1f'%(1-t/39))
        axs[2,num].plot(bioavail_Mn.isel(time=t).to_masked_array().ravel()*1e3,(total_SOC+total_cellulose+total_lignin).isel(time=t).to_masked_array().ravel(),'o',label='Year %d'%(t+1),c='%1.1f'%(1-t/39))

    axs[0,num].set(title='Warming = %d, Ndep = %d\nLitter C stock'%(d['warming'].item(),d['Ndep'].item()),xlabel='Bioavailable soil Mn (mmol m$^{-2}$)',ylabel='C stock (kg C m$^{-2}$)')
    axs[1,num].set(title='Mineral-associated soil C stock',xlabel='Bioavailable soil Mn (mmol m$^{-2}$)',ylabel='C stock (kg C m$^{-2}$)')
    axs[2,num].set(title='Total soil and forest floor C stock',xlabel='Bioavailable soil Mn (mmol m$^{-2}$)',ylabel='C stock (kg C m$^{-2}$)')
    axs[0,num].tick_params(labelleft=True)
    axs[1,num].tick_params(labelleft=True)
    axs[2,num].tick_params(labelleft=True)
    axs[0,num].legend()

for num in range(axs.size):
    letter_label(axs.ravel()[num],num)

save_one_fig(Figure3_expanded)




f,axs=plt.subplots(nrows=len(controldata['birnrate']),ncols=len(controldata['soil_pH']),clear=True,num='Figure S4 Mn bioavailability',figsize=(12,10),sharex=True,sharey=True)
maxval=1.0
minval=1e-5
totalMn=((results_bybirnrate['Total Mn++']+results_bybirnrate['Total Mn+++']+results_bybirnrate['Total chelated_Mn+++'])/1000*results_bybirnrate.Porosity.mean()*results_bybirnrate.saturation + \
            (results_bybirnrate['Birnessite2 VF']*7/molar_volume_birnessite) + \
            results_bybirnrate['Total Sorbed Mn++']/100**3).load().coarsen(time=oneyr,boundary='trim').mean()*100**2*(results_bybirnrate.z_bottom-results_bybirnrate.z_top)
birnessite=results_bybirnrate['Birnessite2 VF'].load().coarsen(time=oneyr,boundary='trim').mean()*7/molar_volume_birnessite*100**2*(results_bybirnrate.z_bottom-results_bybirnrate.z_top)

bioavail_Mn=(totalMn-birnessite).sum(dim='depth')

for ph in range(len(controldata['soil_pH'])):
    for birnrate in range(len(controldata['birnrate'])):
        # Total Mn in mol/m2
        # row=4-birnrate
        row=birnrate
        h=axs[row,ph].pcolormesh(totalMn['time']/365,
                -results_bybirnrate['z_middle'].isel(soil_pH=0,birnrate=0),1-birnessite.isel(soil_pH=ph,birnrate=birnrate)/totalMn.isel(soil_pH=ph,birnrate=birnrate),
                shading='auto',cmap=plt.get_cmap('viridis'),norm=matplotlib.colors.LogNorm(vmin=minval,vmax=maxval))
        axs[row,ph].set_xlabel('Time (years)')
        if (row)==0:
            axs[row,ph].set_title('pH = %1.1f'%controldata['soil_pH'][ph].load().item(),pad=10,fontsize='large',fontweight='bold')
        if ph==0:
            axs[row,ph].text(-0.55,0.5,'Birnessite\nrate scale = %s\n(mmol Mn m$^{-3}$ y$^{-1}$)'%birnrate_labels[birnrate],
            rotation=90,va='center',ha='center',transform=axs[row,ph].transAxes,fontsize='large',fontweight='bold')

        axs[row,ph].set_ylabel('Depth (cm)')
        axs[row,ph].tick_params(labelleft=True,labelbottom=True)
            
cb=f.colorbar(h,ax=axs)
cb.set_label('Bioavailable fraction of Mn',fontsize='large')

for num in range(axs.size):
    letter_label(axs.ravel()[num],num,x=-0.1)

save_one_fig(f)

f,axs=plt.subplots(ncols=len(redoxdata['soil_pH']),nrows=len(redoxdata['anox_lenscales']),clear=True,num='Figure S3 Mn redistribution redox',figsize=(4*len(redoxdata['soil_pH']),8),sharex=True,sharey=True,squeeze=False)
totalMn_redox=(((redoxdata['Total Mn++']+redoxdata['Total Mn+++']+redoxdata['Total chelated_Mn+++'])/1000*redoxdata.Porosity.mean()*redoxdata.saturation + \
            (redoxdata['Birnessite2 VF']*7/molar_volume_birnessite) + \
            redoxdata['Total Sorbed Mn++']/100**3).load().coarsen(time=oneyr,boundary='trim').mean()*100**2*(redoxdata.z_bottom-redoxdata.z_top))
maxval=(totalMn_redox-totalMn_redox.isel(time=0)).max().compute()
minval=(totalMn_redox-totalMn_redox.isel(time=0)).min().compute()
for ph in range(len(redoxdata['soil_pH'])):
    for anox in range(len(redoxdata['anox_lenscales'])):
        row=anox
        # Total Mn in mol/m2
        h=axs[row,ph].pcolormesh(totalMn_redox['time']/365,
                -redoxdata['z_middle'].isel(soil_pH=0,anox_lenscales=0).squeeze(),(totalMn_redox.isel(soil_pH=ph,anox_lenscales=anox)-totalMn_redox.isel(soil_pH=ph,anox_lenscales=anox,time=0)).squeeze(),
                shading='auto',vmin=min(minval,-maxval),vmax=max(maxval,-minval),cmap=plt.get_cmap('RdBu_r'))
        axs[row,ph].set_xlabel('Time (years)')
        # if row==0:
        #     axs[row,ph].set_title('pH = %1.1f'%redoxdata['soil_pH'][ph].load().item(),pad=10,fontsize='large',fontweight='bold')
        # if ph==0:
        #     axs[row,ph].text(-0.6,0.5,f'Drainage time\nscale = {redoxdata["anox_lenscales"][anox].load().item()} days',
        #         rotation=90,va='center',transform=axs[row,ph].transAxes,fontsize='large',fontweight='bold')
        axs[row,ph].set_title(f'Drainage time\nscale = {redoxdata["anox_lenscales"][anox].load().item()} days')
        axs[row,ph].set_ylabel('Depth (cm)')
        axs[row,ph].tick_params(labelleft=True,labelbottom=True)
            
cb=f.colorbar(h,ax=axs)
cb.set_label('Change in Mn stock (mol m$^{-2}$)',fontsize='large')

for num in range(axs.size):
    letter_label(axs.ravel()[num],num,x=-0.1)

save_one_fig(f)


total_litter_redox=((redoxdata['Total Sorbed Cellulose']+redoxdata['Total Sorbed Lignin']).load().coarsen(time=oneyr,boundary='trim').mean()*12e-3*(redoxdata.z_bottom-redoxdata.z_top)/100).sum(dim='depth').load().squeeze()
birnessite_redox=(redoxdata['Birnessite2 VF'].load().coarsen(time=oneyr,boundary='trim').mean()*7/molar_volume_birnessite*100**2*(redoxdata.z_bottom-redoxdata.z_top)).squeeze()
total_SOC_redox=(((redoxdata['Total DOM3']+redoxdata['Total DOM1']+redoxdata['Total DOM2']).load().\
            coarsen(time=oneyr,boundary='trim').mean()*redoxdata.saturation*redoxdata.Porosity.mean())\
                *12e-3/1000*100**3*(redoxdata.z_bottom-redoxdata.z_top)/100).sum(dim='depth').squeeze()
rootuptake_redox=(redoxdata['Total Tracer2']/1000*redoxdata.Porosity.mean()*redoxdata.saturation).coarsen(time=oneyr,boundary='trim').mean()*100**2*(redoxdata.z_bottom-redoxdata.z_top)

f,axs=plt.subplots(nrows=4,ncols=1,num='Redox effects',clear=True,figsize=(4,7.5))

cmap=plt.get_cmap('Reds')
redox=redoxdata["anox_lenscales"]
axs[2].plot(redox,total_litter_redox.isel(time=39),'o')

axs[1].plot(redox,redoxdata['litter_Mn'].isel(litter_year=39).squeeze(),'o')

axs[0].plot(redox,(totalMn_redox.squeeze()-birnessite_redox).isel(depth=0,time=39),'o')

axs[3].plot(redox,total_litter_redox.isel(time=39)+total_SOC_redox.isel(time=39),'o')

axs[2].set(title='Forest floor C stock',xlabel='Drainage time scale (days)',ylabel='C stock (kg C m$^{-2}$)')
axs[1].set(title='Leaf litter Mn concentration',xlabel='Drainage time scale (days)',ylabel='Mn concentration (mmol kg$^{-1}$)')
axs[0].set(title='Forest floor bioavailable Mn',xlabel='Drainage time scale (days)',ylabel='Bioavailable Mn (mmol m$^{-2}$)')
axs[3].set(title='Total Soil C stock',xlabel='Drainage time scale (days)',ylabel='C stock (kg C m$^{-2}$)')

for num in range(axs.size):
    letter_label(axs.ravel()[num],num,x=0,title=True)

save_one_fig(f)



f,axs=plt.subplots(ncols=len(controldata['soil_pH']),nrows=len(controldata['birnrate']),clear=True,num='Figure S3 Mn redistribution',figsize=(12,10),sharex=True,sharey=True)
maxval=(totalMn-totalMn.isel(time=0)).max().compute()
minval=(totalMn-totalMn.isel(time=0)).min().compute()
for ph in range(len(controldata['soil_pH'])):
    for birnrate in range(len(controldata['birnrate'])):
        row=birnrate
        # Total Mn in mol/m2
        h=axs[row,ph].pcolormesh(totalMn['time']/365,
                -results_bybirnrate['z_middle'].isel(soil_pH=0,birnrate=0),totalMn.isel(soil_pH=ph,birnrate=birnrate)-totalMn.isel(soil_pH=ph,birnrate=birnrate,time=0),
                shading='auto',vmin=min(minval,-maxval),vmax=max(maxval,-minval),cmap=plt.get_cmap('RdBu_r'))
        axs[row,ph].set_xlabel('Time (years)')
        if row==0:
            axs[row,ph].set_title('pH = %1.1f'%controldata['soil_pH'][ph].load().item(),pad=10,fontsize='large',fontweight='bold')
        if ph==0:
            axs[row,ph].text(-0.55,0.5,'Birnessite\nrate scale = %s \n(mmol Mn m$^{-3}$ y$^{-1}$)'%birnrate_labels[birnrate],
                rotation=90,va='center',ha='center',transform=axs[row,ph].transAxes,fontsize='large',fontweight='bold')

        axs[row,ph].set_ylabel('Depth (cm)')
        axs[row,ph].tick_params(labelleft=True,labelbottom=True)
            
cb=f.colorbar(h,ax=axs)
cb.set_label('Change in Mn stock (mol m$^{-2}$)',fontsize='large')

for num in range(axs.size):
    letter_label(axs.ravel()[num],num,x=-0.1)

save_one_fig(f)


rootuptake_redox=(redoxdata['Total Tracer2']/1000*redoxdata.Porosity.mean()*redoxdata.saturation).coarsen(time=oneyr,boundary='trim').mean()*100**2*(redoxdata.z_bottom-redoxdata.z_top)
rootuptake=(results_bybirnrate['Total Tracer2']/1000*results_bybirnrate.Porosity.mean()*results_bybirnrate.saturation).coarsen(time=oneyr,boundary='trim').mean()*100**2*(results_bybirnrate.z_bottom-results_bybirnrate.z_top)

f,axs=plt.subplots(ncols=2,num='Mn leaching',figsize=(8,4),clear=True)
(((totalMn+rootuptake).isel(time=0).sum(dim='depth')-(totalMn+rootuptake).isel(time=39).sum(dim='depth'))
        /(totalMn+rootuptake).isel(time=0).sum(dim='depth')*100/4).T.plot(ax=axs[0],vmin=0,
                                            cmap='Blues'#,cbar_kwargs={'label':'Cumulative Mn leached (%/decade)'}
                                            )
axs[0].set(title='Cumulative Mn leached (Soil chem)',xlabel='Soil pH',ylabel='Birnessite dissolution rate scale\n(mmol Mn m$^{-3}$ soil year$^{-1}$)',yscale='log',ylim=(0.2,88))

axs[1].plot(satfrac.mean(dim='depth'),(((totalMn_redox+rootuptake_redox).isel(time=0).sum(dim='depth')
        -(totalMn_redox+rootuptake_redox).isel(time=39).sum(dim='depth'))
                /(totalMn_redox+rootuptake_redox).isel(time=0).sum(dim='depth')*100/4).T.squeeze())
axs[1].set(title='Cumulative Mn leached (Drainage)',xlabel='Mean saturated soil fraction',ylabel='Cumulative Mn leached (%/decade)')

save_one_fig(f)

molar_volume_manganite = 24.45 # cm3/mol
molar_volume_MnOH2am = 22.3600
molar_volume_birnessite = 251.1700
Mn_molarmass=54.94   


def plot_output(output,axs,subsample=24,do_legend=False,**kwargs):
    def fixdata(data):
        return data.coarsen(time=subsample,boundary='trim').mean().dropna(dim='time')
    outdata=output[['Total Sorbed Cellulose','Total Sorbed Lignin','Total DOM1','Total DOM3','Total Mn++',
                    'Total Mn+++','Birnessite2 VF','Total Tracer2','Free H+','Porosity','BD','saturation','CEC H+',
                    'Total Sorbed Al+++','Total Sorbed Ca++','Total Sorbed K+','Total Sorbed Mg++','Total Sorbed Na+','Total Sorbed Mn++','Total Sorbed Mn+++']].load()
    for num in range(len(output.depth)):
        out=outdata.isel(depth=num)
        t=fixdata(out.time/(365))
        porosity=out['Porosity'].mean()
        saturation=out['saturation']
        BD=out['BD']
        axs[num,0].plot(t,fixdata(out['Total Sorbed Cellulose'])*12/100**3,label='Cellulose',c='C0',**kwargs)
        axs[num,0].plot(t,fixdata(out['Total Sorbed Lignin'])*12/100**3,label='Lignin',c='C1',**kwargs)
        axs[num,0].plot(t,fixdata(out['Total DOM1']*12/1000*porosity*saturation),c='C2',label='DOM',**kwargs)
        # axs[num,0].plot(t,out['Total Sorbed DOM1']*12/100**3,label='Sorbed DOM')
        axs[num,0].plot(t,fixdata(out['Total DOM3']*12/1000*porosity*saturation),c='C3',label='Sorbed DOM',**kwargs)
        axs[num,0].set_ylabel('C density\n(g C cm$^{-3}$)')
        
        axs[num,1].plot(t,fixdata(out['Total Mn++']*1e6/1000*porosity*saturation*Mn_molarmass/BD),label='Mn$^{+\!\!+}$',c='C0',**kwargs)
        axs[num,1].plot(t,fixdata(out['Total Mn+++']*1e6/1000*porosity*saturation*Mn_molarmass/BD),label='Mn$^{+\!\!+\!\!+}$',c='C1',**kwargs)
        # axs[num,1].plot(t,layers[num].output_DF['Manganite VF']/molar_volume_manganite*1e6*Mn_molarmass/BD,label='Manganite')
        axs[num,2].plot(t,fixdata(out['Birnessite2 VF']*7/molar_volume_birnessite*1e6*Mn_molarmass/BD),label='Birnessite',c='C0',**kwargs)
        
        # annual_root_uptake=fixdata(out['Total Tracer2']).groupby(floor(t)).max() 
        annual_root_uptake=(out['Total Tracer2']).coarsen(time=oneyr,boundary='trim').max()
        axs[num,1].plot(annual_root_uptake.time/365,annual_root_uptake*1e6/1000*porosity.mean()*saturation*Mn_molarmass/BD,label='Annual Mn$^{+\!\!+}$ root uptake',c='C2',**kwargs)
        
        axs[num,3].plot(t,-np.log10(fixdata(out['Free H+'])),label='pH',c='C0',**kwargs)
        
        # ax=axs[num,1].twinx()
        # axs[num,2].plot(t,layers[num].output_DF['Mn(OH)2(am) VF']/molar_volume_MnOH2am*1e6*Mn_molarmass/BD,label='Mn(OH)$_2$(am)',ls='-')
        axs[num,2].plot(t,fixdata(out['Total Mn++']+out['Total Mn+++'])*1e6/1000*porosity*saturation*Mn_molarmass/BD + 
                            fixdata(out['Birnessite2 VF']*7/molar_volume_birnessite)*1e6*Mn_molarmass/BD,c='k',label='Total Mn',**kwargs)
        
        axs[num,1].set_ylabel('Mn concentration\n($\mu$g g$^{-1}$)')
        axs[num,2].set_ylabel('Mn concentration\n($\mu$g g$^{-1}$)')
        axs[num,3].set_ylabel('pH')
        # axs[num,1].set_ylim(*axs[0,1].get_ylim())

        for n,cation in enumerate(['Al+++','Ca++','K+','Mg++','Na+','Mn++','Mn+++']):
            axs[num,4].plot(t,fixdata(out['Total Sorbed '+cation])/(BD*1e-3*100**3)*1000,c='C'+str(n),label=cation,**kwargs)
        axs[num,4].plot(t,fixdata(out['CEC H+'])/(BD*1e-3*100**3)*1000,label='H+',c='C'+str(n+1),**kwargs)
        
        axs[num,4].set_ylabel('Exch conc\n(mmol/kg)')
        

        
    axs[-1,0].set_xlabel('Time (years)')
    if do_legend:
        axs[0,0].legend()
        axs[0,1].legend()
        axs[0,2].legend()
        axs[0,4].legend()

    axs[-1,1].set_xlabel('Time (years)')
    axs[-1,2].set_xlabel('Time (years)')
    axs[-1,3].set_xlabel('Time (years)')
    axs[-1,4].set_xlabel('Time (years)')
    axs[0,0].set_title('Organic Carbon')
    axs[0,1].set_title('Manganese ions')
    axs[0,2].set_title('Total Mn')
    axs[0,3].set_title('pH')
    axs[0,4].set_title('Exchangeable cations')


data_allNdeps=alldata.sel(warming=0)



data_warmings=alldata.sel(Ndep=0).squeeze()

warmings=data_warmings['warming'].to_masked_array()


f,axs=plt.subplots(nrows=1,ncols=3,num='Figure 4 Warming time series',clear=True,sharey='row',figsize=(12,4))
litter=total_litter_all.sel(Ndep=0).squeeze()

h=axs[0].plot(litter.time/365,litter.isel(warming=0,soil_pH=0,birnrate=birnrate_baseline_i),'b-',label='pH = %1.1f, Warming = %d\u00B0C'%(data_warmings['soil_pH'][0],warmings[0]))
h=axs[0].plot(litter.time/365,litter.isel(warming=0,soil_pH=-1,birnrate=birnrate_baseline_i),'b--',label='pH = %1.1f, Warming = %d\u00B0C'%(data_warmings['soil_pH'][-1],warmings[0]))
h=axs[0].plot(litter.time/365,litter.isel(warming=1,soil_pH=0,birnrate=birnrate_baseline_i),'y-',label='pH = %1.1f, Warming = %d\u00B0C'%(data_warmings['soil_pH'][0],warmings[1]))
h=axs[0].plot(litter.time/365,litter.isel(warming=1,soil_pH=-1,birnrate=birnrate_baseline_i),'y--',label='pH = %1.1f, Warming = %d\u00B0C'%(data_warmings['soil_pH'][-1],warmings[1]))
h=axs[0].plot(litter.time/365,litter.isel(warming=2,soil_pH=0,birnrate=birnrate_baseline_i),'r-',label='pH = %1.1f, Warming = %d\u00B0C'%(data_warmings['soil_pH'][0],warmings[2]))
h=axs[0].plot(litter.time/365,litter.isel(warming=2,soil_pH=-1,birnrate=birnrate_baseline_i),'r--',label='pH = %1.1f, Warming = %d\u00B0C'%(data_warmings['soil_pH'][-1],warmings[2]))

litter=total_litter_all.isel(Ndep=1).squeeze()

h=axs[1].plot(litter.time/365,litter.isel(warming=0,soil_pH=0,birnrate=birnrate_baseline_i),'b-',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][0],warmings[0]))
h=axs[1].plot(litter.time/365,litter.isel(warming=0,soil_pH=-1,birnrate=birnrate_baseline_i),'b--',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][-1],warmings[0]))
h=axs[1].plot(litter.time/365,litter.isel(warming=1,soil_pH=0,birnrate=birnrate_baseline_i),'y-',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][0],warmings[1]))
h=axs[1].plot(litter.time/365,litter.isel(warming=1,soil_pH=-1,birnrate=birnrate_baseline_i),'y--',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][-1],warmings[1]))
h=axs[1].plot(litter.time/365,litter.isel(warming=2,soil_pH=0,birnrate=birnrate_baseline_i),'r-',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][0],warmings[2]))
h=axs[1].plot(litter.time/365,litter.isel(warming=2,soil_pH=-1,birnrate=birnrate_baseline_i),'r--',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][-1],warmings[2]))


litter=total_litter_all.isel(Ndep=2).squeeze()

h=axs[2].plot(litter.time/365,litter.isel(warming=0,soil_pH=0,birnrate=birnrate_baseline_i),'b-',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][0],warmings[0]))
h=axs[2].plot(litter.time/365,litter.isel(warming=0,soil_pH=-1,birnrate=birnrate_baseline_i),'b--',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][-1],warmings[0]))
h=axs[2].plot(litter.time/365,litter.isel(warming=1,soil_pH=0,birnrate=birnrate_baseline_i),'y-',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][0],warmings[1]))
h=axs[2].plot(litter.time/365,litter.isel(warming=1,soil_pH=-1,birnrate=birnrate_baseline_i),'y--',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][-1],warmings[1]))
h=axs[2].plot(litter.time/365,litter.isel(warming=2,soil_pH=0,birnrate=birnrate_baseline_i),'r-',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][0],warmings[2]))
h=axs[2].plot(litter.time/365,litter.isel(warming=2,soil_pH=-1,birnrate=birnrate_baseline_i),'r--',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][-1],warmings[2]))

axs[2].set_ylabel('Total particulate C (kg C m$^{-2}$)')
axs[2].set_xlabel('Year')
axs[2].set_title('N dep = 150 kg N ha$^{-1}$ year$^{-1}$')
axs[2].tick_params(labelleft=True)

axs[1].set_ylabel('Total particulate C (kg C m$^{-2}$)')
axs[1].set_xlabel('Year')
axs[1].set_title('N dep = 50 kg N ha$^{-1}$ year$^{-1}$')
axs[1].tick_params(labelleft=True)

axs[0].set_ylabel('Total particulate C (kg C m$^{-2}$)')
axs[0].set_xlabel('Year')
axs[0].set_title('N dep = 0 kg N ha$^{-1}$ year$^{-1}$')
axs[0].legend(ncol=2)
axs[0].tick_params(labelleft=True)

for num in range(axs.size):
    letter_label(axs.ravel()[num],num)

save_one_fig(f)


totalMn=((alldata['Total Mn++']+alldata['Total Mn+++']+alldata['Total chelated_Mn+++'])/1000*alldata.Porosity.mean()*alldata.saturation + \
            (alldata['Birnessite2 VF']*7/molar_volume_birnessite) + \
            alldata['Total Sorbed Mn++']/100**3).load().coarsen(time=oneyr,boundary='trim').mean()*100**2*(alldata.z_bottom-alldata.z_top)
birnessite=alldata['Birnessite2 VF'].load().coarsen(time=oneyr,boundary='trim').mean()*7/molar_volume_birnessite*100**2*(alldata.z_bottom-alldata.z_top)
bioavail_Mn=(totalMn-birnessite).sum(dim='depth')

f,a=plt.subplots(ncols=3,num='Warming by N dep',clear=True,figsize=(12.5,4))
br=3
maxval=total_litter_all.squeeze().isel(birnrate=br,time=40).max().compute()
minval=total_litter_all.squeeze().isel(birnrate=br,time=40).min().compute()
levels=np.arange(np.round(minval,1)-0.05,np.round(maxval,1)+0.05,0.05)

a[0].contourf(bioavail_Mn.isel(birnrate=br,Ndep=0,time=40,warming=0).squeeze()*1e3,alldata['warming'],total_litter_all.squeeze().isel(birnrate=br,Ndep=0,time=40).T,levels=levels)
a[1].contourf(bioavail_Mn.isel(birnrate=br,Ndep=1,time=40,warming=0).squeeze()*1e3,alldata['warming'],total_litter_all.squeeze().isel(birnrate=br,Ndep=1,time=40).T,levels=levels)
h=a[2].contourf(bioavail_Mn.isel(birnrate=br,Ndep=2,time=40,warming=0).squeeze()*1e3,alldata['warming'],total_litter_all.squeeze().isel(birnrate=br,Ndep=2,time=40).T,levels=levels)
cb=f.colorbar(ax=a[2],mappable=h)
cb.set_label('Forest floor C stock (kg C m$^{-2}$)')
a[0].set(xlabel='Mn bioavailability (mmol m$^{-2}$)',ylabel='Warming (C)',title='N dep = 0 kg N ha$^{-1}$ year$^{-1}$')
a[1].set(xlabel='Mn bioavailability (mmol m$^{-2}$)',ylabel='Warming (C)',title='N dep = 50 kg N ha$^{-1}$ year$^{-1}$')
a[2].set(xlabel='Mn bioavailability (mmol m$^{-2}$)',ylabel='Warming (C)',title='N dep = 150 kg N ha$^{-1}$ year$^{-1}$')

save_one_fig(f)

f,a=plt.subplots(ncols=3,num='N dep by Warming',clear=True,figsize=(12.5,4))
maxval=total_litter_all.squeeze().isel(birnrate=2,time=40).max().compute()
minval=total_litter_all.squeeze().isel(birnrate=2,time=40).min().compute()
levels=np.arange(np.round(minval,1)-0.05,np.round(maxval,1)+0.05,0.05)

a[0].contourf(bioavail_Mn.isel(birnrate=2,Ndep=0,time=40,warming=0).squeeze()*1e3,alldata['Ndep'],total_litter_all.squeeze().isel(birnrate=2,warming=0,time=40).T,levels=levels)
a[1].contourf(bioavail_Mn.isel(birnrate=2,Ndep=0,time=40,warming=1).squeeze()*1e3,alldata['Ndep'],total_litter_all.squeeze().isel(birnrate=2,warming=1,time=40).T,levels=levels)
h=a[2].contourf(bioavail_Mn.isel(birnrate=2,Ndep=0,time=40,warming=2).squeeze()*1e3,alldata['Ndep'],total_litter_all.squeeze().isel(birnrate=2,warming=2,time=40).T,levels=levels)
cb=f.colorbar(ax=a[2],mappable=h)
cb.set_label('Forest floor C stock (kg C m$^{-2}$)')
a[0].set(xlabel='Mn bioavailability (mmol m$^{-2}$)',ylabel='N dep (kg N ha$^{-1}$ year$^{-1}$)',title='Warming = 0 C')
a[1].set(xlabel='Mn bioavailability (mmol m$^{-2}$)',ylabel='N dep (kg N ha$^{-1}$ year$^{-1}$)',title='Warming = 2 C')
a[2].set(xlabel='Mn bioavailability (mmol m$^{-2}$)',ylabel='N dep (kg N ha$^{-1}$ year$^{-1}$)',title='Warming = 5 C')

save_one_fig(f)

f,a=plt.subplots(ncols=2,nrows=3,num='Mn bioavail by treatment',clear=True,figsize=(11,4))
ph=1
brs=[2,4]
cm=plt.get_cmap('Reds')
for br in range(2):
    # h=a[0,br].contourf(bioavail_Mn['time']/365,alldata['warming'],bioavail_Mn.isel(birnrate=brs[br],Ndep=0,soil_pH=ph).squeeze())
    # cb=f.colorbar(h,ax=a[0,br])
    # cb.set_label('Bioavailable Mn (mmol m$^{-2}$)')
    a[0,br].plot(bioavail_Mn['time']/365,bioavail_Mn.isel(birnrate=brs[br],Ndep=0,soil_pH=ph,warming=0).squeeze(),c=cm(0.25))
    a[0,br].plot(bioavail_Mn['time']/365,bioavail_Mn.isel(birnrate=brs[br],Ndep=0,soil_pH=ph,warming=1).squeeze(),c=cm(0.50))
    a[0,br].plot(bioavail_Mn['time']/365,bioavail_Mn.isel(birnrate=brs[br],Ndep=0,soil_pH=ph,warming=2).squeeze(),c=cm(0.75))
    a[0,br].set(ylabel='Warming (C)',xlabel='Time (years)',title='Bioavailable Mn by treatment')
a[0,0].set_title('Low birnessite dissolution\nMn bioavailability')
a[0,1].set_title('High birnessite dissolution\nMn bioavailability')

minval=min(total_litter_all.isel(birnrate=2,Ndep=0,soil_pH=ph).squeeze().min().compute(),total_litter_all.isel(birnrate=3,Ndep=0,soil_pH=ph).squeeze().min().compute())
maxval=max(total_litter_all.isel(birnrate=2,Ndep=0,soil_pH=ph).squeeze().max().compute(),total_litter_all.isel(birnrate=3,Ndep=0,soil_pH=ph).squeeze().max().compute())
for br in range(2):
    # h=a[2,br].contourf(total_litter_all['time'],alldata['warming'],total_litter_all.isel(birnrate=brs[br],Ndep=0,soil_pH=ph).squeeze(),levels=np.linspace(minval,maxval,10))
    # cb=f.colorbar(h,ax=a[2,br])
    # cb.set_label('Particulate organic matter (kg C m$^{-2}$)')
    a[2,br].plot(bioavail_Mn['time']/365,total_litter_all.isel(birnrate=brs[br],Ndep=0,soil_pH=ph,warming=0).squeeze(),c=cm(0.25))
    a[2,br].plot(bioavail_Mn['time']/365,total_litter_all.isel(birnrate=brs[br],Ndep=0,soil_pH=ph,warming=1).squeeze(),c=cm(0.50))
    a[2,br].plot(bioavail_Mn['time']/365,total_litter_all.isel(birnrate=brs[br],Ndep=0,soil_pH=ph,warming=2).squeeze(),c=cm(0.75))
    a[2,br].set(ylabel='Warming (C)',xlabel='Time (years)',title='Particulate organic matter by treatment')

minval=min(alldata['litter_Mn'].isel(birnrate=2,Ndep=0,soil_pH=ph).squeeze().min().compute(),alldata['litter_Mn'].isel(birnrate=3,Ndep=0,soil_pH=ph).squeeze().min().compute())
maxval=max(alldata['litter_Mn'].isel(birnrate=2,Ndep=0,soil_pH=ph).squeeze().max().compute(),alldata['litter_Mn'].isel(birnrate=3,Ndep=0,soil_pH=ph).squeeze().max().compute())
for br in range(2):
    # h=a[1,br].contourf(alldata['litter_year'],alldata['warming'],alldata['litter_Mn'].isel(birnrate=brs[br],Ndep=0,soil_pH=ph).squeeze())
    # cb=f.colorbar(h,ax=a[1,br])
    # cb.set_label('Litter Mn concentration (mmol kg$^{-1}$)')
    a[1,br].plot(alldata['litter_year'],alldata['litter_Mn'].isel(birnrate=brs[br],Ndep=0,soil_pH=ph,warming=0).squeeze(),c=cm(0.25))
    a[1,br].plot(alldata['litter_year'],alldata['litter_Mn'].isel(birnrate=brs[br],Ndep=0,soil_pH=ph,warming=1).squeeze(),c=cm(0.50))
    a[1,br].plot(alldata['litter_year'],alldata['litter_Mn'].isel(birnrate=brs[br],Ndep=0,soil_pH=ph,warming=2).squeeze(),c=cm(0.75))
    a[1,br].set(ylabel='Warming (C)',xlabel='Time (years)',title='Litter Mn concentration by treatment')


save_one_fig(f)

total_litter=total_litter_all.sel(Ndep=0).squeeze()
#z_bottom and z_top are in cm. This calculation results in mol/m2 (first converts mol/L to mol/cm3, then multiplies by depth in cm to get to mol/m2 in each layer)
totalMn=((data_warmings['Total Mn++']+data_warmings['Total Mn+++']+data_warmings['Total chelated_Mn+++'])/1000*data_warmings.Porosity.mean()*data_warmings.saturation + \
            (data_warmings['Birnessite2 VF']*7/molar_volume_birnessite) + \
            data_warmings['Total Sorbed Mn++']/100**3).load().coarsen(time=oneyr,boundary='trim').mean()*100**2*(data_warmings.z_bottom-data_warmings.z_top)
birnessite=data_warmings['Birnessite2 VF'].load().coarsen(time=oneyr,boundary='trim').mean()*7/molar_volume_birnessite*100**2*(data_warmings.z_bottom-data_warmings.z_top)
bioavail_Mn=(totalMn-birnessite).sum(dim='depth')

Figure2,axs=plt.subplots(nrows=2,ncols=1,num='Figure 2',clear=True,figsize=(4,7.5))

cmap=plt.get_cmap('Reds')

# axs[0].plot(bioavail_Mn.isel(warming=0,time=39).to_masked_array().ravel()*1e3,total_litter.isel(warming=0,time=39).to_masked_array().ravel(),'o',c='k',label='+ O\u00B0C')
# axs[0].plot(bioavail_Mn.isel(warming=1,time=39).to_masked_array().ravel(),total_litter.isel(warming=1,time=39).to_masked_array().ravel(),'o',c=cmap(0.5),label='+ 2\u00B0C')
# axs[0].plot(bioavail_Mn.isel(warming=2,time=39).to_masked_array().ravel(),total_litter.isel(warming=2,time=39).to_masked_array().ravel(),'o',c=cmap(0.75),label='+ 5\u00B0C')

# axs[0].legend()

axs[0].plot(bioavail_Mn.isel(warming=0,time=39).to_masked_array().ravel()*1e3,data_warmings['litter_Mn'].isel(warming=0,litter_year=39).to_masked_array().ravel(),'o',c='k',label='+ O\u00B0C')
# axs[1].plot(bioavail_Mn.isel(warming=1,time=39).to_masked_array().ravel(),data_warmings['litter_Mn'].isel(warming=1,litter_year=39).to_masked_array().ravel(),'o',c=cmap(0.5),label='+ 2\u00B0C')
# axs[1].plot(bioavail_Mn.isel(warming=2,time=39).to_masked_array().ravel(),data_warmings['litter_Mn'].isel(warming=2,litter_year=39).to_masked_array().ravel(),'o',c=cmap(0.75),label='+ 5\u00B0C')

axs[1].plot(bioavail_Mn.isel(warming=0,time=39).to_masked_array().ravel()*1e3,(totalMn-birnessite).isel(depth=0,warming=0,time=39).to_masked_array().ravel()*1e3,'o',c='k',label='+ O\u00B0C')
# axs[2].plot(bioavail_Mn.isel(warming=1,time=39).to_masked_array().ravel(),(totalMn-birnessite).isel(depth=0,warming=1,time=39).to_masked_array().ravel(),'o',c=cmap(0.5),label='+ 2\u00B0C')
# axs[2].plot(bioavail_Mn.isel(warming=2,time=39).to_masked_array().ravel(),(totalMn-birnessite).isel(depth=0,warming=2,time=39).to_masked_array().ravel(),'o',c=cmap(0.75),label='+ 5\u00B0C')


# axs[0].set(title='Forest floor C stock',xlabel='Bioavailable soil Mn (mmol m$^{-2}$)',ylabel='C stock (kg C m$^{-2}$)')
axs[0].set(title='Leaf litter Mn concentration',xlabel='Bioavailable soil Mn (mmol m$^{-2}$)',ylabel='Mn concentration (mmol kg$^{-1}$)')
axs[1].set(title='Forest floor bioavailable Mn',xlabel='Bioavailable soil Mn (mmol m$^{-2}$)',ylabel='Bioavailable Mn (mmol m$^{-2}$)')

for num in range(axs.size):
    letter_label(axs.ravel()[num],num,x=0,title=True)

save_one_fig(Figure2)



total_litter=total_litter_all.squeeze()
totalMn=((alldata['Total Mn++']+alldata['Total Mn+++']+alldata['Total chelated_Mn+++'])/1000*alldata.Porosity.mean()*alldata.saturation + \
            (alldata['Birnessite2 VF']*7/molar_volume_birnessite) + \
            alldata['Total Sorbed Mn++']/100**3).load().coarsen(time=oneyr,boundary='trim').mean()*100**2*(alldata.z_bottom-alldata.z_top)
birnessite=alldata['Birnessite2 VF'].load().coarsen(time=oneyr,boundary='trim').mean()*7/molar_volume_birnessite*100**2*(alldata.z_bottom-alldata.z_top)
bioavail_Mn=(totalMn-birnessite).sum(dim='depth').squeeze()

f,axs=plt.subplots(nrows=3,ncols=2,num='Figure 5 treatment bars',clear=True,figsize=(6.5,7.5))
cm=plt.get_cmap('Reds')
cm_N=plt.get_cmap('Blues')

yr=39

axs[2,0].plot(alldata['litter_Mn'].isel(warming=0,litter_year=yr,Ndep=0,soil_pH=1).squeeze().to_masked_array(),total_litter.isel(warming=0,time=yr,Ndep=0,soil_pH=1).to_masked_array().ravel(),'-o',lw=0.5,c=cm(0.25),label='+ 0\u00B0C')
axs[2,0].plot(alldata['litter_Mn'].isel(warming=1,litter_year=yr,Ndep=0,soil_pH=1).squeeze().to_masked_array(),total_litter.isel(warming=1,time=yr,Ndep=0,soil_pH=1).to_masked_array().ravel(),'-o',lw=0.5,c=cm(0.5),label='+ 2\u00B0C')
axs[2,0].plot(alldata['litter_Mn'].isel(warming=2,litter_year=yr,Ndep=0,soil_pH=1).squeeze().to_masked_array(),total_litter.isel(warming=2,time=yr,Ndep=0,soil_pH=1).to_masked_array().ravel(),'-o',lw=0.5,c=cm(0.75),label='+ 5\u00B0C')
for num in range(5):
    axs[2,0].plot(alldata['litter_Mn'].isel(birnrate=num,litter_year=yr,Ndep=0,soil_pH=1).squeeze().to_masked_array(),total_litter.isel(birnrate=num,time=yr,Ndep=0,soil_pH=1).to_masked_array().ravel(),':',c=cm(0.5))

axs[2,1].plot(alldata['litter_Mn'].isel(warming=0,litter_year=yr,Ndep=0,soil_pH=1).squeeze().to_masked_array(),total_litter.isel(warming=0,time=yr,Ndep=0,soil_pH=1).to_masked_array().ravel(),'-o',lw=0.5,c=cm_N(0.25),label='Ndep = 0')
axs[2,1].plot(alldata['litter_Mn'].isel(warming=0,litter_year=yr,Ndep=1,soil_pH=1).squeeze().to_masked_array(),total_litter.isel(warming=0,time=yr,Ndep=1,soil_pH=1).to_masked_array().ravel(),'-o',lw=0.5,c=cm_N(0.5),label='Ndep = 50')
axs[2,1].plot(alldata['litter_Mn'].isel(warming=0,litter_year=yr,Ndep=2,soil_pH=1).squeeze().to_masked_array(),total_litter.isel(warming=0,time=yr,Ndep=2,soil_pH=1).to_masked_array().ravel(),'-o',lw=0.5,c=cm_N(0.75),label='Ndep = 150')
for num in range(5):
    axs[2,1].plot(alldata['litter_Mn'].isel(birnrate=num,warming=0,litter_year=yr,soil_pH=1).squeeze().to_masked_array(),total_litter.isel(birnrate=num,time=yr,warming=0,soil_pH=1).to_masked_array().ravel(),':',c=cm_N(0.5))


# axs[1,1].bar(np.arange(5),total_litter.isel(warming=0,time=yr,Ndep=0,soil_pH=1).to_masked_array(),color=cm(0.25),width=0.25,label='+ 0\u00B0C')
# axs[1,1].bar(np.arange(5)+0.25,total_litter.isel(warming=1,time=yr,Ndep=0,soil_pH=1).to_masked_array(),color=cm(0.5),width=0.25,label='+ 0\u00B0C')
# axs[1,1].bar(np.arange(5)+0.5,total_litter.isel(warming=2,time=yr,Ndep=0,soil_pH=1).to_masked_array(),color=cm(0.75),width=0.25,label='+ 5\u00B0C')

# axs[1,0].bar(np.arange(5),total_litter.isel(warming=0,time=yr,Ndep=0,soil_pH=1).to_masked_array(),color=cm_N(0.25),width=0.25,label='Ndep = 0')
# axs[1,0].bar(np.arange(5)+0.25,total_litter.isel(warming=0,time=yr,Ndep=1,soil_pH=1).to_masked_array(),color=cm_N(0.5),width=0.25,label='Ndep = 50')
# axs[1,0].bar(np.arange(5)+0.5,total_litter.isel(warming=0,time=yr,Ndep=2,soil_pH=1).to_masked_array(),color=cm_N(0.75),width=0.25,label='Ndep = 150')



axs[1,0].bar(np.arange(5),alldata['litter_Mn'].isel(warming=0,litter_year=yr,Ndep=0,soil_pH=1).squeeze().to_masked_array(),color=cm(0.25),width=0.25,label='+ 0\u00B0C')
axs[1,0].bar(np.arange(5)+0.25,alldata['litter_Mn'].isel(warming=1,litter_year=yr,Ndep=0,soil_pH=1).squeeze().to_masked_array(),color=cm(0.5),width=0.25,label='+ 2\u00B0C')
axs[1,0].bar(np.arange(5)+0.5,alldata['litter_Mn'].isel(warming=2,litter_year=yr,Ndep=0,soil_pH=1).squeeze().to_masked_array(),color=cm(0.75),width=0.25,label='+ 5\u00B0C')

axs[1,1].bar(np.arange(5),alldata['litter_Mn'].isel(warming=0,litter_year=yr,Ndep=0,soil_pH=1).squeeze().to_masked_array(),color=cm_N(0.25),width=0.25,label='Ndep = 0')
axs[1,1].bar(np.arange(5)+0.25,alldata['litter_Mn'].isel(warming=0,litter_year=yr,Ndep=1,soil_pH=1).squeeze().to_masked_array(),color=cm_N(0.5),width=0.25,label='Ndep = 50')
axs[1,1].bar(np.arange(5)+0.5,alldata['litter_Mn'].isel(warming=0,litter_year=yr,Ndep=2,soil_pH=1).squeeze().to_masked_array(),color=cm_N(0.75),width=0.25,label='Ndep = 150')

axs[0,0].bar(np.arange(5),bioavail_Mn.isel(warming=0,time=yr,Ndep=0,soil_pH=1).squeeze().to_masked_array(),color=cm(0.25),width=0.25,label='+ 0\u00B0C')
axs[0,0].bar(np.arange(5)+0.25,bioavail_Mn.isel(warming=1,time=yr,Ndep=0,soil_pH=1).squeeze().to_masked_array(),color=cm(0.5),width=0.25,label='+ 2\u00B0C')
axs[0,0].bar(np.arange(5)+0.5,bioavail_Mn.isel(warming=2,time=yr,Ndep=0,soil_pH=1).squeeze().to_masked_array(),color=cm(0.75),width=0.25,label='+ 5\u00B0C')

axs[0,1].bar(np.arange(5),bioavail_Mn.isel(warming=0,time=yr,Ndep=0,soil_pH=1).squeeze().to_masked_array(),color=cm_N(0.25),width=0.25,label='Ndep = 0')
axs[0,1].bar(np.arange(5)+0.25,bioavail_Mn.isel(warming=0,time=yr,Ndep=1,soil_pH=1).squeeze().to_masked_array(),color=cm_N(0.5),width=0.25,label='Ndep = 50')
axs[0,1].bar(np.arange(5)+0.5,bioavail_Mn.isel(warming=0,time=yr,Ndep=2,soil_pH=1).squeeze().to_masked_array(),color=cm_N(0.75),width=0.25,label='Ndep = 150')

axs[0,0].legend()
axs[0,1].legend()

axs[0,1].set(title='Soil bioavailable Mn',xlabel='Birnessite dissolution rate scale\n(mmol Mn m$^{-3}$ soil year$^{-1}$)',ylabel='Bioavailable Mn (mmol m$^{-2}$)',xticks=np.arange(5)+0.25,xticklabels=birnrate_labels,yscale='log')
axs[1,1].set(title='Leaf litter Mn concentration',xlabel='Birnessite dissolution rate scale\n(mmol Mn m$^{-3}$ soil year$^{-1}$)',ylabel='Mn concentration (mmol kg$^{-1}$)',xticks=np.arange(5)+0.25,xticklabels=birnrate_labels,ylim=(0,205))
axs[2,1].set(title='Particulate organic matter stock',xlabel='Litter Mn concentration (mmol kg$^{-1}$)',ylabel='POM stock (kg C m$^{-2}$)',ylim=(0.9,1.6))
axs[0,0].set(title='Soil bioavailable Mn',xlabel='Birnessite dissolution rate scale\n(mmol Mn m$^{-3}$ soil year$^{-1}$)',ylabel='Bioavailable Mn (mmol m$^{-2}$)',xticks=np.arange(5)+0.25,xticklabels=birnrate_labels,yscale='log')
axs[1,0].set(title='Leaf litter Mn concentration',xlabel='Birnessite dissolution rate scale\n(mmol Mn m$^{-3}$ soil year$^{-1}$)',ylabel='Mn concentration (mmol kg$^{-1}$)',xticks=np.arange(5)+0.25,xticklabels=birnrate_labels,ylim=(0,205))
axs[2,0].set(title='Particulate organic matter stock',xlabel='Litter Mn concentration (mmol kg$^{-1}$)',ylabel='POM stock (kg C m$^{-2}$)',ylim=(0.9,1.6))


# axs[1,0].plot(bioavail_Mn.isel(warming=0,Ndep=0,time=39).to_masked_array().ravel(),alldata['litter_Mn'].isel(warming=0,Ndep=0,litter_year=39).to_masked_array().ravel(),'o',c=cmap(0.25),label='+ O\u00B0C')
# axs[1,1].plot(bioavail_Mn.isel(warming=0,Ndep=0,time=39).to_masked_array().ravel(),alldata['litter_Mn'].isel(warming=0,Ndep=0,litter_year=39).to_masked_array().ravel(),'o',c=cmap(0.25),label='+ O\u00B0C')
# axs[1,1].plot(bioavail_Mn.isel(warming=1,Ndep=0,time=39).to_masked_array().ravel(),alldata['litter_Mn'].isel(warming=1,Ndep=0,litter_year=39).to_masked_array().ravel(),'o',c=cmap(0.5),label='+ 2\u00B0C')
# axs[1,1].plot(bioavail_Mn.isel(warming=2,Ndep=0,time=39).to_masked_array().ravel(),alldata['litter_Mn'].isel(warming=2,Ndep=0,litter_year=39).to_masked_array().ravel(),'o',c=cmap(0.75),label='+ 5\u00B0C')

# axs[1,0].plot(bioavail_Mn.isel(warming=0,Ndep=1,time=39).to_masked_array().ravel(),alldata['litter_Mn'].isel(warming=0,Ndep=1,litter_year=39).to_masked_array().ravel(),'^',c=cmap(0.5),label='+ O\u00B0C')
# axs[1].plot(bioavail_Mn.isel(warming=1,Ndep=1,time=39).to_masked_array().ravel(),alldata['litter_Mn'].isel(warming=1,Ndep=1,litter_year=39).to_masked_array().ravel(),'^',c=cmap(0.5),label='+ 2\u00B0C')
# axs[1].plot(bioavail_Mn.isel(warming=2,Ndep=1,time=39).to_masked_array().ravel(),alldata['litter_Mn'].isel(warming=2,Ndep=1,litter_year=39).to_masked_array().ravel(),'^',c=cmap(0.75),label='+ 5\u00B0C')

# axs[1,0].plot(bioavail_Mn.isel(warming=0,Ndep=2,time=39).to_masked_array().ravel(),alldata['litter_Mn'].isel(warming=0,Ndep=2,litter_year=39).to_masked_array().ravel(),'s',c=cmap(0.75),label='+ O\u00B0C')
# axs[1].plot(bioavail_Mn.isel(warming=1,Ndep=2,time=39).to_masked_array().ravel(),alldata['litter_Mn'].isel(warming=1,Ndep=2,litter_year=39).to_masked_array().ravel(),'s',c=cmap(0.5),label='+ 2\u00B0C')
# axs[1].plot(bioavail_Mn.isel(warming=2,Ndep=2,time=39).to_masked_array().ravel(),alldata['litter_Mn'].isel(warming=2,Ndep=2,litter_year=39).to_masked_array().ravel(),'s',c=cmap(0.75),label='+ 5\u00B0C')

# axs[2,0].plot(bioavail_Mn.isel(warming=0,Ndep=0,time=39,soil_pH=1).to_masked_array().ravel(),(totalMn-birnessite).isel(depth=0,warming=0,Ndep=0,time=39,soil_pH=1).to_masked_array().ravel(),'o',c=cmap(0.25),label='+ O\u00B0C')
# axs[2,1].plot(bioavail_Mn.isel(warming=0,Ndep=0,time=39,soil_pH=1).to_masked_array().ravel(),(totalMn-birnessite).isel(depth=0,warming=0,Ndep=0,time=39,soil_pH=1).to_masked_array().ravel(),'o',c=cmap(0.25),label='+ O\u00B0C')
# axs[2,1].plot(bioavail_Mn.isel(warming=1,Ndep=0,time=39,soil_pH=1).to_masked_array().ravel(),(totalMn-birnessite).isel(depth=0,warming=1,Ndep=0,time=39,soil_pH=1).to_masked_array().ravel(),'o',c=cmap(0.5),label='+ 2\u00B0C')
# axs[2,1].plot(bioavail_Mn.isel(warming=2,Ndep=0,time=39,soil_pH=1).to_masked_array().ravel(),(totalMn-birnessite).isel(depth=0,warming=2,Ndep=0,time=39,soil_pH=1).to_masked_array().ravel(),'o',c=cmap(0.75),label='+ 5\u00B0C')

# axs[2,0].plot(bioavail_Mn.isel(warming=0,Ndep=1,time=39,soil_pH=1).to_masked_array().ravel(),(totalMn-birnessite).isel(depth=0,warming=0,Ndep=1,time=39,soil_pH=1).to_masked_array().ravel(),'^',c=cmap(0.5),label='+ O\u00B0C')
# # axs[2].plot(bioavail_Mn.isel(warming=1,Ndep=1,time=39).to_masked_array().ravel(),(totalMn-birnessite).isel(depth=0,warming=1,Ndep=1,time=39).to_masked_array().ravel(),'^',c=cmap(0.5),label='+ 2\u00B0C')
# # axs[2].plot(bioavail_Mn.isel(warming=2,Ndep=1,time=39).to_masked_array().ravel(),(totalMn-birnessite).isel(depth=0,warming=2,Ndep=1,time=39).to_masked_array().ravel(),'^',c=cmap(0.75),label='+ 5\u00B0C')

# axs[2,0].plot(bioavail_Mn.isel(warming=0,Ndep=2,time=39,soil_pH=1).to_masked_array().ravel(),(totalMn-birnessite).isel(depth=0,warming=0,Ndep=2,time=39,soil_pH=1).to_masked_array().ravel(),'s',c=cmap(0.75),label='+ O\u00B0C')
# # axs[2].plot(bioavail_Mn.isel(warming=1,Ndep=2,time=39).to_masked_array().ravel(),(totalMn-birnessite).isel(depth=0,warming=1,Ndep=2,time=39).to_masked_array().ravel(),'s',c=cmap(0.5),label='+ 2\u00B0C')
# # axs[2].plot(bioavail_Mn.isel(warming=2,Ndep=2,time=39).to_masked_array().ravel(),(totalMn-birnessite).isel(depth=0,warming=2,Ndep=2,time=39).to_masked_array().ravel(),'s',c=cmap(0.75),label='+ 5\u00B0C')




for num in range(axs.size):
    letter_label(axs.ravel()[num],num)

save_one_fig(f)



import Manganese_network as Mn
import decomp_network
reaction_network=Mn.make_network(leaf_Mn_mgkg=0.0,Mn2_scale=1e-4,Mn3_scale=1e-12) 

networkfig=plt.figure('Reaction network',clear=True)
pos={'Birnessite2': (644.0, 594.0),
 'Mn+++': (635.0, 522.0),
 'chelated_Mn+++':(650,530),
 'O2(aq)': (865.0, 450.0),
 'Cellulose': (72.0, 378.0),
 'DOM1': (346.0, 234.0),
 'DOM2': (130.0, 378.0),
 'HCO3-': (667.0, 90.0),
 'Mn++': (683.0, 378.0),
 'Tracer2': (552.0, 234.0),
 'Mg++': (204.0, 505.0),
 'Ca++': (284.0, 505.0),
 'Na+': (359.0, 505.0),
 'K+': (431.0, 505.0),
 'Al+++': (120.0, 505.0),
 'DOM3': (203.0, 90.0),
 'Hydrolysis': (147.0, 306.0),
 'DOM aerobic respiration': (884.0, 162.0),
 'DOM1 Mn+++ reduction': (667.0, 162.0),
 'DOM1 Mn+++ abiotic reduction': (422.0, 162.0),
 'DOM sorption': (216.0, 162.0),
 'Mn Peroxidase': (346.0, 306.0),
 'Bacterial Mn++ oxidation': (873.0, 306.0),
 'Cation exchange': (335.0, 420.0),
 'Root uptake of Mn++': (552.0, 306.0),
 'DOM desorption': (164.0, 18.0)}
drawn,pos=decomp_network.draw_network_with_reactions(reaction_network,
        omit=['NH4+','Rock(s)','gas','secondary','H+','>Carboxylate-','Carboxylic_acid','H2O','Root_biomass',
            'Sorption_capacity','Lignin depolymerization','Lignin exposure','Lignin','DOM2 aerobic respiration'],
        font_size='medium',node_size=500,font_color='k',arrowstyle='->',arrowsize=10.0,edge_color='gray',node_alpha=1.0,
        namechanges={'cellulose':'Cellulose','DOM1':'DOM','O2(aq)':'O$_2$(aq)','CH4(aq)':'CH$_4$(aq)','HCO3-':'CO$_2$','DOM2':'Lignin','sorbed_DOM1':'Sorbed DOM',
                     'Fe(OH)2':'Fe(OH)$_2$','Fe(OH)3':'Fe(OH)$_3$','Mn++':r'Mn$^\mathrm{+\!\!+}$','Mn+++':r'Mn$^\mathrm{+\!\!+\!\!\!+}$','Acetate-':'Acetate',
                     'DOM3':'MAOM','Tracer2':'Root uptake','Birnessite2':'Birnessite'},connectionstyle='arc3, rad=0.2')#,pos=pos)

save_one_fig(networkfig)

######
# Plot annual decomp for comparison with litter decomposition studies
#####


incubation_data=xarray.combine_by_coords(incubation_list).squeeze()
incubation_data=controldata.isel(depth=0).squeeze()
import pandas
datadir='/home/b0u/Mn data/'
Berg_massloss1=pandas.read_csv(datadir+'Berg et al 2015/Fig 4A.csv')
Berg_massloss2=pandas.read_csv(datadir+'Berg et al 2015/Fig 4B.csv')

f,axs=plt.subplots(1,1,num='Figure S2 Litter annual mass loss',clear=True,figsize=(6,4.6))
incubation_length=1
addyrs=0
litter_massloss=(1-(incubation_data['Total Sorbed Lignin']+incubation_data['Total Sorbed Cellulose']).isel(time=slice(oneyr+addyrs*oneyr,None,oneyr*incubation_length))/(incubation_data['Total Sorbed Lignin']+incubation_data['Total Sorbed Cellulose']).isel(time=1+addyrs*oneyr)).to_masked_array().ravel()*100
lignin_massloss=(1-(incubation_data['Total Sorbed Lignin']).isel(time=slice(oneyr+addyrs*oneyr,None,oneyr*incubation_length),)/(incubation_data['Total Sorbed Lignin']).isel(time=1+addyrs*oneyr)).to_masked_array().ravel()*100

# Mn_conc=incubation_data['litter_Mn'].to_masked_array().ravel().compressed()[1::incubation_length]*leafC_mass_conv
Mn_conc=(incubation_data['Total Mn++']*(incubation_data['Porosity']*incubation_data['saturation']*incubation_data['dz']/100*1e3)*1e3/(.163/.4)).isel(time=slice(1,None,oneyr*incubation_length)).to_masked_array().ravel()

Mn_conc=(incubation_data['Total Mn++']*(incubation_data['Porosity']*incubation_data['saturation']*incubation_data['dz']/100*1e3)*1e3/(.163/.4)).coarsen(time=oneyr,boundary='trim').max()
lignindep=.163*.25/.05/12*1e3
# lignin_massloss=100-(incubation_data['Total Sorbed Lignin']-incubation_data['Total Sorbed Lignin'].shift(time=oneyr)).isel(time=slice(oneyr*5,oneyr*15,oneyr)).to_masked_array().ravel()/lignindep*100
lignin=incubation_data['Total Sorbed Lignin']+incubation_data['Total DOM2']*1000*0.5
lignin_massloss=(lignin.coarsen(time=oneyr,boundary='trim').max()-lignin.coarsen(time=oneyr,boundary='trim').min()).compute()/lignindep*100

# x=np.argsort(Mn_conc)
# axs[0].plot(Mn_conc[x],litter_massloss[x],'-o',label='Model mass loss')
axs.plot(Mn_conc[:,:,15:30:4].to_masked_array().ravel(),lignin_massloss[:,:,15:30:4].to_masked_array().ravel(),'o',label='Model lignin mass loss')
# axs[0].plot(Berg_massloss1['x'],Berg_massloss1['Massloss'],'+',label='Berg multiple species')
axs.plot(Berg_massloss2['x']/(1e-3*Mn_molarmass),Berg_massloss2['Massloss'],'+',label='Berg et al. (2013) measurements')
# axs[0].plot(Davey_data['Mn mg-g-1 DM'],Davey_data['Limit value']*(1-exp(-Davey_data['k value']*oneyr*1*100/Davey_data['Limit value'])),'x',label='Davey Oak')  
axs.set_xlabel('Litter Mn concentration (mmol kg$^{-1}$)')
axs.set_ylabel('Mass loss (%)')
axs.legend()
axs.set_title('One year mass loss')

save_one_fig(f)

# Davey_data=pandas.read_excel(datadir+'Davey et al 2007/Davey et al 2007 table data.xlsx',sheet_name=0)

# axs[1].set_title('Mass loss over time')
# t=np.linspace(0,1000,20)
# cmap=plt.get_cmap('viridis')
# for site in range(len(Davey_data)):
#     m=Davey_data['Limit value'][site]
#     k=Davey_data['k value'][site]
#     davey_line=axs[1].plot(t,m*(1-np.exp(-k*t/(m/100))),c=cmap(Davey_data['Mn mg-g-1 DM'][site]/(1e-3*Mn_molarmass)/Mn_conc.max()),ls=':')
    
# initial_mass=(incubation_data['Total Sorbed Lignin']+incubation_data['Total Sorbed Cellulose']).isel(time=1)   
# # t=arange(0,oneyr*incubation_length-dt,dt)
# for ph in range(len(np.atleast_1d(incubation_data.soil_pH))):
#     for birnrate in range(len(np.atleast_1d(incubation_data.birnrate))):
#         for rep in range(len(incubation_data['litter_year'])//incubation_length):
#             massloss=(1-(incubation_data['Total Sorbed Lignin']+incubation_data['Total Sorbed Cellulose']).isel(soil_pH=ph,birnrate=birnrate)[1+oneyr*rep*incubation_length:oneyr*(rep+1)*incubation_length]/initial_mass.isel(soil_pH=ph,birnrate=birnrate))*100
#             mod_line=axs[1].plot(massloss['time']-massloss['time'][0],massloss,c=cmap(incubation_data['litter_Mn'].isel(soil_pH=ph,birnrate=birnrate,litter_year=rep*incubation_length)*leafC_mass_conv/Mn_conc.max()))

# axs[1].set_xlabel('Time (days)')
# axs[1].set_ylabel('Mass loss (%)')
# cb=f.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(vmin=0,vmax=Mn_conc.max())),ax=axs[1])
# cb.set_label('Litter Mn concentration (mmol kg$^{-1}$)')
# axs[1].legend(handles=davey_line+mod_line,labels=['Davey et al. (2007) measurements','Model'])

# for num in range(axs.size):
#     letter_label(axs.ravel()[num],num)

# save_one_fig(f)

f,axs=plt.subplots(nrows=6,num='Incubation time series',clear=True,figsize=(4,6))
axs[0].plot(incubation_data['Total Mn++'].isel(birnrate=4,soil_pH=0))
axs[0].plot(incubation_data['Total Mn++'].isel(birnrate=1,soil_pH=0))
axs[0].plot(incubation_data['Total Mn++'].isel(birnrate=4,soil_pH=4))
axs[0].plot(incubation_data['Total Mn++'].isel(birnrate=1,soil_pH=4))
axs[0].set(title='Mn(II) concentration')

axs[1].plot(incubation_data['Total chelated_Mn+++'].isel(birnrate=4,soil_pH=0))
axs[1].plot(incubation_data['Total chelated_Mn+++'].isel(birnrate=1,soil_pH=0))
axs[1].plot(incubation_data['Total chelated_Mn+++'].isel(birnrate=4,soil_pH=4))
axs[1].plot(incubation_data['Total chelated_Mn+++'].isel(birnrate=1,soil_pH=4))
axs[1].set(title='Chelated Mn(III) concentration')

axs[2].plot(incubation_data['Total DOM2'].isel(birnrate=4,soil_pH=0))
axs[2].plot(incubation_data['Total DOM2'].isel(birnrate=1,soil_pH=0))
axs[2].plot(incubation_data['Total DOM2'].isel(birnrate=4,soil_pH=4))
axs[2].plot(incubation_data['Total DOM2'].isel(birnrate=1,soil_pH=4))
axs[2].set(title='DOM2 concentration')

axs[3].plot(incubation_data['Total Sorbed Lignin'].isel(birnrate=4,soil_pH=0))
axs[3].plot(incubation_data['Total Sorbed Lignin'].isel(birnrate=1,soil_pH=0))
axs[3].plot(incubation_data['Total Sorbed Lignin'].isel(birnrate=4,soil_pH=4))
axs[3].plot(incubation_data['Total Sorbed Lignin'].isel(birnrate=1,soil_pH=4))
axs[3].set_title('Lignin')

axs[4].plot(incubation_data['Total Tracer'].isel(birnrate=4,soil_pH=0))
axs[4].plot(incubation_data['Total Tracer'].isel(birnrate=1,soil_pH=0))
axs[4].plot(incubation_data['Total Tracer'].isel(birnrate=4,soil_pH=4))
axs[4].plot(incubation_data['Total Tracer'].isel(birnrate=1,soil_pH=4))
axs[4].set_title('Cumulative CO2')

axs[5].plot(incubation_data['Total DOM1'].isel(birnrate=4,soil_pH=0))
axs[5].plot(incubation_data['Total DOM1'].isel(birnrate=1,soil_pH=0))
axs[5].plot(incubation_data['Total DOM1'].isel(birnrate=4,soil_pH=4))
axs[5].plot(incubation_data['Total DOM1'].isel(birnrate=1,soil_pH=4))
axs[5].set_title('DOM1')

save_one_fig(f)

def save_all_figs(dirname,format='png',**kwargs):
    for fname in plt.get_figlabels():
        fname_fixed=fname.replace('/','-')
        print(fname_fixed)
        plt.figure(fname_fixed).savefig('{dirname:s}/{fname:s}.{format}'.format(dirname=dirname,format=format,fname=fname_fixed),**kwargs)





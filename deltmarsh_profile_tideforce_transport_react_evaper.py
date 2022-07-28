#Created by Ben Sulman for the Manganese soil profile, and revised by Jiaze Wang to simulate delta marsh soil profile,02/18/2021

from run_alquimia import get_alquimiavector,ffi,lib,check_status,init_alquimia,convert_condition_to_alquimia,print_metadata
import decomp_network
import delmar_network_tide_evp as Mar
from matplotlib import pyplot
import numpy
import xarray
import pdb


class layer:
    def __init__(self,volume,saturation=0.9,temperature=20.0,water_density=1000.0,porosity=0.25,pressure=101325.0,BD=1.5,CEC=None,rateconstants={},diffquo={}):
        self.volume=volume # m3
        self.saturation=saturation   # out of 1.0
        self.temperature=temperature # degrees C
        self.water_density=water_density # kg/m3
        self.porosity=porosity # out of 1.0
        self.pressure=pressure # Pa
        self.mineral_specific_surface_area={}
        self.mineral_volume_fraction={}
        self.surface_site_density={}
        self.total_immobile={}
        self.total_mobile={}
        self.free_mobile={}
        self.rateconstants=rateconstants
        self.mineral_rate_cnst={}
        self.aux_ints=[]
        self.aux_doubles=[]
        self.mineral_reaction_rate={}
        self.secondary_free_ion_concentration={}
        self.flow_in={}
        self.flow_out={}
        self.diffquo=diffquo
        self.BD=BD   # g/cm3
        self.CEC=CEC # meq/kg       #cjw-- remove CEC
#cjw above function defines some essential attributes for each instance
        # Some things that I'm not using are skipped

    def copy_from_alquimia(self,data):
        self.volume=data.properties.volume
        self.saturation=data.properties.saturation
        self.temperature=data.state.temperature
        self.water_density=data.state.water_density
        self.porosity=data.state.porosity
        self.pressure=data.state.aqueous_pressure

        self.mineral_names=get_alquimiavector(data.meta_data.mineral_names)
        for num in range(data.meta_data.mineral_names.size):
            self.mineral_rate_cnst[self.mineral_names[num]]=data.properties.mineral_rate_cnst.data[num]
            self.mineral_specific_surface_area[self.mineral_names[num]]=data.state.mineral_specific_surface_area.data[num]
            self.mineral_volume_fraction[self.mineral_names[num]]=data.state.mineral_volume_fraction.data[num]
            self.mineral_reaction_rate[self.mineral_names[num]]=data.aux_output.mineral_reaction_rate.data[num]

        self.primary_names=get_alquimiavector(data.meta_data.primary_names)
        for num in range(len(self.primary_names)):
            self.total_immobile[self.primary_names[num]]=data.state.total_immobile.data[num]
            self.total_mobile[self.primary_names[num]]=data.state.total_mobile.data[num]
            self.free_mobile[self.primary_names[num]]=data.aux_output.primary_free_ion_concentration.data[num]

        for num in range(len(self.secondary_names)):
            self.secondary_free_ion_concentration[self.secondary_names[num]]=data.aux_output.secondary_free_ion_concentration.data[num]

        self.surface_site_names=get_alquimiavector(data.meta_data.surface_site_names)
        for num in range(len(self.surface_site_names)):
            self.surface_site_density[self.surface_site_names[num]]=data.state.surface_site_density.data[num]

        self.aux_ints=get_alquimiavector(data.aux_data.aux_ints)
        self.aux_doubles=get_alquimiavector(data.aux_data.aux_doubles)


    def copy_to_alquimia(self,data):
        data.properties.volume=self.volume
        data.properties.saturation=self.saturation
        data.state.temperature=self.temperature
        data.state.water_density=self.water_density
        data.state.porosity=self.porosity
        data.state.aqueous_pressure=self.pressure


        for num in range(len(self.mineral_names)):
            data.properties.mineral_rate_cnst.data[num]=self.mineral_rate_cnst[self.mineral_names[num]]
            data.state.mineral_specific_surface_area.data[num]=self.mineral_specific_surface_area[self.mineral_names[num]]
            data.state.mineral_volume_fraction.data[num]=self.mineral_volume_fraction[self.mineral_names[num]]

        for num in range(len(self.surface_site_names)):
            data.state.surface_site_density.data[num]=self.surface_site_density[self.surface_site_names[num]]

        for num in range(len(self.primary_names)):
            data.state.total_immobile.data[num]=self.total_immobile[self.primary_names[num]]
            data.state.total_mobile.data[num]=self.total_mobile[self.primary_names[num]]

        for num in range(len(self.aux_ints)):
            data.aux_data.aux_ints.data[num]=self.aux_ints[num]
        for num in range(len(self.aux_doubles)):
            data.aux_data.aux_doubles.data[num]=self.aux_doubles[num]


    def setup_output(self,nsteps,dt):
        self.output={
            'total_mobile':numpy.ma.masked_all((nsteps,len(self.total_mobile)),dtype=float),
            'free':numpy.ma.masked_all((nsteps,len(self.total_mobile)),dtype=float),
            'immobile':numpy.ma.masked_all((nsteps,len(self.total_mobile)),dtype=float),
            'aq_complex':numpy.ma.masked_all((nsteps,len(self.secondary_names)),dtype=float),
            'mineral_VF':numpy.ma.masked_all((nsteps,len(self.mineral_rate_cnst)),dtype=float),
            'mineral_rate':numpy.ma.masked_all((nsteps,len(self.mineral_rate_cnst)),dtype=float),
            'time':numpy.arange(nsteps,dtype=float)*dt,
            'actual_dt':numpy.ma.masked_all(nsteps,dtype=float),'ncuts':numpy.ma.masked_all(nsteps,dtype=int),
            'porosity':numpy.ma.masked_all(nsteps,dtype=float),
            'saturation':numpy.ma.masked_all(nsteps,dtype=float),
#            'CEC H+':numpy.ma.masked_all(nsteps,dtype=float),             #cjw-- remove CEC
            'flow_in':numpy.ma.masked_all((nsteps,len(self.total_mobile)),dtype=float),
            'flow_out':numpy.ma.masked_all((nsteps,len(self.total_mobile)),dtype=float),
        }

    def write_output(self,step,dt,num_cuts=0):
        self.output['total_mobile'][step,:]=numpy.array([self.total_mobile[name] for name in self.primary_names])
        self.output['free'][step,:]=numpy.array([self.free_mobile[name] for name in self.primary_names])
        self.output['aq_complex'][step,:]=numpy.array([self.secondary_free_ion_concentration[name] for name in self.secondary_names])
        self.output['immobile'][step,:]=numpy.array([self.total_immobile[name] for name in self.primary_names])
        self.output['mineral_VF'][step,:]=numpy.array([self.mineral_volume_fraction[name] for name in self.mineral_names])
        self.output['time'][step]=step*dt
        self.output['mineral_rate'][step,:]=numpy.array([self.mineral_reaction_rate[name] for name in self.mineral_names])
        self.output['porosity'][step]=self.porosity
        self.output['saturation'][step]=self.saturation      #cjw update saturation
        self.output['actual_dt'][step]=dt/2**num_cuts
        self.output['ncuts'][step]=num_cuts
        # Keep track of how much H+ is in the carboxylate buffer rather than the CEC site, so we can plot CEC-exchangeable H+ later.
        # Assumes there is one ion exchange site, and that the buffering sorption site is on Rock(s)
        # Aux doubles in this spot stores free surface site density (probably best not to rely on aux_doubles in general though)
#        self.output['CEC H+'][step]=self.total_immobile['H+']-(self.surface_site_density['>Carboxylate-']*self.mineral_volume_fraction['Rock(s)']-self.aux_doubles[len(self.total_mobile)*2+len(self.secondary_free_ion_concentration)+1+self.surface_site_names.index('>Carboxylate-')])  ##cjw--remove CEC
        self.output['flow_in'][step,:]=numpy.array([self.flow_in.get(name,0.0) for name in self.primary_names])
        self.output['flow_out'][step,:]=numpy.array([self.flow_out.get(name,0.0) for name in self.primary_names])

    def convert_output(self):
        import pandas
        output_DF=pandas.DataFrame(index=self.output['time'])
        output_DF['Porosity']=self.output['porosity']
        output_DF['ncuts']=self.output['ncuts']
        output_DF['actual_dt']=self.output['actual_dt']
        output_units={'time':'s','Porosity':'NA','ncuts':'NA','actual_dt':'s'}
        total=pandas.DataFrame(self.output['total_mobile'],columns=self.primary_names,index=self.output['time']).add_prefix('Total ')
        for col in total.columns:
            output_DF[col]=total[col]
        output_units.update([('Total '+s,'M') for s in self.primary_names])
        free=pandas.DataFrame(self.output['free'],columns=self.primary_names,index=self.output['time']).add_prefix('Free ')
        for col in free.columns:
            output_DF[col]=free[col]
        output_units.update([('Free '+s,'M') for s in self.primary_names])
        sorbed=pandas.DataFrame(self.output['immobile'],columns=self.primary_names,index=self.output['time']).add_prefix('Total Sorbed ')
        for col in sorbed.columns:
            output_DF[col]=sorbed[col]
        output_units.update([('Total Sorbed '+s,'mol/m^3') for s in self.primary_names])
        mineral_VF=pandas.DataFrame(self.output['mineral_VF'],columns=self.mineral_names,index=self.output['time']).add_suffix(' VF')
        for col in mineral_VF.columns:
            output_DF[col]=mineral_VF[col]
        output_units.update([(s+' VF','m^3 mnrl/m^3 bulk') for s in self.mineral_names])
        mineral_rate=pandas.DataFrame(self.output['mineral_rate'],columns=self.mineral_names,index=self.output['time']).add_suffix(' Rate')
        for col in mineral_rate.columns:
            output_DF[col]=mineral_rate[col]
        output_units.update([(s+' Rate','mol/m^3/sec') for s in self.mineral_names])
        secondary=pandas.DataFrame(self.output['aq_complex'],columns=self.secondary_names,index=self.output['time'])
        for col in secondary.columns:
            output_DF[col]=secondary[col]
        output_units.update([(s,'M') for s in self.secondary_names])

        flow_in=pandas.DataFrame(self.output['flow_in'],columns=self.primary_names,index=self.output['time']).add_suffix(' inflow')
        for col in flow_in.columns:
            output_DF[col]=flow_in[col]
        output_units.update([(s+' inflow','mol/m2/sec') for s in self.primary_names])
        flow_out=pandas.DataFrame(self.output['flow_out'],columns=self.primary_names,index=self.output['time']).add_suffix(' outflow')
        for col in flow_out.columns:
            output_DF[col]=flow_out[col]
        output_units.update([(s+' outflow','mol/m2/sec') for s in self.primary_names])

#        output_DF['CEC H+']=pandas.DataFrame(self.output['CEC H+'],index=self.output['time'])
#        output_units['CEC H+']='mol/m^3'

        self.output_DF=output_DF.reset_index(drop=True).set_index(output_DF.index/(24*3600))
        self.output_units=output_units


    def run_onestep(self,chem,data,dt,status,min_dt=0.1,num_cuts=0,diffquo={},bc=None,flux_tol=10,truncate_concentration=0.0,rateconstants={},min_cuts=0,TRc={},EvR={},PwR=0.0):
        self.copy_to_alquimia(data)
        bc_tmp=bc
        converged=False
        porosity=data.state.porosity
        saturation=data.properties.saturation
        volume=data.properties.volume
        max_cuts=num_cuts
        actual_dt=dt/2**num_cuts

        for num,reactname in enumerate(get_alquimiavector(data.meta_data.aqueous_kinetic_names)):
            data.properties.aqueous_kinetic_rate_cnst.data[num]=rateconstants[reactname]

        for spec in diffquo.keys():
            pos=get_alquimiavector(data.meta_data.primary_names).index(spec)
            if data.state.total_mobile.data[pos] < truncate_concentration and data.state.total_mobile.data[pos] != 0.0:
                for num,reactname in enumerate(get_alquimiavector(data.meta_data.aqueous_kinetic_names)):
                    if '->' in reactname and spec in reactname.split('->')[0].split():  # reactants
                        data.properties.aqueous_kinetic_rate_cnst.data[num]=0.0

        cut_for_flux=False

        for spec in diffquo.keys():
            pos=get_alquimiavector(data.meta_data.primary_names).index(spec)
<<<<<<< HEAD

            bc_conc=bc.total_mobile.data[pos]
            local_conc=data.state.total_mobile.data[pos]

=======
            
            bc_conc=bc.total_mobile.data[pos]
            local_conc=data.state.total_mobile.data[pos]
            
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629
            if local_conc<0:
                raise RuntimeError('Initial concentration of %s < 0'%spec)

        if cut_for_flux or num_cuts<min_cuts:
            converged=False
        else:
            chem.ReactionStepOperatorSplit(ffi.new('void **',data.engine_state),actual_dt,ffi.addressof(data.properties),ffi.addressof(data.state),
                                            ffi.addressof(data.aux_data),status)
            data.state.porosity=porosity
            #data.properties.saturation=saturation
            converged=status.converged
            # Check for negative concentrations
            # if (numpy.array(get_alquimiavector(data.state.total_mobile))<0).any() or (numpy.array(get_alquimiavector(data.state.total_immobile))<0).any():
            #     converged=False

        if converged:
            check_status(status,False)

            chem.GetAuxiliaryOutput(ffi.new('void **',data.engine_state),ffi.addressof(data.properties),ffi.addressof(data.state),
                                            ffi.addressof(data.aux_data),ffi.addressof(data.aux_output),status)
            check_status(status,False)

        if converged:
            vol_tmp=volume*porosity*saturation
            vol_new=vol_tmp+PwR*actual_dt
            #print('water volume change %1.5f,%1.5f,%1.5f',PwR*actual_dt,vol_new,vol_tmp)
            data.properties.saturation=vol_new/(volume*porosity)
            for spec in diffquo.keys():
                if spec in TRc.keys():
                    pos=get_alquimiavector(data.meta_data.primary_names).index(spec)
                    if spec != 'SOM':
                        data.state.total_mobile.data[pos] = ((data.state.total_mobile.data[pos] + TRc[spec]*actual_dt)*vol_tmp*1000 + EvR[spec]*actual_dt)/(vol_new*1000)
                    else:
                        data.state.total_immobile.data[pos] = data.state.total_immobile.data[pos] + TRc[spec]*actual_dt
                    #if spec=='NO3-':
                    #    print('transport induced change %3.20f,%2.30f,%d'%(TRc[spec]*actual_dt,data.state.total_mobile.data[pos],actual_dt))
                else:
                    pos=get_alquimiavector(data.meta_data.primary_names).index(spec)
                    data.state.total_mobile.data[pos] = data.state.total_mobile.data[pos]
                    data.state.total_immobile.data[pos] = data.state.total_immobile.data[pos]
                #if spec in TRc.keys():
                #pos=get_alquimiavector(data.meta_data.primary_names).index(spec)
                #data.state.total_mobile.data[pos] = ((data.state.total_mobile.data[pos] + TRc[spec]*actual_dt)*vol_tmp*1000 + EvR[spec]*actual_dt)/(vol_new*1000)
                    ##if spec=='DOM1':
                    ##    print('transport induced change %2.30f,%2.30f,%2.30f,%d'%(EvR[spec]*actual_dt,TRc[spec]*actual_dt,data.state.total_mobile.data[pos],actual_dt))
                #else:
                    #pos=get_alquimiavector(data.meta_data.primary_names).index(spec)
                    #data.state.total_mobile.data[pos] = (data.state.total_mobile.data[pos]*vol_tmp*1000)/(vol_new*1000)
            self.copy_from_alquimia(data)
            return max_cuts
        else:

            if actual_dt/2<min_dt:
                if cut_for_flux:
                    raise RuntimeError('Pflotran failed to converge (because of boundary fluxes) after %d cuts to dt = %1.2f s'%(num_cuts,actual_dt))
                else:
                    raise RuntimeError('Pflotran failed to converge after %d cuts to dt = %1.2f s'%(num_cuts,actual_dt))
<<<<<<< HEAD

=======
 
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629
            #flux_tmp=numpy.zeros(data.meta_data.primary_names.size)
            # data will be reset to layer contents at beginning of next call
            # Run it twice, because we cut the time step in half
            ncuts=self.run_onestep(chem,data,dt,status,min_dt,num_cuts=num_cuts+1,diffquo=diffquo,bc=bc_tmp,truncate_concentration=truncate_concentration,rateconstants=rateconstants,TRc=TRc,EvR=EvR,PwR=PwR)
            # If this completes, it means that run_onestep was successful
            self.copy_from_alquimia(data)
            bc_tmp=bc

            if ncuts>max_cuts:
                max_cuts=ncuts

            # This starts from ncuts so it doesn't have to try all the ones that failed again
            for n in range(2**(ncuts-(num_cuts+1))):
                ncuts2=self.run_onestep(chem,data,dt,status,min_dt,num_cuts=ncuts,diffquo=diffquo,bc=bc_tmp,truncate_concentration=truncate_concentration,rateconstants=rateconstants,TRc=TRc,EvR=EvR,PwR=PwR)
<<<<<<< HEAD

=======
                
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629
                self.copy_from_alquimia(data)
                bc_tmp=bc
                if ncuts2>max_cuts:
                    max_cuts=ncuts2
<<<<<<< HEAD

=======
            
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629
            return max_cuts


##cjw remove leaf_Mn
def convert_to_xarray(layers,t0=0.0,drop_nas=True,convert_output=True):
    for l in layers:
        if convert_output or not hasattr(l,'output_DF'):
            l.convert_output()
    data_array = xarray.concat([xarray.Dataset.from_dataframe(layer.output_DF) for layer in layers],dim='layer').rename({'index':'time','layer':'depth'})
    data_array['dz']=xarray.DataArray([layer.volume for layer in layers],dims='depth',attrs={'units':'cm'})*100     #cjw why multiply by 100
    data_array['z_bottom']=data_array['dz'].cumsum()
    data_array['z_top']=data_array['z_bottom']-data_array['dz']
    data_array['z_middle']=data_array['z_bottom']-data_array['dz']/2
    data_array['z_bottom'].attrs['units']='cm'
    data_array['z_middle'].attrs['units']='cm'
    data_array['z_top'].attrs['units']='cm'

    data_array['depth']=data_array['z_middle']

    for var in layers[0].output_units:
        data_array[var].attrs['units']=layers[0].output_units[var]
# cjw: revise needed to construct layers for marsh site.
    # To do: Add other layer properties/attributes and also leaf Mn concentration
    data_array['saturation']=xarray.DataArray([layer.saturation for layer in layers],dims='depth',attrs={'units':'fraction'})
    data_array['BD']=xarray.DataArray([layer.BD for layer in layers],dims='depth',attrs={'units':'g cm-3'})
#    data_array['CEC']=xarray.DataArray([layer.CEC for layer in layers],dims='depth',attrs={'units':'meq kg-1'})

    data_array['time']=data_array['time']+t0
    data_array['time'].attrs['units']='days'

#cjw    if leaf_Mn is not None:
#cjw        data_array['litter_Mn']=xarray.DataArray(leaf_Mn,dims='litter_year',attrs={'units':'mmol/kg dry mass'})

    if drop_nas:
        data_array=data_array.dropna(dim='time')

    return data_array

def copy_to_layers(data_xarray,layers):
    for var in data_xarray.variables:
        for depth in range(len(data_xarray.depth)):
            if var.startswith('Total Sorbed'):
                layers[depth].total_immobile[var[len('Total Sorbed '):]]=data_xarray[var].dropna(dim='time').isel(time=-1,depth=depth).item()
            elif var.startswith('Total '):
                layers[depth].total_mobile[var[len('Total '):]]=data_xarray[var].dropna(dim='time').isel(time=-1,depth=depth).item()
            elif var.startswith('Free '):
                layers[depth].free_mobile[var[len('Free '):]]=data_xarray[var].dropna(dim='time').isel(time=-1,depth=depth).item()
            elif var.endswith('VF'):
                layers[depth].mineral_volume_fraction[var[:-3]]=data_xarray[var].dropna(dim='time').isel(time=-1,depth=depth).item()

rate_scale=2e-10
truncate_conc=1e-30
thresh=truncate_conc*1.01

reaction_network=Mar.make_network() # We will add the Mn along with leaf litter manually instead of generating it through decomposition

rateconstants={
'1.00e+00 DOM1  -> 3.33e-01 Acetate-  + 3.33e-01 HCO3-  + 6.67e-01 H+  + 1.33e+00 H2(aq)  + 3.33e-01 Tracer':rate_scale*1e2,
'1.00e+00 DOM1  + 1.00e+00 O2(aq)  -> 1.00e+00 HCO3-  + 1.00e+00 H+  + 1.00e+00 Tracer':rate_scale*1e2,
'4.00e+00 H2(aq)  + 1.00e+00 HCO3-  + 1.00e+00 H+  -> 1.00e+00 CH4(aq)  + 3.00e+00 H2O':rate_scale*0.31,
'1.00e+00 Acetate-  + 2.00e+00 NO3-  + 1.00e+00 H+  -> 2.00e+00 HCO3-  + 1.00e+00 N2(aq)  + 0.00e+00 N2O(aq)  + 2.00e+00 H2O  + 2.00e+00 Tracer':rate_scale*0.0036,
'1.00e+00 Acetate-  + 2.00e+00 O2(aq)  -> 2.00e+00 HCO3-  + 2.00e+00 H+  + 2.00e+00 Tracer':rate_scale*1e2,
'1.00e+00 CH4(aq)  + 1.00e+00 O2(aq)  -> 1.00e+00 HCO3-  + 1.00e+00 H+  + 1.00e+00 H2O  + 1.00e+00 Tracer':rate_scale*2.31,
'2.00e+00 NH4+  + 4.00e+00 O2(aq)  -> 2.00e+00 NO3-  + 2.00e+00 H2O  + 4.00e+00 H+':rate_scale*10,
#'1.00e+00 Acetate- + 1.00e+00 NO3- + 1.00e+00 H2O + 1.00e+00 H+ -> 1.00e+00 NH4+ + 2.00e+00 HCO3- + 2.00e+00 Tracer':rate_scale,
#'1.00e+00 DOM1 + 5.00e-01 NO3- + 5.00e-01 H2O <-> 5.00e-01 NH4+ + 1.00e+00 HCO3-':rate_scale,
#'1.00e+00 DOM1 + 5.00e-01 NO3- + 5.00e-01 H2O -> 5.00e-01 NH4+ + 1.00e+00 HCO3-':rate_scale,
'1.00e+00 DOM1  + 5.00e-01 NO3-  + 5.00e-01 H2O  -> 5.00e-01 NH4+  + 1.00e+00 HCO3-':rate_scale,
#'1.00e+00 HS- + 1.00e+00 NO3- + 1.00e+00 H+ + 1.00e+00 H2O <-> 1.00e+00 SO4-- + 1.00e+00 NH4+':rate_scale*0.1,
#'1.00e+00 HS- + 1.00e+00 NO3- + 1.00e+00 H+ + 1.00e+00 H2O -> 1.00e+00 SO4-- + 1.00e+00 NH4+':rate_scale*0.1,
'1.00e+00 HS-  + 1.00e+00 NO3-  + 1.00e+00 H+  + 1.00e+00 H2O  -> 1.00e+00 SO4--  + 1.00e+00 NH4+':rate_scale*0.01,
'1.00e+00 Acetate-  + 8.00e+00 Fe+++  -> 2.00e+00 HCO3-  + 8.00e+00 Fe++  + 9.00e+00 H+  + 2.00e+00 Tracer':rate_scale,
'1.00e+00 Fe++  + 2.50e-01 O2(aq)  + 1.00e+00 H+  <-> 1.00e+00 Fe+++  + 5.00e-01 H2O':1.0e0*1.0e1*0,
'1.00e+00 Fe++  + 2.50e-01 O2(aq)  + 1.00e+00 H+  -> 1.00e+00 Fe+++  + 5.00e-01 H2O':rate_scale*100,
'1.00e+00 CH4(aq)  + 8.00e+00 Fe+++  + 3.00e+00 H2O  -> 1.00e+00 HCO3-  + 8.00e+00 Fe++  + 9.00e+00 H+  + 1.00e+00 Tracer':rate_scale*0.01,
'1.00e+00 CH4(aq)  + 1.00e+00 NO3-  -> 1.00e+00 HCO3-  + 1.00e+00 NH4+  + 1.00e+00 Tracer':rate_scale*0.01,
'1.00e+00 Acetate-  + 1.00e+00 SO4--  -> 2.00e+00 HCO3-  + 2.00e+00 HS-  + 2.00e+00 Tracer':rate_scale*0.01,
'1.00e+00 CH4(aq)  + 1.00e+00 SO4--  -> 1.00e+00 HCO3-  + 1.00e+00 HS-  + 1.00e+00 H2O  + 1.00e+00 Tracer':rate_scale,
<<<<<<< HEAD
'1.00e+00 Acetate-  -> 1.00e+00 CH4(aq)  + 1.00e+00 HCO3-  + 1.00e+00 Tracer':rate_scale*10,
=======
'1.00e+00 Acetate-  -> 1.00e+00 CH4(aq)  + 1.00e+00 HCO3-  + 1.00e+00 Tracer':rate_scale*100,
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629
'SOM decay to CO2 (SOMDEC sandbox)': 1e-6,
'SOM decay to DOM1 (SOMDEC sandbox)':1e-7,
}

#cjw
input_file='deltamarsh.in'

decomp_network.PF_network_writer(reaction_network).write_into_input_deck('SOMdecomp_template.txt',input_file,log_formulation=False,truncate_concentration=truncate_conc)
#cjwdecomp_network.PF_network_writer(reaction_network).write_into_input_deck('SOMdecomp_template.txt','deltamarsh.in',length_days=30,log_formulation=False)

#cjw incubation_length=5 # Years of litter decomp

#cjw molar_mass={'N':14.007}
#cjw  molar_mass={'Mg++':24.305,'Al+++':26.982,'K+':39.098,'Ca++':40.078,'Mn':54.938,'Na+':22.99,'N':14.007}
#cjw  init_exch_cations={'Mg++':1.5,'Al+++':7.0,'K+':1.3,'Ca++':5.0,'Na+':0.2,'Mn++':0.3} # mmol/kg. From Jin et al 2010 Table 3

# Read secondary complex names from input file since Alquimia does not provide them

##cjw run_alquimia.run_simulation already include this

with open(input_file,'r') as infile:
    secondary_names=[]
    for line in infile:
        if 'SECONDARY_SPECIES' in line.split('#')[0]:
            break
    for line in infile:
        l=line.strip().split('#')[0]
        if l.startswith('END') or l.startswith('/'):
            break
        if len(l)>0:
            secondary_names.append(l)

import time
starting_time=time.time()

rateconstants_warmed=rateconstants.copy()
#cjw pH for fresh, brackish and salt is from DeLaune, 1983 (average data from 0-50cm).
for pH in numpy.arange(6.3,7.3):

    chem,data,sizes,status=init_alquimia(input_file,hands_off=False)
        #cjw below is the codes for layer setting up.
        # Set up layers
        # Top (organic) layer should be thinner and have lower bulk density though
        # Low bulk density causes simulation to slow or crash though. Actually CEC being too low (<100 combined with BD<1) is the problem
#cjw        layers=[layer(0.05,rateconstants=rateconstants_warmed,BD=0.05,porosity=0.5)]+[layer(0.1,BD=0.25,rateconstants=rateconstants_warmed) for num in range(3)]
    #layers=[layer(0.01,rateconstants=rateconstants_warmed,BD=0.08,porosity=0.95,saturation=0.8)]+[layer(0.05,BD=0.1,rateconstants=rateconstants_warmed,porosity=0.95,saturation=0.85)]+[layer(0.05,BD=0.15,rateconstants=rateconstants_warmed,porosity=0.92,saturation=0.85)]+[layer(0.05,BD=0.2,rateconstants=rateconstants_warmed,porosity=0.9,saturation=0.88)]+[layer(0.05,BD=0.25,rateconstants=rateconstants_warmed,porosity=0.85,saturation=0.88)]+[layer(0.1,BD=0.35,rateconstants=rateconstants_warmed,porosity=0.75,saturation=0.9)]
    #station 2825
<<<<<<< HEAD
    #layers=[layer(0.01,rateconstants=rateconstants_warmed,BD=0.15,porosity=0.95,saturation=0.85)]+[layer(0.05,BD=0.15,rateconstants=rateconstants_warmed,porosity=0.95,saturation=0.85)]+[layer(0.05,BD=0.13,rateconstants=rateconstants_warmed,porosity=0.95,saturation=0.88)]+[layer(0.05,BD=0.1,rateconstants=rateconstants_warmed,porosity=0.96,saturation=0.89)]+[layer(0.05,BD=0.1,rateconstants=rateconstants_warmed,porosity=0.96,saturation=0.91)]+[layer(0.1,BD=0.09,rateconstants=rateconstants_warmed,porosity=0.97,saturation=0.92)]
    ##station 3166
    #layers=[layer(0.01,rateconstants=rateconstants_warmed,BD=0.04,porosity=0.98,saturation=0.93)]+[layer(0.05,BD=0.06,rateconstants=rateconstants_warmed,porosity=0.98,saturation=0.91)]+[layer(0.05,BD=0.07,rateconstants=rateconstants_warmed,porosity=0.97,saturation=0.91)]+[layer(0.05,BD=0.07,rateconstants=rateconstants_warmed,porosity=0.97,saturation=0.92)]+[layer(0.05,BD=0.07,rateconstants=rateconstants_warmed,porosity=0.97,saturation=0.92)]+[layer(0.1,BD=0.06,rateconstants=rateconstants_warmed,porosity=0.98,saturation=0.92)]
    ##station 3169
    layers=[layer(0.01,rateconstants=rateconstants_warmed,BD=0.27,porosity=0.90,saturation=0.74)]+[layer(0.05,BD=0.18,rateconstants=rateconstants_warmed,porosity=0.93,saturation=0.81)]+[layer(0.05,BD=0.17,rateconstants=rateconstants_warmed,porosity=0.93,saturation=0.82)]+[layer(0.05,BD=0.11,rateconstants=rateconstants_warmed,porosity=0.96,saturation=0.88)]+[layer(0.05,BD=0.07,rateconstants=rateconstants_warmed,porosity=0.97,saturation=0.91)]+[layer(0.1,BD=0.06,rateconstants=rateconstants_warmed,porosity=0.97,saturation=0.92)]

=======
    layers=[layer(0.01,rateconstants=rateconstants_warmed,BD=0.15,porosity=0.95,saturation=0.85)]+[layer(0.05,BD=0.15,rateconstants=rateconstants_warmed,porosity=0.95,saturation=0.85)]+[layer(0.05,BD=0.13,rateconstants=rateconstants_warmed,porosity=0.95,saturation=0.88)]+[layer(0.05,BD=0.1,rateconstants=rateconstants_warmed,porosity=0.96,saturation=0.89)]+[layer(0.05,BD=0.1,rateconstants=rateconstants_warmed,porosity=0.96,saturation=0.91)]+[layer(0.1,BD=0.09,rateconstants=rateconstants_warmed,porosity=0.97,saturation=0.92)]
    ##station 3166
    #layers=[layer(0.01,rateconstants=rateconstants_warmed,BD=0.04,porosity=0.98,saturation=0.93)]+[layer(0.05,BD=0.06,rateconstants=rateconstants_warmed,porosity=0.98,saturation=0.91)]+[layer(0.05,BD=0.07,rateconstants=rateconstants_warmed,porosity=0.97,saturation=0.91)]+[layer(0.05,BD=0.07,rateconstants=rateconstants_warmed,porosity=0.97,saturation=0.92)]+[layer(0.05,BD=0.07,rateconstants=rateconstants_warmed,porosity=0.97,saturation=0.92)]+[layer(0.1,BD=0.06,rateconstants=rateconstants_warmed,porosity=0.98,saturation=0.92)]
    ##station 3169
    #layers=[layer(0.01,rateconstants=rateconstants_warmed,BD=0.27,porosity=0.90,saturation=0.74)]+[layer(0.05,BD=0.18,rateconstants=rateconstants_warmed,porosity=0.93,saturation=0.81)]+[layer(0.05,BD=0.17,rateconstants=rateconstants_warmed,porosity=0.93,saturation=0.82)]+[layer(0.05,BD=0.11,rateconstants=rateconstants_warmed,porosity=0.96,saturation=0.88)]+[layer(0.05,BD=0.07,rateconstants=rateconstants_warmed,porosity=0.97,saturation=0.91)]+[layer(0.1,BD=0.06,rateconstants=rateconstants_warmed,porosity=0.97,saturation=0.92)]
    
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629
    suldpth=[0.01,0.06,0.11,0.16,0.21,0.31]   #depth in meters for sulfur inicond profile
    sulidx=0
    for l in layers:
        l.secondary_names=secondary_names
    for l in layers:
#           l.initcond=Mar.pools.copy()
        l.initcond=decomp_network.change_constraints(Mar.pools,{'H+':'%1.1f P'%pH,'Fe++':'1e-2',
                                                        'Ca++':'1e-30','Na+':'1e-30','Cl-':'1e-30','Fe+++':'1e-15','CH4(aq)':'1e-15',
                                                        'SO4--':'1e-30',#})     #'%1.8f'%((0.75-suldpth[sulidx])*1e-15/30)})
                                                        'O2(aq)':'0.0002'})               #cjw: set the sulfate profile so4 decrease with depth (0.75-suldpth[sulidx])/30)
        sulidx=sulidx+1
#cjw    for l in layers[0:3]:
#cjw        l.initcond=decomp_network.change_constraints(Mar.pools,{'Fe++':'%1.8f'%(l.volume/100)})      #cjw Fe++ profile
        # initcond=decomp_network.change_site_density(initcond, '>DOM1', 1e4)
#cjw        sed_air_interface=layers[0].initcond
#cjw        sed_water_interface=layers[0].initcond
    pools_tide=Mar.pools.copy()
    for n,p in enumerate(pools_tide):
        if p['name']=='O2(aq)':
            pools_tide[n]=pools_tide[n].copy()
            pools_tide[n].update(constraints={'initial':'0.2 G O2(g)'})
        #elif p['name']=='SO4--':
        #    pools_tide[n]=pools_tide[n].copy()
        #    pools_tide[n].update(constraints={'initial':'1e-15'})
        #elif p['name']=='Ca++':
        #    pools_tide[n]=pools_tide[n].copy()
        #    pools_tide[n].update(constraints={'initial':'1e-15'})
        #elif p['name']=='CH4(aq)':
        #    pools_tide[n]=pools_tide[n].copy()
        #    pools_tide[n].update(constraints={'initial':'1e-15'})
        #elif p['name']=='Fe+++':
        #    pools_tide[n]=pools_tide[n].copy()
        #    pools_tide[n].update(constraints={'initial':'1e-15'})
    #bc=pools_tide
    bc=layers[0].initcond

#cjw    bc=None
    dt=3600
    repyr=11
    nyears=1*repyr
    nsteps=365*24//(dt//3600)*nyears
        # Set up initial condition
    for l in layers:
            # Initialize state data
        data.properties.volume=l.volume
        data.properties.saturation=l.saturation

        data.state.temperature=l.temperature
        data.state.water_density=l.water_density
        data.state.porosity=l.porosity
        data.state.aqueous_pressure=l.pressure

<<<<<<< HEAD
        #print('layer porosity,BD: %1.4f,%1.4f'%(l.porosity,l.BD))
=======
        #print('layer porosity,BD: %1.4f,%1.4f'%(l.porosity,l.BD)) 
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629


            # Set properties: surface site density and mineral rate constants
            # This is necessary when running in hands-on mode
        if l.initcond is not None:
            for constraint in l.initcond:
                if constraint['kind']=='surf_complex':
                    sitename=constraint['name']
                    sitenum=get_alquimiavector(data.meta_data.surface_site_names).index(sitename)
                    print('Applying site density: {name:s} (position={pos:d}): {dens:1.1f}'.format(name=sitename,pos=sitenum,dens=constraint['site_density']))
                    data.state.surface_site_density.data[sitenum]=constraint['site_density']
                elif constraint['kind']=='mineral':
                    name=constraint['name']
                    num=get_alquimiavector(data.meta_data.mineral_names).index(name)
                    data.properties.mineral_rate_cnst.data[num]=float(constraint['rate'].split()[0].replace('d','e'))

            # Set up boundary condition if applicable
        if bc is not None:
            bc_cond=convert_condition_to_alquimia(bc,'initial')
            bc_state=ffi.new('AlquimiaState *')
            bc_state.temperature=l.temperature #layers[0].temperature           #replace layers[0]
            bc_state.water_density=l.water_density #layers[0].water_density
            bc_state.porosity=l.porosity #layers[0].porosity
            bc_state.aqueous_pressure=l.pressure #layers[0].pressure
            bc_auxdata=ffi.new('AlquimiaAuxiliaryData *')
            lib.AllocateAlquimiaState(sizes,bc_state)
            for num in range(data.state.surface_site_density.size):
                bc_state.surface_site_density.data[num]=data.state.surface_site_density.data[num]
            lib.AllocateAlquimiaAuxiliaryData(sizes,bc_auxdata)
            chem.ProcessCondition(ffi.new('void **',data.engine_state),bc_cond,ffi.addressof(data.properties),bc_state,bc_auxdata,status)
            check_status(status,False)
        else:
            bc_state=None

            # Aqueous kinetic rate constants also need to be specified in hands-off mode
            # Alquimia interface always sets backward rate to zero
        for num,reactname in enumerate(get_alquimiavector(data.meta_data.aqueous_kinetic_names)):
            data.properties.aqueous_kinetic_rate_cnst.data[num]=l.rateconstants[reactname]

        init_cond=convert_condition_to_alquimia(l.initcond,'initial')
        chem.ProcessCondition(ffi.new('void **',data.engine_state),init_cond,ffi.addressof(data.properties),ffi.addressof(data.state),ffi.addressof(data.aux_data),status)
        check_status(status,False)

            # Pflotran sets porosity based on minerals or something? Needs to be reset
        data.state.porosity=l.porosity

        l.copy_from_alquimia(data)
        l.setup_output(nsteps+1,dt)

            ##### At this point, the model should be initialized ##########
        print('''

            *****************************************************
            Successfully initialized alquimia geochemical engine
            *****************************************************

        ''')

    print('''


            *************************************
            * Starting simulation with pH = %1.1f, Ndep = %03d, Warming = %d *
            *************************************


        '''%(pH,0,20))        #cjw remove warming



    #initial_HCO3 = l.total_mobile['HCO3-']
    #initial_O2 = l.total_mobile['O2(aq)']

        # Flow rate cm/s = 10 L/m2/s, positive is downward
#    flow_rate=numpy.linspace(1e-6,1e-7,len(layers)) # Rate declines linearly with depth, assumes removal or accumulation in lower layers
    flow_rate=numpy.zeros(len(layers))*1e-7
#    flow_rate=numpy.linspace(1e-7,1e-8,len(layers))
    min_dt=0.1
    truncate_concentration=1e-20

#cjw update wet-dry condition for dynamic boundary conditions
##if water level > marsh elevation, then bc=sed_water_interface, else bc=sed_air_interface
###cjw compute tide water level below to
    import scipy as sp
    import mpmath as mp
    import numpy as np
    import math
    import pandas as pd
    import operator
    import seaborn as sns

<<<<<<< HEAD
##cjw parameter needed for transport
=======
##cjw parameter needed for transport 
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629
    tide_t=range(0,nsteps)
    Nlyr=range(0,len(layers))
    tdt=3600

<<<<<<< HEAD
    wtl=pd.read_csv('./10WaterLevel_3169.csv')
    salinity=pd.read_csv('./10Salinity_3169.csv')
=======
    wtl=pd.read_csv('./WaterLevel_2825.csv')
    salinity=pd.read_csv('./Salinity_2825.csv')
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629

#c    abmsl=2.59 #units in meter
#c    alpha=[163,154.6,176.1,153.8,37.4,30.8,37.2,19.2]    #phase angle in degrees M2,S2,N2,K2,K1,O1,P1,Q1
#c    aM2=[0.013,0.007,0.005,0.002,0.114,0.114,0.036,0.025]       #amplitude in meters M2,S2,N2,K2,K1,O1,P1,Q1
#c    w=[28.984104,30.0,28.43973,30.082138,15.041069,13.943035,14.958931,13.398661] #angual speed degree/hour
#c    tide_t=range(0,nsteps)        #t in hours
#c    Zs=2.59+0.1   #marsh elevation above mean sea level in meters
    Zs=0.0
    Zt=[]
    sf=[]       #salinity in tide
#c    for i in tide_t:
#c        sz=abmsl
#c        for j in range(0,len(w)):
#c            z=aM2[j]*math.cos(w[j]*math.pi/(180*dt)*i*dt+alpha[j])
#c            sz=sz+z
#c        Zt+=[sz,]
##create salinity time series in tide, this is just an arbitrary pattern for testing
#c    for i in tide_t:
#c        sal = 36*math.cos(w[1]*i+alpha[1])
#c        sf+=[abs(sal),]
    Zt=wtl['WaterLevel'].values.tolist()
    sf=salinity['Salinity'].values.tolist()
<<<<<<< HEAD
##replace NAN value with mean value
    rep_wtl=np.nanmean(Zt)
=======
##replace NAN value with mean value 
    rep_wtl=np.nanmean(Zt)   
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629
    rep_sal=np.nanmean(sf)

    Zt=[rep_wtl if math.isnan(x) else x for x in Zt]
    sf=[rep_sal if math.isnan(x) else x for x in sf]
<<<<<<< HEAD
##subset data to match rainfall data, here subset data to year 2009.
    Zt=[x for x in Zt[0:24*365]]
    sf=[x for x in sf[0:24*365]]

=======
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629
    if len(Zt) < nsteps//repyr:               #run twice of the observation data.
        Zt=Zt+[rep_wtl]*(nsteps//repyr-len(Zt))
        sf=sf+[rep_sal]*(nsteps//repyr-len(sf))
    Zt=Zt*repyr
    sf=sf*repyr

#cjw create a tidal forcing concentration of compounds in tide for water-sediment interface boundary
##tconc_frac is the fraction of each primary species relative to salinity in water.
#    rstcl=0.14   #unit is mg/ml:mg/ml
<<<<<<< HEAD
    salfc=0.00180665
=======
#    salfc=0.00180665
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629
#ppt to sulfate concetration
#    sulfc=salfc/0.14 ## sulfate concetration is mg/L, and 1e-3*mg/L/molar_mass=M
    tconc_fc={
    #'DOM1':0.01,
    #'Acetate-':0.01,
<<<<<<< HEAD
    'HCO3-':0,#0.01,
    'O2(aq)':0,
    'NO3-':0,#1.29046e-10,
    'NH4+':0,#1.29046e-10,
    'SO4--':0.14*1e-3/salfc, #1.29046e-2,
    'Cl-':1e-3/salfc,
    'Na+':0.5769*1e-3/salfc,     ##include potassium 1.13%,sodium 30.6%
    'Ca++':0.024*1e-3/salfc,
    'HS-':1e-6*1e-3/salfc,
    #'Fe+++':1e-15,
    #'Fe++':1e-2,
    'CH4(aq)':0,#0.1,
    'H2(aq)':0,#0.01,
    'N2(aq)':0,#0.01,
    'N2O(aq)':0,#0.001,
    'H+':0,#0.01,         #this could be related to the alklinity in porewater and tide
    'Tracer':0,#0.01,
=======
    'HCO3-':0,
    'O2(aq)':0,
    'NO3-':1.29046e-10,
    'NH4+':1.29046e-10,
    'SO4--':1.80665e-3*0.14, #1.29046e-2, 
    'Cl-':1.80665e-3,
    'Na+':1.80665e-3*0.5769,     ##include potassium 1.13%,sodium 30.6%
    'Ca++':1.80665e-3*0.024,
    'HS-':1e-8,
    #'Fe+++':1e-15,
    #'Fe++':1e-2,
    'CH4(aq)':0,
    'H2(aq)':0,
    'N2(aq)':0,
    'N2O(aq)':0,
    'H+':0,         #this could be related to the alklinity in porewater and tide
    'Tracer':0,
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629
    }

    molar_mass={#'SOM':180.156,
    #'HRimm':59.052,
    #'DOM1':180.156,
    #'Acetate-':59.052,
    'HCO3-':61.0168,
    'O2(aq)':31.998,
    'NO3-':62.0049,
    'NH4+':18.039,
    'SO4--':96.06,
    'Cl-':35.454,
    'Na+':22.98977,
    'Ca++':40.0780,
    'HS-':33.1,
    #'Fe+++':55.845,
    #'Fe++':55.845,
    'CH4(aq)':16.04,
    'H2(aq)':2,
    'N2(aq)':28.0134,
    'N2O(aq)':44.013,
    'H+':1,         #this could be related to the alklinity in porewater and tide
    'Tracer':1,
    }   ##molar mass is g/mol
    tide_conc={'salt':sf}
    for spec in tconc_fc.keys():
        if spec=='O2(aq)':
            tmp=[(i*tconc_fc[spec]+9*1e-3)/molar_mass[spec] for i in sf]
        else:
            tmp=[i*tconc_fc[spec]/molar_mass[spec] for i in sf]
        tide_conc[spec]=tmp
<<<<<<< HEAD

=======
    
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629
    ## evaptranspiration
    #t=range(0,24*365) #t in hours, dt is time step in hours
    tday=range(0,365*1)
    Etd=[]
    Eah=[]
    Rg=[]
    for i in tday:
        Etmp=3.3+2.56*math.sin(2*math.pi*i/365-1.74)  #unit is l/m2/day
        Etd+=[Etmp,]

    for j in range(0,24*365):
        dj=j-24*(math.ceil(j/24)-1)
        day=math.ceil(j/24)-1
        if dj > 6 and dj <= 18:
            Ehtmp=(-(dj**2-24*dj+144)*Etd[day]/288+0.125*Etd[day])*1e-3     ##unit is m3/h per unit area m2
        elif dj < 6 or dj > 18:
<<<<<<< HEAD
            Ehtmp=0.0
=======
            Ehtmp=0.0   
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629
        Eah+=[Ehtmp,]                    ##unit is m3/h per unit area m2.
##compute pulse rain fall
    #N_rain = 24*365*1 # sample count total hours for a year
    #P_rain = 24*3.15  # period one event every 3.15 days
    #D_rain = 2   # width of pulse rain event last for 2
    #sg_rain = np.arange(N_rain) % P_rain < D_rain
    #for i in range(0,24*365):
    #    if sg_rain[i] == True:
    #        rain = 0.001#10*1e-10 #8.6      # unit is m3/h per unit area m2
    #    elif sg_rain[i] == False:
<<<<<<< HEAD
    #        rain = 0.0
    #    Rg+=[rain,]
    prcp=pd.read_csv('./precip_h.csv')
    ptmp=prcp['HPCP'].values.tolist()
    Rg=[x*0.001 for x in ptmp[24*366:24*(366+365)]]

=======
    #        rain = 0.0      
    #    Rg+=[rain,]
    prcp=pd.read_csv('./precip_h.csv')
    ptmp=prcp['HPCP'].values.tolist()
    Rg=[x*0.001 for x in ptmp[24*366:24*(365+366)]]
    
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629
    Eah=Eah*repyr
    Rg=Rg*repyr
#Vp:volume of pore water in 30cm pedon(unit 1/m2)
#Rg: is the empeirically determined rainfall amount in unit of l/m2/h (set rainfall to a constant number of 3.4 l/m2/day=3.4/24 l/m2/h)
<<<<<<< HEAD
#fc: is the field capacity(dimensionless quotient of porewater volume
#    at field capacity to the saturated volume) 0.934+/-0.008 (0.966+/-0.029)
#V0 is the saturated water volume (l/m2)or porosity of the 30-cm pedon
#sd is the specific drainage rate (h-1). sd=K/(n*L), K is the hydraulic conductivities in cm/h,
=======
#fc: is the field capacity(dimensionless quotient of porewater volume 
#    at field capacity to the saturated volume) 0.934+/-0.008 (0.966+/-0.029)
#V0 is the saturated water volume (l/m2)or porosity of the 30-cm pedon
#sd is the specific drainage rate (h-1). sd=K/(n*L), K is the hydraulic conductivities in cm/h, 
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629
#   (1.45cm/h, 0.21cm/h), n is porosities (0.43,67.7) which is dimensionless,L is the depth or length of the cores (30cm)
#    based on the suggesting values at Goat Island and Oyster Landing, sd here iare 0.0103h-1 and 0.1114h-1.
# the following equations are subject to the constraint that Vp can never exceed V0.
    fc=0.9#0.934
    #sd=0.0#0.0103 # h-1
    #Rg=Eah*5  #3.4/24 # l/m2/h
    #V0=layers[num].volume*layers[num].porosity #field capacity of water volum for each layer with 1m2 surface unit is M3/m2. 300*0.677    #15/24 #l/m2, porosity n=V_v/V_t=0.67, V_t=0.3m*1m^2=0.3m^3=300l
    #dt=1 #time step is 1 hour
<<<<<<< HEAD
    #Vp=[]
=======
    #Vp=[]    
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629
    #Rg=[]
#salt excretion
#Se is the rate of salt excretion in unit of g.m-2.h-1
#Sp is the salt content in 30cm soil pedon in unit of g/m2
<<<<<<< HEAD
#k2 is the ion secretion factor
#   and is expressed as a fraction of the theoretical maximum ion uptake
#   that could be taken up passively in the transpiration stream, Here k2 is set to 0.03
#k3 is the quotient of transpiration/evapotranspiration
=======
#k2 is the ion secretion factor 
#   and is expressed as a fraction of the theoretical maximum ion uptake 
#   that could be taken up passively in the transpiration stream, Here k2 is set to 0.03
#k3 is the quotient of transpiration/evapotranspiration 
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629
#   and is set equal to 0.5 for simulation of high marsh sites dominated by dwarf spartina
    k2=0.03
    k3=0.5


#Se[t]=k2*k3*Sp[t]*Eah[t]/Vp[t]

#salt is lost by drainage and secretion
#Sp is salt content contained within 30cm deep soil section
#if Vp[t]>fc*V0
#Sp[t+1]=Sp[t]-(Sp[t]*sd+Se[t])*dt
#else if  Vp[t]<=fc*V0
#Sp[t+1]=Sp[t]-Se[t]*dt
#salt ppt to g/m2 conversion factor 1ppt=1000mg/L=1000g/m3, pedon is 30cm depth
## high tide exchanges
#sf is the salinity of infiltrating flood water in unit of g/L.
<<<<<<< HEAD
#kd is an empirically derived diffusion coeficient which is 0.0725 cm2/h. It was determined by
#    kd=0.0725
#    Lgth=1    #30/10   #core length in cm/10
#    cfc=1
#    Sini=3906         #salt mass content g/m2=ppt*Vp=30ppt


=======
#kd is an empirically derived diffusion coeficient which is 0.0725 cm2/h. It was determined by 
#    kd=0.0725
#    Lgth=1    #30/10   #core length in cm/10 
#    cfc=1
#    Sini=3906         #salt mass content g/m2=ppt*Vp=30ppt

 
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629

    #    soillyr=[l.volume for l in layers]
    zbio=0.1                #maximum bioturbation depth in meter below which bioturbation decrease exponentially
    thick=[l.volume for l in layers]
    mid_thick=[i/2 for i in thick]
    zdpth=np.array([l.volume for l in layers]).cumsum()
    mid_dpth=list(map(operator.sub,zdpth,mid_thick))

    #wInf=(3.3*10**(-0.87478367-0.00043512*zdpth[len(layers)-1]))*0.01/(3600*24*365)   ## advection rate or sedimentation rate in cm/yr -> m/s and zdpth in meter is the depth of the bottom layer
    porInf=0.8               ## porosity at the bottom layer or deepest depth of the column
    por0=0.95
    coeffp=4*0.01                 ## unit in cm, coeffecient for exponential porosity change -> change unit to m
    por=porInf+(por0-porInf)*np.exp([i*(-1/coeffp) for i in mid_dpth])  #porosity within layers
    pori=porInf+(por0-porInf)*np.exp([i*(-1/coeffp) for i in zdpth])  #porosity at interfaces
<<<<<<< HEAD

=======
    
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629
    alpha=np.array([0]+[(thick[i]+thick[i+1])/2 for i in range(len(thick)-1)])   ##[0] create a 0 value for alpha to have an extra column or element ie without [0], alpha will have 3 element, with it will have 4 element

    Temp=0                   ## temperature in degree.
    coeffDb=1*0.01                # coefficient for exponential bioturbation decrease unit in cm ->m
    sig=0.3*1e-3             ## diffusion boundary layer thickness in meter (which means it is 0.3 mm)
<<<<<<< HEAD
    #Db0=15*(wInf*(365*24*3600)/0.01)**0.6*(1e-4/(24*3600))     ## bioturbation diffusion coeff in unit m^2/s (cm^2/yr -> m^2/s) above zbio
=======
    #Db0=15*(wInf*(365*24*3600)/0.01)**0.6*(1e-4/(24*3600))     ## bioturbation diffusion coeff in unit m^2/s (cm^2/yr -> m^2/s) above zbio 
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629

    #ad=0.0336                     ## ion-specific coefficient with unit of cm2/d/degree C for diffusion coeff, ad is different among substrates -> change unit to m2/s/degree C
    adFc=0.0001/(24*3600)
    ad={'DOM1':0.0336,
    'Acetate-':0.0336,
    'HCO3-':0.0336,
    'O2(aq)':0.0386,
    'NO3-':0.0336,
    'NH4+':0.0336,
<<<<<<< HEAD
    'SO4--':0.0336, #1.29046e-2,
=======
    'SO4--':0.0336, #1.29046e-2, 
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629
    'Cl-':0.0336,
    'Na+':0.0336,     ##include potassium 1.13%,sodium 30.6%
    'Ca++':0.0336,
    'HS-':0.0242,
    'Fe+++':0.0336,
    'Fe++':0.0336,
    'CH4(aq)':0.0336,
    'H2(aq)':0.0336,
    'N2(aq)':0.0336,
    'N2O(aq)':0.0336,
    'H+':0.0336,         #this could be related to the alklinity in porewater and tide
    'Tracer':0.0336,
    }
    for nm in ad.keys():
        ad[nm]=ad[nm]*adFc
    #D0=1e-9                  ## zero-degree coefficient m2/s, change among substrates
    DFc=1.1574e-9             ##diffusion convert factor from cm2/day to m2/s
    D0={'DOM1':0.2,
    'Acetate-':0.2,
    'HCO3-':0.415589,
    'O2(aq)':0.955,
    'NO3-':0.845*0.01,
    'NH4+':0.847*0.01,
<<<<<<< HEAD
    'SO4--':0.432, #1.29046e-2,
=======
    'SO4--':0.432, #1.29046e-2, 
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629
    'Cl-':0.432,
    'Na+':0.432,     ##include potassium 1.13%,sodium 30.6%
    'Ca++':0.432,
    'HS-':0.842,
    'Fe+++':0.432,
    'Fe++':0.432,
    'CH4(aq)':0.955,
    'H2(aq)':0.955,
    'N2(aq)':0.955,
    'N2O(aq)':0.955,
    'H+':0.415589,         #this could be related to the alklinity in porewater and tide
    'Tracer':0.415589,
    }
    for nm in ad.keys():
        D0[nm]=D0[nm]*DFc
#    Dt={}
#    Dsed={}
#    Dsedi={}
<<<<<<< HEAD

    #Dt= D0 + ad*Temp                  ## free solution diffusion coefficient at ambient temperature T
    #Dts=Dt*por**2            ## diffusion coefficient within layers.
    #Dtsi=Dt*pori**2          ## diffusion coefficient at interface layers

=======
    
    #Dt= D0 + ad*Temp                  ## free solution diffusion coefficient at ambient temperature T 
    #Dts=Dt*por**2            ## diffusion coefficient within layers.
    #Dtsi=Dt*pori**2          ## diffusion coefficient at interface layers
  
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629
    #Dsedb=[Fir[i]*Dts[i] for i in range(len(Dts))]                   ## enhanced difussion coeff within layer with bioturbation
    #Dsedbi=[Firi[i]*Dtsi[i] for i in range(len(Dtsi))]               ## enhanced difussion coeff at layer interface with bioturbation

    #Dsed=[Fir[i]*Dts[i]+Db0*np.exp(-1*(mid_dpth[i]-zbio)/coeffDb) for i in range(len(Dts))]                    ## enhanced diffusion coeff within layer below bioturbation
    #Dsedi=[Firi[i]*Dtsi[i]+Db0*np.exp(-1*(zdpth[i]-zbio)/coeffDb) for i in range(len(Dtsi))]                ## enhanced diffusion coeff at layer interface below bioturbation depth
<<<<<<< HEAD

    #Db0=15*wInf**0.6*(1e-4/(365*24*3600))     ## bioturbation diffusion coeff in unit m^2/s (cm^2/yr -> m^2/s) above zbio
    #Dbz=[Db0*np.exp(-1*(mid_dpth[i]-zbio)/coeffDb) for i in range(len(mid_dpth))]    #bioturbation coeff within layer below zbio
    #Dbzi=[Db0*np.exp(-1*(zdpth[i]-zbio)/coeffDb) for i in range(len(zdpth))]          #bioturbation coeff at layer interface below zbio

    #Dsblk=Dsedb+Db0      ## diff coeff in model will be sum of Dbio + enhanced difussion coeff for dissolved species, for soild it will be Dbio


=======
    
    #Db0=15*wInf**0.6*(1e-4/(365*24*3600))     ## bioturbation diffusion coeff in unit m^2/s (cm^2/yr -> m^2/s) above zbio 
    #Dbz=[Db0*np.exp(-1*(mid_dpth[i]-zbio)/coeffDb) for i in range(len(mid_dpth))]    #bioturbation coeff within layer below zbio
    #Dbzi=[Db0*np.exp(-1*(zdpth[i]-zbio)/coeffDb) for i in range(len(zdpth))]          #bioturbation coeff at layer interface below zbio
    
    #Dsblk=Dsedb+Db0      ## diff coeff in model will be sum of Dbio + enhanced difussion coeff for dissolved species, for soild it will be Dbio
    
    
>>>>>>> e6c29c9268945aefa3c939cb7f808345c053a629
    #c=np.zeros((5,len(t)),dtype=np.float64)+0.01     #initial
    #c[:,0]=[0.05,0.04,0.03,0.02,0.01]
    #upbc=1  #boundary condition
### cjw Zt is tide water level
##update bc state by different layer property
#cjw    idxbc=0
    #Ini_tidec=tide_conc
    z=numpy.array([0]+[l.volume for l in layers]).cumsum()
    z_mid=(z[:-1]+z[1:])/2

#    for l in range(len(layers)):
#        for spec in layers[l].primary_names:
#            if spec in Dsed.keys():
#                layers[l].diffquo[spec]=Dsed[spec][l]
#        'SO4--':0.00001**((l+1)*0.75)}
#            'H2(aq)':0.001**((l+1)*0.75),
#            'N2(aq)':0.001**((l+1)*0.75),
#            'N2O(aq)':0.001**((l+1)*0.75)
#            }

#    layers[0].surface_site_density['>DOM1']=1e2

    for l in layers:
        l.write_output(0,dt)

    flow_in=numpy.zeros(len(layers)+1,dtype=float)
    flow_out=numpy.zeros(len(layers)+1,dtype=float)


    t0=time.time()
    tprev=t0

    tstart=0

        # Restart from existing state?
#cjw        restart_state='Mn_saved.nc'
    restart_state=None
    if isinstance(restart_state,str):
        restart_state=xarray.open_dataset(restart_state)
    if isinstance(restart_state,xarray.Dataset):
        copy_to_layers(restart_state, layers)
        tstart=restart_state.dropna(dim='time')['time'].isel(time=-1).item()
    elif isinstance(restart_state,list) and isinstance(restart_state[0],layer):
        layers=restart_state

#cjw tmp for bc state to be reset as the way it is at dry condition; set dumpy layers for flux diffusion
    #bc_tmp=bc_state
    #soillyr=[l.volume for l in layers]   #layer thickness
    excret_spec=['Cl-','Na+','SO4--','Ca++']
    gas_species={'O2(aq)':0.00027235027,#0.000282, O2 partial pressure/kH=0.20950/769.23
    'CH4(aq)':2.37998572e-09,#0.00175,#1e-9, partial pressure/ henry constant=1.7*1e-6/714.29 under standard condition T=25, Wania R. et al., 2010
    'H2(aq)':0.0006084006,#0.00078,#4.325e-10,   partial pressure/ henry constant=0.78/1282.05
    'N2(aq)':0.00068,#0.000714,
    'N2O(aq)':16.4*1e-9,}    ## air equilibrium concentration CO2 is 0.00039/29.41=1.32607956e-05
    immobile_species=['SOM']
    Fimb={'SOM':2e-9}    #immobile species flux at sediment water interface unit is mol/m3-bluk
    SOMtranR=numpy.zeros((len(layers),nsteps))
#cjw define diffusion dict
    Porw_vol=numpy.zeros((len(layers),nsteps),dtype=float)
    Porw_vol2=numpy.zeros((len(layers),nsteps),dtype=float)
    Rain_tmp=numpy.zeros((len(layers),nsteps),dtype=float)
    ##flux at sediment water surface
    Jswi=numpy.zeros((len(D0),nsteps),dtype=float)
    #Dt={}
    #Dsed={}
    #Dsedi={}
    ##cjw define dict for flux
#    Rc={}                     #flux between layers
#cjw question about the unit of total_mobile (mol/L?) and immobile species (immobile seems to have mol/m3-bulk)
    success=True
    for step in range(nsteps):
        #cjw define diffusion dict
        Dt={}
        Dsed={}
        Dsedi={}
        ##cjw define dict for flux
        Rc={}                     #flux between layers
        dzw=Zt[step]-Zs
        #dzw=abs(dzw)
#cjw diffusion changes with water depth
        if dzw >0:
            Fir=max(1,15.9*dzw**(-0.43))                     ## depth is in meter, diffusion enhancement factor to represent bio irrigation impact within layers
            wInf=(3.3*10**(-0.87478367-0.00043512*dzw))*0.01/(3600*24*365)   ## advection rate or sedimentation rate in cm/yr -> m/s and zdpth in meter is the depth of the bottom layer
            Db0=5.2*10**(0.76241122-0.00039724*dzw)*(0.0001/(365*24*3600))             ##unit is cm2/yr -> m2/s
        else:
            Fir=1                                            ## diffusion enhancement factor caused by bio irrigation at layer interface
            wInf=3.9e-11
            Db0=5.2*10**(0.76241122-0.00039724*dzw)*(0.0001/(365*24*3600))             ##unit is cm2/yr -> m2/s
        #Db0=15*(wInf*(365*24*3600)/0.01)**0.6/365*(1e-4/(24*3600))     ## bioturbation diffusion coeff in unit m^2/s (cm^2/d -> m^2/s) above zbio
        #Db0=5.2*10**(0.76241122-0.00039724*dzw)*(0.0001/(365*24*3600))             ##unit is cm2/yr -> m2/s
        Dbi=[Db0 if zdpth[i] <= zbio else Db0*np.exp(-1*(zdpth[i]-zbio)/coeffDb) for i in range(len(zdpth))]
        if dzw > 0:
            for spec in D0.keys():
                Dt[spec] = D0[spec]*0.01 + ad[spec]*Temp       ## free solution diffusion coefficient at ambient temperature T
                tmp=[Dt[spec]*por[i]**2*Fir+Db0 if mid_dpth[i] <= zbio else Dt[spec]*por[i]**2*Fir+Db0*np.exp(-1*(mid_dpth[i]-zbio)/coeffDb) for i in range(len(por))]
                tmpi=[Dt[spec]*pori[i]**2*Fir+Db0 if zdpth[i] <= zbio else Dt[spec]*pori[i]**2*Fir+Db0*np.exp(-1*(zdpth[i]-zbio)/coeffDb) for i in range(len(pori))]
                Dsed[spec]=tmp                            ## enhanced difussion coeff within layer with/without bioturbation
                Dsedi[spec]=tmpi                          ## enhanced diffusion coeff at layer interface below bioturbation depth
        else:
            for spec in D0.keys():
                Dt[spec] = D0[spec]*0.01 + ad[spec]*Temp       ## free solution diffusion coefficient at ambient temperature T
                tmp=[Dt[spec]*por[i]**2*Fir+Db0 if mid_dpth[i] <= zbio else Dt[spec]*por[i]**2*Fir+Db0*np.exp(-1*(mid_dpth[i]-zbio)/coeffDb) for i in range(len(por))]
                tmpi=[Dt[spec]*pori[i]**2*Fir+Db0 if zdpth[i] <= zbio else Dt[spec]*pori[i]**2*Fir+Db0*np.exp(-1*(zdpth[i]-zbio)/coeffDb) for i in range(len(pori))]
                Dsed[spec]=tmp                            ## enhanced difussion coeff within layer with/without bioturbation
                Dsedi[spec]=tmpi                          ## enhanced diffusion coeff at layer interface below bioturbation depth
#cjw compute diffusion for each timpstep for each layer
        for l in range(len(layers)):
            for spec in layers[l].primary_names:
                if spec in Dsed.keys():
                    layers[l].diffquo[spec]=Dsed[spec][l]
            #print('Dsedi %2.20f at layer l %2d'%(Dsedi['SO4--'][l],l))
        for num in range(len(layers)):
            if num < len(layers)-1:
                dz=(thick[num+1]+thick[num])/2
            if num > 0:
                dz2=(thick[num]+thick[num-1])/2
            if step > 0:
                layers[num].saturation=Porw_vol[num,step-1]
            #layers[num].porosity=por[num]
##transport rate computation starts here
            if num == 0:
                #layers[0].copy_to_alquimia(data)
                if bc is not None:
                    primarynames=layers[num].primary_names
                if dzw > 0:
                    if Zt[step-1] <= Zs: #should seperate infiltration and diffusion, the following is diffusion process, but the very first moment of flooding is infiltration.
                        for spec in primarynames:
                            if spec in D0.keys():
                                if spec in tconc_fc.keys():            ##the upper boundary of sediment column or SWI is imposed to be the tide_conc.
                                    Cw_tmp=tide_conc[spec][step]              ##concentration in overlaying water supposed is uniform within water column, otherwise this value should be bottomwater concentration
                                    Cswi=Cw_tmp #(layers[num].total_mobile[spec]-Cw_tmp)*sig/(sig+thick[0]/2)+Cw_tmp    ##concentration at sediment-water interface
                                    #Jflx=layers[num].porosity**2*Dsedi[spec][num]*(tide_conc[spec][step]-layers[num].total_mobile[spec])*1e3/(soillyr[0])   #(lyrdpth[0]+lyrdpth[1]) mol/m2/s ; Jflx positive means downward flux or flux into soil and negative means upward flux outside of soil
                                else:
                                    #Jflx=0
                                    Cswi=layers[num].total_mobile[spec]
                                #spec_mol_layer=layers[num].total_mobile[spec]*layers[num].volume*1000*layers[num].porosity*layers[num].saturation
                                #rtmp4=Jflx*1/(layers[num].volume*1000*layers[num].porosity*layers[num].saturation)
                                #rtmp3=0.0
                                #cje compute flux from underlie layer to top layer
                                rtmp1=pori[num]*Dsedi[spec][num]*(layers[num+1].total_mobile[spec]-layers[num].total_mobile[spec])/(layers[num].porosity*dz*thick[num])
                                rtmp2=wInf*porInf*(alpha[num+1]*layers[num].total_mobile[spec]+(1-alpha[num+1])*layers[num+1].total_mobile[spec])/(layers[num].porosity*thick[num])
                                rtmp3=por[num]*Dsed[spec][num]*(layers[num].total_mobile[spec]-Cswi)/(layers[num].porosity*(thick[num]/2+sig)*thick[num])
                                rtmp4=wInf*porInf*(0.5*Cswi+(1-0.5)*layers[num].total_mobile[spec])/(layers[num].porosity*thick[num])
                                Rtmp=(rtmp1-rtmp2-rtmp3+rtmp4)*0
                                ##cjw SWI flux -- first order fick's law
                                Jflx=(rtmp3*1000-rtmp4*1000)*layers[num].porosity*thick[num]/thick[num]                         ##unit is mol/m2/s
                                Jidx=list(D0.keys()).index(spec)
                                Jswi[Jidx,step]=Jflx
                            elif spec in immobile_species:
                                rtmp1=(1-pori[num])*Dbi[num]*(layers[num+1].total_immobile[spec]-layers[num].total_immobile[spec])/((1-layers[num].porosity)*dz*thick[num])
                                rtmp2=wInf*(1-porInf)*(alpha[num+1]*layers[num].total_immobile[spec]+(1-alpha[num+1])*layers[num+1].total_immobile[spec])/((1-layers[num].porosity)*thick[num])
                                rtmp3=Fimb[spec]/((1-layers[num].porosity)*thick[num])    #when flooded, immobile flux reduce to 0.1*Fimb
                                rtmp4=0.0
                                Rtmp=(rtmp3-rtmp4+rtmp1-rtmp2)*0
                            else:
                                Rtmp=0.0
                            if spec not in Rc.keys():
                                Rc[spec]=[Rtmp]     ##suppose the flux happens at one unit m2 surface
                            else:
                                Rc[spec].append(Rtmp)
                    elif Zt[step-1] > Zs:
                        for spec in primarynames:
                            if spec in D0.keys():
                                if spec in tconc_fc.keys():            ##the upper boundary of sediment column or SWI is imposed to be the tide_conc.
                                    Cw_tmp=tide_conc[spec][step]              ##concentration in overlaying water supposed is uniform within water column, otherwise this value should be bottomwater concentration
                                    Cswi=Cw_tmp #(layers[num].total_mobile[spec]-Cw_tmp)*sig/(sig+thick[0]/2)+Cw_tmp    ##concentration at sediment-water interface
                                    #Jflx=layers[num].porosity**2*Dsedi[spec][num]*(tide_conc[spec][step]-layers[num].total_mobile[spec])*1e3/(soillyr[0])   #(lyrdpth[0]+lyrdpth[1]) mol/m2/s ; Jflx positive means downward flux or flux into soil and negative means upward flux outside of soil
                                else:
                                    #Jflx=0
                                    Cswi=layers[num].total_mobile[spec]
                                #spec_mol_layer=layers[num].total_mobile[spec]*layers[num].volume*1000*layers[num].porosity*layers[num].saturation
                                #rtmp4=Jflx*1/(layers[num].volume*1000*layers[num].porosity*layers[num].saturation)
                                #rtmp3=0.0
                                #cje compute flux from underlie layer to top layer
                                rtmp1=pori[num]*Dsedi[spec][num]*(layers[num+1].total_mobile[spec]-layers[num].total_mobile[spec])/(layers[num].porosity*dz*thick[num])
                                rtmp2=wInf*porInf*(alpha[num+1]*layers[num].total_mobile[spec]+(1-alpha[num+1])*layers[num+1].total_mobile[spec])/(layers[num].porosity*thick[num])
                                rtmp3=por[num]*Dsed[spec][num]*(layers[num].total_mobile[spec]-Cswi)/(layers[num].porosity*(thick[num]/2+sig)*thick[num])
                                rtmp4=wInf*porInf*(0.5*Cswi+(1-0.5)*layers[num].total_mobile[spec])/(layers[num].porosity*thick[num])
                                Rtmp=rtmp1-rtmp2-rtmp3+rtmp4
                                ##cjw SWI flux -- first order fick's law
                                Jflx=(rtmp3*1000-rtmp4*1000)*layers[num].porosity*thick[num]/thick[num]                         ##unit is mol/m2/s
                                Jidx=list(D0.keys()).index(spec)
                                Jswi[Jidx,step]=Jflx
                            elif spec in immobile_species:
                                rtmp1=(1-pori[num])*Dbi[num]*(layers[num+1].total_immobile[spec]-layers[num].total_immobile[spec])/((1-layers[num].porosity)*dz*thick[num])
                                rtmp2=wInf*(1-porInf)*(alpha[num+1]*layers[num].total_immobile[spec]+(1-alpha[num+1])*layers[num+1].total_immobile[spec])/((1-layers[num].porosity)*thick[num])
                                rtmp3=Fimb[spec]/((1-layers[num].porosity)*thick[num])    #when flooded, immobile flux reduce to 0.1*Fimb
                                rtmp4=0.0
                                Rtmp=rtmp3-rtmp4+rtmp1-rtmp2
                            else:
                                Rtmp=0.0
                            if spec not in Rc.keys():
                                Rc[spec]=[Rtmp]     ##suppose the flux happens at one unit m2 surface
                            else:
                                Rc[spec].append(Rtmp)
                                #print('Jflx %2.30f'%Rtmp)
                elif dzw <=0 :
                    for spec in primarynames:
                        #pos=get_alquimiavector(data.meta_data.primary_names).index(spec)  ##this should work because data is copied to alquimia#is empty now
                        if spec in D0.keys():
                            if spec in gas_species.keys():
                                Cw_tmp=gas_species[spec]
                                if spec == 'O2(aq)':
                                    #Cswi=0.002
                                    #Cw_tmp=gas_species[spec]#0.00028
                                    Cswi=Cw_tmp#(layers[num].total_mobile[spec]-Cw_tmp)*sig/(sig+thick[0]/2)+Cw_tmp
                                else:
                                    if layers[num].total_mobile[spec] > Cw_tmp:
                                        Cswi=Cw_tmp#(layers[num].total_mobile[spec]-Cw_tmp)*sig/(sig+thick[0]/2)+Cw_tmp
                                    else:
                                        Cswi=0.0#layers[num].total_mobile[spec]  
                            else:
                                Cswi=layers[num].total_mobile[spec]                        
                            rtmp1=pori[num]*Dsedi[spec][num]*(layers[num+1].total_mobile[spec]-layers[num].total_mobile[spec])/(layers[num].porosity*dz*thick[num])
                            rtmp2=wInf*porInf*(alpha[num+1]*layers[num].total_mobile[spec]+(1-alpha[num+1])*layers[num+1].total_mobile[spec])/(layers[num].porosity*thick[num])
                            rtmp3=por[num]*Dsed[spec][num]*100*(layers[num].total_mobile[spec]-Cswi)/(layers[num].porosity*thick[num]*(sig+thick[num]/2))  ##surface o2 diffusion increase 2 order (100)
                            rtmp4=wInf*porInf*(0.5*Cswi+(1-0.5)*layers[num].total_mobile[spec])/(layers[num].porosity*thick[num])

                            Rtmp=rtmp1-rtmp2-rtmp3+rtmp4

                            Jflx=(rtmp3*1000-rtmp4*1000)*layers[num].porosity*thick[num]    ##positive means flux is upward, negative means flux is in soil
                            Jidx=list(D0.keys()).index(spec)
                            Jswi[Jidx,step]=Jflx
                        elif spec in immobile_species:
                            rtmp1=(1-pori[num])*Dbi[num]*(layers[num+1].total_immobile[spec]-layers[num].total_immobile[spec])/((1-layers[num].porosity)*dz*thick[num])
                            rtmp2=wInf*(1-porInf)*(alpha[num+1]*layers[num].total_immobile[spec]+(1-alpha[num+1])*layers[num+1].total_immobile[spec])/((1-layers[num].porosity)*thick[num])
                            rtmp3=Fimb[spec]/((1-layers[num].porosity)*thick[num])
                            rtmp4=0.0
                            Rtmp=rtmp3+rtmp4+rtmp1-rtmp2
                            #print('SOM %1.5f',Rtmp)
                        else:
                            Rtmp=0.0
                        if spec not in Rc.keys():
                            Rc[spec]=[Rtmp]     ##suppose the flux happens at one unit m2 surface
                        else:
                            Rc[spec].append(Rtmp)
            elif num < len(layers)-1 and num >0:
                #layers[num+1].porosity=por[num+1]
                primarynames=layers[num].primary_names
                for spec in primarynames:
                    if spec in D0.keys():
                        #print('Dsedi %2.20f at layer l %2d'%(Dsedi[spec][num],num))
                        rtmp1=pori[num]*Dsedi[spec][num]*(layers[num+1].total_mobile[spec]-layers[num].total_mobile[spec])/(layers[num].porosity*dz*thick[num])
                        rtmp2=wInf*porInf*(alpha[num+1]*layers[num].total_mobile[spec]+(1-alpha[num+1])*layers[num+1].total_mobile[spec])/(layers[num].porosity*thick[num])
                        rtmp3=pori[num-1]*Dsedi[spec][num-1]*(layers[num].total_mobile[spec]-layers[num-1].total_mobile[spec])/(layers[num].porosity*thick[num]*dz2)
                        rtmp4=wInf*porInf*(alpha[num]*layers[num-1].total_mobile[spec]+(1-alpha[num])*layers[num-1].total_mobile[spec])/(layers[num].porosity*thick[num])
                        Rtmp=rtmp1-rtmp2-rtmp3+rtmp4
                        #Rtmp=pori[num+1]*Dsedi[spec][num+1]*(layers[num+1].total_mobile[spec]-layers[num].total_mobile[spec])/(layers[num].porosity*dz*thick[num]) - wInf*porInf*(alpha[num+1]*layers[num].total_mobile[spec]+(1-alpha[num+1])*layers[num+1].total_mobile[spec])/(layers[num].porosity*thick[num]) - pori[num]*Dsedi[spec][num]*(layers[num].total_mobile[spec]-layers[num-1].total_mobile[spec])/(layers[num].porosity*thick[num]*dz2) + wInf*porInf*(alpha[num]*layers[num-1].total_mobile[spec]+(1-alpha[num])*layers[num-1].total_mobile[spec])/(layers[num].porosity*thick[num])
                    elif spec in immobile_species:
                        rtmp1=(1-pori[num])*Dbi[num]*(layers[num+1].total_immobile[spec]-layers[num].total_immobile[spec])/((1-layers[num].porosity)*dz*thick[num])
                        rtmp2=wInf*(1-porInf)*(alpha[num+1]*layers[num].total_immobile[spec]+(1-alpha[num+1])*layers[num+1].total_immobile[spec])/((1-layers[num].porosity)*thick[num])
                        rtmp3=(1-pori[num-1])*Dbi[num-1]*(layers[num].total_immobile[spec]-layers[num-1].total_immobile[spec])/((1-layers[num].porosity)*thick[num]*dz2)
                        rtmp4=wInf*(1-porInf)*(alpha[num]*layers[num-1].total_immobile[spec]+(1-alpha[num])*layers[num-1].total_immobile[spec])/((1-layers[num].porosity)*thick[num])
                        Rtmp=rtmp1-rtmp2-rtmp3+rtmp4
                    else:
                        Rtmp=0.0
                    if spec not in Rc.keys():
                        Rc[spec]=[Rtmp]     ##suppose the flux happens at one unit m2 surface
                    else:
                        Rc[spec].append(Rtmp)
            elif num==len(layers)-1:
                primarynames=layers[num].primary_names
                for spec in primarynames:
                    Cbtm=layers[num].total_mobile[spec]      ##bottom boundary layer for mobile concentration
                    Sbtm=layers[num].total_immobile[spec]    ## bottom boundary for immobile concentration
                    if spec in D0.keys():
                        #Rtmp=0.0
                        rtmp1=0.0
                        rtmp2=wInf*porInf*(alpha[num]*layers[num].total_mobile[spec]+(1-alpha[num])*Cbtm)/(layers[num].porosity*thick[num])
                        rtmp3=pori[num-1]*Dsedi[spec][num-1]*(layers[num].total_mobile[spec]-layers[num-1].total_mobile[spec])/(layers[num].porosity*thick[num]*dz2)
                        rtmp4=wInf*porInf*(alpha[num]*layers[num-1].total_mobile[spec]+(1-alpha[num])*layers[num-1].total_mobile[spec])/(layers[num].porosity*thick[num])
                        Rtmp=rtmp1-rtmp2-rtmp3+rtmp4
                    elif spec in immobile_species:
                        #Rtmp=0.0
                        rtmp1=0.0
                        rtmp2=wInf*(1-porInf)*(alpha[num]*layers[num].total_immobile[spec]+(1-alpha[num])*Cbtm)/((1-layers[num].porosity)*thick[num])
                        rtmp3=(1-pori[num-1])*Dbi[num-1]*(layers[num].total_immobile[spec]-layers[num-1].total_immobile[spec])/((1-layers[num].porosity)*thick[num]*dz2)
                        rtmp4=wInf*(1-porInf)*(alpha[num]*layers[num-1].total_immobile[spec]+(1-alpha[num])*layers[num-1].total_immobile[spec])/((1-layers[num].porosity)*thick[num])
                        Rtmp=rtmp1-rtmp2-rtmp3+rtmp4
                    else:
                        Rtmp=0.0
                    if spec not in Rc.keys():
                        Rc[spec]=[Rtmp]     ##suppose the flux happens at one unit m2 surface
                    else:
                        Rc[spec].append(Rtmp)
        #dCO2=layers[0].total_mobile['HCO3-']-initial_HCO3
        #layers[0].total_mobile['HCO3-']=layers[0].total_mobile['HCO3-']-dCO2
        #layers[0].total_immobile['H+']=layers[0].total_immobile['H+']-dCO2*1000*layers[0].porosity*layers[0].saturation

#cjw below is the evaptranspiration processes, rain , infiltration, drainge and excretion processes
        #if dzw<=0: #Zt[step] <= Zs:               #water level should be the current timestep, unlike the Vp and Sp which are last timestep
        Sp_o={}
        Se_o={}
        Sp_n={}
        Vp=[]
        flow_in=numpy.zeros(len(layers)+1,dtype=float)
        flow_out=numpy.zeros(len(layers)+1,dtype=float)
        Eah[step]=Eah[step]*2#*0.1#Vsat=[]
        if dzw <= 0:
            #Rain_tmp[0,step]=Rg[step]
            cp0=0
            cp1=0
            for num in range(len(layers)):
                cp0+=layers[num].volume*layers[num].porosity
                cp1+=layers[num].volume*layers[num].porosity*layers[num].saturation
            Wlylim=cp0-cp1+Eah[step]    #water limitation for entire soil core in m3
            for num in range(len(layers)):
                #sd is the specific drainage rate (h-1). sd=K/(n*L), K is the hydraulic conductivities in cm/h,
#   (1.45cm/h, 0.21cm/h), n is porosities (0.43,67.7) which is dimensionless,L is the depth or length of the cores (30cm)
                sd=0.2*0.01/(layers[num].porosity*zdpth[num])                            #draign rate h-1 0.21
                V0=layers[num].volume*layers[num].porosity                                # water volume in layer unit is m3
                #print('saturation %1.5f',layers[num].saturation)
                Vtmp_o=layers[num].volume*layers[num].porosity*layers[num].saturation     #initialize pore water volume in m3
                #print('saturation %1.5f,%1.5f',layers[num].saturation,Vtmp_o)
                if Vtmp_o > V0:
                    if Vtmp_o-V0 > 1e-8:
                        print('error1: water content exceed the limitation V0 %1.5f with Vtmp_o %1.5f at step %6d at layer %d',V0,Vtmp_o,step,num)
                    else:
                        Vtmp_o=V0
                elif Vtmp_o < V0:
                    dV=V0-Vtmp_o
                else:
                    #Vtmp_o=V0
                    dV=0.0
                Vlim=dV+Eah[step]*thick[num]/zdpth[-1]  #water capacity in the current layer: the inflow of water shall not pass this limit

                if num ==0:
                    Vtmp1=0
                    V_in=0.0
                    #print('flow in init %1.5f',flow_in[num])
                    if Rg[step] >= Wlylim:
                        flow_in[num]=Wlylim
                        #print('water limit %1.5f,%1.5f,%1.5f',Wlylim,Rg[step],Vlim)
                    else:
                        flow_in[num]=Rg[step]
                    #print('water limit %1.5f,%1.5f',flow_in[num],Vlim)
                    if Vtmp_o >= fc*V0:
                        V_out=Vtmp_o*sd
                    elif Vtmp_o < fc*V0:
                        V_out=0.0
                    Rain_tmp[num,step]=flow_in[num]+V_in-V_out         #water buget in soil
                    if Rain_tmp[num,step] > Vlim :
                        #V_out=V_out+Rain_tmp[num,step]-(V0-Vtmp_o)
                        flow_in[num+1]=Rain_tmp[num,step]-Vlim
                        Rain_tmp[num,step]=Rain_tmp[num,step]-flow_in[num+1]
                        #print('water limit num=1 %1.5f,%1.5f',flow_in[num+1],Vlim)
                elif num >0 :
                    Vtmp1=layers[num-1].volume*layers[num-1].porosity*layers[num-1].saturation
                    Vup0=layers[num-1].volume*layers[num-1].porosity
                    sd1=0.2*0.01/(layers[num-1].porosity*zdpth[num-1])   #0.21
                    #print('conectivity %1.5f,%1.5f',sd1,sd)
                    if Vtmp1 >= fc*Vup0:
                        V_in=Vtmp1*sd1
                    elif Vtmp1 < fc*Vup0:
                        V_in=0.0
                    if Vtmp_o > fc*V0:
                        V_out=Vtmp_o*sd
                    elif Vtmp_o <= fc*V0:
                        V_out=0.0
                    Rain_tmp[num,step]=flow_in[num]+V_in-V_out                    #water buget in soil
                    #print('water limit num>0 %1.5f,%1.5f',flow_in[num],Vlim)
                    if Rain_tmp[num,step] > Vlim:
                        #V_out=V_out+Rain_tmp[num,step]-(V0-Vtmp_o)
                        flow_in[num+1]=Rain_tmp[num,step]-Vlim
                        Rain_tmp[num,step]=Rain_tmp[num,step]-flow_in[num+1]
                if Vtmp_o-V0 > 1e-9:
                    print('error %1.5f',Wlylim)#error1: water content exceed the limitation V0 %1.5f with Vtmp_o %1.5f at step %6d at layer %d',V0,Vtmp_o,step,num)
                elif Vtmp_o <= V0:
                    #print('unsaturated condition')
                    Vtmp_n=Vtmp_o+(Rain_tmp[num,step]-Eah[step]*thick[num]/zdpth[len(layers)-1])  ##per second    ## this is not right because draingae from upper layer should infil lower layer
                #else:
                #    Vtmp_n=V0
## compute mass change of constitutes
                primarynames=layers[num].primary_names
                for spec in primarynames:
                    if num==0:
                        Sp_in=0.0
                    else:
                        Sp_in=layers[num-1].total_mobile[spec]*V_in*1000
                    Sp_out=layers[num].total_mobile[spec]*V_out*1000             ##drainage out
                    #print('after drainage %1.5f'%V_in)
                    Sp_tmp=layers[num].total_mobile[spec]*Vtmp_o*1000            ## total molar per unit square m2
            #    #if spec not in Sp_o.keys():
            #    #    Sp_o[spec]=Sp_tmp    #initialize pore water mol mass in each layer mol per unit area.
            #    #else:
            #    #    Sp_o[spec].append(Sp_tmp)
#cjw compute excretion rate
                    if spec in excret_spec:
                        Se_o[spec]=k2*k3*Sp_tmp*(Eah[step]*thick[num]/zdpth[len(layers)-1])/(Vtmp_o)            ## excretion mol/h
                        #print('spec %s Sp_n %2.50f at layer %2d'%(spec,Se_o[spec],num))
                    else:
                        Se_o[spec]=0.0
                    Sde_tmp=Sp_in-Sp_out-Se_o[spec]


                    if spec in gas_species.keys():
                        Sde_tmp=0.0/3600
                    else:
                        Sde_tmp=Sde_tmp/3600
                    if spec not in Sp_n.keys():
                        Sp_n[spec]=[Sde_tmp]                    ## rate of change/reduce of total mass by evaptranspiration + drainage in mass
                    else:
                        Sp_n[spec].append(Sde_tmp)
                Porw_vol[num,step]=Vtmp_n/V0
                #print('dry water volume change %2.6f',Porw_vol[num,step])
                Vchng=(Vtmp_n-Vtmp_o)/3600
                Vp+=[Vchng,]                           ## porewater volume change rate m3/s
        elif dzw > 0:
            for num in range(len(layers)):
                V0=layers[num].volume*layers[num].porosity                                # water volume in layer unit is m3
                Vtmp_o=layers[num].volume*layers[num].porosity*layers[num].saturation     #initialize pore water volume in m3

                primarynames=layers[num].primary_names
                for spec in primarynames:
                    if Zt[step-1] <= Zs: # the very first moment flood happens, and infiltration happens
                #print('first high tide true',j)
                        if Vtmp_o < V0:
                            Vtmp_n=V0
                            #for spec in primarynames:
                            if spec in tide_conc.keys():
                        #for spec in primarynames:
                                Sde_tmp= ((V0-Vtmp_o)*1000*tide_conc[spec][step])   ## infiltration process at the very first moment of inundation; rate of change/reduce of total mass by evaptranspiration
                            else:
                                Sde_tmp=0.0
                            #if spec not in Sp_n.keys():
                            #    Sp_n[spec]=[Sde_tmp]                    ## rate of change/reduce of total mass by evaptranspiration + drainage in mass
                            #else:
                            #    Sp_n[spec].append(Sde_tmp)
                        else:
                            #Vtmp_o = V0
                            Vtmp_n=V0
                            Sde_tmp=0.0
                    else:     # the remainder of the flood
                #print('remainder high tide true',j)
                        #if Vtmp_o == V0:
                        Vtmp_n=V0
                    #for spec in primarynames:
                        Sde_tmp=0.0
                    if spec in gas_species.keys():
                        Sde_tmp=0.0/3600
                    #elif spec in excret_spec:
                    #    Sde_tmp=Sde_tmp/3600
                    else:
                        Sde_tmp=Sde_tmp/3600
                    if spec not in Sp_n.keys():
                        Sp_n[spec]=[Sde_tmp]                    ## rate of change/reduce of total mass by evaptranspiration + drainage in mass
                    else:
                        Sp_n[spec].append(Sde_tmp)

                Porw_vol[num,step]=Vtmp_n/V0
                Vchng=(Vtmp_n-Vtmp_o)/3600
                #print('water volume change %2.6f',Porw_vol[num,step])
                Vp+=[Vchng,]                           ## porewater volume change rate m3/s

##cjw above is the first try on migrating boundary condition setup into time step loop.

        try:
            for n,l in enumerate(layers):
                l.copy_to_alquimia(data)
                dq={}
                if bc is not None:
                    primarynames=get_alquimiavector(data.meta_data.primary_names)
                    for spec in primarynames:
                        if spec in l.diffquo.keys():
                            if numpy.iterable(l.diffquo[spec]):
                                dq[spec]=l.diffquo[spec][step%len(l.diffquo[spec])]
                            else:
                                dq[spec]=l.diffquo[spec]
                        else:
                            dq[spec]=0.0
            #cjw get Transport Rate or flux for run_onestep
                #print('before porewater saturation %2.30f at layer %d'%(layers[n].saturation ,n))
                PwcR=Vp[n]#*0.001#*0.0    ##get porewater volume change rate m3/s for the current layer
                Rc_lyr={}
                Sde_lyr={}  ## change caused by excretion and drainge and infiltration and evapration
                for spec in layers[0].primary_names:
                    if spec in Rc.keys():
                        #Rc_lyr[spec]=Rc[spec][(step-1)*len(layers)+n]
                        Rc_lyr[spec]=Rc[spec][n]
                        #print('Rc_lyr %2.30f for %s at layer %d'%(Rc_lyr[spec],spec,n))
                        Sde_lyr[spec]=Sp_n[spec][n]
                    else:
                        Rc_lyr[spec]==0.0
                        Sde_lyr[spec]=0.0
                    SOMtranR[n,step]=Rc['SOM'][n]
                num_cuts=l.run_onestep(chem,data,dt,status,min_dt=min_dt,diffquo=dq,bc=bc_state,truncate_concentration=truncate_concentration,rateconstants=l.rateconstants,TRc=Rc_lyr,EvR=Sde_lyr,PwR=PwcR)
                l.copy_from_alquimia(data)
                Porw_vol2[n,step]=l.saturation#=Porw_vol[n,step]       #cjw update saturation level

                #print('after porewater saturation %2.30f at layer %d'%(layers[n].saturation ,n))
                # Write output
                l.write_output(step+1,dt,num_cuts)

        except RuntimeError as err:
            print('ERROR on timestep %d, layer %d: %s'%(step,n,err))
            print('Returning output so far')
            l.write_output(step,dt,num_cuts)
            success=False
            break
        except KeyboardInterrupt as err:
            print('INTERRUPTED on timestep %d'%(step))
            print('Returning output so far')
            l.write_output(step,dt,num_cuts)
            success=False
            break




        if step%100==0 and step>0:
            t1=time.time()
            cuts=[l.output['ncuts'][step-100:step].mean() for l in layers]
            mean_dt=numpy.mean([l.output['actual_dt'][step-100:step].mean() for l in layers])
            print('*** Step {step:d} of {nsteps:d} ({nyears:1.1f} of {totalyears:d} years). Time elapsed: {t:1.1f} min ({tperstep:1.1f} s per {steplength:1.1f} hour timestep). Mean cuts: {meancuts:s} Mean dt: {meandt:1.1f} s ***'.format(
                    step=step,nsteps=nsteps,t=(t1-t0)/60,tperstep=(t1-tprev)/25,meancuts=str(cuts),meandt=mean_dt,steplength=dt/3600,nyears=step*dt/(3600*24*365),totalyears=int(nsteps*dt/(3600*24*365))))
            tprev=t1


#cjw not sure what is in output
    output=convert_to_xarray(layers,t0=tstart)
    if isinstance(restart_state,xarray.Dataset):
        output=xarray.concat([restart_state,output.isel(time=slice(1,None))],dim='time')


    import datetime
    today=datetime.datetime.today()
    output.to_netcdf(path='./Evap_output.nc')
    #LtranR.to_xarray().to_netcdf(path='./TransportRate.nc')
    #Jswi.to_xarray().to_netcdf(path='./Surfaceflux.nc')
    #Jswi.to_csv(path='./Jswi.csv')
    #LtranR.to_csv(path='./LtranR.csv')
    #numpy.save("./Evap_LtranR.npy",LtranR)
    numpy.save("./Evap_Porw.npy",Porw_vol2)
    numpy.save("./Evap_Jswi.npy",Jswi)
#cjw        output.to_netcdf('Mn_output/Mn_pH{ph:1.1f}_Ndep{Ndep:03d}_warming{warming:d}_{year:04d}-{month:02d}-{day:02d}.nc'.format(ph=pH,Ndep=int(Ndep/(1000/molar_mass['N']/100**2/(365*24*3600) )),year=today.year,month=today.month,day=today.day,warming=warming))
#cjw        convert_to_xarray([incubation_layer],leaf_Mn=leaf_Mn_concs).to_netcdf('Mn_output/Mn_incubations_pH{ph:1.1f}_Ndep{Ndep:03d}_warming{warming:d}_{year:04d}-{month:02d}-{day:02d}.nc'.format(ph=pH,Ndep=int(Ndep/(1000/molar_mass['N']/100**2/(365*24*3600) )),year=today.year,month=today.month,day=today.day,warming=warming))


print('\n\n\n Simulation finished. Total time: %1.1f minutes\n'%((time.time()-starting_time)/60))


def plot_output(output,axs,subsample=1,do_legend=True,**kwargs):
    for num in range(len(output.depth)):
        out=output.isel(depth=num,time=slice(None,None,subsample)).dropna(dim='time')
        t=out.time*24/(1)
        #print(len(t),len(out['Total DOM1']))
        #std=(31+28+31+30+31+30+31)*24
        #etd=(31+28+31+30+31+30+31+30)*24
        #std=(31+28+31+30+31+30+31)*24
        #etd=(31+28+31+30+31+30+31+14)*24
        #std=31*24
        #etd=45*24
        #etd=(31+28+31+30)*24
        numyr=1*repyr

        etd=365*24*numyr
        subt=t#[std:etd]
        egnr=24*0+365*24*10
        porosity=out['Porosity']
        saturation=out['saturation']
        BD=out['BD']
        axs[num,1].plot(subt[egnr:etd],out['Total DOM1'][egnr:etd]*1e6,c='C2',label='DOM',**kwargs)
        axs[num,0].plot(subt[egnr:etd],out['Total O2(aq)'][egnr:etd]*1e6,label='O$_2$',c='C0',**kwargs)
        #axs[num,1].plot(subt[egnr:etd],out['Total DOM1'][egnr:etd]*1e6,label='NO$_3$$^-$',c='C5',**kwargs)
        axs[num,2].plot(subt[egnr:etd],out['Total Cl-'][egnr:etd]*1e6,label='Cl$^-$',c='C1',**kwargs)
        #axs[num,3].plot(subt,out['Total Fe+++'][std:etd]/(out['Total Fe+++'][std:etd]+out['Total Fe++'][std:etd]),label='Fe$^3$$^+$',c='C5',**kwargs)
        #axs[num,3].plot(t,out['Total Na+']*1e3,label='Na$^+$',c='C5',**kwargs)
        #axs[num,4].plot(subt[egnr:etd], SOMtranR[num,egnr:etd])
        axs[num,4].plot(subt[egnr:etd],out['Total Sorbed SOM'][egnr:etd]*1e6,label='SOM',c='C1',**kwargs)
        axs[num,3].plot(subt[egnr:etd],out['Total SO4--'][egnr:etd]*1e6,c='C3',label='SO$_4$$^2$$^-$',**kwargs)
        #axs[num,4].plot(saturation,c='C2',label='saturation',**kwargs)
        #axs[num,4].plot(subt[egnr:etd],out['Total SO4--']/out['Total Cl-'][egnr:etd],c='C2',label='SO$_4$$^2$$^-$',**kwargs)
        #axs[num,4].plot(subt[egnr:etd],(numpy.zeros(len(subt))+0.14)[egnr:etd],c='black',label='SO$_4$$^2$$^-$',**kwargs)
        axs[num,5].plot(subt[egnr:etd],out['Total CH4(aq)'][egnr:etd]*1e6,c='C4',label='CH$_4$',**kwargs)
        #axs[num,0].set_ylabel('mM')
        #axs[num,1].set_ylabel('mM')
        axs[num,0].set_ylabel('\u03BCM')

        #axs[num,0].set_yscale('log')
        #axs[num,1].set_yscale('log')
        #axs[num,2].set_yscale('log')
        #axs[num,3].set_yscale('log')
        #axs[num,4].set_yscale('log')
        #axs[num,5].set_yscale('log')
        #axs[num,3].set_ylabel('mM')
        #axs[num,4].set_ylabel('mM')
        ###cjw methane saturation conc in water is 1746 umol/L
        if num == 0:
            #axs[num,0].set_ylim([0,100])
            axs[num,0].set_ylim([0,300])
            #axs[num,2].set_ylim([200,1500])
            #axs[num,3].set_ylim([0,90])
            #axs[num,4].set_ylim([0,0.2])
            #axs[num,5].set_ylim([0,0.2])
        elif num==1:
            axs[num,0].set_ylim([0,2])
            #axs[num,1].set_ylim([0,0.1])
            #axs[num,2].set_ylim([0,120])
            #axs[num,3].set_ylim([0,5])
            #axs[num,4].set_ylim([0,0.2])
            #axs[num,5].set_ylim([0,2])
        elif num==2:
            axs[num,0].set_ylim([0,0.01])
            #axs[num,1].set_ylim([0,1e-5])
            #axs[num,2].set_ylim([0,5])
            #axs[num,3].set_ylim([0,1])
            #axs[num,4].set_ylim([0,0.2])
            #axs[num,5].set_ylim([0,800])
        elif num>=3:
            axs[num,0].set_ylim([0,0.005])
            #axs[num,1].set_ylim([0,1e-5])
            #axs[num,2].set_ylim([0,0.2])
            #axs[num,3].set_ylim([0,1e-5])
            #axs[num,4].set_ylim([0,0.2])
            #axs[num,5].set_ylim([0,800])

        #    axs[num,0].set_ylim([350,700])
        #elif num >= 2 :
        #    axs[num,0].set_ylim([350,1e4])
        #if num >= 3 :
        #    axs[num,0].set_ylim([0,2000])
        #    axs[num,1].set_ylim([0,0.5])
        #else:
        #    axs[num,0].set_ylim([0,25])

        #elif num==1:
        #    axs[num,0].set_ylim([0,0.1])
            #axs[num,1].set_ylim([0,0.0015])
        #else:
        #    axs[num,0].set_ylim([0,2])
        #    axs[num,1].set_ylim([0,0.0015])
        #axs[num,0].set_ylim([0,2000])


        #if num < 4:
        #    axs[num,4].set_ylim([0,100])
        #elif num >=4 :
        #    axs[num,4].set_ylim([0,1200])

#        axs[num,0].set_ylabel('C density\n(g C cm$^{-3}$)')

#        axs[num,1].plot(t,out['Total Mn++']*1e6/1000*porosity*saturation*Mn_molarmass/BD,label='Mn$^{+\!\!+}$',c='C0',**kwargs)
#        axs[num,1].plot(t,out['Total Mn+++']*1e6/1000*porosity*saturation*Mn_molarmass/BD,label='Mn$^{+\!\!+\!\!+}$',c='C1',**kwargs)
        # axs[num,1].plot(t,layers[num].output_DF['Manganite VF']/molar_volume_manganite*1e6*Mn_molarmass/BD,label='Manganite')
#        axs[num,2].plot(t,out['Birnessite2 VF']*7/molar_volume_birnessite*1e6*Mn_molarmass/BD,label='Birnessite',c='C0',**kwargs)

#        annual_root_uptake=out['Total Tracer2'].groupby(numpy.floor(t)).max()
#        axs[num,1].plot(annual_root_uptake.time,annual_root_uptake*1e6/1000*porosity.mean()*saturation*Mn_molarmass/BD,label='Annual Mn$^{+\!\!+}$ root uptake',c='C2',**kwargs)

#        axs[num,3].plot(t,-numpy.log10(out['Free H+']),label='pH',c='C0',**kwargs)

        # ax=axs[num,1].twinx()
        # axs[num,2].plot(t,layers[num].output_DF['Mn(OH)2(am) VF']/molar_volume_MnOH2am*1e6*Mn_molarmass/BD,label='Mn(OH)$_2$(am)',ls='-')
#        axs[num,2].plot(t,(out['Total Mn++']+out['Total Mn+++'])*1e6/1000*porosity*saturation*Mn_molarmass/BD +
#                            (out['Birnessite2 VF']*7/molar_volume_birnessite)*1e6*Mn_molarmass/BD,c='k',label='Total Mn',**kwargs)

#        axs[num,1].set_ylabel('Mn concentration\n($\mu$g g$^{-1}$)')
#        axs[num,2].set_ylabel('Mn concentration\n($\mu$g g$^{-1}$)')
#        axs[num,3].set_ylabel('pH')
        # axs[num,1].set_ylim(*axs[0,1].get_ylim())

#        for n,cation in enumerate(['Al+++','Ca++','K+','Mg++','Na+','Mn++','Mn+++']):
#            axs[num,4].plot(t,output['Total Sorbed '+cation].isel(depth=num,time=slice(None,None,subsample)).dropna(dim='time')/(output.BD.isel(depth=2)*1e-3*100**3)*1000,c='C'+str(n),label=cation,**kwargs)
#        axs[num,4].plot(t,output['CEC H+'].isel(depth=num,time=slice(None,None,subsample)).dropna(dim='time')/(output.BD.isel(depth=2)*1e-3*100**3)*1000,label='H+',c='C'+str(n+1),**kwargs)

#        axs[num,4].set_ylabel('Exch conc\n(mmol/kg)')



    axs[-1,0].set_xlabel('Time (Month)')
    if do_legend:
        axs[0,0].legend()
        axs[0,1].legend()
        axs[0,2].legend()
        axs[2,3].legend()
        axs[2,4].legend()
        #axs[2,4].legend(['porewater','tidal water'])
        axs[0,5].legend()

    axs[-1,1].set_xlabel('Time (Day)')
    axs[-1,2].set_xlabel('Time (Day)')
    axs[-1,3].set_xlabel('Time (Day)')
    axs[-1,4].set_xlabel('Time (Day)')
    axs[-1,5].set_xlabel('Time (Day)')
    axs[0,1].set_title('Total DOM1')
    axs[0,0].set_title('Total Oxygen')
    axs[0,2].set_title('Total Cl$^-$')
    #axs[0,1].set_title('Total NO$_3$$^-$')
    #axs[0,3].set_title('Total Fe$^3$$^+$')
    axs[0,3].set_title('Total SO$_4$$^2$$^-$')
    #axs[0,4].set_title('Ratio SO$_4$$^2$$^-$/Cl$^-$')
    #axs[0,4].set_title('saturation')
    axs[0,5].set_title('CH$_4$')

hydrofig,(ax1,ax2,ax3)=pyplot.subplots(3, 1,figsize=(8,6))
ax1.plot(Zt,linewidth=2)
ax2.plot(tide_conc['Cl-'],linewidth=2)
#ax2.plot(Ini_tidec['Cl-'],linewidth=0.1)
ax3.plot(tide_conc['SO4--'],linewidth=2)
#ax3.plot(Ini_tidec['SO4--'],linewidth=0.1)
pyplot.savefig('hydrofig_evap.pdf')

sat,axs=pyplot.subplots(1, 1,figsize=(6,8))
axs.plot(Porw_vol2[0,:])
axs.plot(Porw_vol2[1,:])
axs.plot(Porw_vol2[2,:])
axs.plot(Porw_vol2[3,:])
axs.plot(Porw_vol2[4,:])
axs.plot(Porw_vol2[5,:])
#difffig,axd=pyplot.subplots(1, 1,figsize=(6,8))
#axd.plot(Dsed['CH4(aq)'],mid_dpth,linewidth=2)
#pyplot.gca().invert_yaxis()
pyplot.savefig('porewater_evap.pdf')

#Ratefig,axd=pyplot.subplots(1, 1,figsize=(6,8))
#pltmp=[i*1e6 for i in Rc['CH4(aq)']]
#axd.plot(pltmp,mid_dpth,linewidth=2)
#pyplot.gca().invert_yaxis()
#pyplot.savefig('Ratefig.pdf')

#Diffbio,axbio=pyplot.subplots(1, 1,figsize=(8,8))
#Diffb=[5.2*10**(0.76241122-0.00039724*i)*(0.0001/(365*24*3600)) for i in Zt]
#axbio.plot(Diffb,linewidth=2)
#pyplot.savefig('Diffbio.pdf')

#porosityfig,(axp1,axp2)=pyplot.subplots(1, 2,figsize=(6,8))
#axp1.plot(por,mid_dpth,linewidth=2)
#axp2.plot(pori,zdpth,linewidth=2)
#axp1.invert_yaxis()
#axp2.invert_yaxis()
#pyplot.savefig('porosityfig.pdf')

f,axs=pyplot.subplots(ncols=6,nrows=len(output.depth),sharex=True,clear=True,num='Simulation results',figsize=(14,9))
plot_output(output,axs)
pyplot.setp(axs, xticks=[(59+365*10)*24,(120+365*10)*24,(181+365*10)*24,(243+365*10)*24,(304+365*10)*24,24*365*11], xticklabels=['02','04','06','08','10','12'])
#pyplot.setp(axs, xticks=[31,34,37,40,44], xticklabels=['31','34','37','40','44'])
#pyplot.setp(axs, xticks=[1,14,31,45,59,73,90,104,120], xticklabels=['1',' ','31',' ','59',' ','90',' ','120'])
#pyplot.setp(axs, xticks=[212,219,226,233,240], xticklabels=['212',' ','226',' ','240'])
#pyplot.setp(axs, xticks=[212,215,218,221,226], xticklabels=['212','215','218','221','226'])
pyplot.tight_layout()

pyplot.savefig('f_evap.pdf')
pyplot.show()
#cjw networkfig=pyplot.figure('Reaction network',clear=True)
#cjw drawn=decomp_network.draw_network_with_reactions(reaction_network,omit=['NH4+','Rock(s)','gas','secondary','H+','>Carboxylate-','Carboxylic_acid'],
#c        font_size='medium',node_size=1500,font_color='k',arrowstyle='->',arrowsize=10.0,edge_color='gray',node_alpha=1.0,
#c        namechanges={'cellulose':'Cellulose','DOM1':'DOM','O2(aq)':'O$_2$(aq)','CH4(aq)':'CH$_4$(aq)','HCO3-':'HCO$_3^-$','DOM2':'Exposed lignin','sorbed_DOM1':'Sorbed DOM',
#c                     'Fe(OH)2':'Fe(OH)$_2$','Fe(OH)3':'Fe(OH)$_3$','Mn++':r'Mn$^\mathrm{+\!\!+}$','Mn+++':r'Mn$^\mathrm{+\!\!+\!\!\!+}$','Acetate-':'Acetate',})
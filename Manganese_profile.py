from run_alquimia import get_alquimiavector,ffi,lib,check_status,init_alquimia,convert_condition_to_alquimia,convert_rateconstants
import decomp_network
import Manganese_network as Mn
import numpy
import xarray
import time

class layer:
    def __init__(self,volume,saturation=1.0,temperature=20.0,water_density=1000.0,porosity=0.25,pressure=101325.0,BD=1.5,CEC=50.0,rateconstants={},diffquo={}):
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
        self.CEC=CEC # meq/kg
        
        # Some things that I'm not using are skipped

    def __str__(self):
        out='Layer with volume = %1.2e m, porosity = %1.2f, CEC = %1.2f'%(self.volume,self.porosity,self.CEC)
        out=out + '\n                    Mobile    Immobile\n'
        for spec in self.total_mobile:
            out=out + spec.ljust(20)+'%1.2e   %1.2e\n'%(self.total_mobile[spec],self.total_immobile[spec])
        out=out+'Rateconstants:\n'
        for key in sorted(self.rateconstants):
            out=out+key.ljust(30)+': %1.2e\n'%self.rateconstants[key]
        
        return out
        
    def copy_from_alquimia(self,data):
        self.volume=data.properties.volume
        self.saturation=data.properties.saturation
        self.temperature=data.state.temperature
        self.water_density=data.state.water_density
        self.porosity=data.state.porosity
        self.pressure=data.state.aqueous_pressure
        
        self.mineral_names=get_alquimiavector(data.meta_data.mineral_names)
        for num in range(data.meta_data.mineral_names.size):
            self.mineral_rate_cnst[self.mineral_names[num]]=data.properties.mineral_rate_cnst.data[num] # Alquimia doesn't actually update this
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
            'CEC H+':numpy.ma.masked_all(nsteps,dtype=float),
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
        self.output['actual_dt'][step]=dt/2**num_cuts
        self.output['ncuts'][step]=num_cuts
        # Keep track of how much H+ is in the carboxylate buffer rather than the CEC site, so we can plot CEC-exchangeable H+ later. 
        # Assumes there is one ion exchange site, and that the buffering sorption site is on Rock(s)
        # Aux doubles in this spot stores free surface site density (probably best not to rely on aux_doubles in general though)
        self.output['CEC H+'][step]=self.total_immobile['H+']-(self.surface_site_density['>Carboxylate-']*self.mineral_volume_fraction['Rock(s)']-self.aux_doubles[len(self.total_mobile)*2+len(self.secondary_free_ion_concentration)+1+self.surface_site_names.index('>Carboxylate-')])
        self.output['flow_in'][step,:]=numpy.array([self.flow_in.get(name,0.0) for name in self.primary_names])
        self.output['flow_out'][step,:]=numpy.array([self.flow_out.get(name,0.0) for name in self.primary_names])
        
    def convert_output(self,downsample=10):
        import pandas
        output_DF=pandas.DataFrame(index=self.output['time'])
        output_DF['time']=self.output['time']
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
        
        output_DF['CEC H+']=pandas.DataFrame(self.output['CEC H+'],index=self.output['time'])
        output_units['CEC H+']='mol/m^3'

        ind=numpy.arange(len(output_DF))
        output_DF=output_DF.groupby(numpy.floor(ind/downsample)*downsample).mean()
        
        self.output_DF=output_DF.reset_index(drop=True).set_index(output_DF['time']/(24*3600))
        self.output_DF['time']=self.output_DF['time']/(24*3600)
        self.output_units=output_units
        
        
    def run_onestep(self,chem,data,dt,status,min_dt=0.1,num_cuts=0,diffquo={},bc=None,flux_tol=0.15,truncate_concentration=0.0,rateconstants={},min_cuts=0):
        self.copy_to_alquimia(data)
        converged=False
        porosity=data.state.porosity
        
        max_cuts=num_cuts
        actual_dt=dt/2**num_cuts
        
        for num,reactname in enumerate(get_alquimiavector(data.meta_data.aqueous_kinetic_names)):
            data.properties.aqueous_kinetic_rate_cnst.data[num]=rateconstants[reactname]
        
        for spec in diffquo.keys():
            pos=get_alquimiavector(data.meta_data.primary_names).index(spec)
            if data.state.total_mobile.data[pos] < truncate_concentration and data.state.total_mobile.data[pos] != 0.0:
                # print('Truncating concentration of {spec:s} from {conc:1.1g}'.format(spec=spec,conc=data.state.total_mobile.data[pos]))
                for num,reactname in enumerate(get_alquimiavector(data.meta_data.aqueous_kinetic_names)):
                    if '->' in reactname and spec in reactname.split('->')[0].split():  # reactants
                        
                        # print('Setting rate of reaction %s to zero'%reactname)
                        data.properties.aqueous_kinetic_rate_cnst.data[num]=0.0
        
        cut_for_flux=False
        flux=numpy.zeros(data.meta_data.primary_names.size)
        for spec in diffquo.keys():
            pos=get_alquimiavector(data.meta_data.primary_names).index(spec)
            bc_conc=bc.total_mobile.data[pos]
            local_conc=data.state.total_mobile.data[pos]
            
            
            # flux[pos] = (bc_conc-local_conc)*diffquo[spec]*actual_dt
            # Run half step of flux first, then other half after chemistry
            flux[pos] = (bc_conc-local_conc)*(1-numpy.exp(-diffquo[spec]*actual_dt/2))
            
            data.state.total_mobile.data[pos] = data.state.total_mobile.data[pos] + flux[pos]
            
            if bc_conc<0:
                raise ValueError('Boundary condition concentration of %s < 0'%spec)
            if local_conc<0:
                raise RuntimeError('Initial concentration of %s < 0'%spec)

            # Here we test if the flux was fast enough to change concentration significantly relative to bc value
            # If the change is more than flux_tol*bc, we cut the time step to resolve changes better
            # if flux[pos] != 0:
            #     print(data.state.total_mobile.data[pos],local_conc,data.state.total_mobile.data[pos] - local_conc,bc_conc)
            if abs(data.state.total_mobile.data[pos] - local_conc) > bc_conc*flux_tol:
                # print('Cutting time step to resolve flux better')
                cut_for_flux=True
                
                
        if cut_for_flux or num_cuts<min_cuts:
            converged=False
        else:
            chem.ReactionStepOperatorSplit(ffi.new('void **',data.engine_state),actual_dt,ffi.addressof(data.properties),ffi.addressof(data.state),
                                            ffi.addressof(data.aux_data),status)
            data.state.porosity=porosity
            converged=status.converged
            # Check for negative concentrations
            # if (numpy.array(get_alquimiavector(data.state.total_mobile))<0).any() or (numpy.array(get_alquimiavector(data.state.total_immobile))<0).any():
            #     converged=False
            
        if converged:
            check_status(status,False)
            
            chem.GetAuxiliaryOutput(ffi.new('void **',data.engine_state),ffi.addressof(data.properties),ffi.addressof(data.state),
                                            ffi.addressof(data.aux_data),ffi.addressof(data.aux_output),status)
            check_status(status,False)
            
            # Second half of flux
            flux=numpy.zeros(data.meta_data.primary_names.size)
            for spec in diffquo.keys():
                pos=get_alquimiavector(data.meta_data.primary_names).index(spec)
                bc_conc=bc.total_mobile.data[pos]
                local_conc=data.state.total_mobile.data[pos]
                
                
                # flux[pos] = (bc_conc-local_conc)*diffquo[spec]*actual_dt
                # Run half step of flux first, then other half after chemistry
                flux[pos] = (bc_conc-local_conc)*(1-numpy.exp(-diffquo[spec]*actual_dt/2))
                
                data.state.total_mobile.data[pos] = data.state.total_mobile.data[pos] + flux[pos]
                
                if bc_conc<0:
                    raise ValueError('Boundary condition concentration of %s < 0'%spec)
                if local_conc<0:
                    raise RuntimeError('Initial concentration of %s < 0'%spec)

                # Here we test if the flux was fast enough to change concentration significantly relative to bc value
                # If the change is more than flux_tol*bc, we cut the time step to resolve changes better
                if abs(data.state.total_mobile.data[pos] - local_conc) > bc_conc*flux_tol:
                    # print('Cutting time step to resolve flux better')
                    cut_for_flux=True
                    converged=False

        if converged:
            self.copy_from_alquimia(data)
            return max_cuts
        else:

            if actual_dt/2<min_dt:
                print(self)
                print('diffquo:')
                for spec in diffquo:
                    print(spec.ljust(14)+'%1.2e'%diffquo[spec])
                if cut_for_flux:
                    raise RuntimeError('Pflotran failed to converge (because of boundary fluxes) after %d cuts to dt = %1.2f s'%(num_cuts,actual_dt))
                else:
                    raise RuntimeError('Pflotran failed to converge after %d cuts to dt = %1.2f s'%(num_cuts,actual_dt))
            
            # data will be reset to layer contents at beginning of next call
            # Run it twice, because we cut the time step in half
            ncuts=self.run_onestep(chem,data,dt,status,min_dt,num_cuts=num_cuts+1,diffquo=diffquo,bc=bc,truncate_concentration=truncate_concentration,rateconstants=rateconstants,flux_tol=flux_tol)
            # If this completes, it means that run_onestep was successful
            self.copy_from_alquimia(data)
            if ncuts>max_cuts:
                max_cuts=ncuts
                
            # This starts from ncuts so it doesn't have to try all the ones that failed again
            for n in range(2**(ncuts-(num_cuts+1))):
                ncuts2=self.run_onestep(chem,data,dt,status,min_dt,num_cuts=ncuts,diffquo=diffquo,bc=bc,truncate_concentration=truncate_concentration,rateconstants=rateconstants,flux_tol=flux_tol)
                self.copy_from_alquimia(data)
                if ncuts2>max_cuts:
                    max_cuts=ncuts2

            return max_cuts
        
        

def convert_to_xarray(layers,t0=0.0,leaf_Mn=None,drop_nas=True,convert_output=True):
    for l in layers:
        if convert_output or not hasattr(l,'output_DF'):
            l.convert_output()
    data_array = xarray.concat([xarray.Dataset.from_dataframe(layer.output_DF) for layer in layers],dim='layer').rename({'layer':'depth'})
    data_array['dz']=xarray.DataArray([layer.volume for layer in layers],dims='depth',attrs={'units':'cm'})*100
    data_array['z_bottom']=data_array['dz'].cumsum()
    data_array['z_top']=data_array['z_bottom']-data_array['dz']
    data_array['z_middle']=data_array['z_bottom']-data_array['dz']/2
    data_array['z_bottom'].attrs['units']='cm'
    data_array['z_middle'].attrs['units']='cm'
    data_array['z_top'].attrs['units']='cm'
    
    data_array['depth']=data_array['z_middle']
    
    for var in layers[0].output_units:
        data_array[var].attrs['units']=layers[0].output_units[var]
        
    # To do: Add other layer properties/attributes and also leaf Mn concentration
    data_array['saturation']=xarray.DataArray([layer.saturation for layer in layers],dims='depth',attrs={'units':'fraction'})
    data_array['BD']=xarray.DataArray([layer.BD for layer in layers],dims='depth',attrs={'units':'g cm-3'})
    data_array['CEC']=xarray.DataArray([layer.CEC for layer in layers],dims='depth',attrs={'units':'meq kg-1'})
    
    data_array['time']=data_array['time']+t0
    data_array['time'].attrs['units']='days'
    
    if leaf_Mn is not None:
        data_array['litter_Mn']=xarray.DataArray(leaf_Mn,dims='litter_year',attrs={'units':'mmol/kg dry mass'})
    
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

reaction_network=Mn.make_network(leaf_Mn_mgkg=0.0,Mn2_scale=1e-4,Mn3_scale=1e-15,NH4_scale=1e-2,DOM_scale=1.0) # We will add the Mn along with leaf litter manually instead of generating it through decomposition

rateconstants={
    'DOM aerobic respiration':1e-7,
    # 'DOM2 aerobic respiration':0.5e-9*0,
    'Hydrolysis':1.5/(365*24*3600),
    'Lignin exposure':1.0/(365*24*3600),
    # 'Lignin depolymerization':0.01/(365*24*3600)*0,
    'Mn Peroxidase':0.2e-8,     # Manganese Peroxidase (produces chelated Mn+++)
    'Mn chelate lignin':0.2e-2, # chelated Mn+++ reacting with lignin/DOM2
    'Mn chelate loss':1.0e-2,   # Chelated Mn+++ spontaneously decomposing
    'Mn-independent lignin degradation':0.4e-8, 
    'Root uptake of Mn++':0.2e-8/100**3,
    # Bandstra et al microbial Mn reduction median rate of 0.0123 mM/hour -> 3.4e-9 M/s
    'DOM1 Mn+++ reduction':2e-9, #1e-7, # microbial Manganese reduction. 
    # Beth says there should be plenty of papers out there about this but maybe not a good synthesis
    'DOM1 Mn+++ abiotic reduction':1e30, # abiotic Mn reduction. Rate constant is multiplied by [Mn+++]^4*[DOM1] so it needs to be very high
    'Bacterial Mn++ oxidation':1e-10, # Bacterial Mn++ oxidation. Rate constant 1e-9 estimated from Fig 8 in Tebo et al 2004 (~35 uM over 10 hours)
    # DOM sorption/desorption
    'DOM desorption':1/(365*24*3600*50), # 1st order rate constant (1/M-s)
    'DOM sorption':1e-10,  # M/(biomass*s)
}

precision=4

incubation_length=5 # Years of litter decomp

molar_mass={'Mg++':24.305,'Al+++':26.982,'K+':39.098,'Ca++':40.078,'Mn':54.938,'Na+':22.99,'N':14.007}
init_exch_cations={'Mg++':1.5,'Al+++':7.0,'K+':1.3,'Ca++':5.0,'Na+':0.2,'Mn++':0.3} # mmol/kg. From Jin et al 2010 Table 3



not_T_sens = [
    'Root uptake of Mn++',
    'DOM1 Mn+++ abiotic reduction',
    'DOM desorption',
    'DOM sorption',
]

def run_sim(Ndep,warming,pH,anox_freq,rateconstants,input_file,Q10=2.0,dt=3600*12,nyears=40,restart_state=None,do_incubation=True,
            fname=None,anox_lenscale=1.0,anox_depthscale=0.125,birnessite_rateconst=2e-11):
    chem,data,sizes,status=init_alquimia(input_file,hands_off=False)
    rateconstants_warmed=rateconstants.copy()
    for react in rateconstants_warmed:
        if react in not_T_sens:
            continue
        rateconstants_warmed[react]=rateconstants_warmed[react]*Q10**(warming/10.0)
    rateconstants_stoich=convert_rateconstants(rateconstants_warmed,reaction_network,precision=precision)
    # Set up layers
    # Top (organic) layer should be thinner and have lower bulk density though
    # Low bulk density causes simulation to slow or crash though. Actually CEC being too low (<100 combined with BD<1) is the problem
    layers=[layer(0.05,rateconstants=rateconstants_stoich,BD=0.425,porosity=0.5,CEC=200.0)]+[layer(0.1,rateconstants=rateconstants_stoich) for num in range(4)]

    # Read secondary complex names from input file since Alquimia does not provide them
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
    
    for l in layers:
        l.secondary_names=secondary_names
    

    # Herndon et al 2014 (BGC): Deep bulk soil Mn concentration ~1500 ug/g = 54 umol/cm3 assuming bulk density=2 g/cm3
    Mn_molarmass=54.94        #g/mol
    molar_volume_birnessite = 251.1700 # Divide by 7 below because birnessite is defined in database as 7 mols of Mn
    Mn_VF=1500e-6/Mn_molarmass*l.BD *molar_volume_birnessite/7
    for l in layers:
        l.initcond=decomp_network.change_constraints(Mn.pools,{'Cellulose':1.0e-8,'Lignin':1.0e-8,'H+':'%1.1f P'%pH,
                                                        # 'Manganite':'1.0d-7  1.d2 m^2/m^3','Mn(OH)2(am)':'%1.2g 1.d2 m^2/m^3'%Mn_VF,
                                                        'Birnessite2':'%1.2g 1.d2 m^2/m^3'%Mn_VF,
                                                        'Mn++':'%1.2g TOTAL_AQ_PLUS_SORB'%(init_exch_cations['Mn++']*1e-6*l.BD*1000/l.porosity), # Jin et al 2020 Table 3: Exch Mn in deeper layers 0.2-0.5 mmol/kg. Convert to mol/L water for AQ_PLUS_SORB (TOTAL_SORB doesn't seem to work)
                                                        'Mg++':'%1.2g TOTAL_AQ_PLUS_SORB'%(init_exch_cations['Mg++']*1e-6*l.BD*1000/l.porosity),
                                                        # 'Ca++':'%1.2g TOTAL_AQ_PLUS_SORB'%(init_exch_cations['Ca++']*1e-6*l.BD*1000/l.porosity),
                                                        'Na+':'%1.2g TOTAL_AQ_PLUS_SORB'%(init_exch_cations['Na+']*1e-6*l.BD*1000/l.porosity),
                                                        'K+':'%1.2g TOTAL_AQ_PLUS_SORB'%(init_exch_cations['K+']*1e-6*l.BD*1000/l.porosity),
                                                        'Al+++':'%1.2g TOTAL_AQ_PLUS_SORB'%(init_exch_cations['Al+++']*1e-6*l.BD*1000/l.porosity),
                                                        'O2(aq)':'0.2 G O2(g)'})

    # initcond=decomp_network.change_site_density(initcond, '>DOM1', 1e4)
    bc=layers[0].initcond
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
            bc_state.temperature=layers[0].temperature
            bc_state.water_density=layers[0].water_density
            bc_state.porosity=layers[0].porosity
            bc_state.aqueous_pressure=layers[0].pressure
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
            
        # CEC also needs to be specified in hands-off mode
        # CEC in PFLOTRAN is in eq/m3, so it must be converted from normal units of meq/kg
        CEC_pf=l.CEC*1e-3*(l.BD*1e-3*100**3)
        print('Applying CEC: %1.2g'%l.CEC)
        data.state.cation_exchange_capacity.data[0]=CEC_pf

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
        * Starting simulation with pH = %1.1f, Ndep = %03d, Warming = %d, anox_freq = %d , anox_len = %1.1f *
        *************************************
        
        
        '''%(pH,int(Ndep/(1000/molar_mass['N']/100**2/(365*24*3600) )),warming,anox_freq,anox_lenscale))


    
    initial_HCO3 = l.total_mobile['HCO3-']
    initial_O2 = l.total_mobile['O2(aq)']

    # Flow rate cm/s = 10 L/m2/s, positive is downward
    # flow_rate=numpy.linspace(1e-7,1e-8,len(layers)) # Rate declines linearly with depth, assumes removal or accumulation in lower layers
    flow_rate=numpy.zeros(len(layers))+1e-7
    min_dt=3 # Seconds
    truncate_concentration=1e-20


    litter_mass=0.163 # kg C/m2/year. Total leaf litter C from Table 1 in Smith et al (2017) Shale Hills C budget paper

    Mn_molarmass=54.94        #g/mol
    C_molarmass=12.01         #g/mol
    litter_ligninfrac=0.25     # Davey et al 2007 Table 2 has oak litter lignin ~250-350 mg/g (.25-.35 g/g)
    litter_Cfrac_mass=0.4       #g/g Davey et al 2007 Table 2 has oak litter C ~ 0.52
    litter_Mn_mg_g_initial=2.0  #mg/g
    
    litter_chem={'Mg++':1120.,'Al+++':35.,'K+':1841.,'Ca++':8721.,'Na+':1e3} # From Beth's data, units of ug/g leaf litter
    


    root_efolding=0.15 # e-folding depth of root biomass in m
    root_biomass_top=0.3 # gC/cm3
    root_Cfrac=0.4

    # Units of z and z_mid are m
    z=numpy.array([0]+[l.volume for l in layers]).cumsum()
    z_mid=(z[:-1]+z[1:])/2 
    for l in range(len(layers)):
        layers[l].total_immobile['Root_biomass']=root_biomass_top*numpy.exp(-z_mid[l]/root_efolding)/(root_Cfrac*12)*100**3
        layers[l].mineral_rate_cnst['Mn(OH)2(am)']=2e-11
        layers[l].mineral_rate_cnst['Birnessite2']=birnessite_rateconst
        # layers[l].surface_site_density['>DOM1']=1e3
        layers[l].total_immobile['Sorption_capacity']=1/12*100**3*0.01 # 1 g/cm3
        # Sinusoid redox state
        # layers[l].diffquo={'O2(aq)':0.001*0.1**(z_mid[l]*10)*(1+numpy.sin(2*numpy.pi*anox_freq/(365*24*3600/dt)*numpy.arange(nsteps)))}
        # Redox state with exponential relaxation, layer-dependent relaxation rate
        if anox_freq>0:
            anox_length=anox_lenscale*numpy.exp(z_mid[l]/anox_depthscale)  # Length of anoxic period in days, by depth
            t_anox=(numpy.arange(nsteps)*dt/(3600*24))%(365//anox_freq) # Time in redox cycle
            # layers[l].diffquo={'O2(aq)':(0.001*(1-numpy.exp(-((numpy.arange(nsteps)*dt/(3600*24))%(365//anox_freq))*0.1**(z_mid[l]*10))))}
            layers[l].diffquo={'O2(aq)':(0.001*numpy.where(t_anox<=anox_length,0.0,1.0))}
        else:
            layers[l].diffquo={'O2(aq)':numpy.atleast_1d(0.001*0.1**(z_mid[l]*10))}
        layers[l].total_mobile['DOM3']=5.3

    # Treat top layer as O horizon with less root biomass and mineral Mn
    layers[0].total_immobile['Root_biomass']=root_biomass_top/(root_Cfrac*12)*100**3*1e-3
    # layers[0].mineral_volume_fraction['Mn(OH)2(am)']=1e-7
    layers[0].mineral_volume_fraction['Birnessite2']=1e-7
    layers[0].surface_site_density['>DOM1']=1e2
    layers[0].total_mobile['DOM3']=1e-10
    layers[0].total_immobile['Sorption_capacity']=1e-10
    

    for l in layers:
        l.write_output(0,dt)
    
    if do_incubation:
        incubation_layer_init=convert_to_xarray([layers[0]])    
        from copy import deepcopy
        incubation_layer=deepcopy(layers[0])
    

    flow_in=numpy.zeros(len(layers),dtype=float)
    flow_out=numpy.zeros(len(layers),dtype=float)
    immobile_specs=['DOM2','Tracer2','Tracer','DOM3','chelated_Mn+++']


    t0=time.time()
    tprev=t0

    leaf_Mn_concs=numpy.ma.masked_all(nyears)

    tstart=0

    # Restart from existing state?
    if isinstance(restart_state,str):
        restart_state=xarray.open_dataset(restart_state)
    if isinstance(restart_state,xarray.Dataset):
        copy_to_layers(restart_state, layers)
        tstart=restart_state.dropna(dim='time')['time'].isel(time=-1).item()
    elif isinstance(restart_state,list) and isinstance(restart_state[0],layer):
        layers=restart_state



    success=True
    for step in range(nsteps):
        # First deposit leaf litter, including last year worth of root Mn uptake
        if (step*dt/3600+tstart*24)%(365*24) == 0: # End of year. Needs dt to be even divisor of one day to work
            layers[0].total_immobile['Cellulose']=layers[0].total_immobile['Cellulose'] + litter_mass*(1-litter_ligninfrac)/layers[0].volume/C_molarmass*1000
            layers[0].total_immobile['Lignin']=layers[0].total_immobile['Lignin'] + litter_mass*litter_ligninfrac/layers[0].volume/C_molarmass*1000
            
            
            # Add up all Mn uptake from previous year to add to leaf litter, and reset those pools
            Mn_uptake_total=0.0 # mol Mn
            for l in layers:
                Mn_uptake_total=Mn_uptake_total+l.total_mobile['Tracer2']*l.volume*1000*l.porosity*l.saturation
                l.total_mobile['Tracer2']=0.0

            layers[0].total_mobile['Mn++']=layers[0].total_mobile['Mn++'] + Mn_uptake_total/(layers[0].porosity*layers[0].saturation*layers[0].volume*1000)
            for spec in litter_chem:
                layers[0].total_mobile[spec]=layers[0].total_mobile[spec] + litter_chem[spec]*1e-6*litter_mass/litter_Cfrac_mass*1000/molar_mass[spec]
        
            print('Leaf litter Mn concentration = %1.1f mmol/kg'%(Mn_uptake_total*1e3/(litter_mass/litter_Cfrac_mass)))
            print('Layer Mn++ concentration: %1.2g M'%layers[0].total_mobile['Mn++'])
            leaf_Mn_concs[int((step*dt)/(365*24*3600))]=Mn_uptake_total*1e3/(litter_mass/litter_Cfrac_mass)
            print('Top layer birnessite: %1.2g ug/g. Bottom layer: %1.2g ug/g'%(layers[0].mineral_volume_fraction['Birnessite2']*7/molar_volume_birnessite*1e6*Mn_molarmass/layers[0].BD,
                            layers[-1].mineral_volume_fraction['Birnessite2']*7/molar_volume_birnessite*1e6*Mn_molarmass/layers[-1].BD))
                            
            # Incubation layer gets reset with just one year of litter biomass
            if do_incubation:
                if (step*dt/3600+tstart*24)%(365*24*incubation_length) == 0:
                    copy_to_layers(incubation_layer_init,[incubation_layer])
                    incubation_layer.total_immobile['Cellulose']=litter_mass*(1-litter_ligninfrac)/incubation_layer.volume/C_molarmass*1000
                    incubation_layer.total_immobile['Lignin']= litter_mass*litter_ligninfrac/incubation_layer.volume/C_molarmass*1000
                    incubation_layer.total_mobile['Mn++']= Mn_uptake_total/(layers[0].porosity*layers[0].saturation*incubation_layer.volume*1000)
        
        # Next advect mobile species. This assumes we can separate advective transport and reactions at this time step length
        # Alternate strategy would be to pass these rates to run_onestep and spread them over the variable timesteps
        for spec in layers[0].primary_names:
            if spec in immobile_specs:
                continue # Skip this one because it's not really mobile, more a workaround for Lignin decomp
            flow_in[:]=0.0
            flow_out[:]=0.0
            for num in range(len(layers)-1):
                # First calculate net flow for each layer
                flow_out[num]=flow_rate[num]*layers[num].total_mobile[spec]*10 # mol/L*cm/s*(10L/m2)/cm -> mol/m2/s
                flow_in[num+1]=flow_rate[num]*layers[num].total_mobile[spec]*10
            # Leaching out of bottom layer (to groundwater/loss)
            flow_out[-1]=flow_rate[-1]*layers[-1].total_mobile[spec]*10
            for num in range(len(layers)):
                # Then calculate change in concentration using amount of water in the layer
                # Units of mobile species are M (mol/L water)
                spec_mol_layer=layers[num].total_mobile[spec]*layers[num].volume*1000*layers[num].porosity*layers[num].saturation # mol of the species in the layer initially
                spec_mol_layer=spec_mol_layer+(flow_in[num]-flow_out[num])*dt
                layers[num].total_mobile[spec]=spec_mol_layer/(layers[num].volume*1000*layers[num].porosity*layers[num].saturation)
                layers[num].flow_in[spec]=flow_in[num]
                layers[num].flow_out[spec]=flow_out[num]

        # Skip flows for incubation layer I guess?
        
        # N deposition as NH4+. Ndep units of mol N/m2/s. Need to convert to mol/L
        layers[0].total_mobile['NH4+'] = layers[0].total_mobile['NH4+']+Ndep*dt/(layers[0].volume*1000*layers[0].porosity*layers[0].saturation)
        if do_incubation:
            incubation_layer.total_mobile['NH4+'] = incubation_layer.total_mobile['NH4+']+Ndep*dt/(incubation_layer.volume*1000*incubation_layer.porosity*incubation_layer.saturation)

        # Just keep HCO3- equilibrated with the atmosphere in top layer under oxic conditions
        # Also remove equal amount of H+ because we are tracking bicarbonate, not CO2(aq). CO2(aq) + H2O = HCO3- + H+
        # H+ gets removed from immobile pool because aqueous H+ concentration might be less than HCO3- concentration
        for num in range(len(layers)):
            if layers[num].diffquo['O2(aq)'][step%len(layers[num].diffquo['O2(aq)'])]>0:
                dCO2=layers[num].total_mobile['HCO3-']-initial_HCO3
                layers[num].total_mobile['HCO3-']=layers[num].total_mobile['HCO3-']-dCO2
                layers[num].total_immobile['H+']=layers[num].total_immobile['H+']-dCO2*1000*layers[num].porosity*layers[num].saturation
        
        if do_incubation:
            dCO2=incubation_layer.total_mobile['HCO3-']-initial_HCO3
            incubation_layer.total_mobile['HCO3-']=incubation_layer.total_mobile['HCO3-']-dCO2
            incubation_layer.total_immobile['H+']=incubation_layer.total_immobile['H+']-dCO2*1000*incubation_layer.porosity*incubation_layer.saturation
        # Equilibrate O2 also
        # layers[0].total_mobile['O2(aq)']=initial_O2
        # Should set up stepper so it can handle equilibrium boundary conditions better      
                
        # Next do chemistry.
        try:
            if do_incubation:
                todo=enumerate(layers+[incubation_layer])
            else:
                todo=enumerate(layers)
            for n,l in todo:
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
                            # print('O2 before: %1.1g'%data.state.total_mobile.data[get_alquimiavector(data.meta_data.primary_names).index('O2(aq)')])
                num_cuts=l.run_onestep(chem,data,dt,status,min_dt=min_dt,diffquo=dq,bc=bc_state,truncate_concentration=truncate_concentration,rateconstants=l.rateconstants,flux_tol=0.8)
                l.copy_from_alquimia(data)
                # Write output
                l.write_output(step+1,dt,num_cuts)
                # print('O2 after: %1.1g'%data.state.total_mobile.data[get_alquimiavector(data.meta_data.primary_names).index('O2(aq)')])
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
            print('*** Step {step:d} of {nsteps:d} ({nyears:1.1f} of {totalyears:d} years). Time elapsed: {t:1.1f} min ({tperstep:1.2f} s per {steplength:1.1f} hour timestep). Mean cuts: {meancuts:s} Mean dt: {meandt:1.1f} s ***'.format(
                    step=step,nsteps=nsteps,t=(t1-t0)/60,tperstep=(t1-tprev)/25,meancuts=str(cuts),meandt=mean_dt,steplength=dt/3600,nyears=step*dt/(3600*24*365),totalyears=int(nsteps*dt/(3600*24*365))),flush=True)
            tprev=t1



    output=convert_to_xarray(layers,t0=tstart,leaf_Mn=leaf_Mn_concs)
    if isinstance(restart_state,xarray.Dataset):
        output=xarray.concat([restart_state,output.isel(time=slice(1,None))],dim='time')
        
        
    import datetime
    today=datetime.datetime.today()
    if fname is None:
        fname='/lustre/or-hydra/cades-ccsi/scratch/b0u/Mn_output/Mn_output_{year:04d}-{month:02d}-{day:02d}.nc'.format(year=today.year,month=today.month,day=today.day)
    gname=make_gname(pH,Ndep,warming,anox_freq,anox_lenscale,birnessite_rateconst)
    
    newdims=['soil_pH','Ndep','warming','redox_cycles','anox_lenscales','birnrate']
    output['soil_pH']=pH
    output['Ndep']=int(Ndep/(1000/molar_mass['N']/100**2/(365*24*3600)))
    output['warming']=warming
    output['redox_cycles']=anox_freq 
    output['anox_lenscales']=anox_lenscale 
    output['warming']=warming
    output['birnrate']=round(birnessite_rateconst*1e12,5)
    output=output.expand_dims(newdims).set_coords(newdims)

    import os
    if not os.path.exists(fname):
        output.to_netcdf(fname,mode='w',group=gname)
    else:
        output.to_netcdf(fname,mode='a',group=gname)
    if do_incubation:
        convert_to_xarray([incubation_layer],leaf_Mn=leaf_Mn_concs).to_netcdf(fname,mode='a',group=gname+'_incubation')

    return success

def make_gname(pH,Ndep,warming,anox_freq,anox_lenscale,birnessite_rateconst):
    gname='pH{ph:1.1f}_Ndep{Ndep:03d}_warming{warming:d}_anox_freq{redox:d}_anox_len{anoxlen:1.1f}_birnrate{birnrate:1.1f}'.format(ph=pH,Ndep=int(Ndep/(1000/molar_mass['N']/100**2/(365*24*3600) )),warming=int(warming),redox=int(anox_freq),anoxlen=anox_lenscale,birnrate=birnessite_rateconst*1e12)
    return gname

def setup_sims():
    import pandas
    # Whalen et al 2018: Background N dep is 8-10 kg N/ha/year
    # Treatments were +50 kg N/ha/year and +150 kgN/ha/year (as NH4NO3)
    # 1 kg N/ha/year
    Ndeps=   [0,50,150,0,50,150,0,50,150]
    warmings=[0, 0,  0,2,2,  2, 5,5 ,5] # Degrees C

    pHs=numpy.arange(4.0,6.5,0.5)
    # pHs=[4.5,6.0]
    # anox_freqs=[12,8,4,1] # anoxic events per year
    anox_freqs=[50] # ~ weekly anoxic periods
    anox_lengths=[0.25,0.5,2]
    birn_rates=1e-12*2.0**numpy.arange(-2,8,2)

    Ndep_sims=[]
    pH_sims=[]
    anox_freq_sims=[]
    anox_len_sims=[]
    warming_sims=[]

    # Consider rearranging this so main thing is warming x Ndep x pH x birn_rate
    # and redox is a separate set?
    birnrate_sims=[]

    # Build one long list of all the sim params first
    for sim in range(len(Ndeps)):
    # for Ndep in numpy.array([0,50,150])*1000/molar_mass['N']/100**2/(365*24*3600): # Converted to mol N/m2/s
        
        warming=warmings[sim]

        for pH in pHs:
            for birn_rate in birn_rates:
                Ndep_sims.append(Ndeps[sim])
                warming_sims.append(warming)
                anox_freq_sims.append(50)
                anox_len_sims.append(1.0)
                pH_sims.append(pH)
                birnrate_sims.append(birn_rate)

    for anox_len in anox_lengths:
        Ndep_sims.append(0)
        warming_sims.append(0)
        anox_freq_sims.append(50)
        anox_len_sims.append(anox_len)
        pH_sims.append(pHs[2])
        birnrate_sims.append(birn_rates[2])

    return pandas.DataFrame({'Ndep':Ndep_sims,'warming':warming_sims,'anox_freq':anox_freq_sims,'anox_len':anox_len_sims,'pH':pH_sims,'birn_rate':birnrate_sims})

if __name__ == '__main__':

    import time,datetime
    import sys
    starting_time=time.time()


    from argparse import ArgumentParser
    
    parser = ArgumentParser()
    parser.add_argument('-f',dest='fname',help='Output file name',default='')
    parser.add_argument('-n',dest='jobnum',help='Job number',default=1)
    parser.add_argument('-N',dest='totaljobs',help='Total number of jobs',default=1)
    options = parser.parse_args()

    if options.fname != '':
        fname=options.fname
    else:
        today=datetime.datetime.today()
        fname='/lustre/or-scratch/cades-ccsi/b0u/Mn_output/Mn_output_{year:04d}-{month:02d}-{day:02d}.nc'.format(year=today.year,month=today.month,day=today.day)
   
    jobnum=int(options.jobnum)-1
    totaljobs=int(options.totaljobs)
    # Set up for parallel jobs
    if jobnum+1>totaljobs:
        raise ValueError('jobnum + 1 > totaljobs')
    if totaljobs>1:
        fname=fname[:-3]+'_%02d.nc'%jobnum

    input_file='/lustre/or-scratch/cades-ccsi/b0u/Mn_output/manganese_%02d.in'%jobnum

    decomp_network.PF_network_writer(reaction_network,precision=precision).write_into_input_deck('SOMdecomp_template.txt',input_file,log_formulation=True,
            CO2name='Tracer',truncate_concentration=1e-25,database='/home/b0u/models/PFLOTRAN/REDOX-PFLOTRAN/hanford.dat')


    allsims=setup_sims()
    # Then just run sims for this job
    sims=list(range(jobnum,len(allsims),totaljobs))

    print('Total number of sims: %d'%len(allsims))
    print('This job: ',sims)
    
    for simnum in sims:
        this_sim=allsims.loc[simnum]
        Ndep=this_sim['Ndep']*1000/molar_mass['N']/100**2/(365*24*3600)
        pH=this_sim['pH']
        warming=this_sim['warming']
        anox_freq=this_sim['anox_freq']
        anox_len=this_sim['anox_len']
        birn_rate=this_sim['birn_rate']
        run_sim(Ndep,warming,pH,anox_freq,rateconstants,input_file,Q10=2.0,dt=3600*6,fname=fname,anox_lenscale=anox_len,birnessite_rateconst=birn_rate,
            do_incubation=False,nyears=80)

    print('\n\n\n Simulations finished. Total time: %1.1f hours\n'%((time.time()-starting_time)/3600),flush=True)

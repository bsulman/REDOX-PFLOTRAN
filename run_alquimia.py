from .decomp_network import get_stoich_from_name
from ._alquimia import ffi,lib
import numpy

import sys

def check_status(status,printmessage=True):
    if status.error != 0:
        message=ffi.string(status.message)
        lib.FreeAlquimiaEngineStatus(status)
        raise RuntimeError(message)
    else:
        if printmessage:
            print(str(ffi.string(status.message),encoding='utf8'))

def print_metadata(data):
    for datatype in ['aqueous_kinetic_names','ion_exchange_names','isotherm_species_names','mineral_names','primary_names','surface_site_names']:
        field=getattr(data.meta_data,datatype)
        print('%d %s:'%(field.size,datatype))
        for num in range(field.size):
            print('   '+str(ffi.string(field.data[num]),encoding='utf8'))

def get_alquimiavector(vec):
    s=vec.size
    if s==0:
        return []
    elif isinstance(vec.data[0],(float,int)):
        return [vec.data[n] for n in range(s)]
    else: # Assume it's a string
        return [str(ffi.string(vec.data[n]),encoding='utf8') for n in range(s)]

def print_alquimia_object(obj,indent=0):
    if hasattr(obj,'size'):
        # print(' '*indent+str(get_alquimiavector(obj)))
        for element in range(obj.size):
            print(' '*indent+'---')
            print_alquimia_object(obj.data[element],indent+2)
            print(' '*indent+'---')
    elif isinstance(obj,(int,float)):
        print(' '*indent+str(obj))
    elif ffi.typeof(obj) == ffi.typeof(ffi.new('char *')):
        print(' '*indent+str(ffi.string(obj),encoding='utf8'))
    elif isinstance(obj,ffi.CData):
        for field in dir(obj):
            print(' '*indent+field)
            print_alquimia_object(getattr(obj,field),indent+2)
    else:
        print(' '*indent+obj)

def convert_rateconstants(rateconstants,reactions,precision=2):
    import decomp_network
    return dict([(decomp_network.get_stoich_from_name(name,reactions,precision),rateconstants[name]) for name in rateconstants])

def convert_output(output,meta_data,secondary_names):
    import pandas
    output_DF=pandas.DataFrame(index=output['time'])
    output_DF['Porosity']=output['porosity']
    output_DF['ncuts']=output['ncuts']
    output_DF['actual_dt']=output['actual_dt']
    output_units={'time':'days','Porosity':'NA','ncuts':'NA','actual_dt':'s'}
    total=pandas.DataFrame(output['total_mobile'],columns=get_alquimiavector(meta_data.primary_names),index=output['time']).add_prefix('Total ')
    for col in total.columns:
        output_DF[col]=total[col]
    output_units.update([('Total '+s,'M') for s in get_alquimiavector(meta_data.primary_names)])
    free=pandas.DataFrame(output['free'],columns=get_alquimiavector(meta_data.primary_names),index=output['time']).add_prefix('Free ')
    for col in free.columns:
        output_DF[col]=free[col]
    output_units.update([('Free '+s,'M') for s in get_alquimiavector(meta_data.primary_names)])
    sorbed=pandas.DataFrame(output['immobile'],columns=get_alquimiavector(meta_data.primary_names),index=output['time']).add_prefix('Total Sorbed ')
    for col in sorbed.columns:
        output_DF[col]=sorbed[col]
    output_units.update([('Total Sorbed '+s,'mol/m^3') for s in get_alquimiavector(meta_data.primary_names)])
    mineral_VF=pandas.DataFrame(output['mineral_VF'],columns=get_alquimiavector(meta_data.mineral_names),index=output['time']).add_suffix(' VF')
    for col in mineral_VF.columns:
        output_DF[col]=mineral_VF[col]
    output_units.update([(s+' VF','m^3 mnrl/m^3 bulk') for s in get_alquimiavector(meta_data.mineral_names)])
    mineral_rate=pandas.DataFrame(output['mineral_rate'],columns=get_alquimiavector(meta_data.mineral_names),index=output['time']).add_suffix(' Rate')
    for col in mineral_rate.columns:
        output_DF[col]=mineral_rate[col]
    output_units.update([(s+' Rate','mol/m^3/sec') for s in get_alquimiavector(meta_data.mineral_names)])
    secondary=pandas.DataFrame(output['aq_complex'],columns=secondary_names,index=output['time'])
    for col in secondary.columns:
        output_DF[col]=secondary[col]
    output_units.update([(s,'M') for s in secondary_names])
    
    return output_DF.reset_index(drop=True).set_index(output_DF.index/(24*3600)),output_units

def run_onestep(chem,data,dt,status,min_dt=0.1,num_cuts=0,diffquo={},bc=None,flux_tol=0.15,truncate_concentration=0.0,rateconstants={},min_cuts=0,hands_off=False):
    converged=False
    porosity=data.state.porosity
    
    max_cuts=num_cuts
    actual_dt=dt/2**num_cuts
    
    if not hands_off:
        for num,reactname in enumerate(get_alquimiavector(data.meta_data.aqueous_kinetic_names)):
            data.properties.aqueous_kinetic_rate_cnst.data[num]=rateconstants[reactname]
        
        for spec in diffquo.keys():
            pos=get_alquimiavector(data.meta_data.primary_names).index(spec)
            if data.state.total_mobile.data[pos] < truncate_concentration and data.state.total_mobile.data[pos] != 0.0:
                print('Truncating concentration of {spec:s} from {conc:1.1g}'.format(spec=spec,conc=data.state.total_mobile.data[pos]))
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
                        
        if cut_for_flux:
            raise RuntimeError('Cut for flux after chemistry')
        
        return max_cuts
    else:
        # Undo influx because we will do it again in the shorter run
        for key in diffquo.keys():
            pos=get_alquimiavector(data.meta_data.primary_names).index(key)
            data.state.total_mobile.data[pos] = data.state.total_mobile.data[pos] - flux[pos]
            
        if actual_dt/2<min_dt:
            if cut_for_flux:
                raise RuntimeError('Pflotran failed to converge (because of boundary fluxes) after %d cuts to dt = %1.2f s'%(num_cuts,actual_dt))
            else:
                raise RuntimeError('Pflotran failed to converge after %d cuts to dt = %1.2f s'%(num_cuts,actual_dt))
            
        # Run it twice, because we cut the time step in half
        ncuts=run_onestep(chem,data,dt,status,min_dt,num_cuts=num_cuts+1,diffquo=diffquo,bc=bc,truncate_concentration=truncate_concentration,rateconstants=rateconstants,hands_off=hands_off)
        if ncuts>max_cuts:
            max_cuts=ncuts
            
        # This starts from ncuts so it doesn't have to try all the ones that failed again
        for n in range(2**(ncuts-(num_cuts+1))):
            ncuts2=run_onestep(chem,data,dt,status,min_dt,num_cuts=ncuts,diffquo=diffquo,bc=bc,truncate_concentration=truncate_concentration,rateconstants=rateconstants,hands_off=hands_off)
            if ncuts2>max_cuts:
                max_cuts=ncuts2

        return max_cuts
    

def convert_condition_to_alquimia(cond,name):
    if cond is None:
        condition=ffi.new("AlquimiaGeochemicalCondition *")
        lib.AllocateAlquimiaGeochemicalCondition(len(name), 0, 0, condition);
        condition.name[0:len(name)]=name.encode()
        return condition


    primary=[]
    immobile=[]
    mineral=[]
    surf_complex=[]
    
    for constraint in cond:
        if constraint['kind']=='primary':
            primary.append(constraint)
        if constraint['kind']=='immobile':
            if 'initCN' in constraint:
                # For flexible CN pools, need to add separate C and N pool constraints
                from .decomp_network import decomp_pool
                immobile.append(decomp_pool(name=constraint['name']+'C',kind='immobile',constraints={name:constraint['constraints'][name]}))
                immobile.append(decomp_pool(name=constraint['name']+'N',kind='immobile',constraints={name:constraint['constraints'][name]/constraint['initCN']}))
            else:
                immobile.append(constraint)
        if constraint['kind']=='mineral':
            mineral.append(constraint)
        if constraint['kind']=='surf_complex':
            surf_complex.append(constraint)
    
    condition=ffi.new('AlquimiaGeochemicalCondition *')
    lib.AllocateAlquimiaGeochemicalCondition(len(name),len(primary)+len(immobile),len(mineral),condition)
    condition.name[0:len(name)]=name.encode()
    
    for num in range(len(primary)):
        constraint=primary[num]['constraints'][name]
        lib.AllocateAlquimiaAqueousConstraint(ffi.addressof(condition.aqueous_constraints.data[num]))
        condition.aqueous_constraints.data[num].primary_species_name[0:len(primary[num]['name'])]=primary[num]['name'].encode()
        if isinstance(constraint,float):
            condition.aqueous_constraints.data[num].constraint_type[0:len('total_aqueous')]='total_aqueous'.encode()
            condition.aqueous_constraints.data[num].value=constraint
        elif len(constraint.split())==1:
            condition.aqueous_constraints.data[num].constraint_type[0:len('total_aqueous')]='total_aqueous'.encode()
            condition.aqueous_constraints.data[num].value=float(constraint.replace('d','e'))
        elif constraint.split()[1]=='P':
            condition.aqueous_constraints.data[num].constraint_type[0:2]='pH'.encode()
            condition.aqueous_constraints.data[num].value=float(constraint.split()[0].replace('d','e'))
        elif constraint.split()[1]=='G':
            condition.aqueous_constraints.data[num].constraint_type[0:3]='gas'.encode()
            condition.aqueous_constraints.data[num].value=float(constraint.split()[0].replace('d','e'))
            condition.aqueous_constraints.data[num].associated_species[0:len(constraint.split()[2])]=constraint.split()[2].encode()
        elif constraint.split()[1]=='M':
            condition.aqueous_constraints.data[num].constraint_type[0:len('mineral')]='mineral'.encode()
            condition.aqueous_constraints.data[num].value=float(constraint.split()[0].replace('d','e'))
            condition.aqueous_constraints.data[num].associated_species[0:len(constraint.split()[2])]=constraint.split()[2].encode()
        elif constraint.split()[1]=='TOTAL_SORB':
            condition.aqueous_constraints.data[num].constraint_type[0:len('total_sorbed')]='total_sorbed'.encode()
            condition.aqueous_constraints.data[num].value=float(constraint.split()[0].replace('d','e'))
        elif constraint.split()[1]=='TOTAL_AQ_PLUS_SORB':
            condition.aqueous_constraints.data[num].constraint_type[0:len('total_aqueous_plus_sorbed')]='total_aqueous_plus_sorbed'.encode()
            condition.aqueous_constraints.data[num].value=float(constraint.split()[0].replace('d','e'))
        elif constraint.split()[1]=='Z':
            condition.aqueous_constraints.data[num].constraint_type[0:len('charge')]='charge'.encode()
            condition.aqueous_constraints.data[num].value=float(constraint.split()[0].replace('d','e'))
        else:
            raise ValueError('Unrecognized constraint type: %s'%constraint)

    for num in range(len(immobile)):
        constraint=immobile[num]['constraints'][name]
        lib.AllocateAlquimiaAqueousConstraint(ffi.addressof(condition.aqueous_constraints.data[num+len(primary)]))
        condition.aqueous_constraints.data[num+len(primary)].primary_species_name[0:len(immobile[num]['name'])]=immobile[num]['name'].encode()
        # Need to check within alquimia to make this an immobile constraint, I guess
        condition.aqueous_constraints.data[num+len(primary)].constraint_type[0:len('total_sorbed')]='total_sorbed'.encode()
        condition.aqueous_constraints.data[num+len(primary)].value=constraint
        
    for num in range(len(mineral)):
        constraint=mineral[num]['constraints'][name]
        lib.AllocateAlquimiaMineralConstraint(ffi.addressof(condition.mineral_constraints.data[num]))
        condition.mineral_constraints.data[num].mineral_name[0:len(mineral[num]['name'])]=mineral[num]['name'].encode()
        condition.mineral_constraints.data[num].volume_fraction=float(constraint.split()[0].replace('d','e'))
        condition.mineral_constraints.data[num].specific_surface_area=float(constraint.split()[1].replace('d','e'))
        
        
    return condition


def init_alquimia(input_file,hands_off=True):
    engine_name='pflotran'
    status=ffi.new('AlquimiaEngineStatus *') 
    lib.AllocateAlquimiaEngineStatus(status)

    chem=ffi.new('AlquimiaInterface *')

    lib.CreateAlquimiaInterface(engine_name.encode(),chem,status)
    check_status(status)


    data=ffi.new('AlquimiaData *')
    sizes=ffi.addressof(data.sizes)
    functionality=ffi.addressof(data.functionality)
    engine_state=ffi.new('void **',data.engine_state)

    # Initialize Petsc first
    from mpi4py import MPI
    initialized = ffi.new('PetscBool *')
    lib.PetscInitialized(initialized)   
    if not initialized[0]:
        err = lib.PetscInitializeNoArguments()
        assert err==0, 'Error in Petsc initialization'

    # input_file='../../alquimia-dev/benchmarks/batch_chem/ABCD_microbial.in'
    chem.Setup(input_file.encode(),ffi.cast('_Bool',hands_off),engine_state,sizes,functionality,status)
    check_status(status)
    data.engine_state=engine_state[0]

    # Allocates memory for all components of data structure
    lib.AllocateAlquimiaData(data)

    # Get metadata
    chem.GetProblemMetaData(engine_state,ffi.addressof(data.meta_data),status)
    check_status(status,False)

    print_metadata(data)
    
    return (chem,data,sizes,status)

def run_simulation(input_file,simlength_days,dt=3600*12,min_dt=0.1,volume=1.0,saturation=1.0,hands_off=True,
                    water_density=1000.0,temperature=25.0,porosity=0.25,pressure=101325.0,initcond=None,bc=None,
                    diffquo={},rateconstants={},truncate_concentration=0.0,CEC=None,printstep=50,run_name=''):
    
    import time
    t0=time.time()
    nsteps=int(simlength_days/(dt/(3600*24)))

    print('Starting initialization: %s'%run_name)

    engine_name='pflotran'
    status=ffi.new('AlquimiaEngineStatus *') 
    lib.AllocateAlquimiaEngineStatus(status)

    chem=ffi.new('AlquimiaInterface *')

    lib.CreateAlquimiaInterface(engine_name.encode(),chem,status)
    check_status(status)


    data=ffi.new('AlquimiaData *')
    sizes=ffi.addressof(data.sizes)
    functionality=ffi.addressof(data.functionality)
    engine_state=ffi.new('void **',data.engine_state)
    
    # Initialize Petsc first
    from mpi4py import MPI # This initializes MPI. Necessary for it to work on CADES due to MPI installation
    initialized = ffi.new('PetscBool *')
    lib.PetscInitialized(initialized)   
    if not initialized[0]:
        err = lib.PetscInitializeNoArguments()
        assert err==0, 'Error in Petsc initialization'

    # input_file='../../alquimia-dev/benchmarks/batch_chem/ABCD_microbial.in'
    chem.Setup(input_file.encode(),ffi.cast('_Bool',hands_off),engine_state,sizes,functionality,status)
    check_status(status)
    data.engine_state=engine_state[0]

    # Allocates memory for all components of data structure
    lib.AllocateAlquimiaData(data)

    # Get metadata
    chem.GetProblemMetaData(engine_state,ffi.addressof(data.meta_data),status)
    check_status(status,False)

    print_metadata(data)
            

    # Set up initial condition
    if porosity <= 0:
        raise ValueError('Porosity must be > 0')
    if water_density <= 0:
        raise ValueError('Water density must be > 0')
    if pressure <= 0:
        raise ValueError('pressure must be > 0')

    # Initialize state data
    data.properties.volume=volume
    data.properties.saturation=saturation

    data.state.temperature=temperature
    data.state.water_density=water_density
    data.state.porosity=porosity
    data.state.aqueous_pressure=pressure
    
    # Set properties: surface site density and mineral rate constants
    # This is necessary when running in hands-on mode
    if initcond is not None:
        for constraint in initcond:
            if constraint['kind']=='surf_complex':
                sitename=constraint['name']
                sitenum=get_alquimiavector(data.meta_data.surface_site_names).index(sitename)
                print('Applying site density: {name:s} (position={pos:d}): {dens:1.1f}'.format(name=sitename,pos=sitenum,dens=constraint['site_density']))
                data.state.surface_site_density.data[sitenum]=constraint['site_density']
            elif constraint['kind']=='mineral':
                name=constraint['name']
                num=get_alquimiavector(data.meta_data.mineral_names).index(name)
                data.properties.mineral_rate_cnst.data[num]=float(constraint['rate'].split()[0].replace('d','e'))
    if CEC is not None:
        print('Applying CEC: %1.2e'%CEC)
        data.state.cation_exchange_capacity.data[0]=CEC
                
    # Aqueous kinetic rate constants also need to be specified in hands-off mode
    # Alquimia interface always sets backward rate to zero
    if not hands_off:
        for num,reactname in enumerate(get_alquimiavector(data.meta_data.aqueous_kinetic_names)):
            data.properties.aqueous_kinetic_rate_cnst.data[num]=rateconstants[reactname]

    init_cond=convert_condition_to_alquimia(initcond,'initial')
    chem.ProcessCondition(ffi.new('void **',data.engine_state),init_cond,ffi.addressof(data.properties),ffi.addressof(data.state),ffi.addressof(data.aux_data),status)
    check_status(status,False)
    
    # uimia_object(data.state)
    print_alquimia_object(data.properties)

    # Pflotran sets porosity based on minerals or something? Needs to be reset
    data.state.porosity=porosity

    # Set up boundary condition if applicable
    if bc is not None:
        bc_cond=convert_condition_to_alquimia(bc,'initial')
        bc_state=ffi.new('AlquimiaState *')
        bc_state.temperature=temperature
        bc_state.water_density=water_density
        bc_state.porosity=porosity
        bc_state.aqueous_pressure=pressure
        bc_auxdata=ffi.new('AlquimiaAuxiliaryData *')
        lib.AllocateAlquimiaState(sizes,bc_state)
        for num in range(data.state.surface_site_density.size):
            bc_state.surface_site_density.data[num]=data.state.surface_site_density.data[num]
        lib.AllocateAlquimiaAuxiliaryData(sizes,bc_auxdata)
        chem.ProcessCondition(ffi.new('void **',data.engine_state),bc_cond,ffi.addressof(data.properties),bc_state,bc_auxdata,status)
        check_status(status,False)
    else:
        bc_state=None

    ##### At this point, the model should be initialized ##########
    print('''

    *****************************************************
    Successfully initialized alquimia geochemical engine
    Simulation: %s
    *****************************************************

    '''%run_name)

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

    c=cell(volume=volume,rateconstants=rateconstants,diffquo=diffquo,secondary_names=secondary_names)
    c.copy_from_alquimia(data)

    c.setup_output(nsteps+1,dt)
    c.write_output(0,dt)


    tprev=t0
    num_cuts=0 # Keep track of cuts from previous step, so it doesn't start from zero each time
    success=True
    for step in range(1,nsteps):
        try:
            c.copy_to_alquimia(data)
            dq={}
            if bc is not None:
                primarynames=get_alquimiavector(data.meta_data.primary_names)
                for spec in primarynames:
                    if spec in c.diffquo.keys():
                        if numpy.iterable(c.diffquo[spec]):
                            dq[spec]=c.diffquo[spec][step%len(c.diffquo[spec])]
                        else:
                            dq[spec]=c.diffquo[spec]
                    else:
                        dq[spec]=0.0
                        # print('O2 before: %1.1g'%data.state.total_mobile.data[get_alquimiavector(data.meta_data.primary_names).index('O2(aq)')])
            num_cuts=c.run_onestep(chem,data,dt,status,min_dt=min_dt,diffquo=dq,bc=bc_state,truncate_concentration=truncate_concentration,rateconstants=rateconstants,hands_off=hands_off)
            c.copy_from_alquimia(data)
            # Write output
            c.write_output(step,dt,num_cuts)
            # print('O2 after: %1.1g'%data.state.total_mobile.data[get_alquimiavector(data.meta_data.primary_names).index('O2(aq)')])
        except RuntimeError as err:
            print('ERROR on timestep %d: %s'%(step,err))
            print('Returning output so far')
            c.write_output(step,dt,num_cuts)
            success=False
            break
        except KeyboardInterrupt as err:
            print('INTERRUPTED on timestep %d'%(step))
            print('Returning output so far')
            c.write_output(step,dt,num_cuts)
            success=False
            break
                
        
        if step%printstep==0:
            t1=time.time()
            print('*** Step {step:d} of {nsteps:d}. Time elapsed: {t:d} s ({tperstep:1.1g} s per {steplength:1.1g} hour timestep). Mean cuts: {meancuts:1.1f} Mean dt: {meandt:1.1f} s ***'.format(
                    step=step,nsteps=nsteps,t=int(t1-t0),tperstep=(t1-tprev)/25,meancuts=c.output['ncuts'][step-10:step].mean(),meandt=c.output['actual_dt'][step-10:step].mean(),steplength=dt/3600),flush=True)
            tprev=t1
        
    

    # Need to free all the Alquimia arrays?
    lib.FreeAlquimiaData(data)
    
    c.convert_output()
    
    if success:
        print('''

        **********************************************************
        *           Successfully finished run %s in %1.1f minutes            
        **********************************************************

        '''%(run_name,(time.time()-t0)/60))
    return c.output_DF,c.output_units


class cell:
    def __init__(self,volume,saturation=1.0,temperature=20.0,water_density=1000.0,porosity=0.25,pressure=101325.0,BD=1.5,CEC=50.0,rateconstants={},diffquo={},secondary_names=[]):
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
        self.diffquo=diffquo
        self.BD=BD   # g/cm3
        self.CEC=CEC # meq/kg
        self.secondary_names=secondary_names
        
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
            'CEC H+':numpy.ma.masked_all(nsteps,dtype=float)
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
        if '>Carboxylate-' in self.surface_site_names:
            self.output['CEC H+'][step]=self.total_immobile['H+']-(self.surface_site_density['>Carboxylate-']*self.mineral_volume_fraction['Rock(s)']-self.aux_doubles[-len(self.surface_site_names)+self.surface_site_names.index('>Carboxylate-')])
        else:
            self.output['CEC H+'][step]=self.total_immobile['H+']

    def convert_output(self):
        import pandas
        output_DF=pandas.DataFrame(index=self.output['time'])
        output_DF['Porosity']=self.output['porosity']
        output_DF['ncuts']=self.output['ncuts']
        output_DF['actual_dt']=self.output['actual_dt']
        output_units={'time':'s','Porosity':'NA','ncuts':'NA','actual_dt':'s'}
        total=pandas.DataFrame(self.output['total_mobile'],columns=self.primary_names,index=self.output['time']).add_prefix('Total ')
        # for col in total.columns:
        #     output_DF[col]=total[col]
        output_units.update([('Total '+s,'M') for s in self.primary_names])
        free=pandas.DataFrame(self.output['free'],columns=self.primary_names,index=self.output['time']).add_prefix('Free ')
        # for col in free.columns:
        #     output_DF[col]=free[col]
        output_units.update([('Free '+s,'M') for s in self.primary_names])
        sorbed=pandas.DataFrame(self.output['immobile'],columns=self.primary_names,index=self.output['time']).add_prefix('Total Sorbed ')
        # for col in sorbed.columns:
        #     output_DF[col]=sorbed[col]
        output_units.update([('Total Sorbed '+s,'mol/m^3') for s in self.primary_names])
        mineral_VF=pandas.DataFrame(self.output['mineral_VF'],columns=self.mineral_names,index=self.output['time']).add_suffix(' VF')
        # for col in mineral_VF.columns:
        #     output_DF[col]=mineral_VF[col]
        output_units.update([(s+' VF','m^3 mnrl/m^3 bulk') for s in self.mineral_names])
        mineral_rate=pandas.DataFrame(self.output['mineral_rate'],columns=self.mineral_names,index=self.output['time']).add_suffix(' Rate')
        # for col in mineral_rate.columns:
        #     output_DF[col]=mineral_rate[col]
        output_units.update([(s+' Rate','mol/m^3/sec') for s in self.mineral_names])
        secondary=pandas.DataFrame(self.output['aq_complex'],columns=self.secondary_names,index=self.output['time'])
        # for col in secondary.columns:
        #     output_DF[col]=secondary[col]
        output_units.update([(s,'M') for s in self.secondary_names])

        output_DF=pandas.concat([output_DF,total,free,sorbed,mineral_VF,mineral_rate,secondary],axis=1)
        
        output_DF['CEC H+']=pandas.DataFrame(self.output['CEC H+'],index=self.output['time'])
        output_units['CEC H+']='mol/m^3'
        
        self.output_DF=output_DF.reset_index(drop=True).set_index(output_DF.index/(24*3600))
        self.output_units=output_units
        
        
    def run_onestep(self,chem,data,dt,status,min_dt=0.1,num_cuts=0,diffquo={},bc=None,flux_tol=0.15,truncate_concentration=0.0,rateconstants={},min_cuts=0,hands_off=False):
        self.copy_to_alquimia(data)
        converged=False
        porosity=data.state.porosity
        
        max_cuts=num_cuts
        actual_dt=dt/2**num_cuts
        
        if not hands_off:
            for num,reactname in enumerate(get_alquimiavector(data.meta_data.aqueous_kinetic_names)):
                data.properties.aqueous_kinetic_rate_cnst.data[num]=rateconstants[reactname]
        
            names=get_alquimiavector(data.meta_data.primary_names)
            for spec in names:
                pos=names.index(spec)
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
                print('Mobile: ')
                print(self.total_mobile)
                print('Free:')
                print(self.free_mobile)
                print('Immobile:')
                print(self.total_immobile)
                print('Mineral VF:')
                print(self.mineral_volume_fraction)
                if cut_for_flux:
                    raise RuntimeError('Pflotran failed to converge (because of boundary fluxes) after %d cuts to dt = %1.2f s'%(num_cuts,actual_dt))
                else:
                    raise RuntimeError('Pflotran failed to converge after %d cuts to dt = %1.2f s'%(num_cuts,actual_dt))
            
            # data will be reset to layer contents at beginning of next call
            # Run it twice, because we cut the time step in half
            ncuts=self.run_onestep(chem,data,dt,status,min_dt,num_cuts=num_cuts+1,diffquo=diffquo,bc=bc,truncate_concentration=truncate_concentration,rateconstants=rateconstants,hands_off=hands_off)
            # If this completes, it means that run_onestep was successful
            self.copy_from_alquimia(data)
            if ncuts>max_cuts:
                max_cuts=ncuts
                
            # This starts from ncuts so it doesn't have to try all the ones that failed again
            for n in range(2**(ncuts-(num_cuts+1))):
                ncuts2=self.run_onestep(chem,data,dt,status,min_dt,num_cuts=ncuts,diffquo=diffquo,bc=bc,truncate_concentration=truncate_concentration,rateconstants=rateconstants,hands_off=hands_off)
                self.copy_from_alquimia(data)
                if ncuts2>max_cuts:
                    max_cuts=ncuts2

            return max_cuts
        




    
def plot_result(result,SOM_ax=None,pH_ax=None,Fe_ax=None,gasflux_ax=None,porewater_ax=None,do_legend=False,gdrywt=False,BD=None,SOC_pct=None,cellulose_SOC_frac=1.0):
    if gdrywt:
        if SOC_pct is None and BD is None:
            raise TypeError('SOC_pct or BD must be a number if gdrywt is True')
        SOC_mol_m3=result['Total Sorbed cellulose'].iloc[0]/cellulose_SOC_frac # mol SOC/m3
        SOC_gC_cm3=SOC_mol_m3*12/100**3
        SOC_gC_gdwt=SOC_pct/100
        cm3_to_dwt=SOC_gC_gdwt/SOC_gC_cm3 # Conversion from /cm3 to /gdwt. This is 1/ bulk density in g/cm3
        if BD is not None:
            cm3_to_dwt=1.0/BD
        print('cm3_to_dwt = %1.3f'%cm3_to_dwt)
    if SOM_ax is not None:
        if gdrywt:
            # Plot in %, i.e. gC/gdrywt *100.
            l=SOM_ax.plot(result['Total Sorbed cellulose']/SOC_mol_m3*SOC_pct+SOC_pct*(1-cellulose_SOC_frac),label='SOM')[0]
            SOM_ax.set_ylabel('Concentration\n(SOC %)')
        else:
            l=SOM_ax.plot(result['Total Sorbed cellulose']*1e-3,label='SOM')[0]
            SOM_ax.set_ylabel('Concentration\n(mmol C/cm$^{-3}$)')

        SOM_ax.set_title('SOM remaining')
        SOM_ax.set_xlabel('Time (days)')

    if pH_ax is not None:
        pH_ax.plot(-numpy.log10(result['Free H+']))

        pH_ax.set_title('pH')
        pH_ax.set_ylabel('pH')
        pH_ax.set_xlabel('Time (days)')
        
    if Fe_ax is not None:
        molar_volume=34.3600 # From database. cm3/mol
        molar_weight = 106.8690
        if gdrywt:
            # Assume we can use SOC % to convert from volume to dry weight
            l=Fe_ax.plot(result['Fe(OH)3 VF']/molar_volume*1e6*cm3_to_dwt   ,label='Fe(OH)3')[0]
            Fe_ax.set_ylabel('Concentration\n($\mu$mol/g dwt)')
            l=Fe_ax.plot(result['Total Fe+++']*result['Porosity']*1e3*cm3_to_dwt   ,label='Fe+++',ls='--')[0]
            
            l=Fe_ax.plot(result['Total Fe++']*result['Porosity']*1e3*cm3_to_dwt ,ls=':'  ,label='Fe++')[0]
        else:
            l=Fe_ax.plot(result['Fe(OH)3 VF']/molar_volume*1e6   ,label='Fe(OH)3')[0]
            Fe_ax.set_ylabel('Concentration\n($\mu$mol cm$^{-3}$)')
        
            # M/L to umol/cm3: 1e6/1e3=1e3
            l=Fe_ax.plot(result['Total Fe+++']*result['Porosity']*1e3   ,label='Fe+++',ls='--')[0]
            
            l=Fe_ax.plot(result['Total Fe++']*result['Porosity']*1e3 ,ls=':'  ,label='Fe++')[0]
        
        Fe_ax.set_title('Fe species')
        
        Fe_ax.set_xlabel('Time (days)')
        if do_legend:
            Fe_ax.legend(fontsize='small')
    
    if gasflux_ax is not None:
        gasflux_ax.set_yscale('log')
        if gdrywt:
            l=gasflux_ax.plot(result.index.values[:-1],numpy.diff(result['Total CH4(aq)']*result['Porosity'])/numpy.diff(result.index.values)*1e3*cm3_to_dwt,label='CH4')[0]
            
            l=gasflux_ax.plot(result.index.values[:-1],numpy.diff(result['Total Tracer']*result['Porosity'])/numpy.diff(result.index.values)*1e3*cm3_to_dwt,label='CO2',ls='--',c='C5')[0]
            gasflux_ax.set_ylabel('Flux rate\n($\mu$mol g dwt$^{-1}$ day$^{-1}$)')
        else:
            l=gasflux_ax.plot(result.index.values[:-1],numpy.diff(result['Total CH4(aq)']*result['Porosity'])/numpy.diff(result.index.values)*1e3,label='CH4')[0]
            
            l=gasflux_ax.plot(result.index.values[:-1],numpy.diff(result['Total Tracer']*result['Porosity'])/numpy.diff(result.index.values)*1e3,label='CO2',ls='--',c='C5')[0]
            gasflux_ax.set_ylabel('Flux rate\n($\mu$mol cm$^{-3}$ day$^{-1}$)')

        gasflux_ax.set_title('Gas fluxes')
        
        gasflux_ax.set_xlabel('Time (days)')
        if do_legend:
            gasflux_ax.legend(fontsize='small')
        
    if porewater_ax is not None:
        porewater_ax.set_yscale('log')
        porewater_ax.plot(result['Total DOM1'],label='DOM')
        porewater_ax.plot(result['Total Acetate-'],label='Acetate',c='C3')
        porewater_ax.plot(result['Total O2(aq)'],'--',label='O2',c='C4')
        porewater_ax.plot(result['Total Fe+++'],'--',label='Fe+++',c='C1')
        porewater_ax.plot(result['Free Fe+++'],':',label='Fe+++',c='C1')
        porewater_ax.plot(result['Total Fe++'],':',label='Fe++',c='C2')
        
        porewater_ax.set_title('Porewater concentrations')
        porewater_ax.set_ylabel('Concentration (M)')
        porewater_ax.set_xlabel('Time (days)')
        
        if do_legend:
            porewater_ax.legend(fontsize='small')
    
    # tight_layout()

if __name__ == '__main__':
    
    
    
    import decomp_network
    pools = [
    decomp_network.decomp_pool(name='cellulose',CN=50,constraints={'initial':1e2},kind='immobile'),
    # decomp_network.decomp_pool(name='HRimm',constraints={'initial':1e-20},kind='immobile'),

    decomp_network.decomp_pool(name='DOM1',CN=50,constraints={'initial':1e-15},kind='primary'),
    decomp_network.decomp_pool(name='H+',kind='primary',constraints={'initial':'4.0 P'}),
    decomp_network.decomp_pool(name='O2(aq)',kind='primary',constraints={'initial':1e-12}),
    decomp_network.decomp_pool(name='HCO3-',kind='primary',constraints={'initial':'400e-6 G CO2(g)'}),
    decomp_network.decomp_pool(name='Fe+++',kind='primary',constraints={'initial':'.37e-3 M Fe(OH)3'}),
    decomp_network.decomp_pool(name='Fe++',kind='primary',constraints={'initial':'0.37e-30'}),
    decomp_network.decomp_pool(name='NH4+',kind='primary',constraints={'initial':1e-15}), # SOMDecomp sandbox requires this
    decomp_network.decomp_pool(name='Tracer',kind='primary',constraints={'initial':1e-15}), # Just to accumulate CO2 loss
    decomp_network.decomp_pool(name='CH4(aq)',kind='primary',constraints={'initial':1e-15}),
    decomp_network.decomp_pool(name='Acetate-',kind='primary',constraints={'initial':1e-15}),

    decomp_network.decomp_pool(name='CO2(g)',kind='gas'),
    decomp_network.decomp_pool(name='O2(g)',kind='gas'),

    decomp_network.decomp_pool(name='CO2(aq)',kind='secondary'),
    decomp_network.decomp_pool(name='OH-',kind='secondary'),
    decomp_network.decomp_pool(name='FeCO3+',kind='secondary'),
    decomp_network.decomp_pool(name='Fe(OH)4-',kind='secondary'),
    decomp_network.decomp_pool(name='Acetic_acid(aq)',kind='secondary'),
    decomp_network.decomp_pool(name='FeCH3COO+',kind='secondary'),
    decomp_network.decomp_pool(name='FeIIIDOM1(aq)',kind='secondary'),
    decomp_network.decomp_pool(name='FeIIIAcetate(aq)',kind='secondary'),
    # decomp_network.decomp_pool(name='FeCO3(aq)',kind='secondary'),

    decomp_network.decomp_pool(name='Fe(OH)3',rate='1.d-5 mol/m^2-sec',constraints={'initial':'1.75d-5  1.d2 m^2/m^3'},kind='mineral'),
    # decomp_network.decomp_pool(name='Goethite',rate='1.d-5 mol/m^2-sec',constraints={'initial':'1.75d-2  1.d1 m^2/m^3'},kind='mineral'),
    # decomp_network.decomp_pool(name='Fe',rate='1.d-7 mol/m^2-sec',constraints={'initial':'1.0e-6  1. m^2/m^3'},kind='mineral'),
    decomp_network.decomp_pool(name='Fe(OH)2',rate='1.d-7 mol/m^2-sec',constraints={'initial':'0.0e-20  1.e-10 m^2/m^3'},kind='mineral'),
    decomp_network.decomp_pool(name='Rock(s)',rate='0.0 mol/m^2-sec',constraints={'initial':'0.5  5.0e3 m^2/m^3'},kind='mineral'),

    decomp_network.decomp_pool(name='>Carboxylate-',kind='surf_complex',mineral='Rock(s)',site_density=10.0,complexes=['>Carboxylic_acid']),
    ]
    

            
            
    pools_atmoO2=pools.copy()
    for n,p in enumerate(pools_atmoO2):
        if p['name']=='O2(aq)':
            pools_atmoO2[n]=pools_atmoO2[n].copy()
            pools_atmoO2[n].update(constraints={'initial':'0.2 G O2(g)'})
            
    
    pools_lowFe=pools.copy()
    for n,p in enumerate(pools_lowFe):
        if p['name']=='Fe(OH)3':
            pools_lowFe[n]=pools_lowFe[n].copy()
            pools_lowFe[n].update(constraints={'initial':'0.0  1.e2 m^2/m^3'})
            
    
    # Using names assigned in alquimia. Not the best solution...
    rateconstants={
           '1.0 Fe++ + 0.25 O2(aq) + 1.0 H+ <-> 1.0 Fe+++ + 0.5 H2O'                                             : 2.0e0*1.0,
           "1.0e+00 DOM1  -> 3.3e-01 Acetate-  + 3.3e-01 HCO3-  + 2.0e+00 H+  + 3.3e-01 Tracer"                  : 1.0e-9*1.0,
           "1.0e+00 DOM1  + 1.0e+00 O2(aq)  -> 1.0e+00 HCO3-  + 1.0e+00 H+  + 1.0e+00 Tracer"                    : 1.0e-8*1.0,
           "1.0e+00 Acetate-  + 2.0e+00 O2(aq)  -> 2.0e+00 HCO3-  + 2.0e+00 H+  + 2.0e+00 Tracer"                : 1.0e-8*1.0,
           "1.0e+00 Acetate-  + 8.0e+00 Fe+++  -> 2.0e+00 HCO3-  + 8.0e+00 Fe++  + 9.0e+00 H+  + 2.0e+00 Tracer" : 5.0e-10*1.0,
           "1.0e+00 Acetate-  -> 1.0e+00 CH4(aq)  + 1.0e+00 HCO3-  + 1.0e+00 Tracer"                             : 1.0e-11*1.0,
           "cellulose decay to DOM1 (SOMDEC sandbox)"                                                                    : 1.0/(365*24*3600),
           "cellulose decay to CO2 (SOMDEC sandbox)"                                                                    : 1.0/(365*24*3600)
    }
    
    
    from numpy import zeros,diff
    dq=0.01 # Diffusion coefficient when aerobic
    O2_const=zeros(365*24)+dq

    O2_initial=zeros(365*24)
    O2_initial[:100*24]=dq
    result_highO2,output_units=run_simulation('fermentation.in',365,3600,initcond=pools_atmoO2,bc=pools_atmoO2,diffquo={'O2(aq)':O2_const},hands_off=True,rateconstants=rateconstants)
    result,output_units=run_simulation('fermentation.in',365,3600,initcond=pools,hands_off=False,rateconstants=rateconstants,bc=pools_atmoO2,diffquo={'O2(aq)':O2_initial})    
    result_lowFe,output_units=run_simulation('fermentation.in',365,3600,initcond=pools_lowFe,hands_off=False,rateconstants=rateconstants,bc=pools_atmoO2,diffquo={'O2(aq)':O2_initial})
    
    O2_periodic=zeros(365*24)
    O2_periodic[24*180:]=dq
    result_periodicO2,output_units=run_simulation('fermentation.in',365*3,3600,initcond=pools_atmoO2,bc=pools_atmoO2,diffquo={'O2(aq)':O2_periodic},hands_off=False,rateconstants=rateconstants)

    # err = lib.PetscFinalize()
    
    import decomp_network
    import plot_pf_output
    from pylab import *
    # result,units=plot_pf_output.convert_units(output,output_units,'M')
    # result_lowFe,units=plot_pf_output.convert_units(output_lowFe,output_units,'M')
    # result_highO2,units=plot_pf_output.convert_units(output_highO2,output_units,'M')
    # result_periodicO2,units=plot_pf_output.convert_units(output_periodicO2,output_units,'M')
    # figure('Network diagram',figsize=(11.8,4.8));clf()
    # ax=subplot(121)
    # decomp_network.draw_network(fermentation_network_Fe_CH4,omit=['secondary','surf_complex','NH4+','Rock(s)'],arrowstyle='-|>')
    # title('Decomposition network diagram (without complexes)')
    # 
    # ax=subplot(122)
    # decomp_network.draw_network(fermentation_network_Fe_CH4,omit=['NH4+','Rock(s)'],arrowstyle='-|>')
    # title('Decomposition network diagram (with aqueous complexes)')
    
        

    
    
    colors={'Anaerobic':'C1','Periodic':'C2','Low Fe':'C3','Aerobic':'C0'}

    fig=figure('Results');clf()
    fig,axes=subplots(5,4,num='Results')
    plot_result(result_highO2,SOM_ax=axes[0,0],pH_ax=axes[1,0],Fe_ax=axes[2,0],porewater_ax=axes[3,0],gasflux_ax=axes[4,0],do_legend=True)
    plot_result(result,SOM_ax=axes[0,1],pH_ax=axes[1,1],Fe_ax=axes[2,1],porewater_ax=axes[3,1],gasflux_ax=axes[4,1])
    plot_result(result_lowFe,SOM_ax=axes[0,2],pH_ax=axes[1,2],Fe_ax=axes[2,2],porewater_ax=axes[3,2],gasflux_ax=axes[4,2])
    plot_result(result_periodicO2,SOM_ax=axes[0,3],pH_ax=axes[1,3],Fe_ax=axes[2,3],porewater_ax=axes[3,3],gasflux_ax=axes[4,3])
    
    axes[0,0].set_title('Aerobic:\n'+axes[0,0].get_title())
    axes[0,1].set_title('Anaerobic:\n'+axes[0,1].get_title())
    axes[0,2].set_title('Low Fe:\n'+axes[0,2].get_title())
    axes[0,3].set_title('Periodically aerobic:\n'+axes[0,3].get_title())
    
    for ax in axes[1,:]:
        ax.set_ylim(3,5.8)
    for ax in axes[0,:]:
        ax.set_ylim(0,0.11)
    for ax in axes[4,:]:
        ax.set_ylim(bottom=1e-11,top=4e-3)
    for ax in axes[3,:]:
        ax.set_ylim(bottom=1e-11,top=0.4)
    for ax in axes[2,:]:
        ax.set_ylim(-0.01,0.53)
    
    tight_layout()
    
    
    # Figures for talk
    fig=figure('Just oxygen',figsize=(6,7.4));clf()
    fig,axes=subplots(3,1,num='Just oxygen')
    plot_result(result_lowFe,SOM_ax=axes[0])
    axes[0].set_ylim(bottom=-0.01)
    
    axes[1].plot(result_lowFe.index.values[:-1],diff(result_lowFe['Total Tracer']*result_lowFe['Porosity'])/diff(result_lowFe.index.values)*1e3,label='CO2',ls='--')
    # axes[2].plot(result_lowFe.index.values[:-1],diff(result_lowFe['Total CH4(aq)']*result_lowFe['Porosity'])/diff(result_lowFe.index.values)*1e6,label='CH4')
    methane_simple=zeros(len(result_lowFe.index.values[:-1]))
    methane_simple[24*100:]=0.02
    axes[2].plot(result_lowFe.index.values[:-1],methane_simple)

    axes[1].set_title('CO$_2$ flux')
    axes[2].set_title('CH$_4$ flux')
    axes[1].set_ylabel('Flux rate\n($\mu$mol cm$^{-3}$ day$^{-1}$)')
    axes[2].set_ylabel('Flux rate\n(nmol cm$^{-3}$ day$^{-1}$)')
    axes[2].set_xlabel('Time (days)')
    axes[2].set_ylim(top=0.031)
    
    xmax=axes[0].get_xlim()[1]
    for ax in axes:
        ax.axvspan(100,xmax,color='b',alpha=0.1)
        ax.set_xlim(right=xmax)
        ax.text(90,ax.get_ylim()[1]*0.9,'Aerobic',ha='right')
        ax.text(110,ax.get_ylim()[1]*0.9,'Inundated',ha='left')
    
    fig.tight_layout()
    
    # With Fe reduction
    fig=figure('With Fe reduction',figsize=(6,8.4));clf()
    fig,axes=subplots(4,1,num='With Fe reduction')
    plot_result(result,SOM_ax=axes[0],Fe_ax=axes[2],gasflux_ax=axes[1],pH_ax=axes[3])
    axes[0].set_ylim(bottom=-0.01)
    axes[2].legend(loc='right',labels=['Iron oxide (solid)','Fe$^{3+}$ (dissolved)','Fe$^{2+}$ (dissolved)'])

    # Turn off acetate to make it cleaner
    axes[1].set_ylim(bottom=0.9e-11)
    axes[0].text(90,0.0,'Aerobic',ha='right')
    axes[0].text(110,0.0,'Inundated',ha='left')
    axes[1].legend(labels=['CH$_4$','CO$_2$'])

    
    xmax=axes[0].get_xlim()[1]
    for ax in axes:
        ax.axvspan(100,xmax,color='b',alpha=0.1)
        ax.set_xlim(right=xmax)

    
    fig.tight_layout()
    
    # Periodic inundation
    fig=figure('Periodic inundation',figsize=(6,8.4));clf()
    fig,axes=subplots(4,1,num='Periodic inundation')
    plot_result(result_periodicO2,SOM_ax=axes[0],Fe_ax=axes[2],gasflux_ax=axes[1],pH_ax=axes[3])
    axes[0].set_ylim(bottom=-0.01)
    axes[2].legend(labels=['Fe(OH)$_3$','Fe$^{3+}$','Fe$^{2+}$'],loc=(0.83,0.2),fontsize='small')

    # Turn off acetate to make it cleaner
    axes[1].set_ylim(bottom=0.9e-11)
    axes[1].legend(labels=['CH$_4$','CO$_2$'])

    
    xmax=axes[0].get_xlim()[1]
    for ax in axes:
        for num in range(3):
            ax.axvspan(num*365,num*365+nonzero(diff(O2_periodic))[0]/24,color='b',alpha=0.1)
        ax.set_xlim(right=xmax)

    
    fig.tight_layout()
    

    # 
    # ax=subplot(234)
    # ax.set_yscale('log')
    # for pool in ['Free DOM1','Free Acetate-','Free Fe+++','Free Fe++','Free O2(aq)']:
    #     l=plot(result[pool],label=pool)[0]
    #     plot(result_lowFe[pool],ls='--',c=l.get_color())
    #     plot(result_periodicO2[pool],ls=':',c=l.get_color())
    #     plot(result_highO2[pool],ls='-.',c=l.get_color())
    # 
    # ylim(1e-15,1e-1)
    # 
    # title('Log concentrations')
    # ylabel('Concentration (M)')
    # xlabel('Time (days)')
    # legend(fontsize='small')
    # 
    # ax=subplot(235)
    # ax.set_yscale('log')
    # for pool in ['Total DOM1','Total Acetate-','Total Fe+++','Total Fe++','Total O2(aq)']:
    #     l=plot(result[pool],label=pool)[0]
    #     plot(result_lowFe[pool],ls='--',c=l.get_color())
    #     plot(result_periodicO2[pool],ls=':',c=l.get_color())
    #     plot(result_highO2[pool],ls='-.',c=l.get_color())
    # 
    # ylim(1e-15,1e-1)
    # 
    # title('Log concentrations')
    # ylabel('Concentration (M)')
    # xlabel('Time (days)')
    # legend(fontsize='small')
    # 
    # 
    # tight_layout()
    show()


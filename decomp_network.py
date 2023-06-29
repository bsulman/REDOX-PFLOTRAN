import networkx as nx
import os

nonCN_immobile_pools=['HRimm','Nmin','Nimp','Nimm','NGASmin','Root_biomass','Sorption_capacity','Plant_NH4_demand','Plant_NO3_demand']

class decomp_pool(dict):
    pass

def pools_list_to_dict(pools):
    if isinstance(pools,dict):
        return pools
    else:
        return {p['name']:p for p in pools}
def pools_dict_to_list(pools):
    if isinstance(pools,dict):
        return list(pools.values())
    else:
        return pools

def get_stoich_from_name(name,reactions,precision=2):
    if isinstance(reactions,list):
        react=pools_list_to_dict(reactions)[name]
    elif isinstance(reactions,dict):
        react=reactions[name]
    elif isinstance(reactions,decomp_network):
        found_react=False
        for r in reactions.edges:
            if reactions.edges[r]['name']==name:
                react=reactions.edges[r]['reaction']
                found_react=True
                break
        if not found_react:
            raise ValueError('Reaction "%s" not found in reaction network'%name)
    else:
        raise ValueError('reactions not in implemented format')
    
    if react['reactiontype']=='MICROBIAL':
        stoich=PF_microbial_reaction_writer(reactions,precision=precision).write_reaction_stoich(react['reactant_pools'],react['product_pools'])
    elif react['reactiontype']=='GENERAL':
        stoich=PF_general_reaction_writer(reactions,precision=precision).write_reaction_stoich(react['reactant_pools'],react['product_pools'])
    elif react['reactiontype']=='SOMDECOMP':
        prods=list(react['product_pools'].keys())
        if len(prods)==1 and (prods[0]=='HCO3-' or 'CO2' in prods[0]):
            stoich='%s decay to %s (SOMDEC sandbox)'%(list(react['reactant_pools'])[0],'CO2')
        else:
            stoich='%s decay to %s (SOMDEC sandbox)'%(list(react['reactant_pools'])[0],list(react['product_pools'])[0])
    else:
        raise ValueError("Reaction type '%s' not implemented"%react['reactiontype'])
    return stoich.strip()

def change_constraint(pools,poolname,newval,constraint='initial',inplace=False):
    'Change a constraint value for a pool in a list of pools. Returns a copy unless inplace is True'
    if not inplace:
        import copy
        out=copy.deepcopy(pools)
    else:
        out=pools
        
    for n,p in enumerate(out):
        if p['name']==poolname:
            out[n].update(constraints={constraint:newval})
            
    if not inplace:
        return out
        
def change_constraints(pools,changes,constraint='initial',inplace=False):
    'Change a constraint value for a pool in a list of pools. Returns a copy unless inplace is True'
    if not inplace:
        import copy
        out=copy.deepcopy(pools)
    else:
        out=pools
        
    for n,p in enumerate(out):
        if p['name'] in changes:
            out[n].update(constraints={constraint:changes[p['name']]})
            
    if not inplace:
        return out

def change_site_density(pools,sitename,newdensity,inplace=False):
    'Change a complexation site density in a list of pools. Returns a copy unless inplace is True'
    if not inplace:
        import copy
        out=copy.deepcopy(pools)
    else:
        out=pools
        
    for n,p in enumerate(out):
        if p['name']==sitename:
            out[n].update(site_density=newdensity)
            
    if not inplace:
        return out

class inhibition(dict):
    pass
    
class monod(dict):
    pass
    
class sorption_isotherm(dict):
    def __init__(self,name,sorbed_species,mineral,k,langmuir_b,sorbed_name):
        self['name']=name
        self['reactant_pools']={sorbed_species:1.0,mineral:1.0}
        self['product_pools']={sorbed_name:1.0}
        self['reactiontype']='sorption_isotherm'
        self['mineral']=mineral
        self['sorbed_species']=sorbed_species
        self['k']=k
        self['langmuir_b']=langmuir_b
        
class ion_exchange(dict):
    def __init__(self,name,cations,mineral,CEC):
        self['name']=name
        self['cations']=cations
        self['mineral']=mineral
        self['reactant_pools']=cations.copy()
        if mineral is not None:
            self['reactant_pools'][mineral]=1.0
            self['product_pools']={mineral:1.0}
        else:
            self['product_pools']=cations.copy()
        self['CEC']=CEC
        self['reactiontype']='ion_exchange'
    
class reaction(dict):
    def __init__(self,name,rate_constant,reactiontype,reactant_pools=None,product_pools=None,stoich=None,rate_units='y',**kwargs):
        if stoich is not None:
            if reactant_pools is not None or product_pools is not None:
                raise ValueError('Specify reactant and product pools, or stoich, not both')
            else:
                reactant_pools,product_pools=self.read_reaction_stoich(stoich)
        if not isinstance(reactant_pools,dict):
            raise TypeError('stoich must be a dictionary of pool names and stoichiometries')
        if not isinstance(product_pools,dict):
            raise TypeError('product_pools must be a dictionary of pool names and stoichiometries')
        self['name']=name
        self['reactant_pools']=reactant_pools
        assert len(reactant_pools)>0, 'There must be at least one reactant (%s)'%name
        self['product_pools']=product_pools
        assert len(product_pools)>0, 'There must be at least one product (%s)'%name
        self['rate_units']=rate_units
        if reactiontype not in ['MICROBIAL','SOMDECOMP','GENERAL']:
            raise ValueError('Only MICROBIAL, GENERAL, and SOMDECOMP reactions are currently implemented')
        self['reactiontype']=reactiontype
        self['rate_constant']=rate_constant
        self.update(kwargs)
    def read_reaction_stoich(self,stoich):
        #REACTION 1.0e+00 DOM1  -> 3.3e-01 Acetate-  + 3.3e-01 HCO3-  + 2.0e+00 H+  + 3.3e-01 Tracer 
        if '->' not in stoich:
            raise ValueError('stoich must be a string in the format "reactants -> products"')
        if 'REACTION' in stoich:
            stoich=stoich[stoich.find('REACTION')+len('REACTION'):]
        stoich=stoich.strip()
        if '<->' in stoich:
            reactant_side,product_side=stoich.split('<->')
        else:
            reactant_side,product_side=stoich.split('->')
        r=reactant_side.split(' + ')
        p=product_side.split(' + ')
        reactants_dict={}
        products_dict={}
        for comp in r:
            if len(comp.split())==1:
                num,chemical=1,comp.strip()
            else:
                num,chemical=comp.split()
            reactants_dict[chemical]=float(num)
        for comp in p:
            if len(comp.split())==1:
                num,chemical=1,comp.strip()
            else:
                num,chemical=comp.split()
            products_dict[chemical]=float(num)
        return reactants_dict,products_dict


# Class for writing network out into Pflotran reaction sandbox format
class PF_writer:
    def __init__(self,network,indent_spaces=2,base_indent=0,precision=2):
        self.level=[]
        self.output=''
        self.indent_spaces=indent_spaces
        self.base_indent=base_indent
        self.network=network
        self.precision=precision
        return 
        
    def reset(self):
        self.__init__(self.network)

    def current_indent(self):
        return self.base_indent+len(self.level)
    def increase_level(self,levelname):
        self.add_line(levelname)
        self.level.append(levelname)
    def decrease_level(self):
        self.level.pop()
        self.output = self.output + ' '*(self.base_indent+len(self.level)*self.indent_spaces) + '/\n'
    def add_line(self,text):
        self.output = self.output + ' '*(self.base_indent+len(self.level)*self.indent_spaces) + text + '\n'
    def write(self,*args,**kwargs):
        raise NotImplementedError('Must use subclass of PF_reaction_writer to write out reactions')
        
class PF_network_writer(PF_writer):
    def write_all_reactions(self,base_indent=0,indent_spaces=2,CO2name='HCO3-',SOMdecomp_Q10=None,moisturefunc=None):
        self.indent_spaces=indent_spaces
        self.base_indent=base_indent

        # Microbial reactions
        already_done = []
        for (ins,outs,reaction) in self.network.edges(data='reaction'):
            if reaction['name'] in already_done:
                continue
            already_done.append(reaction['name'])
            if reaction['reactiontype'] == 'MICROBIAL':
                self.output = self.output + PF_microbial_reaction_writer(reaction,base_indent=self.current_indent(),precision=self.precision).write_reaction()
            elif reaction['reactiontype'] == 'GENERAL':
                self.output = self.output + PF_general_reaction_writer(reaction,base_indent=self.current_indent(),precision=self.precision).write_reaction()
        
        # List all sandbox reactions
        sandbox_reacts=self.network.edge_subgraph([(ins,outs,ks) for ins,outs,ks,r in self.network.edges(data='reaction',keys=True) if r['reactiontype']=='SOMDECOMP'])
        # First write out all pools and their CN ratios
        if len(sandbox_reacts)>0:
            self.increase_level('REACTION_SANDBOX')
            self.increase_level('SOMDECOMP')
            self.increase_level('POOLS')
            for pool in sandbox_reacts.nodes:
                if 'CO2' in pool or 'HCO3-' in pool:
                    continue
                if 'CN' in sandbox_reacts.nodes[pool]:
                    self.add_line(pool.ljust(20)+'{CN:1.{prec}f}'.format(CN=sandbox_reacts.nodes[pool]['CN'],prec=self.precision))
                elif pool.endswith('C') and pool[:-1]+'N' in self.network.nodes:
                    self.add_line(pool[:-1].ljust(20)+'# Variable C:N pool')
                else:
                    self.add_line(pool.ljust(20)+'# Variable C:N pool')
            self.decrease_level()

            # Note: Universal abiotic factors must come earlier in input deck than reactions for them to propagate properly
            if SOMdecomp_Q10 is not None or moisturefunc is not None:
                self.increase_level('Abiotic_Factors')
                if SOMdecomp_Q10 is not None:
                    self.increase_level('TEMPERATURE_RESPONSE_FUNCTION')
                    self.add_line(f'Q10 {SOMdecomp_Q10:1.{self.precision}f}')
                    self.decrease_level()

                if moisturefunc is not None:
                    self.increase_level('MOISTURE_RESPONSE_FUNCTION')
                    self.add_line(moisturefunc)
                    self.decrease_level()
                self.decrease_level()
        
            # Write out all the reactions
            already_done=[]
            for ins,outs,reaction in sandbox_reacts.edges(data='reaction'):
                if reaction['name'] in already_done:
                    continue
                already_done.append(reaction['name'])
                self.output = self.output + PF_sandbox_reaction_writer(reaction,base_indent=self.current_indent(),precision=self.precision).write_reaction()
                        

            self.add_line('CO2_SPECIES_NAME '+CO2name)
            if 'O2(aq)' in self.network.nodes:
                self.add_line('O2_SPECIES_NAME O2(aq)')

            
        for lev in range(len(self.level)):
            self.decrease_level()
        return self.output
            


    def write_into_input_deck(self,templatefile_name,outputfile_name,
            indent_spaces=2,length_days=None,obs_time_hrs=None,log_formulation=True,truncate_concentration=1e-80,database='./hanford.dat',chem_args={},verbose=True,**kwargs):
        base_indent=0
        with open(templatefile_name,'r') as templatefile:
            template_lines=templatefile.readlines()
            
        fmt='%%1.%de'%self.precision
    
        for line in template_lines:
            if 'CHEMISTRY' in line and not line.strip().startswith('#'):
                self.output = self.output + line
                # Start of chemistry block
                # Primary species
                self.increase_level('PRIMARY_SPECIES')
                self.add_line('#### NOTE: Beginning of auto-inserted primary species ####')
                for pool in self.network.nodes:
                    if self.network.nodes[pool]['kind']=='primary':
                        self.add_line( pool )
                self.add_line( '#### NOTE: End of auto-inserted primary species ####')
                self.decrease_level()
                
                # Assume that all primary species should be put into decoupled eq reactions just in case
                self.increase_level('DECOUPLED_EQUILIBRIUM_REACTIONS')
                self.add_line( '#### NOTE: Beginning of auto-inserted primary species ####')
                for pool in self.network.nodes:
                    if self.network.nodes[pool]['kind']=='primary':
                        self.add_line( pool )
                self.add_line( '#### NOTE: End of auto-inserted primary species ####')
                self.decrease_level()
                
                # Secondary species
                self.increase_level('SECONDARY_SPECIES')
                self.add_line('#### NOTE: Beginning of auto-inserted secondary species ####')
                for pool in self.network.nodes:
                    if self.network.nodes[pool]['kind']=='secondary':
                        self.add_line( pool )
                self.add_line( '#### NOTE: End of auto-inserted secondary species ####')
                self.decrease_level()
                
                # Minerals
                self.increase_level('MINERALS')
                self.add_line( '#### NOTE: Beginning of auto-inserted mineral species ####')
                for pool in self.network.nodes:
                    if self.network.nodes[pool]['kind']=='mineral':
                        self.add_line( pool)
                self.add_line( '#### NOTE: End of auto-inserted mineral species ####')
                self.decrease_level()
                
                # Mineral kinetics
                self.increase_level('MINERAL_KINETICS')
                self.add_line( '#### NOTE: Beginning of auto-inserted mineral kinetics ####')
                for pool in self.network.nodes:
                    if self.network.nodes[pool]['kind']=='mineral':
                        self.increase_level(pool)
                        self.add_line('RATE_CONSTANT  ' + self.network.nodes[pool]['rate'])
                        self.decrease_level()
                self.add_line( '#### NOTE: End of auto-inserted mineral kinetics ####')
                self.decrease_level()
                
                # Immobile species
                self.increase_level('IMMOBILE_SPECIES')
                self.add_line( '#### NOTE: Beginning of auto-inserted immobile species ####')
                for pool in self.network.nodes:
                    if self.network.nodes[pool]['kind']=='immobile':
                        if 'CN' in self.network.nodes[pool] or pool in nonCN_immobile_pools: # This should be improved
                            self.add_line( pool)
                        else:
                            self.add_line( pool + 'C')
                            self.add_line( pool + 'N')
                self.add_line( '#### NOTE: End of auto-inserted immobile species ####')
                self.decrease_level()
                
                # Gas species
                # self.increase_level('GAS_SPECIES')
                self.increase_level('PASSIVE_GAS_SPECIES')
                self.add_line( '#### NOTE: Beginning of auto-inserted gas species ####')
                for pool in self.network.nodes:
                    if self.network.nodes[pool]['kind']=='gas':
                        self.add_line( pool)
                self.add_line( '#### NOTE: End of auto-inserted gas species ####')
                self.decrease_level()
                
                # Sorption reactions
                self.increase_level('SORPTION')
                self.add_line( '#### NOTE: Beginning of auto-inserted sorption sites ####')
                for pool in self.network.nodes:
                    if self.network.nodes[pool]['kind']=='surf_complex':
                        self.increase_level('SURFACE_COMPLEXATION_RXN')
                        if 'rates' in self.network.nodes[pool]:
                            self.add_line('MULTIRATE_KINETIC')
                            if 'site_fractions' in self.network.nodes[pool]:
                                self.add_line('SITE_FRACTION '+' '.join([fmt%frac for frac in self.network.nodes[pool]['site_fractions']]))
                            else:
                                self.add_line('SITE_FRACTION '+' '.join([fmt%(1/len(self.network.nodes[pool]['rates'])) for rate in self.network.nodes[pool]['rates']]))
                            self.add_line('RATES '+' '.join([fmt%frac+' ' for frac in self.network.nodes[pool]['rates']]))
                        elif 'rate' in self.network.nodes[pool]:
                            self.add_line('KINETIC')
                        else:
                            self.add_line('EQUILIBRIUM')
                        self.add_line('MINERAL {mineral:s}'.format(mineral=self.network.nodes[pool]['mineral']))
                        self.add_line('SITE {sitename:s} {density:{fmt}}'.format(sitename=pool,density=self.network.nodes[pool]['site_density'],fmt=fmt[1:]))
                        self.increase_level('COMPLEXES')
                        for complex in self.network.nodes[pool]['complexes']:
                            self.add_line(complex)
                        self.decrease_level()
                        self.decrease_level()
                        
                already_done=[]
                for ins,outs,reaction in self.network.edges(data='reaction'):
                    if reaction['reactiontype']=='ion_exchange' and reaction['name'] not in already_done:
                        already_done.append(reaction['name'])
                        self.increase_level('ION_EXCHANGE_RXN')
                        if reaction['mineral'] is not None:
                            self.add_line('MINERAL '+reaction['mineral'])
                        self.add_line('CEC '+fmt%reaction['CEC'])
                        self.increase_level('CATIONS')
                        for cation in reaction['cations']:
                            if reaction['cations'][cation]==1.0:
                                self.add_line('{cation:s} {val:{fmt}} REFERENCE'.format(cation=cation,val=reaction['cations'][cation],fmt=fmt[1:]))
                            else:
                                self.add_line('{cation:s} {val:{fmt}}'.format(cation=cation,val=reaction['cations'][cation],fmt=fmt[1:]))
                        self.decrease_level()
                        self.decrease_level()
                
                isotherm_reacts=self.network.edge_subgraph([(ins,outs,ks) for ins,outs,ks,r in self.network.edges(data='reaction',keys=True) if r['reactiontype']=='sorption_isotherm'])
                if len(isotherm_reacts)>0:
                    self.increase_level('ISOTHERM_REACTIONS')
                    already_done = []
                    for (ins,outs,reaction) in isotherm_reacts.edges(data='reaction'):
                        if reaction['name'] in already_done:
                            continue
                        already_done.append(reaction['name'])
                        self.increase_level(reaction['sorbed_species'])
                        self.add_line('DISTRIBUTION_COEFFICIENT {k:{fmt}}'.format(k=reaction['k']),fmt=fmt[1:])
                        self.add_line('KD_MINERAL_NAME {mineral:s}'.format(mineral=reaction['mineral']))
                        self.add_line('LANGMUIR_B {b:{fmt}}'.format(b=reaction['langmuir_b']),fmt=fmt[1:])
                        self.decrease_level()
                    self.decrease_level()
                self.decrease_level()
                self.add_line( '#### NOTE: End of auto-inserted sorption sites ####')
                
                # Reactions
                self.add_line( '#### NOTE: Beginning of auto-inserted reactions ####')
                self.write_all_reactions(base_indent=base_indent+indent_spaces,indent_spaces=indent_spaces,**kwargs)
                self.add_line( '#### NOTE: End of auto-inserted reactions ####')
                
                if log_formulation:
                    self.add_line('LOG_FORMULATION')
                    
                if truncate_concentration is not None:
                    self.add_line('TRUNCATE_CONCENTRATION {conc:{fmt}}'.format(conc=truncate_concentration,fmt=fmt[1:]))
                    
                self.add_line('DATABASE %s'%os.path.abspath(database))
                
                for arg in chem_args:
                    if isinstance(chem_args[arg],float):
                        self.add_line('{arg:s} {val:{fmt}}'.format(arg=arg,val=chem_args[arg],fmt=fmt[1:]))
                    elif isinstance(chem_args[arg],int):
                        self.add_line('{arg:s} {val:d}'.format(arg=arg,val=chem_args[arg]))
                    elif isinstance(chem_args[arg],str):
                        self.add_line('{arg:s} {val:s}'.format(arg=arg,val=chem_args[arg]))
                    else:
                        raise TypeError('chem_args must be float, int, or str (arg: %s)'%arg)
            
            elif len(line.split())>=2 and line.split()[0]=='CONSTRAINT':
                constraintname=line.split()[1]
                self.output = self.output + line
                self.increase_level('IMMOBILE')
                self.add_line( '#### NOTE: Beginning of auto-inserted immobile species ####')
                for pool in self.network.nodes:
                    if not self.network.nodes[pool]['kind']=='immobile':
                        continue
                    constraint=self.network.nodes[pool]['constraints'].get(constraintname,1e-20)
                    if 'CN' in self.network.nodes[pool] or pool in nonCN_immobile_pools:
                        self.add_line( pool.ljust(20) + ' {const:{fmt}}'.format(const=constraint,fmt=fmt[1:]))
                    else:
                        if 'initCN' not in self.network.nodes[pool]:
                            raise ValueError('initCN must be provided for flexible CN pools: pool %s'%pool)
                        self.add_line( (pool+'C').ljust(20) + '{const:{fmt}}'.format(const=constraint,fmt=fmt[1:]))
                        self.add_line( (pool+'N').ljust(20) + '{const:{fmt}}'.format(const=constraint/self.network.nodes[pool]['initCN'],fmt=fmt[1:]))
                self.add_line( '#### NOTE: End of auto-inserted immobile species ####')
                self.decrease_level()
                
                self.increase_level('CONCENTRATIONS')
                self.add_line( '#### NOTE: Beginning of auto-inserted concentration constraints ####')
                for pool in self.network.nodes:
                    if self.network.nodes[pool]['kind']=='primary':
                        constraint=self.network.nodes[pool]['constraints'].get(constraintname,1e-20)
                        if isinstance(constraint,str):
                            self.add_line( (pool).ljust(20) + constraint)
                        else:
                            self.add_line( (pool).ljust(20) + '{const:{fmt}}'.format(const=constraint,fmt=fmt[1:]))
                            
                self.add_line( '#### NOTE: End of auto-inserted concentration constraints ####')
                self.decrease_level()
                
                self.increase_level('MINERALS')
                self.add_line( '#### NOTE: Beginning of auto-inserted mineral constraints ####')
                for pool in self.network.nodes:
                    if self.network.nodes[pool]['kind']=='mineral':
                        constraint=self.network.nodes[pool]['constraints'].get(constraintname,1e-20)
                        if isinstance(constraint,str):
                            self.add_line( (pool).ljust(20) + self.network.nodes[pool]['constraints'].get(constraintname,1e-20))
                        else:
                            self.add_line( (pool).ljust(20) + '{const:{fmt}} 1.0'.format(const=self.network.nodes[pool]['constraints'].get(constraintname,1e-20),fmt=fmt[1:]))
                            
                self.add_line( '#### NOTE: End of auto-inserted mineral constraints ####')
                self.decrease_level()
                
            elif 'FINAL_TIME' in line and length_days is not None:
                self.add_line( 'FINAL_TIME {ndays:1.{prec}e} d\n'.format(ndays=length_days,prec=self.precision))
            elif 'PERIODIC_OBSERVATION TIME' in line and obs_time_hrs is not None:
                self.add_line('PERIODIC_OBSERVATION TIME %1.1d h'%obs_time_hrs)
            elif len(line.strip())>0 and line.strip().split()[0]=='OUTPUT':
                self.output = self.output + line
                if not verbose:
                    self.add_line('SCREEN OFF')
            else:
                if not line.isspace():
                    base_indent=len(line)-len(line.lstrip())
                self.output = self.output + line
                    
        with open(outputfile_name,'w') as outputfile:
            outputfile.write(self.output)
        return

        
        
    def run_simulation(self,template_file,simulation_name,pflotran_exe,output_suffix='-obs-0.pft',print_output=False,length_days=None,obs_time_hrs=None,
                        log_formulation=True,truncate_concentration=1e-80,database='./hanford.dat',CO2name='HCO3-',**kwargs):
        inputdeck=simulation_name+'_generated.in'
        print('Setting up input deck in %s'%inputdeck)
        self.write_into_input_deck(template_file,inputdeck,length_days=length_days,obs_time_hrs=obs_time_hrs,
                                    log_formulation=log_formulation,database=database,truncate_concentration=truncate_concentration,CO2name=CO2name,**kwargs)
        import subprocess
        cmd='{pflotran_exe:s} -pflotranin {simname:s}_generated.in'.format(pflotran_exe=pflotran_exe,simname=simulation_name)
        print('Running cmd: %s'%cmd)
        status,output = subprocess.getstatusoutput(cmd)
        print('Simulation finished with status %d'%status)
        if print_output:
            print(output)
        if 'ERROR' in output or 'Stopping' in output:
            print(output)
            raise RuntimeError('Pflotran simulation failed')
        import plot_pf_output
        outputfile=simulation_name + '_generated' + output_suffix
        # Set up for more flexibility in output formats
        print('Reading output from file %s'%outputfile)
        output_data,units=plot_pf_output.read_tecfile(outputfile)
        return output_data,units
        
class PF_sandbox_reaction_writer(PF_writer):
    def write_reaction(self,base_indent=None,indent_spaces=None):
        fmt='1.%de'%self.precision
        reaction_data=self.network.copy()
        assert len(reaction_data['reactant_pools']) == 1,'Only one reactant pool is allowed for sandbox reactions'
        assert list(reaction_data['reactant_pools'].values())[0] == 1.0,'Stoichiometry of upstream pool in sandbox reaction must be 1.0'
        if base_indent is not None:
            self.base_indent=base_indent
        self.add_line('# {name:s}'.format(name=reaction_data.pop('name')))
        if 'comment' in reaction_data.keys():
            self.add_line('# {comment:s}'.format(comment=reaction_data.pop('comment')))
        self.increase_level('REACTION')
        upstream_pool=list(reaction_data.pop('reactant_pools').keys())[0]
        self.add_line('UPSTREAM_POOL'.ljust(20)+upstream_pool)
        downstream_pools = reaction_data.pop('product_pools')
        for pool in downstream_pools.keys():
            if 'CO2' not in pool and 'HCO3-' not in pool: # CO2 is included implicitly as the rest of CUE
                self.add_line('DOWNSTREAM_POOL'.ljust(20)+pool.ljust(20)+'{CUE:{fmt}}'.format(CUE=downstream_pools[pool],fmt=fmt))
        turnover_name=reaction_data.pop('turnover_name','TURNOVER_TIME')
        self.add_line(turnover_name.ljust(20)+'{time:{fmt}} {units:s}'.format(time=reaction_data.pop('rate_constant'),units=reaction_data.pop('rate_units'),fmt=fmt))
        for monod in reaction_data.pop('monod_terms',[]):
                self.increase_level('MONOD')
                self.add_line('SPECIES_NAME'.ljust(20)+monod['species'])
                self.add_line('HALF_SATURATION_CONSTANT'+' {const:{fmt}}'.format(const=monod['k'],fmt=fmt))
                if 'threshold' in monod:
                    self.add_line('THRESHOLD_CONCENTRATION {const:{fmt}}'.format(const=monod['threshold'],fmt=fmt))
                if monod.get('pool_normalized',False):
                    self.add_line('POOL_NORMALIZED')
                self.decrease_level()
        for inhib in reaction_data.pop('inhibition_terms',[]):
                self.increase_level('INHIBITION')
                self.add_line('SPECIES_NAME'.ljust(20)+inhib['species'])
                self.add_line('TYPE {inhtype:s}'.format(inhtype=inhib['type']))
                self.add_line('INHIBITION_CONSTANT'+' {const:{fmt}}'.format(const=inhib['k'],fmt=fmt))
                self.decrease_level()
        # Write out the rest of the reaction attributes
        reaction_data.pop('reactiontype',None)
        for param in reaction_data.keys():
            val=reaction_data[param]
            if isinstance(val,str):
                self.add_line(param.ljust(20)+val)
            elif isinstance(val,float):
                self.add_line(param.ljust(20)+f' {val:{fmt}}')
            else:
                raise ValueError('Parameter {param:s} was not a str or float'.format(param=param))
        self.decrease_level()
    
        return self.output
        

class PF_general_reaction_writer(PF_writer):
    def write_reaction_stoich(self,reactants,products):
        line=''
        first = True
        for chem in reactants:
            if first:
                first = False
            else:
                line = line + ' + '
            line = line + '{stoich:1.{prec}e} {chem:s} '.format(stoich=reactants[chem],chem=chem,prec=self.precision)
        line = line + ' <-> '
        first = True
        for chem in products:
            if first:
                first = False
            else:
                line = line + ' + '
            line = line + '{stoich:1.{prec}e} {chem:s} '.format(stoich=products[chem],chem=chem,prec=self.precision)
        return line
    def write_reaction(self,base_indent=None,indent_spaces=None):
        reaction_data=self.network.copy()
        
        if base_indent is not None:
            self.base_indent=base_indent
        self.add_line('# {name:s}'.format(name=reaction_data.pop('name')))
        if 'comment' in reaction_data.keys():
            self.add_line('# {comment:s}'.format(comment=reaction_data.pop('comment')))
        self.increase_level('GENERAL_REACTION')
        line = 'REACTION '
        reactants=reaction_data.pop('reactant_pools')
        products =reaction_data.pop('product_pools')
        line=line+self.write_reaction_stoich(reactants,products)
        self.add_line(line)

        self.add_line('FORWARD_RATE'.ljust(20)+'{time:1.{prec}e}'.format(time=reaction_data.pop('rate_constant'),prec=self.precision))
        self.add_line('BACKWARD_RATE'.ljust(20)+'{time:1.{prec}e}'.format(time=reaction_data.pop('backward_rate_constant',0.0),prec=self.precision))
        # Write out the rest of the reaction attributes
        # What to do with biomass? Maybe better not to use it
        self.decrease_level()
    
        return self.output
        
class PF_microbial_reaction_writer(PF_writer):
    def write_reaction_stoich(self,reactants,products,):
        line=''
        first = True
        for chem in reactants:
            if first:
                first = False
            else:
                line = line + ' + '
            line = line + '{stoich:1.{prec}e} {chem:s} '.format(stoich=reactants[chem],chem=chem,prec=self.precision)
        line = line + ' -> '
        first = True
        for chem in products:
            if first:
                first = False
            else:
                line = line + ' + '
            line = line + '{stoich:1.{prec}e} {chem:s} '.format(stoich=products[chem],chem=chem,prec=self.precision)
        return line
    def write_reaction(self,base_indent=None,indent_spaces=None,):
        reaction_data=self.network.copy()
        
        if base_indent is not None:
            self.base_indent=base_indent
        self.add_line('# {name:s}'.format(name=reaction_data.pop('name')))
        if 'comment' in reaction_data.keys():
            self.add_line('# {comment:s}'.format(comment=reaction_data.pop('comment')))
        self.increase_level('MICROBIAL_REACTION')
        line = 'REACTION '
        reactants=reaction_data.pop('reactant_pools')
        products =reaction_data.pop('product_pools')
        line=line+self.write_reaction_stoich(reactants,products)
        self.add_line(line)

        self.add_line('RATE_CONSTANT'.ljust(20)+'{time:1.{prec}e}'.format(time=reaction_data.pop('rate_constant'),prec=self.precision))
        for monod in reaction_data.pop('monod_terms',[]):
                self.increase_level('MONOD')
                self.add_line('SPECIES_NAME'.ljust(20)+monod['species'])
                self.add_line('HALF_SATURATION_CONSTANT'+' {const:1.{prec}e}'.format(const=monod['k'],prec=self.precision))
                if 'threshold' in monod:
                    self.add_line('THRESHOLD_CONCENTRATION {const:1.{prec}e}'.format(const=monod['threshold'],prec=self.precision))
                self.decrease_level()
        for inhib in reaction_data.pop('inhibition_terms',[]):
                self.increase_level('INHIBITION')
                self.add_line('SPECIES_NAME'.ljust(20)+inhib['species'])
                self.add_line('TYPE {inhtype:s}'.format(inhtype=inhib['type']))
                self.add_line('INHIBITION_CONSTANT'+' {const:1.{prec}e}'.format(const=inhib['k'],prec=self.precision))
                self.decrease_level()
        # Write out the rest of the reaction attributes
        # What to do with biomass? Maybe better not to use it
        if 'biomass' in reaction_data:
            self.increase_level('BIOMASS')
            self.add_line('SPECIES_NAME '+reaction_data['biomass'])
            self.add_line('YIELD {yld:1.{prec}f}'.format(prec=self.precision,yld=reaction_data.pop('biomass_yield',0)))
            self.decrease_level()
        if 'activation_energy' in reaction_data:
            self.add_line('ACTIVATION_ENERGY {Eact:1.{prec}f}'.format(prec=self.precision,Eact=reaction_data.pop('activation_energy')))
        self.decrease_level()
    
        return self.output
        

# Define decomposition network as a subclass of directed graph (DiGraph)
class decomp_network(nx.MultiDiGraph):
    def __init__(self,pools=[],reactions=[]):
        # First initialize the network object
        super().__init__()
        for pool in pools:
            self.add_pool(pool)
        for reaction in reactions:
            self.add_reaction(reaction)
    def add_pool(self,pool):
        pool_data=pool.copy()
        CN=pool_data.pop('CN',None)
        name=pool_data.pop('name')
        kind=pool_data.pop('kind')
        constraints=pool_data.pop('constraints',{'initial':1e-15})
        self.add_node(name,kind=kind,**pool_data)
        if CN is not None:
            self.nodes[name]['CN']=CN
        self.nodes[name]['constraints']=constraints
    def add_reaction(self,reaction):
        reaction_data=reaction.copy()
        reactant_pools=reaction_data.pop('reactant_pools')
        product_pools=reaction_data.pop('product_pools')
        for upstream_pool in reactant_pools.keys():
            assert upstream_pool in self.nodes,'Pool %s must be added before reaction is added'%upstream_pool
            for downstream_pool in product_pools.keys():
                assert downstream_pool in self.nodes or downstream_pool+'C' in self.nodes,'Pool %s must be added before reaction is added'%downstream_pool
                monod_terms=reaction_data.pop('monod_terms',[])
                inhibition_terms=reaction_data.pop('inhibition_terms',[])
                if downstream_pool not in self.nodes:
                    Dadd='C'
                else:
                    Dadd=''
                k=self.add_edge(upstream_pool,downstream_pool+Dadd,name=reaction_data['name'],reaction=reaction)
                for term in monod_terms:
                    for reqitem in ['species','k']:
                        assert reqitem in term.keys(),'Monod terms must include "%s"'%reqitem
                self.edges[(upstream_pool,downstream_pool+Dadd,k)]['monod_terms']=monod_terms
                for term in inhibition_terms:
                    for reqitem in ['species','k','type']:
                        assert reqitem in term.keys(),'Inhibition terms must include "%s"'%reqitem
                self.edges[(upstream_pool,downstream_pool+Dadd,k)]['inhibition_terms']=inhibition_terms
            
        return
        


def categorize_nodes(nodes):
    colors=[]
    categories=[]
    kinds=nodes.nodes('kind')
    for node in nodes:
        if node in node_colors.keys():
            categories.append(node)
        elif kinds[node] in node_colors.keys():
            categories.append(kinds[node])
        elif 'MICROBES' in node:
            categories.append('Microbe')
        elif 'MAOM' in node or 'SOIL' in node:
            categories.append('MAOM')
        elif 'DOM' in node or 'Acetate' in node:
            categories.append('DOM')
        elif 'LITR' in node:
            categories.append('Litter')
        elif 'CWD' in node:
            categories.append('CWD')
        elif 'kind' in nodes[node] and nodes[node]['kind']=='secondary':
            categories.append('Secondary aqueous')
        elif node in ['HCO3-','HS-'] or '(aq)' in node:
            categories.append('Gas')

        elif kinds[node]=='mineral':
            categories.append('Mineral')
        # elif 'Fe' in node or 'Mn' in node:
        #     categories.append('Metal ion')
        elif kinds[node]=='primary':
            categories.append('Primary aqueous')

        elif kinds[node]=='gas':
            categories.append('Gas')
        elif kinds[node]=='surf_complex':
            categories.append('Surface complex')
        elif kinds[node]=='immobile':
            categories.append('POM')
        else:
             categories.append('Unknown')

    return categories
    
def get_reaction_from_database(name,kind,filename='hanford.dat'):
    with open(filename,'r') as dbase:
        database_sections=['primary','secondary','gas','mineral','surf_complex']
        current_section=0
        for line in dbase:
            if line.startswith("'null'"):
                current_section += 1
            if line.startswith("'%s'"%name) and current_section==database_sections.index(kind):
                out=nx.MultiDiGraph()
                out.add_node(name,kind=kind)
                lsplit=line.split()
                if kind=='secondary':
                    offset=1
                elif kind=='mineral':
                    offset=2
                elif kind=='gas':
                    offset=2
                elif kind=='surf_complex':
                    offset=1
                else:
                    raise TypeError('Kind must be secondary, mineral, gas, or surf_complex')
                nspecies=int(lsplit[offset])
                for specnum in range(nspecies):
                    out.add_edge(name.strip('>'),lsplit[offset+2+specnum*2].strip("'"),reactiontype='equilibrium')
                    out.add_edge(lsplit[offset+2+specnum*2].strip("'"),name,reactiontype='equilibrium')
                return out
        if kind == 'surf_complex':
            raise ValueError('Species {species:s} of kind {kind:s} not found in database. For surf_complex, try searching for the complex instead of the site'.format(species=name,kind=kind))
        else:
            raise ValueError('Species {species:s} of kind {kind:s} not found in database'.format(species=name,kind=kind))

node_colors={
    'POM':'C0',
    'Microbe':'C1',
    'DOM':'C2',
    'MAOM':'C3',
    'Litter':'g',
    'CWD':'brown',
    'Mineral':'C4',
    'Primary aqueous':'C5',
    'Secondary aqueous':'C7',
    'Gas':'C9',
    'Metal ion':'orange',
    'Surface complex':'C8',
    'Unknown':'C5',
    'General Reaction':'#ebde34',
    'Microbial Reaction':'#c2f73b',
    'SOMdecomp Reaction':'#c7a82c',
}


def draw_network(network,omit=[],arrowsize=15,font_size='small',font_weight=None,font_color='k',arrowstyle='->',database_file='hanford.dat',do_legend=True,pos=None,node_colors=node_colors,node_alpha=0.8,label_edges=False,namechanges={},**kwargs):
    to_draw=network.copy()
    for p in network.nodes:
        if network.nodes[p]['kind'] in omit or p in omit or p in ['HRimm','Tracer']:
            to_draw.remove_node(p)
        elif network.nodes[p]['kind'] == 'surf_complex':
            for cplx in network.nodes[p]['complexes']:
                to_draw=nx.compose(get_reaction_from_database(cplx,'surf_complex',filename=database_file),to_draw)
        elif network.nodes[p]['kind'] not in ['primary','immobile','implicit','sorbed']:
            to_draw=nx.compose(get_reaction_from_database(p,network.nodes[p]['kind'],filename=database_file),to_draw)
        
    # Get rid of nodes added from database reactions that we want removed
    to_draw.remove_nodes_from([node for node in to_draw.nodes if to_draw.nodes('kind')[node] is None])
            
    if pos is None:
        pos=nx.drawing.nx_agraph.graphviz_layout(to_draw,prog='dot')
    nodecats=categorize_nodes(to_draw)  
    nodecolors=[node_colors[nodecat] for nodecat in nodecats]
    
    nx.draw_networkx_nodes(to_draw,pos=pos,node_color=nodecolors,alpha=node_alpha,**kwargs)
    nx.draw_networkx_edges(to_draw,pos=pos,arrowsize=arrowsize,arrowstyle=arrowstyle,**kwargs)
    nx.draw_networkx_labels(to_draw,pos=pos,labels={n:namechanges.get(n,n) for n in to_draw.nodes},font_size=font_size,font_weight=font_weight,**kwargs)
    
    if label_edges:            
        nx.draw_networkx_edge_labels(to_draw,pos=pos,edge_labels=dict([((e[0],e[1]),e[2]) for e in to_draw.edges(data='name')]),font_size=font_size,fontstyle='italic')
                
    if do_legend:
        from matplotlib.pyplot import Circle,legend
        legend_handles=[]
        legend_labels=[]
        
        for num,node in enumerate(to_draw.nodes):
            if nodecats[num] not in legend_labels:
                legend_labels.append(namechanges.get(nodecats[num],nodecats[num]))
                legend_handles.append(Circle((0,0),radius=5,facecolor=nodecolors[num]))
                
        legend(handles=legend_handles,labels=legend_labels,fontsize='medium',title='Component types')
    
    return to_draw
    
def draw_network_with_reactions(network,omit=[],arrowsize=15,font_size='small',font_weight=None,arrowstyle='->',database_file='hanford.dat',do_legend=True,
            node_colors=node_colors,namechanges={},font_color='k',node_alpha=0.8,node_size=None,edge_color=None,markers={'Reaction':'*'},pos=None,
            width=None,connectionstyle=None,**kwargs):
    to_draw=network.copy()
    
    for p in network.nodes:
        if network.nodes[p]['kind'] in omit or p in omit or p in ['HRimm','Tracer']:
            to_draw.remove_node(p)
        elif network.nodes[p]['kind'] == 'surf_complex':
            for cplx in network.nodes[p]['complexes']:
                # to_draw=nx.compose(get_reaction_from_database(cplx,'surf_complex',filename=database_file),to_draw)
                to_draw.add_edge(p.strip('>'),cplx.strip('>'))
                to_draw.nodes[cplx.strip('>')]['kind']='Sorption Reaction'
                to_draw.add_edge(network.nodes[p]['mineral'],cplx.strip('>'))
                to_draw.remove_node(p)
        elif network.nodes[p]['kind'] not in ['primary','immobile','implicit','sorbed']:
            to_draw=nx.compose(get_reaction_from_database(p,network.nodes[p]['kind'],filename=database_file),to_draw)
    
    for react in network.edges:
        if network.edges[react]['name'] not in to_draw.nodes and network.edges[react]['name'] not in omit:
            e=network.edges[react]['reaction']
            for species in e['reactant_pools']:
                to_draw.add_edge(species,e['name'])
            for species in e['product_pools']:
                to_draw.add_edge(e['name'],species)
            to_draw.nodes[e['name']]['kind']=e['reactiontype'].capitalize().replace('Som','SOM') + ' Reaction'
            to_draw.nodes[e['name']]['reactiontype']=e['reactiontype']
            
    # Get rid of nodes added from database reactions that we want removed
    to_draw.remove_nodes_from([node for node in to_draw.nodes if to_draw.nodes('kind')[node] is None])
    # Get rid of all the original reactions
    to_draw.remove_edges_from(network.edges)

            
    if pos is None:
        pos=nx.drawing.nx_agraph.graphviz_layout(to_draw,prog='dot')
    from numpy import array
    nodecats=array(categorize_nodes(to_draw)  )
    nodecolors=array([node_colors[nodecat] for nodecat in nodecats])
    
    # Non-reactions:
    nonreactions=array(['Reaction' not in to_draw.nodes[n]['kind'] for n in to_draw.nodes])
    minerals=array(['mineral' in to_draw.nodes[n]['kind'] for n in to_draw.nodes])
    nx.draw_networkx_nodes(to_draw,pos=pos,nodelist=array(to_draw.nodes())[nonreactions&~minerals].tolist(),node_color=nodecolors[nonreactions&~minerals],node_size=node_size,node_shape='o',alpha=node_alpha,**kwargs)
    nx.draw_networkx_nodes(to_draw,pos=pos,nodelist=array(to_draw.nodes())[minerals].tolist(),node_color=nodecolors[minerals],node_shape=markers.get('mineral','o'),alpha=node_alpha,node_size=node_size,**kwargs)
        
    nx.draw_networkx_labels(to_draw,pos=pos,labels={n:namechanges.get(n,n) for n in to_draw.nodes},font_size=font_size,font_color=font_color,font_weight=font_weight,**kwargs)
    
    reactionnodes=array(to_draw.nodes())[~nonreactions].tolist()
    nx.draw_networkx_nodes(to_draw,pos=pos,nodelist=reactionnodes,node_shape=markers.get('Reaction','*'),
                node_color=nodecolors[~nonreactions],alpha=node_alpha,node_size=node_size,**kwargs)
    
    nx.draw_networkx_edges(to_draw,pos=pos,connectionstyle=connectionstyle,arrowsize=arrowsize,arrowstyle=arrowstyle,edge_color=edge_color,width=width,**kwargs)
        
    
    
    if do_legend:
        from matplotlib.pyplot import legend,Line2D
        legend_handles=[]
        legend_labels=[]
        
        for num,node in enumerate(to_draw.nodes):
            if namechanges.get(nodecats[num],nodecats[num]) not in legend_labels:
                legend_labels.append(namechanges.get(nodecats[num],nodecats[num]))
                if 'Reaction' in to_draw.nodes[node]['kind']:
                    if markers.get('Reaction','*') != 'None':
                        legend_handles.append(Line2D([0],[0],ls='None',marker=markers.get('Reaction','*'),ms=15.0,color=nodecolors[num]))
                elif 'mineral' in to_draw.nodes[node]['kind']:
                    legend_handles.append(Line2D([0],[0],ls='None',marker=markers.get('mineral','o'),ms=15.0,color=nodecolors[num]))
                else:
                    legend_handles.append(Line2D([0],[0],ls='None',marker='o',ms=15.0,color=nodecolors[num]))
                
        legend(handles=legend_handles,labels=legend_labels,fontsize='large',title='Component types',title_fontsize='large',labelspacing=1.0,ncol=2)
    
    return to_draw,pos
    

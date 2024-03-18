import decomp_network #,plot_pf_output
from numpy import *
from matplotlib import pyplot
import matplotlib

pools = [
decomp_network.decomp_pool(name='Cellulose',CN=50,constraints={'initial':5e3*0.6},kind='immobile'),
decomp_network.decomp_pool(name='Lignin',CN=50,constraints={'initial':5e3*0.4},kind='immobile'),
decomp_network.decomp_pool(name='HRimm',constraints={'initial':1e-20},kind='immobile'),
decomp_network.decomp_pool(name='Root_biomass',constraints={'initial':1e-20},kind='immobile'),

decomp_network.decomp_pool(name='DOM1',CN=50,constraints={'initial':1e-30},kind='primary'),
decomp_network.decomp_pool(name='DOM2',CN=50,constraints={'initial':1e-30},kind='primary'),
decomp_network.decomp_pool(name='H+',kind='primary',constraints={'initial':'5.0 P'}),
decomp_network.decomp_pool(name='O2(aq)',kind='primary',constraints={'initial':1e4}),
decomp_network.decomp_pool(name='HCO3-',kind='primary',constraints={'initial':'400e-6 G CO2(g)'}),
decomp_network.decomp_pool(name='Mn+++',kind='primary',constraints={'initial':'1.0e-30 M Birnessite2'}), # Should I just use MnO4-- as the primary instead? That would allow using original Bernessite def
decomp_network.decomp_pool(name='Mn++',kind='primary',constraints={'initial':'1.0e-30'}), # Change this to reflect observed exchangeable Mn
decomp_network.decomp_pool(name='chelated_Mn+++',kind='primary',constraints={'initial':'1.0e-30'}),
decomp_network.decomp_pool(name='NH4+',kind='primary',constraints={'initial':1e-15}), # SOMDecomp sandbox requires this
decomp_network.decomp_pool(name='Tracer',kind='primary',constraints={'initial':1e-15}), # Just to accumulate CO2 loss
decomp_network.decomp_pool(name='Tracer2',kind='primary',constraints={'initial':1e-15}), # To accumulate root Mn++ uptake
# decomp_network.decomp_pool(name='CH4(aq)',kind='primary',constraints={'initial':1e-15}),
# decomp_network.decomp_pool(name='Acetate-',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='Mg++',kind='primary',constraints={'initial':0.5e-3}),
decomp_network.decomp_pool(name='Ca++',kind='primary',constraints={'initial':0.5e-3}),
decomp_network.decomp_pool(name='Na+',kind='primary',constraints={'initial':2e-3}),
decomp_network.decomp_pool(name='K+',kind='primary',constraints={'initial':2e-3}),
decomp_network.decomp_pool(name='Al+++',kind='primary',constraints={'initial':0.5e-3}),
# decomp_network.decomp_pool(name='Fe+++',kind='primary',constraints={'initial':'0.37e-5'}),
# decomp_network.decomp_pool(name='Fe++',kind='primary',constraints={'initial':'0.37e-2'}),

decomp_network.decomp_pool(name='CO2(g)',kind='gas'),
decomp_network.decomp_pool(name='O2(g)',kind='gas'),
decomp_network.decomp_pool(name='H2O',kind='implicit'),

decomp_network.decomp_pool(name='CO2(aq)',kind='secondary'),
decomp_network.decomp_pool(name='OH-',kind='secondary'),
decomp_network.decomp_pool(name='MnO4--',kind='secondary'),
# decomp_network.decomp_pool(name='Acetic_acid(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='MnIIIDOM2(aq)',kind='secondary'),
# Should add MnIII complex with DOM1
decomp_network.decomp_pool(name='MnIIIDOM1(aq)',kind='secondary'),

# Hui/Beth say pyrolusite probably not precipitating in soils. Look for delta-MnO2? Maybe try 'Mn(OH)2(am)'?
# See Roberts Earth Sci Rev article for some discussion of minerals? More Fe focused though
# July 2020: Beth suggests using Birnessite or other Mn(IV) oxide for precipitation/dissolution source/sink. 
#            Soil and water chemistry (Essington) textbook: Birnessite common in young soils subject to redox fluctuations
#              Precipitate during oxidation of Mn++. Poorly crystalline. Permanent structural charge satisfied by exchangeable interlayer cations
#   Birnessite: (Na, Ca, MnII)(MnIII,MnIV)7O14*2.8H2O
#   Parc et al 1989: Birnessite solubility products: 5 MnIV + 2 Mn3 + 13 O + 5 H2O -> 7 log[Mn++]/[H+]^2 + 3 log(fO2) = 0.48-9.24 (uncertain value)
# Original Hanford database: 'Birnessite' 251.1700 4 -4.0000 'H+' 3.0000 'MnO4--' 5.0000 'Mn++' 7.0000 'H2O' 500.0000 -85.5463 500.0000 500.0000 500.0000 500.0000 500.0000 500.0000 753.5724
# How about we start with 7 MnIII and implicitly oxidize 5 of them?
# 7 Mn+++ + 5 H+ + 1.25 O2  <-> 2 Mn+++ + 5 Mn4+ + 2.5 H2O
# 2 Mn+++ + 5 Mn4+ + 18 H2O <->  Mn7O13*5H2O + 26 H+
# Combined: 7 Mn+++ + 5 H+ + 1.25 O2 + 15.5 H2O <-> Mn7O13*5H2O + 26 H+
#           7 Mn+++ + 1.25 O2 + 15.5 H2O -> Mn7O13*5H2O + 21 H+
# Molar weight should be 215.71 g/mol (http://www.webmineral.com/data/Birnessite.shtml) But database has 753.57
# Density = 3 g/cm3 = 71.9 cm3/mol (https://www.mindat.org/min-680.html) (compared to 251 for database)
# Looks like database molar numbers are multiplied by 3.5 relative to these. I think this is actually correct since official formula is (Mn++++,Mn+++)2O4•1.5(H2O) which must be multiplied by 3.5 to get our 7 total Mn
# Suggests that actual logK is also scaled by 3.5, so original logK would be -24 instead of -85. Still pretty low...
# Database format: 'Birnessite' 251.1700 4 -21.0 'H+' 7.0 'Mn+++' 1.25 O2(aq) 15.5 H2O ... 753.5724
decomp_network.decomp_pool(name='Birnessite2',rate='1.d-15 mol/m^2-sec',constraints={'initial':'0.00122  1.d2 m^2/m^3'},kind='mineral'),
# decomp_network.decomp_pool(name='Mn(OH)2(am)',rate='1.d-16 mol/m^2-sec',constraints={'initial':'0.0d-5  1.d2 m^2/m^3'},kind='mineral'),
# decomp_network.decomp_pool(name='Manganite',rate='0.d-12 mol/m^2-sec',constraints={'initial':'0.0d-5  1.d2 m^2/m^3'},kind='mineral'),

decomp_network.decomp_pool(name='Rock(s)',rate='0.0 mol/m^2-sec',constraints={'initial':'0.5  5.0e3 m^2/m^3'},kind='mineral'),

decomp_network.decomp_pool(name='>Carboxylate-',kind='surf_complex',mineral='Rock(s)',site_density=1.0e4,complexes=['>Carboxylic_acid']),

# This should hopefully work for sorption on Mn minerals. BUT, probably need to run with alquimia to allow site density to change over time
# since site density depending on mineral surface area is not implemented in PFLOTRAN currently
#    Update: This is wrong, documentation out of date. Site density should update as mineral volume fraction changes
# decomp_network.decomp_pool(name='>DOM1',kind='surf_complex',mineral='Rock(s)',site_density=1.0e4,complexes=['>sorbed_DOM1']),
# decomp_network.decomp_pool(name='Sorbed DOM1',kind='isotherm',mineral='Manganite',site_density=0.0e3,complexes=['>sorbed_DOM1']),

# Alternative attempt: Treat sorption as a set of two reactions. To avoid editing pflotran code, use a microbial and a general reaction
# This requires us to manually treat the sorbed stuff as immobile in the code we are using to run the model
decomp_network.decomp_pool(name='DOM3',constraints={'initial':1e-30},kind='primary'),
# Sorption capacity is treated as an immobile biomass term in the sorption equation, which makes eq sorbed DOM proportional to that term
decomp_network.decomp_pool(name='Sorption_capacity',constraints={'initial':1e-20},kind='immobile'),
# sorption rate is microbial DOM1 -> DOM3 (d/dt = V1*sorption_capacity*DOM1/(DOM1+k))
# Desorption rate is general DOM3 -> DOM1 (d/dt = V2*DOM3)
# Equilibrium sorption = V1/V2*sorption_capacity*DOM/(k+DOM) so actual maximum sorption is V1/V2*sorption_capacity
]

# Fraction of Mn+++ that is not recycled in Mn-Peroxidase enzyme loop
def make_network(leaf_Mn_mgkg=25.0,change_constraints={},Mn2_scale=1e-5,Mn3_scale=1e-5,NH4_scale=1e-2,DOM_scale=1.0,CEC=40.0,change_rate={},DOM_sorption_k=1.0):
    Mn_molarmass=54.94 #g/mol
    C_molarmass=12.01
    leaf_Cfrac_mass=0.4
    # Converting from mol Mn/mol leaf C to mg Mn/kg leaf dry mass
    leaf_Mn_frac=leaf_Mn_mgkg / ((Mn_molarmass*1e3)/(C_molarmass*1e-3/leaf_Cfrac_mass))
    reactions = [
            decomp_network.reaction(name='Hydrolysis',reactant_pools={'Cellulose':1.0},product_pools={'DOM1':1.0},
                                            rate_constant=2e-1,rate_units='y', 
                                        # inhibition_terms=[decomp_network.inhibition(species='DOM1',type='MONOD',k=1e-1)],
                                        reactiontype='SOMDECOMP'),
                                        
    # Need a way to make a soluble reactive compound from lignin. Treating it like hydrolysis, but allowing only a small amount to be generated
    # Model runs faster if Mn limitation is on this step (source limit) instead of limiting decomposition rate of DOM2 and building up inhibiting DOM2 concentrations
            decomp_network.reaction(name='Lignin exposure',reactant_pools={'Lignin':1.0},product_pools={'DOM2':1.0},
                                            rate_constant=2e-1,rate_units='y', 
                                        inhibition_terms=[decomp_network.inhibition(species='DOM2',type='MONOD',k=2e-3,pool_normalized=True),
                                        decomp_network.inhibition(species='Cellulose',type='MONOD',k=40), # Suggested by Sun et al 2019, which found that MnP activity only increased late in decomposition
                                        # decomp_network.inhibition(species='Mn++',type='INVERSE_MONOD',k=Mn2_scale),decomp_network.inhibition(species='NH4+',k=NH4_scale,type='MONOD')
                                        ],
                                        reactiontype='SOMDECOMP'),
                                        
            # Direct depolymerization of lignin into DOM1. Basically to allow some decomposition to occur in absence of Mn++ so lignin doesn't build up indefinitely
            # Still deciding whether to do this way or by DOM2 oxidation? Or leaching of DOM2? Or not at all?
            # decomp_network.reaction(name='Lignin depolymerization',reactant_pools={'Lignin':1.0},product_pools={'DOM1':1.0},
            #                                 rate_constant=2e-1,rate_units='y', 
            #                             inhibition_terms=[decomp_network.inhibition(species='Mn++',type='MONOD',k=Mn2_scale)],
            #                                                 reactiontype='SOMDECOMP'),
                                        
    # Mn reduction turns DOM2 (lignin-based) into DOM1 (odixized, decomposable)
            # decomp_network.reaction(name='Lignin Mn reduction',reactant_pools={'DOM2':1.0,'Mn+++':1.0},product_pools={'Mn++':1.0,'H+':1.0,'DOM1':1.0},
            #                                 monod_terms=[decomp_network.monod(species='DOM2',k=2e-11,threshold=1.1e-20),decomp_network.monod(species='Mn+++',k=1.3e-18,threshold=1.1e-25)],
            #                                 # inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=6.25e-8,type='MONOD')],
            #                                 rate_constant=2e-10,reactiontype='MICROBIAL'),
                                            
    # Mn oxidation (Mn peroxidase?)
            # decomp_network.reaction(name='Mn oxidation',reactant_pools={'O2(aq)':1.0,'H+':4.0,'Mn++':4.0},product_pools={'Mn+++':4.0},
            #                                 monod_terms=[decomp_network.monod(species='O2(aq)',k=1e-5,threshold=1.1e-15),decomp_network.monod(species='Mn++',k=1.3e-12,threshold=1.1e-15)],
            #                                 # inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=6.25e-8,type='MONOD')],
            #                                 rate_constant=8e-10,reactiontype='MICROBIAL'),
                                        
    # Mn Peroxidase: Quasi-closed reaction loop for Mn2/3 with some leakage.
        # 1 DOM2 + 1 Mn++ + x H+ -> 1 DOM1 + (1-x) Mn++ + x Mn+++
            # decomp_network.reaction(name='Mn Peroxidase',reactant_pools={'DOM2':1.0,'Mn++':Mn_peroxidase_Mn3_leakage,'H+':Mn_peroxidase_Mn3_leakage},
            #                                              product_pools={'DOM1':1.0,'Mn+++':Mn_peroxidase_Mn3_leakage},
            #                                              monod_terms=[decomp_network.monod(species='DOM2',k=1e-1,threshold=1.1e-20),decomp_network.monod(species='Mn++',k=Mn2_scale*1e1,threshold=1.1e-20)],
            #                                              inhibition_terms=[decomp_network.inhibition(species='NH4+',k=NH4_scale,type='MONOD')],
            #                                              rate_constant=5e-6,reactiontype='MICROBIAL'),

    # Mn Peroxidase: Two step version
        # Mn++ + H+ -> chelated_Mn+++
        # chelated_Mn+++ + DOM2 -> Mn++ + DOM1 + H+
        # 2 chelated_Mn+++ -> Mn++ + Mn4+ ~ 2 chelated_Mn+++ + 0.5 H2O -> Mn++ + Mn+++ + 0.25 O2 + H+
            decomp_network.reaction(name='Mn Peroxidase',stoich='1.0 Mn++ + H+ -> chelated_Mn+++',
                                                         monod_terms=[decomp_network.monod(species='Mn++',k=Mn2_scale,threshold=1.1e-20)],
                                                         inhibition_terms=[decomp_network.inhibition(species='NH4+',k=NH4_scale,type='MONOD'),
                                                                           decomp_network.inhibition(species='DOM2',k=1e-3,type='INVERSE_MONOD')],
                                                         rate_constant=5e-6,reactiontype='MICROBIAL'),  

            decomp_network.reaction(name='Mn chelate lignin',stoich='1.0 chelated_Mn+++ + 1.0 DOM2 -> 1.0 Mn++ + 1.0 DOM1 + 1.0 H+',
                                                        #  monod_terms=[decomp_network.monod(species='chelated_Mn+++',k=Mn3_scale,threshold=1.1e-20),decomp_network.monod(species='DOM2',k=1e-1,threshold=1.1e-20)],
                                                        #  inhibition_terms=[decomp_network.inhibition(species='NH4+',k=NH4_scale,type='MONOD')],
                                                         rate_constant=5e-6,reactiontype='GENERAL'),  
            decomp_network.reaction(name='Mn chelate loss',stoich='2.0 chelated_Mn+++ + 0.5 H2O -> 1.0 Mn++ + 1.0 Mn+++ + 0.25 O2(aq) + 1.0 H+',
                                                        #  monod_terms=[decomp_network.monod(species='chelated_Mn+++',k=Mn2_scale*1e1,threshold=1.1e-20),decomp_network.monod(species='DOM2',k=1e-1,threshold=1.1e-20)],
                                                        #  inhibition_terms=[decomp_network.inhibition(species='NH4+',k=NH4_scale,type='MONOD')],
                                                         rate_constant=5e-6,reactiontype='GENERAL'),                       
                                                             

    # CH2O + H2O -> CO2 + 4H+ + 4 e-
    # O2   + 4H+ + 4 e- -> 2H2O
            decomp_network.reaction(name='DOM aerobic respiration',reactant_pools={'DOM1':1.0,'O2(aq)':1.0},product_pools={'HCO3-':1.0,'H+':1.0,'Tracer':1.0,'Mn++':leaf_Mn_frac},
                                            monod_terms=[decomp_network.monod(species='O2(aq)',k=1e-4,threshold=1.1e-12),decomp_network.monod(species='DOM1',k=DOM_scale,threshold=1.1e-14)],
                                        rate_constant=1.0e-5,reactiontype='MICROBIAL'),
                                        
    # Mn-independent lignin decomposition pathway. Should be slower than Mn peroxidase
            decomp_network.reaction(name='Mn-independent lignin degradation',stoich='1.0 DOM2 -> 1.0 DOM1',
                                            monod_terms=[decomp_network.monod(species='DOM2',k=0.05)],
                                            # inhibition_terms=[decomp_network.inhibition(species='NH4+',k=NH4_scale,type='MONOD')],
                                        rate_constant=1.0e-5,reactiontype='GENERAL'),

    # Manganese reduction reaction
    # CH2O + 2 H2O -> HCO3- + 5 H+ + 4 e-
    # 4 Mn+++ + 4 e- -> 4 Mn++
    # Microbial-mediated manganese reduction, happens under anoxic conditions only
            decomp_network.reaction(name='DOM1 Mn+++ reduction',stoich='1.0 DOM1 + 4.0 Mn+++ -> 1.0 HCO3- + 4.0 Mn++ + 5.0 H+ + 1.0 Tracer',
                                        monod_terms=[decomp_network.monod(species='DOM1',k=DOM_scale/10),decomp_network.monod(species='Mn+++',k=Mn3_scale)],
                                        inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=1e-8,type='MONOD')],
                                        rate_constant=2e-10,reactiontype='MICROBIAL'),

    # Abiotic Mn reduction, happens under oxic conditions
            decomp_network.reaction(name='DOM1 Mn+++ abiotic reduction',stoich='1.0 DOM1 + 4.0 Mn+++ -> 1.0 HCO3- + 4.0 Mn++ + 5.0 H+ + 1.0 Tracer',
                                        rate_constant=2e-5,reactiontype='GENERAL'),

    # Microbial oxidation of Mn(II) under oxic conditions. See Tebo et al 2004. Mostly done by bacteria. Reduces O2 to water. Function unknown but might be chemoautotrophic
    # Mn++ + 0.5 O2 + H2O -> Mn(+4)O2 + 2 H+ (but this includes an intermediate Mn3O4 step)
    # or Mn++ + 0.25 O2 + 1.5 H2O -> Mn(+3)OOH + 2 H+
    # or 3 Mn++ + 0.5 O2 + 3 H2O -> Mn(+3)3O4 + 6 H+
    # Simplified version assuming further oxidation happens in Birnessite precipitation step:
    # Mn++ + H+ + 0.25 O2 -> Mn+++ + 0.5 H2O
    # then (implicitly) 7 Mn+++ + 1.25 O2 + 15.5 H2O -> Mn7O13*5H2O + 21 H+
            decomp_network.reaction(name='Bacterial Mn++ oxidation',stoich='Mn++ + H+ + 0.25 O2(aq) -> Mn+++ + 0.5 H2O',
                                    monod_terms=[decomp_network.monod(species='Mn++',k=Mn2_scale),decomp_network.monod(species='O2(aq)',k=1e-4,threshold=1.1e-12)],
                                    rate_constant=8e-8,reactiontype='MICROBIAL'), # Rate constant estimated from Fig 8 in Tebo et al 2004

    # C2H3O2- + 2 H2O -> 2 CO2 + 7 H+ + 8 e-
    # 2 O2    + 8 H+ + 8 e- -> 4 H2O
            # decomp_network.reaction(name='Acetate aerobic respiration',reactant_pools={'Acetate-':1.0,'O2(aq)':2.0},product_pools={'HCO3-':2.0,'H+':2.0,'Tracer':2.0},
            #                                 monod_terms=[decomp_network.monod(species='O2(aq)',k=1e-5,threshold=1.1e-12),decomp_network.monod(species='Acetate-',k=1e-8,threshold=1.1e-14)],
            #                             rate_constant=1.0e-9,reactiontype='MICROBIAL'),
            # 
            
            # PFLOTRAN will not let you do isotherm sorption if there are any secondary species :-( - reaction_database.F90 line 3424
            # decomp_network.sorption_isotherm(name='DOM sorption',mineral='Rock(s)',sorbed_species='DOM1',k=0.1,langmuir_b=1.0,sorbed_name='Sorbed DOM1'),
            
    # Root uptake of Mn++. Does this need a charge balancing release of anions?
    # According to Haynes (1990), imbalance in root anion vs cation uptake is usually balanced by H+ or OH- release. I guess in this case we need to assume that all other ions are balanced except for Mn?
            decomp_network.reaction(name='Root uptake of Mn++',stoich='1.0 Mn++ -> 1.0 Tracer2 + 2 H+',monod_terms=[decomp_network.monod(species='Mn++',k=1e-7,threshold=1.1e-20)],
                                    biomass='Root_biomass',biomass_yield=0.0,
                                    rate_constant=1.0e-5,reactiontype='MICROBIAL'),
            
            # Cation exchange
            decomp_network.ion_exchange(name='Cation exchange',cations={'Mn++':5.0,'Mn+++':0.6,'Al+++':10.0,'Mg++':1.1,'Ca++':1.0,'Na+':0.5,'K+':0.1,'H+':1.1},CEC=CEC,mineral=None),
            
            # sorption rate is microbial DOM1 -> DOM3 (d/dt = V1*sorption_capacity*DOM1/(DOM1+k))
            # Desorption rate is general DOM3 -> DOM1 (d/dt = V2*DOM3)
            # Equilibrium sorption = V1/V2*sorption_capacity*DOM/(k+DOM) so actual maximum sorption is V1/V2*sorption_capacity
            decomp_network.reaction(name='DOM sorption',stoich='DOM1 -> DOM3',monod_terms=[decomp_network.monod(species='DOM1',k=DOM_sorption_k)],
                                biomass='Sorption_capacity',biomass_yield=0.0,rate_constant=1.0e-10,reactiontype='MICROBIAL'),
            decomp_network.reaction(name='DOM desorption',stoich='DOM3 -> DOM1',reactiontype='GENERAL',rate_constant=1.0e-10),
            
    ] # End of reactions list

    import copy
    pools_copy=copy.deepcopy(pools)
    for n,p in enumerate(pools_copy):
        if p['name'] in change_rate:
            pools_copy[n]['rate']=change_rate[p['name']]
    return decomp_network.decomp_network(pools=decomp_network.change_constraints(pools_copy, change_constraints),reactions=reactions)


def plot_result(result,SOM_ax=None,pH_ax=None,mineral_ax=None,gasflux_ax=None,porewater_ax=None,do_legend=False):

    if SOM_ax is not None:
        l=SOM_ax.plot(result['Cellulose']*1e-3,label='Cellulose')[0]
        SOM_ax.plot(result['Lignin']*1e-3,label='Lignin')[0]

        SOM_ax.set_title('SOM remaining')
        SOM_ax.set_ylabel('Concentration\n(mmol C/cm$^{-3}$)')
        SOM_ax.set_xlabel('Time (days)')
        if do_legend:
            SOM_ax.legend()

    if pH_ax is not None:
        pH_ax.plot(-log10(result['Free H+']))

        pH_ax.set_title('pH')
        pH_ax.set_ylabel('pH')
        pH_ax.set_xlabel('Time (days)')
        
    if mineral_ax is not None:
        molar_volume_birnessite=251.1700 # From database. cm3/mol (first number)
        molar_weight_birnessite = 753.5724 # (last number)
        molar_volume_manganite = 24.45
        molar_volume_pyrolusite = 500.0 # From database. Seems fishy
        l=mineral_ax.plot(result['Birnessite2 VF']/molar_volume_birnessite*1e6   ,label='Birnessite')[0]
        # l=mineral_ax.plot(result['Manganite VF']/molar_volume_manganite*1e6   ,label='Manganite (Mn+++)')[0]
        
        # l=mineral_ax.plot(result['Total Mn+++']*result['Porosity']*1e3   ,label='Mn+++',ls='--')[0]
        
        # l=mineral_ax.plot(result['Total Mn++']*result['Porosity']*1e3 ,ls=':'  ,label='Mn++')[0]
        
        mineral_ax.set_title('Mn minerals')
        mineral_ax.set_ylabel('Concentration\n($\mu$mol/cm$^{-3}$)')
        mineral_ax.set_xlabel('Time (days)')
        if do_legend:
            mineral_ax.legend()
    
    if gasflux_ax is not None:
        gasflux_ax.set_yscale('log')
        
        # l=gasflux_ax.plot(result.index.values[:-1],diff(result['Total CH4(aq)']*result['Porosity'])/diff(result.index.values)*1e3,label='CH4')[0]
        
        l=gasflux_ax.plot(result.index.values[:-1],diff(result['Total Tracer']*result['Porosity'])/diff(result.index.values)*1e3,label='CO2',ls='-',c='C5')[0]

        gasflux_ax.set_title('Gas fluxes')
        gasflux_ax.set_ylabel('Flux rate\n($\mu$mol cm$^{-3}$ day$^{-1}$)')
        gasflux_ax.set_xlabel('Time (days)')
        # if do_legend:
        #     gasflux_ax.legend(fontsize='small')
        
    if porewater_ax is not None:
        porewater_ax.set_yscale('log')
        porewater_ax.plot(result['Total DOM1'],label='DOM (oxidized)')
        porewater_ax.plot(result['Total DOM2'],label='DOM (lignin)')
        # porewater_ax.plot(result['Total Acetate-'],label='Acetate',c='C3')
        porewater_ax.plot(result['Total O2(aq)'],'--',label='O2',c='C4')
        porewater_ax.plot(result['Total Mn+++'],'--',label='Mn+++',c='C1')
        porewater_ax.plot(result['Total Mn++'],':',label='Mn++',c='C2')
        
        porewater_ax.set_title('Porewater concentrations')
        porewater_ax.set_ylabel('Concentration (M)')
        porewater_ax.set_xlabel('Time (days)')
        
        if do_legend:
            porewater_ax.legend(fontsize='small',ncol=2)
    
# result,units=plot_pf_output.convert_units(result,units,'mol/m^3')
if __name__ == '__main__':
    reaction_network=make_network()
    networkfig=pyplot.figure('Reaction network');networkfig.clf()
    drawn=decomp_network.draw_network_with_reactions(reaction_network,omit=['NH4+','Rock(s)','gas','secondary','H+','>Carboxylate-','Carboxylic_acid'],
            font_size='medium',node_size=1500,font_color='k',arrowstyle='->',arrowsize=10.0,edge_color='gray',node_alpha=1.0,
            namechanges={'cellulose':'Cellulose','DOM1':'DOM','O2(aq)':'O$_2$(aq)','CH4(aq)':'CH$_4$(aq)','HCO3-':'HCO$_3^-$','DOM2':'Exposed lignin','sorbed_DOM1':'Sorbed DOM',
                         'Fe(OH)2':'Fe(OH)$_2$','Fe(OH)3':'Fe(OH)$_3$','Mn++':r'Mn$^\mathrm{+\!\!+}$','Mn+++':r'Mn$^\mathrm{+\!\!+\!\!\!+}$','Acetate-':'Acetate',})


    pflotran_exe='../pflotran-interface/src/pflotran/pflotran'
    simlength=60
    simlength=41*30 #Sun et al was 41 months

    # Example results
    results_all,units=decomp_network.PF_network_writer(make_network(leaf_Mn_mgkg=500.0,change_constraints={'Birnessite':'1.0d-5 1.d2 m^2/m^3','Manganite':'0.0973d-5 1.d2 m^2/m^3'})).run_simulation('SOMdecomp_template.txt','manganese',pflotran_exe,print_output=False,length_days=simlength)
    resultsfig,axes=pyplot.subplots(4,1,num='Results',clear=True,figsize=(6.3,8),squeeze=False)
    plot_result(results_all,SOM_ax=axes[0,0],pH_ax=axes[1,0],mineral_ax=axes[2,0],porewater_ax=axes[3,0],do_legend=True)


    cm=pyplot.get_cmap('coolwarm')

    import pandas

    fig,axes=pyplot.subplots(3,1,num='Leaf Mn concentrations',clear=True)
    axes[0].set_title('Lignin remaining')
    axes[0].set_ylabel('Fraction of initial')
    axes[1].set_title('Mn++ concentration')
    axes[1].set_ylabel('Concentration (M)')
    axes[2].set_title('Manganite concentration')
    axes[2].set_ylabel('Concentration\n($\mu$mol/cm$^{-3}$)')
    axes[2].set_xlabel('Time (days)')
    # Range of leaf Mn concentrations
    norm=matplotlib.colors.Normalize(1,5)
    results_all_leafMn={}
    for leafMn in arange(1.0,5.1) :
        result,units=decomp_network.PF_network_writer(make_network(leaf_Mn_mgkg=10.0**leafMn)).run_simulation('SOMdecomp_template.txt','manganese',pflotran_exe,print_output=False,length_days=simlength)
        axes[0].plot(result['Lignin']/result['Lignin'].iloc[0],label='Leaf Mn concentration = 10$^{%d}$ mg/kg'%int(leafMn),c=cm(norm(leafMn)))
        axes[1].plot(result['Total Mn++'],c=cm(norm(leafMn)))
        # axes[2].plot(result['Manganite VF']/24.45*1e6,c=cm(norm(leafMn)))
        birnessite=axes[2].plot(result['Birnessite2 VF']/251.17*1e6,c=cm(norm(leafMn)),ls='--',label='Birnessite')[0]
        results_all_leafMn[leafMn]=result

    axes[0].legend()


    nyears=40
    results_long_leafMn=[]
    nums=linspace(1,35,10)
    nums=logspace(-1,1.5,10)

    cellulose_eq=zeros(len(nums))
    lignin_eq=zeros(len(nums))
    for num in range(len(nums)):
        result,units=decomp_network.PF_network_writer(make_network(leaf_Mn_mgkg=num,Mn2_scale=0.25e-2)).run_simulation('SOMdecomp_template.txt','manganese',pflotran_exe,print_output=False,length_days=nyears*365)
        results_long_leafMn.append(result)

    norm=matplotlib.colors.LogNorm(nums.min(),nums.max())
    fig,axes=pyplot.subplots(ncols=3,num='Multiple years',clear=True)
    volume_factor=2/100*12e-3 # Converting from mol/m3 to kgC/m2 assuming 1 cm depth
    for num in range(len(nums)):
        cellulose=zeros(nyears*365)
        lignin=zeros(nyears*365)
        for yr in range(nyears):
            cellulose[365*yr:]+=interp(arange(nyears*365)/365,results_long_leafMn[num].index/365,results_long_leafMn[num]['Cellulose'])[:365*(nyears-yr)]
            lignin[365*yr:]+=interp(arange(nyears*365)/365,results_long_leafMn[num].index/365,results_long_leafMn[num]['Lignin'])[:365*(nyears-yr)]
        
        axes[1].plot(arange(nyears*365)/365,lignin*volume_factor,ls='--',c=cm(norm(nums[num]))) # Assumes 1 cm thick organic horizon
        cellulose_eq[num]=cellulose[-365:].mean()*volume_factor
        lignin_eq[num]=lignin[-365:].mean()*volume_factor
        
        axes[0].plot(results_long_leafMn[num].index/365,results_long_leafMn[num]['Lignin']*volume_factor,ls='--',c=cm(norm(nums[num])),label='Lignin: %1.1f mg/g Mn'%nums[num])
        axes[2].plot(nums[num],cellulose_eq[num]+lignin_eq[num],'o',c=cm(norm(nums[num])),ms=8.0)
        
    axes[2].plot(nums,cellulose_eq+lignin_eq,'k-')
    axes[0].plot(results_long_leafMn[num].index/365,results_long_leafMn[-1]['Cellulose']*volume_factor,'k-',label='Cellulose')
    axes[1].plot(arange(nyears*365)/365,cellulose*volume_factor,'k-')
    axes[0].legend(loc='upper right')

    axes[0].set(title='One leaf litter cohort',xlabel='Time (years)',ylabel='C stock (kg m$^{-2}$)')
    axes[1].set(title='Cumulative litter layer',xlabel='Time (years)',ylabel='C stock (kg m$^{-2}$)')
    axes[2].set(title='Litter layer vs Mn concentration',xlabel='Leaf Mn concentration (mg g$^{-1}$)',ylabel='C stock (kg m$^{-2}$)')



    fig,axes=pyplot.subplots(3,1,num='Mn leakage amounts',clear=True)
    axes[0].set_title('Lignin remaining')
    axes[0].set_ylabel('Fraction of initial')
    axes[1].set_title('Mn++ concentration')
    axes[1].set_ylabel('Concentration (M)')
    axes[2].set_title('Manganite concentration')
    axes[2].set_ylabel('Concentration\n($\mu$mol/cm$^{-3}$)')
    axes[2].set_xlabel('Time (days)')
    # Range of leaf Mn concentrations
    norm=matplotlib.colors.Normalize(-8,-3)
    for Mnleakage in arange(-8.0,-3.0) :
        result,units=decomp_network.PF_network_writer(make_network(Mn_peroxidase_Mn3_leakage=10.0**Mnleakage,leaf_Mn_mgkg=10.0)).run_simulation('SOMdecomp_template.txt','manganese',pflotran_exe,print_output=False,length_days=simlength)
        axes[0].plot(result['Lignin']/result['Lignin'].iloc[0],label='Mn+++ leakage = 10$^{%d}$'%int(Mnleakage),c=cm(norm(Mnleakage)))
        axes[1].plot(result['Total Mn++'],c=cm(norm(Mnleakage)))
        axes[2].plot(result['Manganite VF']/24.45*1e6,c=cm(norm(Mnleakage)))
        results_all=pandas.concat((results_all,result))

    axes[0].legend()


    fig,axes=pyplot.subplots(3,1,num='Birnessite rates',clear=True)
    axes[0].set_title('Lignin remaining')
    axes[0].set_ylabel('Fraction of initial')
    axes[1].set_title('Mn++ concentration')
    axes[1].set_ylabel('Concentration (M)')
    axes[2].set_title('Manganite concentration')
    axes[2].set_ylabel('Concentration\n($\mu$mol/cm$^{-3}$)')
    axes[2].set_xlabel('Time (days)')
    # Range of leaf Mn concentrations
    norm=matplotlib.colors.Normalize(-16,-10)
    for rate in arange(-16.0,-9+0.1,2) :
        network=make_network(leaf_Mn_mgkg=1.0,Mn_peroxidase_Mn3_leakage=1e-3,change_rate={'Birnessite':'1.d%d mol/m^2-sec'%rate},change_constraints={'Birnessite':'1.0d-5 1.d2 m^2/m^3','Manganite':'0.0973d-5 1.d2 m^2/m^3'})
        result,units=decomp_network.PF_network_writer(network).run_simulation('SOMdecomp_template.txt','manganese',pflotran_exe,print_output=False,length_days=simlength)
        axes[0].plot(result['Lignin']/result['Lignin'].iloc[0],label='Dissolution rate = 10$^{%d}$ mol m$^{-2}$ s$^{-1}$'%int(rate),c=cm(norm(rate)))
        axes[1].plot(result['Total Mn++'],c=cm(norm(rate)))
        manganite=axes[2].plot(result['Manganite VF']/24.45*1e6,c=cm(norm(rate)),label='Manganite')[0]
        # birnessite=axes[2].plot(result['Birnessite VF']/251.17*1e6,c=cm(norm(rate)),ls='--',label='Birnessite')[0]
        results_all=pandas.concat((results_all,result))

    axes[0].legend()
    axes[2].legend(handles=[manganite])

    # fig,axes=pyplot.subplots(1,1,num='Mn2 effect',clear=True)
    # dlignin=results_all['Lignin'].diff()/(results_all.index[1]-results_all.index[0])
    # dlignin=dlignin.mask(dlignin>0)
    # axes.plot(results_all['Free Mn++'],-dlignin,'o')
    # axes.set_title('Lignin degradation rate vs Mn++')
    # axes.set_ylabel('Lignin degradataion rate (day$^{-1}$)')
    # axes.set_xlabel('Mn++ concentration (M)')
    # axes.set_xscale('linear')

    pyplot.show()

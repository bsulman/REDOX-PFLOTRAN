import run_alquimia
#I'm not sure why but if you import run_alquimia after decomp_network it gets messed up with some compiler library
import decomp_network
# Importing matplotlib fails if you don't module load anaconda3 first
from matplotlib import pyplot
import numpy

pools = [
decomp_network.decomp_pool(name='SOM',CN=50,constraints={'initial':1e1},kind='immobile'),
# decomp_network.decomp_pool(name='Lignin',CN=50,constraints={'initial':1e2},kind='immobile'),
decomp_network.decomp_pool(name='HRimm',constraints={'initial':1e-20},kind='immobile'),

decomp_network.decomp_pool(name='DOM1',CN=50,constraints={'initial':1e-15},kind='primary'),
# decomp_network.decomp_pool(name='DOM2',CN=50,constraints={'initial':1e-30},kind='primary'),
decomp_network.decomp_pool(name='H+',kind='primary',constraints={'initial':'6.0 P'}),
decomp_network.decomp_pool(name='O2(aq)',kind='primary',constraints={'initial':'0.2 G O2(g)'}),
decomp_network.decomp_pool(name='HCO3-',kind='primary',constraints={'initial':'400e-6 G CO2(g)'}),
# decomp_network.decomp_pool(name='Mn+++',kind='primary',constraints={'initial':'1.0e-30 M Manganite'}),
# decomp_network.decomp_pool(name='Mn++',kind='primary',constraints={'initial':'1.0e-30'}),
decomp_network.decomp_pool(name='Fe+++',kind='primary',constraints={'initial':'.37e-10 M Fe(OH)3'}),
decomp_network.decomp_pool(name='Fe++',kind='primary',constraints={'initial':'0.37e-15'}),
decomp_network.decomp_pool(name='NH4+',kind='primary',constraints={'initial':1e-15}), # SOMDecomp sandbox requires this
decomp_network.decomp_pool(name='NO3-',kind='primary',constraints={'initial':1e-5}), 
decomp_network.decomp_pool(name='SO4--',kind='primary',constraints={'initial':1e-5}), 
decomp_network.decomp_pool(name='Tracer',kind='primary',constraints={'initial':1e-15}), # Just to accumulate CO2 loss
decomp_network.decomp_pool(name='CH4(aq)',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='H2S(aq)',kind='secondary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='Acetate-',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='H2(aq)',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='N2(aq)',kind='primary',constraints={'initial':1e-15}),

decomp_network.decomp_pool(name='Na+',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='Cl-',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='Ca++',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='HS-',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='N2O(aq)',kind='primary',constraints={'initial':1e-15}),

decomp_network.decomp_pool(name='CO2(g)',kind='gas'),
decomp_network.decomp_pool(name='O2(g)',kind='gas'),
decomp_network.decomp_pool(name='N2(g)*',kind='gas'),
# decomp_network.decomp_pool(name='N2O(g)',kind='gas'),

decomp_network.decomp_pool(name='CO2(aq)',kind='secondary'),
decomp_network.decomp_pool(name='OH-',kind='secondary'),
# decomp_network.decomp_pool(name='MnO4--',kind='secondary'),
decomp_network.decomp_pool(name='Acetic_acid(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='MnIIIDOM2(aq)',kind='secondary'),
decomp_network.decomp_pool(name='FeIIIDOM1(aq)',kind='secondary'),
decomp_network.decomp_pool(name='FeIIIAcetate(aq)',kind='secondary'),

# decomp_network.decomp_pool(name='NH3(aq)',kind='secondary'),
decomp_network.decomp_pool(name='CO3--',kind='secondary'),
# decomp_network.decomp_pool(name='NH4SO4-',kind='secondary'),
# decomp_network.decomp_pool(name='Urea(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='HSO4-',kind='secondary'),
# decomp_network.decomp_pool(name='H2SO4(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='HNO3(aq)    ',kind='secondary'),
# decomp_network.decomp_pool(name='NaNO3(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='NaCl(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='NaSO4-',kind='secondary'),
# decomp_network.decomp_pool(name='NaCO3-',kind='secondary'),
# decomp_network.decomp_pool(name='NaHCO3(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='HCl(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='CaCO3(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='CaCl+',kind='secondary'),
# decomp_network.decomp_pool(name='CaCl2(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='CaHCO3+',kind='secondary'),
# decomp_network.decomp_pool(name='CaSO4(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='CO(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='CO2(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='Acetic_acid(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='S--',kind='secondary'),

# Hui/Beth say pyrolusite probably not precipitating in soils. Look for delta-MnO2? Maybe try 'Mn(OH)2(am)'?
# See Roberts Earth Sci Rev article for some discussion of minerals? More Fe focused though
# decomp_network.decomp_pool(name='Birnessite',rate='0.d-16 mol/m^2-sec',constraints={'initial':'0.0d-5  1.d2 m^2/m^3'},kind='mineral'),
# decomp_network.decomp_pool(name='Mn(OH)2(am)',rate='0.d-16 mol/m^2-sec',constraints={'initial':'0.0d-5  1.d2 m^2/m^3'},kind='mineral'),
# decomp_network.decomp_pool(name='Manganite',rate='1.d-12 mol/m^2-sec',constraints={'initial':'0.0d-5  1.d2 m^2/m^3'},kind='mineral'),
decomp_network.decomp_pool(name='Fe(OH)3',rate='0.d-3 mol/m^2-sec',constraints={'initial':'0.1d-3  1.d2 m^2/m^3'},kind='mineral'),
decomp_network.decomp_pool(name='Pyrite',rate='0.d-3 mol/m^2-sec',constraints={'initial':'0.0d-3  1.d2 m^2/m^3'},kind='mineral'),
decomp_network.decomp_pool(name='Calcite',rate='0.d-3 mol/m^2-sec',constraints={'initial':'0.875d-3  1.d2 m^2/m^3'},kind='mineral'),
decomp_network.decomp_pool(name='Pyrrhotite',rate='0.d-3 mol/m^2-sec',constraints={'initial':'0.0d-3  1.d2 m^2/m^3'},kind='mineral'),

# This is a workaround to allow pH buffering associated with solid substrate (Rock is a placeholder)
# decomp_network.decomp_pool(name='Rock(s)',rate='0.0 mol/m^2-sec',constraints={'initial':'0.5  5.0e3 m^2/m^3'},kind='mineral'),
# decomp_network.decomp_pool(name='>Carboxylate-',kind='surf_complex',mineral='Rock(s)',site_density=0.0e4,complexes=['>Carboxylic_acid']),

decomp_network.decomp_pool(name='H2O',kind='implicit'),
]



conc_scales={
    'DOM1':1e-1,
    'Acetate-':1e-1,
    'Fe+++':1e-10,
    'NH4+':1e-6,
    'HCO3-':1e-6,
    'NO3-':1e-6,
    'SO4--':1e-6,
    'CH4(aq)':1e-6,
    'H2(aq)':1e-6,
    'O2(aq)':1e-4,
}

rate_scale=10e-8
truncate_conc=1e-30
thresh=truncate_conc*0.0
reactions = [
    # decomp_network.reaction(name='Aerobic decomposition',reactant_pools={'SOM':1.0},product_pools={'HCO3-':1.0},reactiontype='SOMDECOMP',
    #                                         rate_constant=1e-8,rate_units='1/sec',turnover_name='RATE_CONSTANT', 
    #                                         # inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=6.25e-8,type='THRESHOLD 1.0d20')]),
    #                                     monod_terms=[decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=thresh)]),
                                        
                                        
    decomp_network.reaction(name='Hydrolysis',stoich='1.0 SOM -> 1.0 DOM1',reactiontype='SOMDECOMP',turnover_name='RATE_CONSTANT',
                                            rate_constant=rate_scale,rate_units='1/sec', #  Jianqiu Zheng et al., 2019: One third of fermented C is converted to CO2
                                        inhibition_terms=[decomp_network.inhibition(species='DOM1',type='MONOD',k=conc_scales['DOM1']),
                                                          # decomp_network.inhibition(species='O2(aq)',type='MONOD',k=1e-11),
                                                          ]),
    
    # Calculating these as per unit carbon, so dividing by the 6 carbons in a glucose
    # C6H12O6 + 4 H2O -> 2 CH3COO- + 2 HCO3- + 4 H+ + 4 H2
    # Should it be inhibited by H2?
    decomp_network.reaction(name='fermentation',reactant_pools={'DOM1':6/6},product_pools={'Acetate-':2/6,'HCO3-':2/6,'H+':4/6,'H2(aq)':4*2/6,'Tracer':2/6}, # balancing pH of FeIII release requires an extra 5.5 H+ to be released here
                                            rate_constant=rate_scale,reactiontype='MICROBIAL', #  Jianqiu Zheng et al., 2019: One third of fermented C is converted to CO2
                                        inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD'),decomp_network.inhibition(species='Acetate-',k=conc_scales['Acetate-'],type='MONOD')],
                                        monod_terms=[decomp_network.monod(species='DOM1',k=conc_scales['DOM1'],threshold=thresh)]),

    # CH2O + H2O -> CO2 + 4H+ + 4 e-
    # O2   + 4H+ + 4 e- -> 2H2O
    decomp_network.reaction(name='DOM oxidation (O2)',stoich='1.0 DOM1 + 1.0 O2(aq) -> 1.0 HCO3- + 1.0 H+ + 1.0 Tracer',
                                            monod_terms=[decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=thresh),decomp_network.monod(species='DOM1',k=conc_scales['DOM1'],threshold=thresh)],
                                        rate_constant=rate_scale,reactiontype='MICROBIAL'),

    # C2H3O2- + 2 H2O -> 2 CO2 + 7 H+ + 8 e-
    # 2 O2    + 8 H+ + 8 e- -> 4 H2O
    decomp_network.reaction(name='Acetate oxidation (O2)',stoich='1.0 Acetate-  + 2.0 O2(aq)  -> 2.0 HCO3-  + 2.0 H+  + 2.0 Tracer ',
                                            monod_terms=[decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=thresh),decomp_network.monod(species='Acetate-',k=conc_scales['Acetate-'],threshold=thresh)],
                                        rate_constant=rate_scale,reactiontype='MICROBIAL'),


    # C2H3O2- + 2 H2O -> 2 CO2 + 7 H+ + 8 e-
    # 8 Fe+++ + 8 e- -> 8 Fe++ 
    decomp_network.reaction(name='Fe(III) reduction',stoich='1.0 Acetate- + 8.0 Fe+++ -> 2.0 HCO3- + 8.0 Fe++ + 9.0 H+ + 2.0 Tracer',
                                            monod_terms=[decomp_network.monod(species='Acetate-',k=conc_scales['Acetate-'],threshold=thresh),decomp_network.monod(species='Fe+++',k=conc_scales['Fe+++'],threshold=thresh)],
                                            inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD'),decomp_network.inhibition(species='NO3-',k=conc_scales['NO3-'],type='MONOD')],
                                            rate_constant=rate_scale*0.5,reactiontype='MICROBIAL'),
                                            
    # Oxidation of Fe++
    # Currently assuming backward reaction rate is zero
    # decomp_network.reaction(name='Fe(II) oxidation',stoich='1.0 Fe++ + 0.25 O2(aq) + 1.0 H+ <-> 1.0 Fe+++ + 0.5 H2O',
    #                                         rate_constant=1e-2,backward_rate_constant=0.0,reactiontype='GENERAL'),

    # Acetoclastic methanogenesis
    # C2H3O2- + H+ -> CH4 + HCO3- + H+
    decomp_network.reaction(name='Acetoclastic methanogenesis',stoich='1.0 Acetate- -> 1.0 CH4(aq) + 1.0 HCO3- + 1.0 Tracer',
                                            monod_terms=[decomp_network.monod(species='Acetate-',k=conc_scales['Acetate-'],threshold=thresh),],
                                            inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD'),
                                            decomp_network.inhibition(species='Fe+++',k=conc_scales['Fe+++'],type='MONOD'),
                                            decomp_network.inhibition(species='NO3-',k=conc_scales['NO3-'],type='MONOD'),
                                            decomp_network.inhibition(species='SO4--',k=conc_scales['SO4--'],type='MONOD')],
                                            rate_constant=rate_scale*0.1,reactiontype='MICROBIAL'),
                                            
    # Hydrogenotrophic methanogenesis
    decomp_network.reaction(name='Hydrogenotrophic methanogenesis',stoich='4.0 H2(aq) + 1.0 HCO3- + 1.0 H+ -> 1.0 CH4(aq) + 3.0 H2O',
                                            monod_terms=[decomp_network.monod(species='H2(aq)',k=conc_scales['H2(aq)'],threshold=thresh),decomp_network.monod(species='HCO3-',k=conc_scales['HCO3-'],threshold=thresh)],
                                            # inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD'),decomp_network.inhibition(species='Fe+++',k=conc_scales['Fe+++'],type='MONOD')],
                                            inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD'),
                                            decomp_network.inhibition(species='Fe+++',k=conc_scales['Fe+++'],type='MONOD'),
                                            decomp_network.inhibition(species='NO3-',k=conc_scales['NO3-'],type='MONOD'),
                                            decomp_network.inhibition(species='SO4--',k=conc_scales['SO4--'],type='MONOD')],
                                            rate_constant=rate_scale*0.1,reactiontype='MICROBIAL'),
    
    # H2 oxidation if oxygen available

    # Nitrification
    decomp_network.reaction(name='Nitrification',stoich='2.0 NH4+ + 4.0 O2(aq) -> 2.0 NO3- + 2.0 H2O + 4.0 H+',
                        monod_terms=[decomp_network.monod(species='NH4+',k=conc_scales['NH4+'],threshold=thresh),decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=thresh)],
                        rate_constant=rate_scale,reactiontype='MICROBIAL'),
                        
    # Denitrification C2H3O2- + 2 NO3- + H+ -> 2 HCO3- + N2 + 2 H2O
    decomp_network.reaction(name='Denitrification',stoich='1.0 Acetate- + 2.0 NO3- + 1.0 H+ -> 2.0 HCO3- + 1.0 N2(aq) + 0.0 N2O(aq) + 2.0 H2O + 2.0 Tracer',
                                            monod_terms=[decomp_network.monod(species='Acetate-',k=conc_scales['Acetate-'],threshold=thresh),decomp_network.monod(species='NO3-',k=conc_scales['NO3-'],threshold=thresh)],
                                            inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD')],
                                            rate_constant=rate_scale*0.8,reactiontype='MICROBIAL'),
                                            
    # Sulfate reduction C2H3O2- + SO4--  -> 2 HCO3- + HS-
    decomp_network.reaction(name='Sulfate reduction',stoich='1.0 Acetate- + 1.0 SO4--  -> 2.0 HCO3- + 2.0 HS- +  2.0 Tracer',
                                            monod_terms=[decomp_network.monod(species='Acetate-',k=conc_scales['Acetate-'],threshold=thresh),decomp_network.monod(species='SO4--',k=conc_scales['SO4--'],threshold=thresh)],
                                            inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD'),
                                                                decomp_network.inhibition(species='NO3-',k=conc_scales['NO3-'],type='MONOD'),
                                                                decomp_network.inhibition(species='Fe+++',k=conc_scales['Fe+++'],type='MONOD')],
                                            rate_constant=rate_scale*0.3,reactiontype='MICROBIAL'),
                                            
                                            
    # Methane oxidation (O2)
    decomp_network.reaction(name='Methane oxidation (O2)',stoich='1.0 CH4(aq)  + 1.0 O2(aq)  -> 1.0 HCO3-  + 1.0 H+ + 1.0 H2O + 1.0 Tracer',
                                            monod_terms=[decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=thresh),decomp_network.monod(species='CH4(aq)',k=conc_scales['CH4(aq)'],threshold=thresh)],
                                        rate_constant=rate_scale,reactiontype='MICROBIAL'),
    
    # Methane oxidation (NO3)
    decomp_network.reaction(name='Methane oxidation (NO3)',stoich='1.0 CH4(aq)  + 1.0 NO3-  -> 1.0 HCO3-  + 1.0 NH4+ + 1.0 Tracer ',
                                            monod_terms=[decomp_network.monod(species='NO3-',k=conc_scales['NO3-'],threshold=thresh),decomp_network.monod(species='CH4(aq)',k=conc_scales['CH4(aq)'],threshold=thresh)],
                                        rate_constant=rate_scale,reactiontype='MICROBIAL'),
    
    # Methane oxidation (SO4)
    decomp_network.reaction(name='Methane oxidation (SO4)',stoich='1.0 CH4(aq)  + 1.0 SO4-- -> 1.0 HCO3-  + 1.0 HS- + 1.0 H2O + 1.0 Tracer ',
                                            monod_terms=[decomp_network.monod(species='SO4--',k=conc_scales['SO4--'],threshold=thresh),decomp_network.monod(species='CH4(aq)',k=conc_scales['CH4(aq)'],threshold=thresh)],
                                        rate_constant=rate_scale,reactiontype='MICROBIAL'),
    
    # Methane oxidation (Fe)
    decomp_network.reaction(name='Methane oxidation (Fe)',stoich='1.0 CH4(aq)  + 8.0 Fe+++ + 3.0 H2O -> 1.0 HCO3-  + 8.0 Fe++ + 9.0 H+ + 1.0 Tracer ',
                                            monod_terms=[decomp_network.monod(species='Fe+++',k=conc_scales['Fe+++'],threshold=thresh),decomp_network.monod(species='CH4(aq)',k=conc_scales['CH4(aq)'],threshold=thresh)],
                                        rate_constant=rate_scale,reactiontype='MICROBIAL'),

]

def plot_result(result,SOM_ax=None,pH_ax=None,Fe_ax=None,gasflux_ax=None,porewater_ax=None,do_legend=False):
    t=result.index
    if SOM_ax is not None:
        l=SOM_ax.plot(t,result['SOM']*1e-3,label='SOM')[0]

        SOM_ax.set_title('SOM remaining')
        SOM_ax.set_ylabel('Concentration\n(mmol C/cm$^{-3}$)')
        SOM_ax.set_xlabel('Time (days)')

    if pH_ax is not None:
        pH_ax.plot(t,-numpy.log10(result['Free H+']))

        pH_ax.set_title('pH')
        pH_ax.set_ylabel('pH')
        pH_ax.set_xlabel('Time (days)')
        
    if Fe_ax is not None:
        molar_volume=34.3600 # From database. cm3/mol
        molar_weight = 106.8690
        l=Fe_ax.plot(t,result['Fe(OH)3 VF']/molar_volume*1e6   ,label='Fe(OH)3')[0]
        
        # M/L to umol/cm3: 1e6/1e3=1e3
        l=Fe_ax.plot(t,result['Total Fe+++']*result['Porosity']*1e3   ,label='Fe+++',ls='--')[0]
        
        l=Fe_ax.plot(t,result['Total Fe++']*result['Porosity']*1e3 ,ls=':'  ,label='Fe++')[0]
        
        Fe_ax.set_title('Fe species')
        Fe_ax.set_ylabel('Concentration\n($\mu$mol/cm$^{-3}$)')
        Fe_ax.set_xlabel('Time (days)')
        if do_legend:
            Fe_ax.legend(fontsize='small')
    
    if gasflux_ax is not None:
        gasflux_ax.set_yscale('log')
        
        l=gasflux_ax.plot(t.values[:-1],numpy.diff(result['Total CH4(aq)']*result['Porosity'])/numpy.diff(result.index.values)*1e3,label='CH4')[0]
        
        l=gasflux_ax.plot(t.values[:-1],numpy.diff(result['Total Tracer']*result['Porosity'])/numpy.diff(result.index.values)*1e3,label='CO2',ls='--',c='C5')[0]

        gasflux_ax.set_title('Gas fluxes')
        gasflux_ax.set_ylabel('Flux rate\n($\mu$mol cm$^{-3}$ day$^{-1}$)')
        gasflux_ax.set_xlabel('Time (days)')
        if do_legend:
            gasflux_ax.legend(fontsize='small')
        
    if porewater_ax is not None:
        porewater_ax.set_yscale('log')
        # porewater_ax.plot(result['Total DOM1'],label='DOM')
        # porewater_ax.plot(result['Total Acetate-'],label='Acetate',c='C3')
        porewater_ax.plot(t,result['Total O2(aq)'],'-',label='O2',c='C0')
        porewater_ax.plot(t,result['Total NO3-'],'-',c='C1',label='NO3-')
        porewater_ax.plot(t,result['Total N2(aq)'],'--',c='C1',label='N2')
        porewater_ax.plot(t,result['Total Fe+++'],'-',label='Fe+++',c='C2')
        # porewater_ax.plot(result['Free Fe+++'],':',label='Fe+++',c='C1')
        porewater_ax.plot(t,result['Total Fe++'],'--',label='Fe++',c='C2')

        porewater_ax.plot(t,result['Total SO4--'],'-',c='C3',label='SO4--')
        porewater_ax.plot(t,result['Total HS-'],'--',c='C3',label='H2S')
        porewater_ax.plot(t,result['Total CH4(aq)'],'--',c='C4',label='CH4')
        
        porewater_ax.set_title('Porewater concentrations')
        porewater_ax.set_ylabel('Concentration (M)')
        porewater_ax.set_xlabel('Time (days)')
        
        if do_legend:
            porewater_ax.legend(fontsize='small')


# All of this stuff is for the plot of the network diagram
mineral_col=265.0
TEA_col=175.0
reduction_col=150.0
oxidation_col=200.0
substrateox_col=reduction_col
substrate_col=130.0
gas_col=TEA_col
redox_seq={
    'O2':1250,
    'NO3-':900,
    # 'Mn4+':700,
    'Fe3+':600,
    'SO4--':400,
    'CO2':200,
    'CH4':100,
    'SOM':1400,
    'DOM':1200,
    'Acetate':900,
    'H2':600
    }
pos={'Pyrrhotite': (mineral_col, (redox_seq['Fe3+']+redox_seq['SO4--'])/2),
 'Pyrite': (mineral_col, (redox_seq['Fe3+']+redox_seq['SO4--']+redox_seq['SO4--']-100)/3),
 'SO4--': (TEA_col, redox_seq['SO4--']),
 'HS-': (TEA_col, redox_seq['SO4--']-75),
 'H2S(aq)': (gas_col, redox_seq['SO4--']),
 'Sulfate reduction': (reduction_col, (redox_seq['SO4--'])),
 'Methane oxidation (SO4)': (oxidation_col, redox_seq['Fe3+']),
 
 'Fe(OH)3': (mineral_col, redox_seq['Fe3+']),
 'Fe+++': (TEA_col, redox_seq['Fe3+']),
 'Fe++': (TEA_col, redox_seq['Fe3+']-75),
 'Fe(II) oxidation': (oxidation_col, (redox_seq['Fe3+']+redox_seq['O2'])/2),
 'Fe(III) reduction': (reduction_col, redox_seq['Fe3+']),
 'Methane oxidation (Fe)': (oxidation_col, redox_seq['Fe3+']),
 
 'NH4+': (TEA_col, redox_seq['NO3-']+100),
 'NO3-': (TEA_col, redox_seq['NO3-']),
 'Denitrification': (reduction_col, redox_seq['NO3-']-40),
 'Nitrification': (substrateox_col, redox_seq['NO3-']+75/2),
 'N2(aq)': (gas_col, redox_seq['NO3-']-150),
 'N2O(aq)': (gas_col, redox_seq['NO3-']-75),
 'Methane oxidation (NO3)': (oxidation_col, (redox_seq['Fe3+'])),
 
 'SOM': (substrate_col, redox_seq['SOM']),
 'Hydrolysis': (reduction_col, (redox_seq['SOM'])),
 'DOM1': (substrate_col, redox_seq['DOM']),
 'fermentation': (reduction_col, (redox_seq['DOM']+redox_seq['Acetate'])/2),
 'Acetate-': (substrate_col, redox_seq['Acetate']),
 'H2(aq)': (substrate_col, redox_seq['H2']),
 
 'O2(aq)': (TEA_col, redox_seq['O2']),
 'HCO3-': (175, redox_seq['CO2']),

 'CH4(aq)': (TEA_col, redox_seq['CH4']),

 
 'DOM oxidation (O2)': (substrateox_col, (redox_seq['O2']+redox_seq['Acetate'])/2+155),
 'Aerobic decomposition': (substrateox_col, (redox_seq['O2']+redox_seq['Acetate'])/2+155),
 'Acetate oxidation (O2)': (substrateox_col, (redox_seq['O2']+redox_seq['Acetate'])/2+155),
 'Methane oxidation (O2)': (oxidation_col, (redox_seq['Fe3+'])),

 'Hydrogenotrophic methanogenesis': (reduction_col, (redox_seq['CO2'])-100),
 'Acetoclastic methanogenesis': (reduction_col, (redox_seq['CO2']+100-100))}

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
    'General Reaction':'#825708',
    'Microbial Reaction':'#e6990b',
    'SOMdecomp Reaction':'#c7a82c',
}


reaction_network=decomp_network.decomp_network(pools=pools,reactions=reactions)

pflotran_exe='../pflotran-interface/src/pflotran/pflotran'
simlength=30

# Run using standard PFLOTRAN and generated input file
precision=3
result,units=decomp_network.PF_network_writer(reaction_network,precision=precision).run_simulation('SOMdecomp_template.txt','saline',pflotran_exe,print_output=False,length_days=simlength,obs_time_hrs=1.0,truncate_concentration=truncate_conc,log_formulation=True)

# Run using python interface to alquimia
#rate_scale=0.0
# decomp_network.PF_network_writer(reaction_network).write_into_input_deck('SOMdecomp_template.txt','saline_alquimia.in',length_days=simlength,log_formulation=True,truncate_concentration=truncate_conc)
rateconstants={
'fermentation':rate_scale,
'DOM oxidation (O2)':rate_scale,
'Hydrogenotrophic methanogenesis':rate_scale*0.1,
'Denitrification':rate_scale*0.8,
'Acetate oxidation (O2)':rate_scale,
'Methane oxidation (O2)':rate_scale,
'Nitrification':rate_scale,
'Fe(III) reduction':rate_scale*0.5,
'Methane oxidation (Fe)':rate_scale,
'Methane oxidation (NO3)':rate_scale,
'Sulfate reduction':rate_scale*0.3,
'Methane oxidation (SO4)':rate_scale,
'Acetoclastic methanogenesis':rate_scale*0.1,
# 'Aerobic decomposition': 1e-6,
'Hydrolysis':rate_scale,
}

rateconstants_stoich=run_alquimia.convert_rateconstants(rateconstants,reactions,precision=precision)
result_alquimia,units_alq=run_alquimia.run_simulation('saline_generated.in',hands_off=False,simlength_days=simlength,dt=3600,initcond=pools,rateconstants=rateconstants_stoich,truncate_concentration=truncate_conc)


networkfig=pyplot.figure('Reaction network',clear=True,figsize=(7.1,8.9))

drawn,pos=decomp_network.draw_network_with_reactions(reaction_network,
        omit=['gas','secondary','H+','>Carboxylate-','Carboxylic_acid','H2O','Na+','Cl-','Calcite','Ca++','Nitrification','NH4+','mineral','HCO3-'],
        font_size='xx-large',node_size=1800,font_color='w',arrowstyle='-|>',arrowsize=20.0,width=1.5,edge_color='k',node_alpha=0.9,node_colors=node_colors,markers={'Reaction':'h','mineral':'8'},
        namechanges={'SOM':'SOM','DOM1':'DOM','O2(aq)':'O$_2$','CH4(aq)':'CH$_4$','HCO3-':'CO$_2$','DOM2':'Exposed lignin','sorbed_DOM1':'Sorbed DOM',
                     'Fe(OH)2':'Fe(OH)$_2$','Fe(OH)3':'Fe(OH)$_3$','Mn++':r'Mn$^\mathrm{+\!\!+}$','Mn+++':r'Mn$^\mathrm{+\!\!+\!\!\!+}$','Acetate-':'Acetate',
                     'H2(aq)':'H$_2$','DOM oxidation (O2)':'Aerobic\nresp','Acetate oxidation (O2)':'Aerobic\nresp','Aerobic decomposition':'Aerobic\nresp',
                     'NH4+':'NH$_4^+$','NO3-':'NO$_3^-$','N2O(aq)':'N$_2$O','N2(aq)':'N$_2$','Fe++':r'Fe$^\mathrm{+\!\!+}$','Fe+++':r'Fe$^\mathrm{+\!\!+\!\!\!+}$',
                     'SO4--':'SO$_4^{--}$','HS-':'H$_2$S','Methane oxidation (O2)':'Methane\noxidation','Methane oxidation (NO3)':'Methane\noxidation',
                     'Methane oxidation (Fe)':'Methane\noxidation','Methane oxidation (SO4)':'Methane\noxidation','Fe(III) reduction':'Fe(III)\nreduction','Sulfate reduction':'Sulfate\nreduction',
                     'Hydrogenotrophic methanogenesis':'Hydrogenotrophic\nmethanogenesis','Acetoclastic methanogenesis':'Acetoclastic\nmethanogenesis',
                     'SOMdecomp Reaction':'SOM Reaction','General Reaction':'Abiotic Reaction','fermentation':'Fermentation','Primary aqueous':'Dissolved ion','Gas':'Dissolved gas'},pos=pos,connectionstyle='arc3, rad=0.2')


networkfig.axes[0].set_facecolor('#8a9ebf')

for p in networkfig.axes[0].patches:
    ex=p.get_path().get_extents()
    if (ex.max[0]<155 or ex.max[0]>188) and ex.min[1]>300:
        p.set_connectionstyle('arc3,rad=-0.1')

resultsfig,axs=pyplot.subplots(num='Results',clear=True,squeeze=False,nrows=3)
# plot_result(result,porewater_ax=axs[1,0],SOM_ax=axs[0,0],do_legend=True)
# axs[1,0].set_ylim(bottom=1e-13)
# axs[0,0].legend(fontsize='medium',ncol=2)
t=result.index
t_alquimia=result_alquimia.index
axs[0,0].plot(t,result['Total O2(aq)']/result['Total O2(aq)'].iloc[0],label='O$_2$')
axs[0,0].plot(t,result['Total NO3-']/(result['Total NO3-']+result['Total N2(aq)']),label='NO$_3^-$')
axs[0,0].plot(t,result['Total Fe+++']/(result['Total Fe+++']+result['Total Fe++']),label='Fe(III)')
axs[0,0].plot(t,result['Total SO4--']/(result['Total SO4--']+result['Total HS-']),label='SO$_4^{--}$')
axs[0,0].plot(t,result['Total CH4(aq)']/result['Total CH4(aq)'].iloc[-1],label='CH$_4$')


axs[0,0].plot(t_alquimia,result_alquimia['Total O2(aq)']/result_alquimia['Total O2(aq)'].iloc[0],ls='--',c='C0')
axs[0,0].plot(t_alquimia,result_alquimia['Total NO3-']/(result_alquimia['Total NO3-']+result_alquimia['Total N2(aq)']),ls='--',c='C1')
axs[0,0].plot(t_alquimia,result_alquimia['Total Fe+++']/(result_alquimia['Total Fe+++']+result_alquimia['Total Fe++']),ls='--',c='C2')
axs[0,0].plot(t_alquimia,result_alquimia['Total SO4--']/(result_alquimia['Total SO4--']+result_alquimia['Total HS-']),ls='--',c='C3')
axs[0,0].plot(t_alquimia,result_alquimia['Total CH4(aq)']/result_alquimia['Total CH4(aq)'].iloc[-1],'--',c='C4')

axs[0,0].set_xlabel('Time (days)',fontsize='medium')
axs[0,0].set_ylabel('Relative concentration',fontsize='medium')
axs[0,0].set_title('Simulated terminal electron acceptors',fontsize='large')
axs[0,0].legend(fontsize='medium',ncol=2)

axs[1,0].plot(t,(result_alquimia['Total O2(aq)']).values-(result['Total O2(aq)']).values,label='O$_2$')
axs[1,0].plot(t,(result_alquimia['Total NO3-']).values-(result['Total NO3-']).values,label='NO$_3^-$')
axs[1,0].plot(t,(result_alquimia['Total Fe+++']).values-(result['Total Fe+++']).values,label='Fe(III)')
axs[1,0].plot(t,(result_alquimia['Total SO4--']).values-(result['Total SO4--']).values,label='SO$_4^{--}$')
axs[1,0].plot(t,(result_alquimia['Total CH4(aq)']).values-(result['Total CH4(aq)']).values,label='CH$_4$')


axs[1,0].set_xlabel('Time (days)',fontsize='medium')
axs[1,0].set_ylabel('Concentration (M)',fontsize='medium')
axs[1,0].set_title('Alquimia minus PFLOTRAN',fontsize='large')
# axs[1,0].legend(fontsize='medium',ncol=2)

axs[2,0].plot(t,result['Total O2(aq)'],label='O$_2$')
axs[2,0].plot(t,result['Total NO3-'],label='NO$_3^-$')
axs[2,0].plot(t,result['Total Fe+++'],label='Fe(III)')
axs[2,0].plot(t,result['Total SO4--'],label='SO$_4^{--}$')
axs[2,0].plot(t,result['Total CH4(aq)'],label='CH$_4$')


axs[2,0].plot(t_alquimia,result_alquimia['Total O2(aq)'],ls='--',c='C0')
axs[2,0].plot(t_alquimia,result_alquimia['Total NO3-'],ls='--',c='C1')
axs[2,0].plot(t_alquimia,result_alquimia['Total Fe+++'],ls='--',c='C2')
axs[2,0].plot(t_alquimia,result_alquimia['Total SO4--'],ls='--',c='C3')
axs[2,0].plot(t_alquimia,result_alquimia['Total CH4(aq)'],'--',c='C4')

axs[2,0].set_yscale('log')
axs[2,0].legend()
axs[2,0].set_xlabel('Time (days)',fontsize='medium')
axs[2,0].set_ylabel('Concentration (M)',fontsize='medium')
axs[2,0].set_title('Simulated terminal electron acceptors',fontsize='large')

pyplot.show()

import run_alquimia
import decomp_network
import numpy

pools = [
decomp_network.decomp_pool(name='cellulose',CN=50,constraints={'initial':8e3},kind='immobile'),
decomp_network.decomp_pool(name='HRimm',constraints={'initial':1e-20},kind='immobile'),

decomp_network.decomp_pool(name='DOM1',CN=50,constraints={'initial':0.5e-1},kind='primary'),
decomp_network.decomp_pool(name='H+',kind='primary',constraints={'initial':'5.0 P'}),
decomp_network.decomp_pool(name='O2(aq)',kind='primary',constraints={'initial':'0.2 G O2(g)'}),
decomp_network.decomp_pool(name='HCO3-',kind='primary',constraints={'initial':'400e-6 G CO2(g)'}),
decomp_network.decomp_pool(name='Fe+++',kind='primary',constraints={'initial':'.37e-10 M Fe(OH)3'}),
decomp_network.decomp_pool(name='Fe++',kind='primary',constraints={'initial':'0.37e-3'}),
decomp_network.decomp_pool(name='NH4+',kind='primary',constraints={'initial':1e-15}), # SOMDecomp sandbox requires this
decomp_network.decomp_pool(name='Tracer',kind='primary',constraints={'initial':1e-15}), # Just to accumulate CO2 loss
decomp_network.decomp_pool(name='CH4(aq)',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='Acetate-',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='H2(aq)',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='Mg++',kind='primary',constraints={'initial':0.5e-3}),
decomp_network.decomp_pool(name='Ca++',kind='primary',constraints={'initial':0.5e-3}),
decomp_network.decomp_pool(name='Na+',kind='primary',constraints={'initial':2e-3}),
decomp_network.decomp_pool(name='K+',kind='primary',constraints={'initial':2e-5}),

decomp_network.decomp_pool(name='CO2(g)',kind='gas'),
decomp_network.decomp_pool(name='O2(g)',kind='gas'),

decomp_network.decomp_pool(name='CO2(aq)',kind='secondary'),
decomp_network.decomp_pool(name='OH-',kind='secondary'),
decomp_network.decomp_pool(name='FeCO3+',kind='secondary'),
decomp_network.decomp_pool(name='Fe(OH)4-',kind='secondary'),
decomp_network.decomp_pool(name='Acetic_acid(aq)',kind='secondary'),
decomp_network.decomp_pool(name='FeCH3COO+',kind='secondary'),
decomp_network.decomp_pool(name='FeIIIDOM1(aq)',kind='secondary'),
decomp_network.decomp_pool(name='FeIIDOM1(aq)',kind='secondary'),
decomp_network.decomp_pool(name='FeIIIAcetate(aq)',kind='secondary'),
decomp_network.decomp_pool(name='Fe(OH)2(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='Fe(OH)2+',kind='secondary'),
decomp_network.decomp_pool(name='FeCO3(aq)',kind='secondary'),
decomp_network.decomp_pool(name='CO3--',kind='secondary'),
decomp_network.decomp_pool(name='CaHCO3+',kind='secondary'),

# See Roberts Earth Sci Rev article for discussion of iron minerals, including some time scales
decomp_network.decomp_pool(name='Fe(OH)3',rate='1.d-6 mol/m^2-sec',constraints={'initial':'0.875d-3  1.d2 m^2/m^3'},kind='mineral'),
# decomp_network.decomp_pool(name='Goethite',rate='1.d-5 mol/m^2-sec',constraints={'initial':'1.75d-2  1.d1 m^2/m^3'},kind='mineral'),
# decomp_network.decomp_pool(name='Fe',rate='1.d-7 mol/m^2-sec',constraints={'initial':'1.0e-6  1. m^2/m^3'},kind='mineral'),
decomp_network.decomp_pool(name='Fe(OH)2',rate='1.d-7 mol/m^2-sec',constraints={'initial':'0.0e-20  1.d2 m^2/m^3'},kind='mineral'),
decomp_network.decomp_pool(name='Rock(s)',rate='0.0 mol/m^2-sec',constraints={'initial':'0.5  5.0e3 m^2/m^3'},kind='mineral'),
# decomp_network.decomp_pool(name='Calcite',rate='1.d-3 mol/m^2-sec',constraints={'initial':'0.0d-3  1.d2 m^2/m^3'},kind='mineral'),

# Maybe this should be a cation exchange reaction instead?
# Probably makes sense to include both
# There should be some Fe(II) sorption or complexation with OM
decomp_network.decomp_pool(name='>Carboxylate-',kind='surf_complex',mineral='Rock(s)',site_density=2000.0*5.0,complexes=['>Carboxylic_acid']), # site density in mol/m3

decomp_network.decomp_pool(name='H2O',kind='implicit'),
]
             
# Herndon et al 2015, BGC: Fe(III) average concentration 0.37 mmol/L (same as mM). SO4- was 0.07, NO3- was 0.03.
# Fe(III) was 60% of dissolved Fe
# ** What form was it in? Does it just dissolve or was it released/complexed by microbial activity?
# DOC 5-15 mmol/L
# Fe probably present as iron oxide precipitates or Fe-OM complexes, not primary minerals
#
# Herndon et al 2015, JGR-B: Organic acids 12% of WEOC in organic layer. formate 48.6 umol/ g SOC; acetate 66, propionate 0.76
# Each mol of acetate produces 1 mol of methane or 8 mol of Fe(II)
# Ferrihydrite specific surface area ~250 m2/g (Kaiser and Guggenberger 2008 https://doi.org/10.1046/j.1365-2389.2003.00544.x)
# 

# Beth: There should be a lot of solid OM functional groups, probably carboxylate/carboxylic acid that exchange protons and buffer pH in these soils
# Fe(II) is probably sorbing onto OM, or staying dissolved
# Fe(III) should precipitate in this system, so if it's dissolved it would be complexed with something

conc_scales={
    'DOM1':1e-1,
    'Acetate-':0.4e-1,
    'Fe+++':1e-10,
    'O2(aq)':1e-4,
    'Fe++':1e-1,
    'H2(aq)':1e-1,
    'HCO3-':1e-1,
}
anox_inhib_conc = 1e-5

truncate_conc=1e-15

def make_reactions(conc_scales=conc_scales,anox_inhib_conc=anox_inhib_conc,Fe_inhibit_CH4=True,rate_scale=1e-8):

    # Turn off inhibition of methanogenesis by Fe+++ by setting inhibition constant very high
    if Fe_inhibit_CH4==True:
        Fe_CH4_inhib=1.0
    else:
        Fe_CH4_inhib=1e10

    reactions = [
    # decomp_network.reaction(name='Aerobic decomposition',reactant_pools={'cellulose':1.0},product_pools={'HCO3-':1.0},reactiontype='SOMDECOMP',
    #                                         rate_constant=0.0,rate_units='1/d',turnover_name='RATE_DECOMPOSITION', 
    #                                         # inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=6.25e-8,type='THRESHOLD 1.0d20')]),
    #                                     monod_terms=[decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=1.1e-9)]),
    
    # From Dave Graham: SOM-C + H2O -> DOM-C
    decomp_network.reaction(name='Hydrolysis',stoich='1.0 cellulose -> 1.0 DOM1',reactiontype='SOMDECOMP',
                                            rate_constant=rate_scale,rate_units='y', #  Jianqiu Zheng et al., 2019: One third of fermented C is converted to CO2
                                        inhibition_terms=[decomp_network.inhibition(species='DOM1',type='MONOD',k=conc_scales['DOM1'],pool_normalized=True),
                                                          # decomp_network.inhibition(species='O2(aq)',k=6.25e-11,type='THRESHOLD -1.0d10')
                                                        #   decomp_network.inhibition(species='O2(aq)',type='MONOD',k=1e-11),
                                                          ]),
    
    # Calculating these as per unit carbon, so dividing by the 6 carbons in a glucose
    # C6H12O6 + 4 H2O -> 2 CH3COO- + 2 HCO3- + 4 H+ + 4 H2
    # Dave Graham: DOM-C + 0.67 H2O -> 0.33 CH3COO- + 0.33 HCO3- + 0.67 H+ + 0.67 H2
    # Should it be inhibited by H2?
    decomp_network.reaction(name='fermentation',reactant_pools={'DOM1':1.0,'H2O':2/3},product_pools={'Acetate-':1/3,'HCO3-':1/3,'H+':2/3,'H2(aq)':2/3,'Tracer':1/3}, # balancing pH of FeIII release requires an extra 5.5 H+ to be released here
                                            rate_constant=rate_scale*40,reactiontype='MICROBIAL', #  Jianqiu Zheng et al., 2019: One third of fermented C is converted to CO2
                                        inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=anox_inhib_conc,type='MONOD'),
                                                        decomp_network.inhibition(species='Acetate-',k=conc_scales['Acetate-']*0.5,type='MONOD')],
                                        monod_terms=[decomp_network.monod(species='DOM1',k=conc_scales['DOM1']*0.1,threshold=1.1e-15)]),

    # CH2O + H2O -> CO2 + 4H+ + 4 e-
    # O2   + 4H+ + 4 e- -> 2H2O
    decomp_network.reaction(name='DOM aerobic respiration',stoich='1.0 DOM1 + 1.0 O2(aq) + 1.0 H2O -> 1.0 HCO3- + 1.0 H+ + 1.0 Tracer',
                                            monod_terms=[decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=0.0),decomp_network.monod(species='DOM1',k=conc_scales['DOM1'],threshold=1.1e-16)],
                                        rate_constant=rate_scale*50,reactiontype='MICROBIAL'),

    # C2H3O2- + 2 H2O -> 2 CO2 + 7 H+ + 8 e-
    # 2 O2    + 8 H+ + 8 e- -> 4 H2O
    decomp_network.reaction(name='Acetate aerobic respiration',stoich='1.0 Acetate-  + 2.0 O2(aq)  -> 2.0 HCO3-  + 2.0 H+  + 2.0 Tracer ',
                                            monod_terms=[decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=0.0),decomp_network.monod(species='Acetate-',k=conc_scales['Acetate-'],threshold=1.1e-16)],
                                        rate_constant=rate_scale*50,reactiontype='MICROBIAL'),


    # C2H3O2- + 2 H2O -> 2 CO2 + 7 H+ + 8 e-
    # 8 Fe+++ + 8 e- -> 8 Fe++ 
    decomp_network.reaction(name='Fe(III) reduction',stoich='1.0 Acetate- + 8.0 Fe+++ + 4.0 H2O -> 2.0 HCO3- + 8.0 Fe++ + 9.0 H+ + 2.0 Tracer',
                                            monod_terms=[decomp_network.monod(species='Acetate-',k=conc_scales['Acetate-'],threshold=1.1e-15),decomp_network.monod(species='Fe+++',k=conc_scales['Fe+++'],threshold=1.1e-15)],
                                            inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=anox_inhib_conc,type='MONOD')],
                                            rate_constant=rate_scale*3,reactiontype='MICROBIAL'),
                                            
    # Oxidation of Fe++
    # Currently assuming backward reaction rate is zero
    decomp_network.reaction(name='Fe(II) abiotic oxidation',stoich='1.0 Fe++ + 0.25 O2(aq) + 1.0 H+ <-> 1.0 Fe+++ + 0.5 H2O',
                                            rate_constant=0.0,backward_rate_constant=0.0,reactiontype='GENERAL'),
                                            
    decomp_network.reaction(name='Fe(II) microbial oxidation',stoich='1.0 Fe++ + 0.25 O2(aq) + 1.0 H+ -> 1.0 Fe+++ + 0.5 H2O',
                                            monod_terms=[decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=0.0),decomp_network.monod(species='Fe++',k=conc_scales['Fe++'],threshold=1.1e-15)],
                                            rate_constant=rate_scale*200,reactiontype='MICROBIAL'),

    # Acetoclastic methanogenesis
    # C2H3O2- + H+ -> CH4 + HCO3- + H+
    # pH dependence: Dunfield et al 1992 has some bell curves. Kotsyurbenko et al 2007 has info on different pathways at different pH. Le Mer and Roger 2001 is a big review with a bit on pH thresholds
    # See methane_ph.py: Optimization against Kotsyurbenko data yields K_M=5.54 and K_I=5.54 for acetaclastic, and k_M=6.75 for hydrogenotrophic
    #                    Also gives rate constant of hydro = 0.415 that of acetaclastic
    decomp_network.reaction(name='Acetaclastic methanogenesis',stoich='1.0 Acetate- + H2O -> 1.0 CH4(aq) + 1.0 HCO3- + 1.0 Tracer',
                                            monod_terms=[decomp_network.monod(species='Acetate-',k=conc_scales['Acetate-'],threshold=1.1e-15),
                                                         ],
                                            inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=anox_inhib_conc,type='MONOD'),
                                                              decomp_network.inhibition(species='Fe+++',k=conc_scales['Fe+++']*Fe_CH4_inhib,type='MONOD'),
                                                              decomp_network.inhibition(species='H+',k=10**-5.54,type='MONOD'),
                                                              decomp_network.inhibition(species='H+',k=10**-5.54,type='INVERSE_MONOD')],
                                            rate_constant=rate_scale*2,reactiontype='MICROBIAL'),

    # Hydrogenotrophic methanogenesis
    decomp_network.reaction(name='Hydrogenotrophic methanogenesis',stoich='4.0 H2(aq) + 1.0 HCO3- + 1.0 H+ -> 1.0 CH4(aq) + 3.0 H2O',
                                        monod_terms=[decomp_network.monod(species='H2(aq)',k=conc_scales['H2(aq)'],threshold=1.1e-15),decomp_network.monod(species='HCO3-',k=conc_scales['HCO3-'],threshold=1.1e-15)],
                                        # inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD'),decomp_network.inhibition(species='Fe+++',k=conc_scales['Fe+++'],type='MONOD')],
                                        inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=anox_inhib_conc,type='MONOD'),
                                        decomp_network.inhibition(species='Fe+++',k=conc_scales['Fe+++']*Fe_CH4_inhib,type='MONOD'),
                                        decomp_network.inhibition(species='H+',k=10**-6.75,type='MONOD'),
                                        # decomp_network.inhibition(species='NO3-',k=conc_scales['NO3-'],type='MONOD'),
                                        # decomp_network.inhibition(species='SO4--',k=conc_scales['SO4--'],type='MONOD')
                                        ],
                                        rate_constant=rate_scale*3*0.32,reactiontype='MICROBIAL'),
    
    # H2 oxidation if oxygen available
    decomp_network.reaction(name='Hydrogen oxidation',stoich='2.0 H2(aq) + 1.0 O2(aq) -> 2.0 H2O',
                                        monod_terms=[decomp_network.monod(species='H2(aq)',k=conc_scales['H2(aq)'],threshold=1.1e-15),
                                                     decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=0.0)],
                                        rate_constant=rate_scale*200,reactiontype='MICROBIAL'),
    
    # Cation exchange
    decomp_network.ion_exchange(name='Cation exchange',cations={'Fe++':0.1,'Fe+++':0.3,'Mg++':1.1,'Ca++':4.1,'Na+':1.0,'K+':0.9,'H+':1.1},CEC=2000.0,mineral=None)

    ]
    return reactions


# Incubation data
import pandas,glob
datadir='/home/b0u/NGEE-Arctic-redox/Data'
Barrow_synthesis_data=pandas.concat([pandas.read_csv(fname,na_values=-9999,header=8,skiprows=[9]) for fname in glob.glob(datadir+'/Barrow_soil_geochem_synthesis/*.csv')])  
units=pandas.read_csv(glob.glob(datadir+'/Barrow_soil_geochem_synthesis/*.csv')[0],na_values=-9999,header=8).iloc[0]
SOC_layermean=Barrow_synthesis_data.groupby(Barrow_synthesis_data['Soil_layer'].str.capitalize())['SOC'].mean()
SOC_layer_microtopo_mean=Barrow_synthesis_data.groupby([Barrow_synthesis_data['Soil_layer'].str.capitalize(),Barrow_synthesis_data['Microtopography'].str.capitalize()=='Trough'])['SOC'].mean()


# SOC units (%) need to be converted to mol C/m3
# Bulk density reported for some cores. Best bet maybe to find relationship between SOC% and BD and assume it's constant accross these incubations
# Although, might be worth checking if they have some bulk density that I haven't found
corephysdata=pandas.read_csv(datadir+'/Barrow_soilcores2/core_physical_chem_data_20180321.csv',header=5,na_values=-9999,skiprows=[6])
corephysunits=pandas.read_csv(datadir+'/Barrow_soilcores2/core_physical_chem_data_20180321.csv',header=5,na_values=-9999).iloc[0]
from scipy.stats import linregress
BD_SOM_fit=linregress(corephysdata[['Organic_Matter_Content','Dry_Bulk_Density']].mask(corephysdata['Organic_Matter_Content']<6.0).dropna().to_numpy()) 
BD_layerest=BD_SOM_fit.intercept + BD_SOM_fit.slope*SOC_layermean/0.6 # Guessing at conversion from SOC to SOM weight
# Bockheim et al 2003, looks like a broader dataset:
#  They use a fit C[%] = -9.7872 ln(BD[g/cm3])+8.2432
# BD_layerest2=numpy.exp((8.2432-SOC_layermean)/9.7872)
BD_layerest2_trough=numpy.exp((8.2432-SOC_layer_microtopo_mean)/9.7872)

BD_porosity_fit = linregress(corephysdata[['Porosity','Dry_Bulk_Density']].dropna())
BD_porosity_fit_organic = linregress(corephysdata[['Dry_Bulk_Density','Porosity']][corephysdata['Organic_Matter_Content']>20].dropna())
BD_porosity_fit_mineral = linregress(corephysdata[['Dry_Bulk_Density','Porosity']][corephysdata['Organic_Matter_Content']<=20].dropna())

waterchemistry=pandas.read_csv(datadir+'/Barrow_porewater_chem/BGC_BarrowWaterChemistry_2013_2014_v1.csv',header=8,skiprows=[9,10],na_values=-9999)
waterchemistry_units=pandas.read_csv(datadir+'/Barrow_porewater_chem/BGC_BarrowWaterChemistry_2013_2014_v1.csv',header=8,na_values=-9999).iloc[0]

cellulosefrac=0.25
molar_volume_FeOH3=34.3600 #cm3/mol
porosity=0.5

def get_layer(Core_ID,layer,minT=4.0):
    if isinstance(Core_ID,int):
        Core_ID='NGADG%04d'%Core_ID
    data_layer=Barrow_synthesis_data[(Barrow_synthesis_data['Soil_layer'].str.capitalize()==layer.capitalize())&(Barrow_synthesis_data['Core_ID']==Core_ID)&(Barrow_synthesis_data['Incubation_Temperature']>minT)]
    if len(data_layer)==0:
        raise ValueError('No data found for %s layer in core %s'%(layer,Core_ID))
    return data_layer


# CEC numbers estimated from http://www.soilquality.org.au/factsheets/cation-exchange-capacity
def make_initcond(Core_ID,layer,cellulosefrac=0.25,porosity=None,minT=4.0,FeII_FeVF_factor=1.5,CEC_OM=200,CEC_mineral=25,oxic=False,Fe_SSA=1e2,carboxylate_conc=0.1,BD_factor=1.0,otherconstraints={},BD_method='Bockheim'):
    if isinstance(Core_ID,int):
        Core_ID='NGADG%04d'%Core_ID
    data_layer=Barrow_synthesis_data[(Barrow_synthesis_data['Soil_layer'].str.capitalize()==layer.capitalize())&(Barrow_synthesis_data['Core_ID']==Core_ID)&(Barrow_synthesis_data['Incubation_Temperature']>minT)]
    if len(data_layer)==0:
        raise ValueError('No data found for %s layer in core %s'%(layer,Core_ID))
    
    SOC_layer=data_layer['SOC'].mean()/100 # In fraction, not percent
    
    # Based on Bockheim et al 2003 fit. Units are g/cm3
    # From eyeballing the graph in that paper, spread seems to be around +/- 50%
    # And BD seems to be understimated at OC>30%
    if BD_method == 'Bockheim' or oxic:
        BD_layer=max(numpy.exp((8.2432-SOC_layer*100)/9.7872),0.125)
    elif BD_method =='porosity':
        # Alternate method: From corephysdata, relationship between porosity and BD is almost linear. So if we can estimate porosity from water content, we can infer BD
        if layer=='Mineral':
            BD_layer = BD_porosity_fit_mineral.intercept/100/(data_layer['Moisture'].max() - BD_porosity_fit_mineral.slope/100)
        elif layer=='Organic':
            BD_layer = BD_porosity_fit_organic.intercept/100/(data_layer['Moisture'].max() - BD_porosity_fit_organic.slope/100)
        else:
            BD_layer = BD_porosity_fit.intercept/100/(data_layer['Moisture'].max() - BD_porosity_fit.slope/100)
    else:
        raise ValueError('BD_method must be "Bockheim" or "porosity"')

    BD_layer = BD_layer*BD_factor
    
    if porosity is None:
        if oxic:
            raise ValueError('Cannot calculate porosity for oxic incubation')
        porosity=data_layer['Moisture'].max()*BD_layer

    
    # Set initial iron concentration based on maximum dissolved Fe(II), multiplied by a correction factor to account for Fe(II) sorption and/or incomplete Fe(III) reduction by end of incubation
    Fe_VF=(data_layer['Fe_II'].max()-data_layer['Fe_II'].min())*FeII_FeVF_factor*1e-6*BD_layer*molar_volume_FeOH3
    FeII_t0=data_layer['Fe_II'].dropna().iloc[0]*1e-6*BD_layer/porosity*1e3
    
    # CEC affects sorption of cations as well as buffering capacity. Calculate as weighted average of CEC for organic and mineral components
    # Units of input CEC are meq/100g which is standard in soil science measurements. Convert to eq/m3 for PFLOTRAN
    CEC=(CEC_OM*SOC_layer + CEC_mineral*(1-SOC_layer))*1e-2*(BD_layer*1e3)
    
    pH_t0=data_layer['pH'].dropna().iloc[0]
    
    newconstraints=otherconstraints.copy()
    newconstraints['cellulose']=SOC_layer/12*BD_layer*100**3*cellulosefrac
    newconstraints['Fe(OH)3']='%1.4e  %1.2e m^2/m^3'%(Fe_VF,Fe_SSA)
    newconstraints['Fe++']=FeII_t0
    newconstraints['H+']='%1.2f P'%pH_t0
    
    # Cations. Using mean for now, could change to be depth-aware
    newconstraints['Mg++']=waterchemistry['Magnesium'].mean()*1e-3
    newconstraints['Ca++']=waterchemistry['Calcium'].mean()*1e-3
    newconstraints['Na+']=waterchemistry['Sodium'].mean()*1e-3
    newconstraints['K+']=waterchemistry['Potassium'].mean()*1e-3

    newconstraints['DOM1']=data_layer['WEOC'].dropna().iloc[0]*1e-6/BD_layer*1000/porosity
    newconstraints['Acetate-']=data_layer['TOAC'].dropna().iloc[0]*1e-6/BD_layer*1000/porosity
                      
    if oxic:
        newconstraints['O2(aq)']='0.2 G O2(g)'
    
    new_initcond=decomp_network.change_constraints(pools,newconstraints)
    # if oxic:
    #     new_initcond=decomp_network.change_site_density(new_initcond, '>Carboxylate-', SOC_layer*100*0.5e1*BD_layer)
    # else:
    new_initcond=decomp_network.change_site_density(new_initcond, '>Carboxylate-', SOC_layer*carboxylate_conc*BD_layer*100**3)
    
    return new_initcond,BD_layer,SOC_layer,CEC,porosity
    

cores_trough=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Anoxic')&(Barrow_synthesis_data['Incubation_Temperature']>4)&(Barrow_synthesis_data['Microtopography'].str.lower()=='trough')]['Core_ID'].unique()
cores_rim=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Anoxic')&(Barrow_synthesis_data['Incubation_Temperature']>4)&(Barrow_synthesis_data['Microtopography'].str.lower()=='rim')]['Core_ID'].unique()
cores_center=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Anoxic')&(Barrow_synthesis_data['Incubation_Temperature']>4)&(Barrow_synthesis_data['Microtopography'].str.lower()=='center')]['Core_ID'].unique()
cores_oxic=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Oxic')&(Barrow_synthesis_data['Incubation_Temperature']>4)]['Core_ID'].unique()

good_cores_trough_anoxic=9
good_cores_rim_anoxic=5
good_cores_center_anoxic=[17,3]
good_cores_oxic=3 # Note, Core 3 is oxic in organic horizon and anoxic in mineral. No crossover.


# To do: Move these into loop so we can test BD uncertainty range
BD_factor = 1.0
pools_atmoO2_organic,BD_atmoO2_organic,SOC_atmoO2_organic,CEC_atmoO2_organic,porosity_atmoO2_organic=make_initcond(3,'Organic',cellulosefrac=.05,oxic=True,porosity=porosity,BD_factor=BD_factor)
pools_organic_trough,BD_organic_trough,SOC_organic_trough,CEC_organic_trough,porosity_organic_trough=make_initcond(9,'Organic',BD_factor=BD_factor)   
pools_organic_nottrough,BD_organic_nottrough,SOC_organic_nottrough,CEC_organic_nottrough,porosity_organic_nottrough=make_initcond(5,'Organic',BD_factor=BD_factor)   
# result_lowFe_organic,output_units=run_alquimia.run_simulation('Arctic_redox_generated.in',simlength,timestep,initcond=pools_lowFe_organic,hands_off=False,rateconstants=rateconstants,bc=pools_atmoO2_organic,diffquo={'O2(aq)':O2_initial},truncate_concentration=truncate_conc)

# No oxic mineral horizon
# pools_atmoO2_mineral,BD_atmoO2_mineral,SOC_atmoO2_mineral,CEC_atmoO2_mineral=make_initcond(3,'Mineral',cellulosefrac=.05*.04,oxic=True)
# result_highO2_mineral,output_units=run_alquimia.run_simulation('Arctic_redox_generated.in',simlength,timestep,initcond=pools_atmoO2_mineral,bc=pools_atmoO2_mineral,diffquo={'O2(aq)':O2_const},hands_off=False,rateconstants=rateconstants,truncate_concentration=truncate_conc,CEC=CEC_atmoO2_mineral)
pools_mineral_trough,BD_mineral_trough,SOC_mineral_trough,CEC_mineral_trough,porosity_mineral_trough=make_initcond(9,'Mineral',BD_factor=BD_factor)
pools_mineral_nottrough,BD_mineral_nottrough,SOC_mineral_nottrough,CEC_mineral_nottrough,porosity_mineral_nottrough=make_initcond(5,'Mineral',BD_factor=BD_factor)   

simlength=150
initfrac=0.0
timestep=3600
dq=0.01 # Diffusion coefficient when aerobic
oxicfrac=0.1

pH_obs=float(decomp_network.pools_list_to_dict(pools_organic_trough)['H+']['constraints']['initial'].split()[0]) 

def convert_to_xarray(df,units):
    ds=df.rename_axis(index='time').to_xarray()
    # ds['time'].units='days'
    for v in units:
        ds[v].attrs['units']=units[v]
    return ds

if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-f',dest='fname',help='Output file name',default='')
    parser.add_argument('-n',dest='jobnum',help='Job number',default=1)
    parser.add_argument('-N',dest='totaljobs',help='Total number of jobs',default=1)
    options = parser.parse_args()

    jobnum=int(options.jobnum)-1
    totaljobs=int(options.totaljobs)

    import datetime
    if options.fname != '':
        fname=options.fname
    else:
        today=datetime.datetime.today()
        filename='Arctic_Fe_output/output_{year:04d}-{month:02d}-{day:02d}.nc'.format(year=today.year,month=today.month,day=today.day)

    if jobnum+1>totaljobs:
        raise ValueError('jobnum + 1 > totaljobs')
    if totaljobs>1:
        filename=filename[:-3]+'_%02d.nc'%jobnum
        deckname='Arctic_Fe_output/Arctic_redox_generated_%02d.in'%jobnum
    else:
        deckname='Arctic_Fe_output/Arctic_redox_generated.in'
    # import os
    # if os.path.exists(filename):
    #     x=2
    #     while True:
    #         if not os.path.exists(filename[:-3]):
    #             filename=filename[:-3]+'_%d.nc'%x
    #             break
    #         x+=1

        
    O2_const=numpy.zeros(365*24)+dq
    O2_initial=numpy.zeros(simlength*24)
    O2_initial[:int(simlength*24*initfrac)]=dq

    incubations=[
        'highO2_organic',
        'organic_trough',
        'organic_nottrough',
        'mineral_trough',
        'mineral_nottrough'
    ]

    simtypes=[]
    nperiods=[]
    pHsims=[]
    BDsims=[]
    inhibs=[]
    Fescales=[]

    for inhib in [True,False,'noFe']:
        # For comparison with data
        for sim in incubations:
            for BD_method in ['Bockheim']:
                for Fescale in [0.1,0.5,1.0,2.0]:
                    simtypes.append(sim)
                    nperiods.append(None)
                    pHsims.append(None)
                    BDsims.append(BD_method)
                    inhibs.append(inhib)
                    Fescales.append(Fescale)

        # Organic horizon oxic-anoxic
        for ndryperiods in [1,2,3,4,5]:
            for pH in [4.5,5.0,5.5]:
                for Fescale in [0.1,0.5,1.0,2.0]:
                    simtypes.append('organic_trough')
                    nperiods.append(ndryperiods)
                    pHsims.append(pH)
                    BDsims.append('Bockheim')
                    inhibs.append(inhib)
                    Fescales.append(Fescale)

        # Mineral horizon oxic-anoxic
        for ndryperiods in [1,2,3,4,5]:
            for pH in [4.5,5.0,5.5]:
                for Fescale in [0.1,0.5,1.0,2.0]:
                    simtypes.append('mineral_trough')
                    nperiods.append(ndryperiods)
                    pHsims.append(pH)
                    BDsims.append('Bockheim')
                    inhibs.append(inhib)
                    Fescales.append(Fescale)

    sims_thisjob = list(range(jobnum,len(simtypes),totaljobs))
    print('Total number of sims: %d'%len(simtypes))
    print('Number of sims in this job: %d'%len(sims_thisjob))
    print('This job: ',sims_thisjob)


    for simnum in sims_thisjob:
        BD_method=BDsims[simnum]
        if simtypes[simnum] == 'highO2_organic':
            initcond,BD,SOC,CEC,porosity=make_initcond(3,'Organic',cellulosefrac=.25,oxic=True,porosity=0.6,BD_method=BD_method)
        elif simtypes[simnum] == 'organic_trough':
            initcond,BD,SOC,CEC,porosity=make_initcond(9,'Organic',BD_method=BD_method)   
        elif simtypes[simnum] == 'organic_nottrough':
            initcond,BD,SOC,CEC,porosity=make_initcond(5,'Organic',BD_method=BD_method)   
        elif simtypes[simnum] == 'mineral_trough':
            initcond,BD,SOC,CEC,porosity=make_initcond(9,'Mineral',BD_method=BD_method)
        elif simtypes[simnum] == 'mineral_nottrough':
            initcond,BD,SOC,CEC,porosity=make_initcond(5,'Mineral',BD_method=BD_method)   
        else:
            raise ValueError(f'Unknown simulation type: {simtypes[simnum]}')

        run_name=simtypes[simnum]
        if nperiods[simnum] is not None:
            run_name = run_name + '_nperiods_%d'%nperiods[simnum]

            spinup_len=0 #15
            spinup=numpy.zeros(24*spinup_len)
            # spinup[:5*24]=dq
            O2_periodic=numpy.zeros(int((simlength-spinup_len)/nperiods[simnum]*24))
            # O2_periodic[:int((simlength-spinup_len)/nperiods[simnum]*24*oxicfrac)]=dq
            O2_periodic[-int((simlength-spinup_len)/nperiods[simnum]*24*oxicfrac):]=dq
            O2_periodic=numpy.resize(O2_periodic,(simlength-spinup_len)*24)
            diffquo={'O2(aq)':numpy.concatenate((spinup,O2_periodic))}
        else:
            run_name = run_name + '_BD_%s'%BD_method
            if simtypes[simnum].startswith('highO2'):
                diffquo={'O2(aq)':O2_const}
            else:
                diffquo={'O2(aq)':O2_initial}
        if pHsims[simnum] is not None:
            initcond=decomp_network.change_constraint(initcond,'H+','%1.2f P'%pHsims[simnum])
            run_name = run_name + '_pH_%1.1f'%pHsims[simnum]
        if inhibs[simnum] == False:
            run_name = run_name + '_noFeCH4Inhibition'
        if inhibs[simnum] == 'noFe':
            run_name = run_name + '_noFe'
            initcond=decomp_network.change_constraint(initcond,'Fe(OH)3','0.0d-3  1.d0 m^2/m^3')
            initcond=decomp_network.change_constraint(initcond,'Fe++',1e-10)
        if Fescales[simnum]!=1.0:
            run_name = run_name + '_Fescale_%1.1f'%Fescales[simnum]
        
        # Generate PFLOTRAN input file with correct reactions
        rate_scale=0.75e-8
        reactions=make_reactions(Fe_inhibit_CH4=inhibs[simnum],rate_scale=rate_scale)
        rateconstants_named={
            'Fe(II) abiotic oxidation'              : 0.0,
            'Fe(II) microbial oxidation'            : rate_scale*200,
            'Hydrogen oxidation'                    : rate_scale*200,
            'fermentation'                          : rate_scale*4*Fescales[simnum],
            "DOM aerobic respiration"               : rate_scale*25,
            "Acetate aerobic respiration"           : rate_scale*40,
            "Fe(III) reduction"                     : rate_scale*3,
            "Acetaclastic methanogenesis"           : rate_scale*2, 
            "Hydrogenotrophic methanogenesis"       : rate_scale*3*0.32, # Scaling factor between acetaclastic and hydrogenotrophic from optimization of Kotsyurbenko et al data 
            # "Aerobic decomposition"                 : rate_scale*10.0*0, #1.0/(365*24*3600)*1.0
            "Hydrolysis"                            : rate_scale*5.0 #1.0/(365*24*3600)*1.0
        }

        if inhibs[simnum] == 'noFe':
            rateconstants_named['Fe(III) reduction']=0.0

        rateconstants=run_alquimia.convert_rateconstants(rateconstants_named,reactions=reactions)

        reaction_network =  decomp_network.decomp_network(pools=pools,reactions=reactions)

        decomp_network.PF_network_writer(reaction_network).write_into_input_deck('SOMdecomp_template.txt',deckname,log_formulation=True,CO2name='Tracer',truncate_concentration=1e-25,database='/home/b0u/models/PFLOTRAN/REDOX-PFLOTRAN/hanford.dat')
            
        
        result,units=run_alquimia.run_simulation(deckname,simlength,timestep,run_name=run_name,initcond=initcond,
                        hands_off=False,rateconstants=rateconstants,bc=pools_atmoO2_organic,
                        diffquo=diffquo,truncate_concentration=truncate_conc,CEC=CEC,porosity=porosity) 

        if simnum==sims_thisjob[0]:
            mode='w'
        else:
            mode='a'
        convert_to_xarray(result,units).to_netcdf(filename,group=run_name,mode=mode)

    # result_highO2_organic,output_units=run_alquimia.run_simulation(deckname,simlength,timestep,run_name='highO2_organic',initcond=pools_atmoO2_organic,bc=pools_atmoO2_organic,diffquo={'O2(aq)':O2_const},hands_off=False,rateconstants=rateconstants,truncate_concentration=truncate_conc,CEC=CEC_atmoO2_organic,porosity=porosity_atmoO2_organic)
    # result_organic_trough,output_units=run_alquimia.run_simulation(deckname,simlength,timestep,run_name='organic_trough',initcond=pools_organic_trough,hands_off=False,rateconstants=rateconstants,bc=pools_atmoO2_organic,diffquo={'O2(aq)':O2_initial},truncate_concentration=truncate_conc,CEC=CEC_organic_trough,porosity=porosity_organic_trough) 
    # result_organic_nottrough,output_units=run_alquimia.run_simulation(deckname,simlength,timestep,run_name='organic_nottrough',initcond=pools_organic_nottrough,hands_off=False,rateconstants=rateconstants,bc=pools_atmoO2_organic,diffquo={'O2(aq)':O2_initial},truncate_concentration=truncate_conc,CEC=CEC_organic_nottrough,porosity=porosity_organic_nottrough) 
    # result_mineral_trough,output_units=run_alquimia.run_simulation(deckname,simlength,timestep,run_name='mineral_trough',initcond=pools_mineral_trough,hands_off=False,rateconstants=rateconstants,bc=pools_atmoO2_organic,diffquo={'O2(aq)':O2_initial},truncate_concentration=truncate_conc,CEC=CEC_mineral_trough,porosity=porosity_mineral_trough)    
    # result_mineral_nottrough,output_units=run_alquimia.run_simulation(deckname,simlength,timestep,run_name='mineral_nottrough',initcond=pools_mineral_nottrough,hands_off=False,rateconstants=rateconstants,bc=pools_atmoO2_organic,diffquo={'O2(aq)':O2_initial},truncate_concentration=truncate_conc,CEC=CEC_mineral_nottrough,porosity=porosity_mineral_nottrough) 

    

    # periodic_out={}
    
    # for ndryperiods in range(1,5):
    #     periodic_out[ndryperiods]={}
    #     for pH in pHs:

    #         O2_periodic=numpy.zeros(int(simlength/ndryperiods)*24)
    #         O2_periodic[:int(simlength/ndryperiods*24*oxicfrac)]=dq
    #         result_periodicO2,output_units=run_alquimia.run_simulation('Arctic_redox_generated.in',simlength,timestep,run_name='Periodic cycles=%d, pH=%1.1f'%(ndryperiods,pH),
    #             initcond=decomp_network.change_constraint(pools_organic_trough,'H+','%1.2f P'%pH),bc=pools_atmoO2_organic,
    #             diffquo={'O2(aq)':O2_periodic},hands_off=False,rateconstants=rateconstants,truncate_concentration=truncate_conc,CEC=CEC_organic_trough)
    #         periodic_out[ndryperiods][pH]=result_periodicO2


    # convert_to_xarray(result_highO2_organic,output_units).to_netcdf(filename,group='highO2_organic')
    # convert_to_xarray(result_organic_trough,output_units).to_netcdf(filename,group='organic_trough',mode='a')
    # convert_to_xarray(result_organic_nottrough,output_units).to_netcdf(filename,group='organic_nottrough',mode='a')
    # convert_to_xarray(result_mineral_trough,output_units).to_netcdf(filename,group='mineral_trough',mode='a')
    # convert_to_xarray(result_mineral_nottrough,output_units).to_netcdf(filename,group='mineral_nottrough',mode='a')

    # import xarray
    # periodic_x=xarray.Dataset()
    # for ndryperiods in periodic_out.keys():
    #     pH_sims=xarray.concat([convert_to_xarray(periodic_out[ndryperiods][pH],output_units) for pH in pHs],dim='pH',join='outer')
    #     periodic_x=xarray.concat([periodic_x,pH_sims],dim='ndryperiods')

    # periodic_x['pH']=pHs
    # periodic_x['ndryperiods']=periodic_x['ndryperiods']+1

    # periodic_x.to_netcdf(filename,group='periodic_sims',mode='a')
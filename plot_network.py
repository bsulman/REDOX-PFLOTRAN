import network_for_ELM,decomp_network
import matplotlib.pyplot as plt

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


# All of this stuff is for the plot of the network diagram
mineral_col=180.0
substrate_col=130.0
reduction_col=140.0
TEA_col=155.0
oxidation_col=165.0
substrateox_col=reduction_col

gas_col=TEA_col
redox_seq={
    'O2':1300,
    'NO3-':800,
    # 'Mn4+':700,
    'Fe3+':400,
    'SO4--':650,
    'CO2':200,
    'CH4':100,
    'SOM':1600,
    'DOM':1200,
    'Acetate':900,
    'H2':600
    }
dy_SOM=120
pos={'Pyrrhotite': (mineral_col, (redox_seq['Fe3+']+redox_seq['SO4--'])/2),
 'Pyrite': (mineral_col, (redox_seq['Fe3+']+redox_seq['SO4--']+redox_seq['SO4--']-100)/3-100),
 'SO4--': (TEA_col, redox_seq['SO4--']),
 'HS-': (TEA_col, redox_seq['SO4--']-140),
 'H2S(aq)': (gas_col, redox_seq['SO4--']),
 'Sulfate reduction': (reduction_col, (redox_seq['SO4--'])),
 'Methane oxidation (SO4)': (oxidation_col, redox_seq['SO4--']-200),
 'Sulfide oxidation': (oxidation_col, redox_seq['SO4--']),
 'Hydrogen oxidation': (oxidation_col, redox_seq['H2']-400),
 
 'Fe(OH)3': (mineral_col, redox_seq['Fe3+']),
 'Goethite': (mineral_col+10, redox_seq['Fe3+']),
 'Fe+++': (TEA_col, redox_seq['Fe3+']),
 'Fe++': (TEA_col, redox_seq['Fe3+']-75),
 'Fe(II) microbial oxidation': (oxidation_col, (redox_seq['Fe3+']+redox_seq['O2'])/2),
 'Fe(III) reduction': (reduction_col, redox_seq['Fe3+']),
 'Methane oxidation (Fe)': (oxidation_col, redox_seq['Fe3+']),
 
 'NH4+': (TEA_col, redox_seq['NO3-']+100),
 'NO3-': (TEA_col, redox_seq['NO3-']),
 'Denitrification': (reduction_col, redox_seq['NO3-']-40),
 'Plant NH4 uptake': (oxidation_col-5, redox_seq['NO3-']+100),
 'Plant NO3 uptake': (oxidation_col-10, redox_seq['NO3-']),
 'Nitrification': (substrateox_col, redox_seq['NO3-']+75/2),
 'N2(aq)': (gas_col, redox_seq['NO3-']-150),
 'N2O(aq)': (gas_col, redox_seq['NO3-']-75),
 'Methane oxidation (NO3)': (oxidation_col, (redox_seq['Fe3+'])),

 'CWD': (substrate_col-4, redox_seq['SOM']+dy_SOM*7),
 'LITR1': (substrate_col-4, redox_seq['SOM']+dy_SOM*6),
 'LITR2': (substrate_col-4, redox_seq['SOM']+dy_SOM*5),
 'LITR3': (substrate_col-4, redox_seq['SOM']+dy_SOM*4),
 'SOIL1': (substrate_col-2, redox_seq['SOM']+dy_SOM*3),
 'SOIL2': (substrate_col-2, redox_seq['SOM']+dy_SOM*2),
 'SOIL3': (substrate_col-2, redox_seq['SOM']+dy_SOM*1),
 'SOIL4': (substrate_col-2, redox_seq['SOM']+0),
 'CWD fragmentation': (reduction_col-4, redox_seq['SOM']+dy_SOM*7),
 'LITR1 decomp': (reduction_col-3, redox_seq['SOM']+dy_SOM*6),
 'LITR2 decomp': (reduction_col-2, redox_seq['SOM']+dy_SOM*5),
 'LITR3 decomp': (reduction_col-1, redox_seq['SOM']+dy_SOM*4),
 'SOIL1 decomp': (reduction_col, (redox_seq['SOM']+dy_SOM*3)),
 'SOIL2 decomp': (reduction_col+1, (redox_seq['SOM']+dy_SOM*2)),
 'SOIL3 decomp': (reduction_col+2, (redox_seq['SOM']+dy_SOM)),
 'SOIL4 decomp': (reduction_col+3, (redox_seq['SOM']+0)),

 'DOM1': (substrate_col, redox_seq['DOM']+100),
 'Fermentation': (reduction_col, (redox_seq['DOM']+redox_seq['Acetate'])/2-30),
 'Acetate-': (substrate_col, redox_seq['Acetate']),
 'H2(aq)': (substrate_col, redox_seq['H2']),
 'H+': (170,1100),
 
 'O2(aq)': (TEA_col, redox_seq['O2']),
 'CO2(aq)': (175, redox_seq['CO2']-150),

 'CH4(aq)': (TEA_col, redox_seq['CH4']),

 
 'DOM1 respiration': (substrateox_col, (redox_seq['O2']+redox_seq['Acetate'])/2+155),
 'Aerobic decomposition': (substrateox_col, (redox_seq['O2']+redox_seq['Acetate'])/2+155),
 'Acetate aerobic respiration': (substrateox_col, (redox_seq['O2']+redox_seq['Acetate'])/2-200),
 'Methane oxidation': (oxidation_col, (1000)),

 'Hydrogenotrophic methanogenesis': (reduction_col, (redox_seq['CO2'])-100),
 'Acetaclastic methanogenesis': (reduction_col, (redox_seq['CO2']+100-50))}


opts=dict(
omit=['secondary','H2O','Na+','Cl-','Halite','Calcite','Ca++','Tracer','Tracer2',
        'HRimm','Nimm','Nimp','Nmin','Plant_NH4_demand','Plant_NO3_demand',
        'Fe(OH)2','CO2(g)*','O2(g)','CH4(g)','gas',#'LITR2','LITR3','CWD',
        'LITR2 decomp anox','LITR3 decomp anox','LITR1 decomp anox',
        'SOIL1 decomp anox','SOIL2 decomp anox','SOIL3 decomp anox','SOIL4 decomp anox',
        'Nitrification','NO3-','NH4+','Plant NO3 uptake','Plant NH4 uptake'],

pos=pos,

font_size='medium',font_weight='normal',node_size=800,font_color='w',arrowstyle='-|>',arrowsize=15.0,width=1.0,edge_color='k',node_alpha=0.9,node_colors=node_colors,markers={'Reaction':'None','mineral':'8'},
namechanges={'MAOM':'SOM','DOM1':'DOM','LITR1':'Litter','O2(aq)':'O$_2$','CH4(aq)':'CH$_4$','HCO3-':'CO$_2$','DOM2':'Exposed lignin','sorbed_DOM1':'Sorbed DOM',
        'Fe(OH)2':'Fe(OH)$_2$','Fe(OH)3':'Fe(OH)$_3$','Mn++':r'Mn$^\mathrm{+\!\!+}$','Mn+++':r'Mn$^\mathrm{+\!\!+\!\!\!+}$','Acetate-':'Acetate',
        'H2(aq)':'H$_2$','DOM oxidation (O2)':'Aerobic\nresp','Acetate oxidation (O2)':'Aerobic\nresp','Aerobic decomposition':'Aerobic\nresp','Acetate aerobic respiration':'Acetate\naerobic resp.',
        'NH4+':'NH$_4^+$','NO3-':'NO$_3^-$','N2O(aq)':'N$_2$O','N2(aq)':'N$_2$','Fe++':r'Fe$^\mathrm{+\!\!+}$','Fe+++':r'Fe$^\mathrm{+\!\!+\!\!\!+}$',
        'SO4--':'SO$_4^{--}$','HS-':'H$_2$S','Methane oxidation (O2)':'Methane\noxidation','Methane oxidation (NO3)':'Methane\noxidation',
        'Methane oxidation (Fe)':'Methane\noxidation','Methane oxidation (SO4)':'Methane\noxidation','Fe(III) reduction':'Fe(III)\nreduction','Sulfate reduction':'Sulfate\nreduction',
        'Hydrogenotrophic methanogenesis':'Hydrogenotrophic\nmethanogenesis','Acetaclastic methanogenesis':'Acetoclastic\nmethanogenesis',
        'SOMdecomp Reaction':'SOM Reaction','General Reaction':'Abiotic Reaction','fermentation':'Fermentation','Primary aqueous':'Dissolved ion','Gas':'Dissolved gas',
        'CO2(aq)':'CO$_2$','CWD fragmentation':'CWD frag','Plant NH4 uptake':'Plant NH$_4$\nuptake','Plant NO3 uptake':'Plant NO$_3$ uptake',
        'DOM1 respiration':'DOM\nrespiration','Fe(II) microbial oxidation':'Fe(II)\noxidation','Pyrrhotite':"FeS"},connectionstyle='arc3, rad=0.1')

f,a=plt.subplots(num='Reaction network',clear=True,figsize=(9.2,6.4))
to_draw,pos=decomp_network.draw_network_with_reactions(network_for_ELM.make_aqueous_network(Fe=True,sulfate=True,methane=True),do_legend=True,ax=a,**opts)

opts.pop('pos',None)

for p in a.patches:
    if p._posA_posB[0][0]>p._posA_posB[1][0]:
        p.set_linestyle(':')
        p.set_connectionstyle('arc3,rad=-0.1')
a.set(facecolor='#8a9ebf',title='ELM-PFLOTRAN decomposition and redox reaction network')

f,a=plt.subplots(num='Reaction network 2',clear=True,figsize=(9.2,8.4))
to_draw,pos=decomp_network.draw_network_with_reactions(network_for_ELM.make_aqueous_network(Fe=True,sulfate=True,methane=True),do_legend=True,ax=a,**opts)
a.set(facecolor='#8a9ebf',title='ELM-PFLOTRAN decomposition and redox reaction network')

for p in a.patches:

    if p._posA_posB[0][0]>p._posA_posB[1][0]:
        p.set_linestyle(':')
        p.set_connectionstyle('arc3,rad=-0.1')




plt.show()
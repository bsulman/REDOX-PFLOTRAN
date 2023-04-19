import decomp_network
import plot_pf_output

pflotran_exe='../pflotran-interface/src/pflotran/pflotran'

# Simple decomposition network 1: 1 POM pool, one DOM pool, first order
POM=decomp_network.decomp_pool(name='POM',CN=15,initval=1e3,kind='immobile')
DOM=decomp_network.decomp_pool(name='DOM1',CN=15,initval=1e-3,kind='primary')
CO2=decomp_network.decomp_pool(name='CO2(aq)',kind='primary',initval=1e-15)
POM_dissolve=decomp_network.reaction(name='POM dissolution (no inhibition)',reactant_pools={'POM':1.0},product_pools={'DOM1':1.0},rate_constant=0.1,reactiontype='SOMDECOMP')
simple_network_DOM1 = decomp_network.decomp_network(pools=[POM,DOM,CO2],reactions=[POM_dissolve])


DOM1_result,DOM1_units=decomp_network.PF_network_writer(simple_network_DOM1).run_simulation('SOMdecomp_template.txt','DOM1',pflotran_exe,CO2name='CO2(aq)')

# Simple decomposition network 2: 1 POM pool, one DOM pool, inhibited by DOM
POM=decomp_network.decomp_pool(name='POM',CN=15,initval=1e3,kind='immobile')
DOM=decomp_network.decomp_pool(name='DOM1',CN=15,initval=1e-3,kind='primary')
POM_dissolve_inhib=decomp_network.reaction(name='POM dissolution (inhibition)',reactant_pools={'POM':1.0},product_pools={'DOM1':1.0},rate_constant=0.1,
                    inhibition_terms=[{'species':'DOM1','type':'MONOD','k':1e-2}],reactiontype='SOMDECOMP')
simple_network_DOM2 = decomp_network.decomp_network(pools=[POM,DOM,CO2],reactions=[POM_dissolve_inhib])


DOM2_result,DOM2_units=decomp_network.PF_network_writer(simple_network_DOM2).run_simulation('SOMdecomp_template.txt','DOM2',pflotran_exe,CO2name='CO2(aq)')

# Microbial decomposition network with priming effects
microbes=decomp_network.decomp_pool(name='MICROBES',CN=10,initval=1e-3,kind='immobile')
SOM1=decomp_network.decomp_pool(name='SOM1',CN=10,initval=1.1e3,kind='immobile')


priming_network = decomp_network.decomp_network(pools=[POM,microbes,SOM1,CO2])
priming_network.add_reaction(decomp_network.reaction(name='POM decomposition',reactant_pools={'POM':1.0},product_pools={'MICROBES':0.6,CO2['name']:0.4},rate_constant=0.2,monod_terms=[{'species':'MICROBES','k':1e0}],reactiontype='SOMDECOMP'))
priming_network.add_reaction(decomp_network.reaction(name='SOM1 decomposition',reactant_pools={'SOM1':1.0},product_pools={'MICROBES':0.1,CO2['name']:0.9},rate_constant=10,monod_terms=[{'species':'MICROBES','k':1e0}],reactiontype='SOMDECOMP'))
priming_network.add_reaction(decomp_network.reaction(name='Microbial turnover',reactant_pools={'MICROBES':1.0},product_pools={CO2['name']:0.5,'SOM1':0.5},rate_constant=0.05,inhibition_terms=[{'species':'MICROBES','k':1e-5,'type':'INVERSE_MONOD'}],reactiontype='SOMDECOMP'))
priming_result,priming_units=decomp_network.PF_network_writer(priming_network).run_simulation('SOMdecomp_template.txt','priming',pflotran_exe,CO2name='CO2(aq)',length_days=100)

# A more complex decomposition network approximating the structure of the CORPSE model

decomp_network_CORPSE=decomp_network.decomp_network()
decomp_network_CORPSE.add_pool(decomp_network.decomp_pool(name='LABILE_POM',CN=    15.0,initval=0.3,kind='immobile'))
decomp_network_CORPSE.add_pool(decomp_network.decomp_pool(name='RESISTANT_POM',CN= 15.0,initval=0.7,kind='immobile'))
decomp_network_CORPSE.add_pool(decomp_network.decomp_pool(name='DOM1',CN=    15.0,initval=1e-5,kind='primary'))
decomp_network_CORPSE.add_pool(decomp_network.decomp_pool(name='DOM2',CN= 15.0,initval=1e-5,kind='primary'))
decomp_network_CORPSE.add_pool(decomp_network.decomp_pool(name='SOIL_MICROBES',CN= 10.0,initval=1e-4,kind='immobile'))
decomp_network_CORPSE.add_pool(decomp_network.decomp_pool(name='NECROMASS',CN=     10.0,initval=1e-5,kind='immobile'))
decomp_network_CORPSE.add_pool(decomp_network.decomp_pool(name='DOM3',CN= 10.0,initval=1e-5,kind='primary'))
decomp_network_CORPSE.add_pool(decomp_network.decomp_pool(name='MAOM',CN=          10.0,initval=0.5,kind='immobile'))
decomp_network_CORPSE.add_pool(decomp_network.decomp_pool(name='HCO3-',kind='primary'))

# Microbial decomposition of unprotected SOM
decomp_network_CORPSE.add_reaction(decomp_network.reaction(reactant_pools={'LABILE_POM':1.0},product_pools={'SOIL_MICROBES':0.6,'HCO3-':0.4},
                rate_constant=0.1,rate_units='y',name='Labile POM microbial decomposition',reactiontype='SOMDECOMP',
                                monod_terms=[{'species':'SOIL_MICROBES','k':0.1e-3}]))
decomp_network_CORPSE.add_reaction(decomp_network.reaction(reactant_pools={'RESISTANT_POM':1.0},product_pools={'SOIL_MICROBES':0.1,'HCO3-':1-0.1},
                    rate_constant=0.5,rate_units='y',name='Resistant POM microbial decomposition',reactiontype='SOMDECOMP',
                                monod_terms=[{'species':'SOIL_MICROBES','k':0.1e-3}]))
decomp_network_CORPSE.add_reaction(decomp_network.reaction(reactant_pools={'NECROMASS':1.0},product_pools={'SOIL_MICROBES':0.6,'HCO3-':1-0.6},
                    rate_constant=0.1,rate_units='y',name='Necromass microbial decomposition',reactiontype='SOMDECOMP',
                                monod_terms=[{'species':'SOIL_MICROBES','k':0.1e-3}]))
                                
# Microbial biomass turnover. Here, CUE would determine maintenance respiration loss as a fraction of turnover
decomp_network_CORPSE.add_reaction(decomp_network.reaction(reactant_pools={'SOIL_MICROBES':1.0},product_pools={'NECROMASS':0.6,'HCO3-':1-0.6},
                    rate_constant=0.1,rate_units='y',name='Microbial biomass turnover',reactiontype='SOMDECOMP',
                                inhibition_terms=[{'species':'SOIL_MICROBES','k':0.1e-11,'type':'INVERSE_MONOD'}]))
                                
# Physical/mineral protection of necromass. CUE=1.0 because it is an abiotic conservative reaction
decomp_network_CORPSE.add_reaction(decomp_network.reaction(reactant_pools={'NECROMASS':1.0},product_pools={'MAOM':1.0},
                                rate_constant=0.01,rate_units='y',name='Physical/mineral protection of necromass',reactiontype='SOMDECOMP',
                                inhibition_terms=[{'species':'MAOM','k':1e-7,'type':'MONOD'}]))
# First order turnover of MAOM back to necromass
decomp_network_CORPSE.add_reaction(decomp_network.reaction(reactant_pools={'MAOM':1.0},product_pools={'NECROMASS':1.0},reactiontype='SOMDECOMP',
                                rate_constant=50.0,rate_units='y',name='First order turnover of MAOM back to necromass')) 

# Dissolution into DOM
decomp_network_CORPSE.add_reaction(decomp_network.reaction(reactant_pools={'LABILE_POM':1.0},product_pools={'DOM1':1.0},reactiontype='SOMDECOMP',
                                rate_constant=0.1,rate_units='y',name='Dissolution of labile POM into DOM',
                                inhibition_terms=[{'species':'DOM1','k':1e-12,'type':'MONOD'}]))
                                
decomp_network_CORPSE.add_reaction(decomp_network.reaction(reactant_pools={'RESISTANT_POM':1.0},product_pools={'DOM2':1.0},
                                rate_constant=0.1,rate_units='y',name='Dissolution of resistant POM into DOM',reactiontype='SOMDECOMP',
                                inhibition_terms=[{'species':'DOM1','k':1e-12,'type':'MONOD'}]))
decomp_network_CORPSE.add_reaction(decomp_network.reaction(reactant_pools={'NECROMASS':1.0},product_pools={'DOM3':1.0},
                                rate_constant=0.1,rate_units='y',name='Dissolution of necromass into DOM',reactiontype='SOMDECOMP',
                                inhibition_terms=[{'species':'DOM1','k':1e-12,'type':'MONOD'}]) )     



CORPSE_result,CORPSE_units=decomp_network.PF_network_writer(decomp_network_CORPSE).run_simulation('SOMdecomp_template.txt','CORPSE',pflotran_exe,CO2name='HCO3-')


# CTC decomposition network
decomp_network_CTC=decomp_network.decomp_network()

decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='SOIL1',CN=  12.*14.007/12.011 ,constraints={'initial':1e-10/12.011},kind='immobile') )
decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='SOIL2',CN=  12.*14.007/12.011 ,constraints={'initial':1e-10/12.011},kind='immobile'))
decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='SOIL3',CN=  10.*14.007/12.011 ,constraints={'initial':1e-10/12.011},kind='immobile'))
decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='SOIL4',CN=  10.*14.007/12.011 ,constraints={'initial':1e-10/12.011},kind='immobile'))
decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='LITR1',constraints={'initial':1e3/12.011},initCN=20*14.007/12.011 ,kind='immobile'))
decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='LITR2',constraints={'initial':1e-10/12.011},initCN=20*14.007/12.011 ,kind='immobile')    )
decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='LITR3',constraints={'initial':1e-10/12.011},initCN=20*14.007/12.011 ,kind='immobile'))
decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='CWD',constraints={'initial':1e-10/12.011},initCN=20*14.007/12.011 ,kind='immobile'))
decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='CO2(aq)',kind='primary'))
decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='NH4+',kind='primary'))
decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='HRimm',constraints={'initial':1e-10},kind='immobile'))

# CWD decomposition to  litter
decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'CWD':1.0},product_pools={'LITR2':0.76,'LITR3':0.24},
                            rate_constant=0.00010,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',reactiontype='SOMDECOMP',
                            name='CWD fragmentation'))

# Litter decomposition
decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'LITR1':1.0},product_pools={'SOIL1':0.61,'CO2(aq)':1-0.61},
                rate_constant=1.204,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='LITR1 decomposition',reactiontype='SOMDECOMP'))
decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'LITR2':1.0},product_pools={'SOIL2':0.45,'CO2(aq)':1-0.45},
                rate_constant=0.0726,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='LITR2 decomposition',reactiontype='SOMDECOMP'))
decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'LITR3':1.0},product_pools={'SOIL3':0.71,'CO2(aq)':1-0.71},
                rate_constant=0.0141,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='LITR3 decomposition',reactiontype='SOMDECOMP'))

# SOM decomposition
decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'SOIL1':1.0},product_pools={'SOIL2':0.72,'CO2(aq)':1-0.72},
            rate_constant=0.0726,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='SOIL1 decomp',reactiontype='SOMDECOMP'))
decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'SOIL2':1.0},product_pools={'SOIL3':0.54,'CO2(aq)':1-0.54},
            rate_constant=0.0141,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='SOIL2 decomp',reactiontype='SOMDECOMP'))
decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'SOIL3':1.0},product_pools={'SOIL4':0.45,'CO2(aq)':1-0.45},
            rate_constant=0.00141,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='SOIL3 decomp',reactiontype='SOMDECOMP'))
decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'SOIL4':1.0},product_pools={'CO2(aq)':1.0},
            rate_constant=0.0001,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='SOIL4 decomp',reactiontype='SOMDECOMP'))

CTC_result,CTC_units=decomp_network.PF_network_writer(decomp_network_CTC).run_simulation('SOMdecomp_template.txt','CTC',pflotran_exe,CO2name='CO2(aq)',length_days=1)

from pylab import *
figure('CTC pools',clear=True)
for pool in ['LITR1C','SOIL1','SOIL2','SOIL3','SOIL4']:
    plot(CTC_result[pool]*12/100**2,label=pool.replace('SOIL','Soil').replace('LITR1C','Litter').replace('CWDC','CWD'),lw=2.0)
title('CTC SOM pools')
ylabel('Concentration (g C cm$^{-3})$',fontsize='large')
legend(fontsize='large',ncol=1)
xlabel('Time (days)',fontsize='large')
xticks(fontsize='large')
yticks(fontsize='large')

f,ax=subplots(num='CTC diagram',clear=True)
layout=decomp_network.nx.drawing.nx_agraph.graphviz_layout
pos=layout(decomp_network_CTC,prog='dot')
# pos['CO2(aq)']=(80,18)
decomp_network.draw_network(decomp_network_CTC,omit=['secondary','Rock(s)','NH4+'],arrowstyle='-|>',font_size='x-large',node_size=1800,pos=pos,node_alpha=0.95,font_weight='bold',
    namechanges={'CO2(aq)':'CO$_2$','LITR1':'Litter1','LITR2':'Litter2','LITR3':'Litter3'},do_legend=False,font_color='w',connectionstyle='arc3, rad=0.1')


import pandas
figure('Simulation results',figsize=(10,8));clf()
ax=subplot(241)
DOM1_result,DOM1_units=plot_pf_output.convert_units(DOM1_result,DOM1_units,'M')
DOM2_result,DOM2_units=plot_pf_output.convert_units(DOM2_result,DOM2_units,'M')
plot(DOM1_result['Total DOM1'],c='C0',label='DOM')
plot(DOM1_result['POM'],c='C1',label='POM')
plot(DOM2_result['Total DOM1'],c='C0',label='DOM (with inhibition k=%1.1g)'%simple_network_DOM2.edges[('POM','DOM1',0)]['inhibition_terms'][0]['k'],ls='--')
plot(DOM2_result['POM'],c='C1',label='POM (with inhibition k=%1.1g)'%simple_network_DOM2.edges[('POM','DOM1',0)]['inhibition_terms'][0]['k'],ls='--')
title('DOM simulations')
ylabel('Concentration (M)')
legend()

ax=subplot(242)
plot(priming_result['POM'],c='C1')
plot(priming_result['SOM1'],c='C2')
plot(priming_result['MICROBES'],c='C3')
axins=ax.inset_axes([0.5,0.2,0.4,0.3])
axins.semilogy(priming_result['MICROBES'],c='C3')
title('Priming simulations')
ylabel('Concentration (mol/m3)')
legend()

ax=subplot(243)
for pool in ['SOIL1','SOIL2','SOIL3','SOIL4','LITR1C','LITR2C','LITR3C','CWDC']:
    plot(CTC_result[pool],label=pool)
title('CTC simulations')
ylabel('Concentration (mol C/m3)')
legend()

ax=subplot(244)
for pool in ['LABILE_POM','RESISTANT_POM','SOIL_MICROBES','NECROMASS','MAOM']:
    plot(CORPSE_result[pool],label=pool)
title('CORPSE simulations')
ylabel('Concentration (mol C/m3)')
legend()
    
layout=decomp_network.nx.drawing.nx_agraph.graphviz_layout
ax=subplot(245)
pos=layout(simple_network_DOM1,prog='dot')
decomp_network.draw_network_with_reactions(simple_network_DOM1,ax=ax,arrowsize=20)

ax=subplot(246)
pos=layout(priming_network,prog='dot')
decomp_network.nx.draw_networkx(priming_network,pos=pos,with_labels=True,ax=ax,nodes=priming_network.nodes,node_color=decomp_network.make_nodecolors(priming_network.nodes),arrowsize=20)

ax=subplot(247)
pos=layout(decomp_network_CTC,prog='dot')
decomp_network.nx.draw_networkx(decomp_network_CTC,pos=pos,with_labels=True,ax=ax,nodes=decomp_network_CTC.nodes,node_color=decomp_network.make_nodecolors(decomp_network_CTC.nodes),arrowsize=20)

ax=subplot(248)
pos=layout(decomp_network_CORPSE,prog='dot')
decomp_network.nx.draw_networkx(decomp_network_CORPSE,pos=pos,with_labels=True,ax=ax,nodes=decomp_network_CORPSE.nodes,node_color=decomp_network.make_nodecolors(decomp_network_CORPSE.nodes),arrowsize=20)


tight_layout()
show()



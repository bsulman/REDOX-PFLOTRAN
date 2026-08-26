from Arctic_Fe_CH4 import *
import sys
import datetime
       
O2_const=numpy.zeros(365*24)+dq
O2_initial=numpy.zeros(simlength*24)
O2_initial[:int(simlength*24*initfrac)]=dq



# Set this up so it can be run in parallel and pick which sim or iteration to do from number in command line arg
# Todo: allow param value iterations to be run in parallel as well to further speed things up (looks like a CADES node could run up to 32 threads at once)

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-f',dest='fname',help='Output file name',default='')
parser.add_argument('-n',dest='jobnum',help='Job number',default=0)
parser.add_argument('-N',dest='totaljobs',help='Total number of jobs',default=1)
options = parser.parse_args()

jobnum=int(options.jobnum)
totaljobs=int(options.totaljobs)+1

if options.fname is not '':
    fname=options.fname
else:
    today=datetime.datetime.today()
    filename='Arctic_Fe_output/paramtest_output_{year:04d}-{month:02d}-{day:02d}.nc'.format(year=today.year,month=today.month,day=today.day)

if jobnum+1>totaljobs:
    raise ValueError('jobnum + 1 > totaljobs')
if totaljobs>1:
    filename=filename[:-3]+'_%02d.nc'%jobnum
    deckname='Arctic_redox_generated_%02d.in'%jobnum
else:
    deckname='Arctic_redox_generated.in'

incubations=[
    'highO2_organic',
    'organic_trough',
    'organic_nottrough',
    'mineral_trough',
    'mineral_nottrough'
]

filemode='w'

def scale_carbconc(pools,scale):
    oldval=decomp_network.pools_list_to_dict(pools)['>Carboxylate-']['site_density']
    return decomp_network.change_site_density(pools,'>Carboxylate-',oldval*scale)

simtypes=[]
ascales=[]
frates=[]
carbconcs=[]

# List of all simulations
for simtype in incubations:
    for acetate_scale in 10**numpy.linspace(-2,1,5):
        for fermentation_rate in numpy.linspace(2,100,5):
            for carb_conc in [0.5,1.0,1.5,2.0]:
                simtypes.append(simtype)
                ascales.append(acetate_scale)
                frates.append(fermentation_rate)
                carbconcs.append(carb_conc)

sims_thisjob = list(range(jobnum,len(simtypes),totaljobs))
print('Total number of sims: %d'%len(simtypes))
print('This job: ',sims_thisjob)

for simnum in sims_thisjob:
    simtype=simtypes[simnum]
    acetate_scale=ascales[simnum]
    fermentation_rate=frates[simnum]
    carb_conc=carbconcs[simnum]

    suffix = '_acetatescale_%1.2f_fermentationrate_%1.2f_carbconc=%1.2f'%(acetate_scale,fermentation_rate,carb_conc)
    scales=conc_scales.copy()
    scales['Acetate-']=acetate_scale
    reactions=make_reactions(scales)
    rates=rateconstants_named.copy()
    rates['fermentation']=rate_scale*fermentation_rate

    rateconstants=run_alquimia.convert_rateconstants(rates,reactions=reactions)

    reaction_network =  decomp_network.decomp_network(pools=pools,reactions=reactions)

    decomp_network.PF_network_writer(reaction_network).write_into_input_deck('SOMdecomp_template.txt',deckname,log_formulation=True,CO2name='Tracer',truncate_concentration=1e-25)

    if 'highO2_organic' == simtype:
        result_highO2_organic,output_units=run_alquimia.run_simulation(deckname,simlength,timestep,run_name='highO2_organic'+suffix,
                    initcond=scale_carbconc(pools_atmoO2_organic,carb_conc),bc=pools_atmoO2_organic,diffquo={'O2(aq)':O2_const},
                    hands_off=False,rateconstants=rateconstants,truncate_concentration=truncate_conc,CEC=CEC_atmoO2_organic,porosity=porosity_atmoO2_organic)
        convert_to_xarray(result_highO2_organic,output_units).to_netcdf(filename,group='highO2_organic'+suffix,mode=filemode)
        filemode='a'
    
    if 'organic_trough' == simtype:
        result_organic_trough,output_units=run_alquimia.run_simulation(deckname,simlength,timestep,run_name='organic_trough'+suffix,
                    initcond=scale_carbconc(pools_organic_trough,carb_conc),hands_off=False,rateconstants=rateconstants,
                    bc=pools_atmoO2_organic,diffquo={'O2(aq)':O2_initial},truncate_concentration=truncate_conc,
                    CEC=CEC_organic_trough,porosity=porosity_organic_trough) 
        convert_to_xarray(result_organic_trough,output_units).to_netcdf(filename,group='organic_trough'+suffix,mode=filemode)
        filemode='a'
    
    if 'organic_nottrough' == simtype:
        result_organic_nottrough,output_units=run_alquimia.run_simulation(deckname,simlength,timestep,run_name='organic_nottrough'+suffix,
                    initcond=scale_carbconc(pools_organic_nottrough,carb_conc),hands_off=False,rateconstants=rateconstants,
                    bc=pools_atmoO2_organic,diffquo={'O2(aq)':O2_initial},truncate_concentration=truncate_conc,
                    CEC=CEC_organic_nottrough,porosity=porosity_organic_nottrough) 
        convert_to_xarray(result_organic_nottrough,output_units).to_netcdf(filename,group='organic_nottrough'+suffix,mode=filemode)
        filemode='a'
    
    if 'mineral_trough' == simtype:
        result_mineral_trough,output_units=run_alquimia.run_simulation(deckname,simlength,timestep,run_name='mineral_trough'+suffix,
                    initcond=scale_carbconc(pools_mineral_trough,carb_conc),hands_off=False,rateconstants=rateconstants,
                    bc=pools_atmoO2_organic,diffquo={'O2(aq)':O2_initial},truncate_concentration=truncate_conc,
                    CEC=CEC_mineral_trough,porosity=porosity_mineral_trough)    
        convert_to_xarray(result_mineral_trough,output_units).to_netcdf(filename,group='mineral_trough'+suffix,mode=filemode)
        filemode='a'
    
    if 'mineral_nottrough' == simtype:
        result_mineral_nottrough,output_units=run_alquimia.run_simulation(deckname,simlength,timestep,run_name='mineral_nottrough'+suffix,
                    initcond=scale_carbconc(pools_mineral_nottrough,carb_conc),hands_off=False,rateconstants=rateconstants,
                    bc=pools_atmoO2_organic,diffquo={'O2(aq)':O2_initial},truncate_concentration=truncate_conc,
                    CEC=CEC_mineral_nottrough,porosity=porosity_mineral_nottrough) 
        convert_to_xarray(result_mineral_nottrough,output_units).to_netcdf(filename,group='mineral_nottrough'+suffix,mode=filemode)
        filemode='a'


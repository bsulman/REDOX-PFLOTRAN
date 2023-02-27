# Guide to ELM-alquimia-PFLOTRAN

Benjamin Sulman

_Updated 2022-01-28_

## Quick start guide:
Set up a directory to store all this stuff:

        cd $HOME
        mkdir ELM-alquimia
        export BASEDIR=$HOME/ELM-alquimia

You can of course organize these codes however you want, as long as you make sure the paths in commands point to the right directories.

1.	Install software packages
    
    *	Follow directions at https://github.com/bsulman/REDOX-PFLOTRAN/blob/master/README_INSTALL. This will install and compile PFLOTRAN and alquimia and clone all the python scripts in REDOX-PFLOTRAN
    
    *	Install Offline Model Testbed (OLMT) and check out alquimia branch:

            cd $BASEDIR
            git clone https://github.com/dmricciuto/OLMT.git
            cd OLMT
            git checkout bsulman/coastal_main

        Note: Your python installation will need to have the netCDF4 package installed to run OLMT. netCDF4 should be installed with the python version you get on cades (`module load python`) but if you are using an anaconda python environment you may need to install the netCDF4 package.
    
    *	Clone/checkout correct ELM code:

            cd $BASEDIR
            git clone -b coastal_main_E3SMv2 git@github.com:bsulman/E3SM.git
            cd E3SM
            git submodule update --init --recursive

    *	To set up environment for ELM compilation, add the following to ~/.bashrc (assuming PFLOTRAN is installed in $HOME/ELM-alquimia):

            # Modules on CADES
            module load PE-gnu/1.0
            module load boost/1.61.0
            module load cmake/3.12.0
            module load zlib/1.2.8
            module load nco
            module load python
            module load openmpi/1.10.3
            module load git
            module load perl

            export CCSI_USERTOOLS=/software/user_tools/current/cades-ccsi
            
            # parallel-enabled hdf5 with openmpi-1.10.2/gcc5.3.0
            export HDF5_PATH=/software/dev_tools/swtree/cs400_centos7.2_pe2016-08/hdf5-parallel/1.8.17/centos7.2_gnu5.3.0
            export PATH=$HDF5_PATH/bin:$PATH
            export LD_LIBRARY_PATH=$HDF5_PATH/lib:$LD_LIBRARY_PATH
            
            
            module load mkl/2017
            BLASLAPACK_LIBDIR=/software/dev_tools/swtree/cs400_centos7.2_pe2016-08/mkl/2017/centos7.2_gnu5.3.0/lib

            # Make this your actual PFLOTRAN dir
            export PFLOTRAN_DIR=$HOME/ELM-alquimia/pflotran-elm-interface/src/pflotran
            export CLM_PFLOTRAN_SOURCE_DIR=$PFLOTRAN_DIR

            export PETSC_DIR=/software/user_tools/current/cades-ccsi/petsc-x
            export PETSC_ARCH=openmpi-1.10-gcc-5.3
            export PETSC_PATH=$PETSC_DIR/$PETSC_ARCH

2.	Generate PFLOTRAN input decks:

        cd $BASEDIR/REDOX-PFLOTRAN
        mkdir -p ELM_decks
        python network_for_ELM.py

    This command generates four PFLOTRAN input decks:
    * `CTC_alquimia_forELM.in`: Approximate recreation of normal ELM soil organic matter and inorganic N pools in PFLOTRAN without additional chemistry
    * `CTC_alquimia_forELM_adspinup.in`: Normal ELM pools with rate constants set for accelerated decomposition spinup
    * `CTC_alquimia_forELM_O2consuming.in`: ELM pools with oxygen consumption, dissolved CO2, DOM, and Fe cycling
    * `CTC_alquimia_forELM_O2consuming_adspinup.in`: ELM pools with oxygen consumption, dissolved CO2, DOM, and Fe cycling with modified rate constants for accelerated decomposition spinup

    The input deck generation uses a set of python scripts for inserting reactions and chemical species into a PFLOTRAN input deck template. The template it uses in this case is `SOMdecomp_template.txt`. The pools and reactions are all specified in the `network_for_ELM.py` script and can be changed by editing the script.

    Note: The C:N ratio of DOM is currently hard-coded into the ELM code, and reactions involving DOM are balanced in the input deck so changing either of those may cause N conservation errors and crash the model.

3.	Run simulation with OLMT:

        export BASEDIR=$HOME/ELM-alquimia
        cd $BASEDIR/OLMT
        mkdir -p ~/cases
        
        python site_fullrun.py --site US-PHM --caseidprefix test_alquimia  --nyears_ad_spinup 50 --nyears_final_spinup 50 --tstep 1 --cpl_bypass --machine cades --no_dynroot --spinup_vars --sitegroup Wetland --gswp3 --nyears_transient 51 --nofire --model_root $BASEDIR/E3SM --nopftdyn --ccsm_input /nfs/data/ccsi/proj-shared/E3SM/inputdata --caseroot ~/cases --runroot /lustre/or-scratch/cades-ccsi/$USER/  --mpilib openmpi --pio_version 2 --hist_nhtfrq_trans -24 --hist_mfilt_trans 365 --hist_mfilt_spinup 12 --hist_nhtfrq_spinup 0 --cn_only --alquimia $BASEDIR/REDOX-PFLOTRAN/ELM_decks/CTC_alquimia_forELM_O2consuming.in --alquimia_ad $BASEDIR/REDOX-PFLOTRAN/ELM_decks/CTC_alquimia_forELM_O2consuming_adspinup.in --trans_varlist "TOTVEGC,TOTSOMC,TOTLITC,soil_O2,HR,GPP,NEE,SMINN,SMINN_TO_PLANT,DIC_vr,SIC_vr,H2OSOI,watsat,SOIL1C_vr,SOIL2C_vr,SOIL3C_vr,SOIL4C_vr,LITR1C_vr,LITR2C_vr,LITR3C_vr,DOC_vr,soil_Fe2,soil_FeOxide,soil_pH,soil_sulfate,soil_sulfide,CH4_vr,chem_dt,SOILLIQ,SOILICE,QFLX_ADV"

    This should compile ELM and run a simulation with coupler bypass and alquimia turned on, using the expanded reaction network with oxygen, DOM, and iron, going through accelerated spinup, normal spinup, and historical simulations. It will run significantly slower than a normal ELM simulation.
    
    Notes:
    * `--cn_only`: PFLOTRAN setup does not yet support P cycling so we run with only C and N
    * `--alquimia`: This flag turns on alquimia compilation in OLMT and is followed by the path to the input deck that PFLOTRAN should use. In this case, the one for the expanded decomposition network that we just generated. You can run the other reaction network produced by `network_for_ELM.py` by switching this flag to the other input deck that was generated.
    * `--alquimia_ad`: The path after this specifies the input deck to use for accelerated decomposition spinup. If not specified, it uses the same deck for both. This is necessary because the ELM-alquimia-PFLOTRAN connection uses the decomposition rate constants from the input deck and cannot automatically change them for different spinup stages.
    * `--trans_varlist`: Turns on specific outputs for the transient (historical) simulation. New variables specific to alquimia include DOC, dissolved Fe(II), dissolved O2, pH, and the actual time stepping length that alquimia uses (`chem_dt`) after using variable time stepping to ensure a valid chemistry solution. Shorter `chem_dt` due to chemistry nonconvergence is the main reason the simulation might run significantly more slowly.
    * `--pio_version 2`: Parallel input-output version. The branch of this code off E3SM v1 only works with PIO version 1 on CADES, and the current branch off E3SM v2 used in these directions only works with PIO version 2.

4. Once simulation is finished, plot ten years of output (in this case, 1880-1889). Depends on an anaconda environment "myanaconda3" having been set up following the REDOX-PFLOTRAN installation instructions:

        cd $BASEDIR/REDOX-PFLOTRAN
        module load anaconda3
        conda activate myanaconda3
        python plot_ELM_alquimia_result.py /lustre/or-scratch/cades-ccsi/$USER/test_alquimia_US-PHM_ICB20TRCNRDCTCBC/run/test_alquimia*.h0.188?-02-01-00000.nc

    This script plots several figures on the screen so you will need to be logged into CADES with X-11 forwarding turnedon (`ssh -X` ...).
    Figures:
    * **Carbon time series**: Upper panel shows total vegetation, litter, and soil organic matter pools. Lower panel shows soil surface CO<sub>2</sub> flux which is calculated from the actual equilibration of surface soil CO<sub>2</sub> concentration with the atmosphere boundary condition.
    * **Time step**: This figure shows the actual time step lengths (in seconds) that the chemistry solver used. The chemistry solver starts at the ELM time step (3600 s) and cuts the length in half if the chemistry fails to reach a valid solution or the gas transport seems too fast for the time step length. Shorter chemistry time steps make the model run slower because it is solving the chemistry more times per ELM time step.
    * **Water and oxygen**: Heat map time series (time vs. depth) and profiles of soil water content, oxygen concentration, DOC concentration, and DIC concentration. The solid, dashed, and dotted lines on the heat maps correspond with the individual profiles on the profile plots.
    * **Redox**: Heat map time series and profiles for dissolved Fe(II), Fe oxide minerals, sulfate, sulfide, and pH.
    * **Carbon**: Heat map time series and profiles for total soil carbon, DOC, and methane with depth.



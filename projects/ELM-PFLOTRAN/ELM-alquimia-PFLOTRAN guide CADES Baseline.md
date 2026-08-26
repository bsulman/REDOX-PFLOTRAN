# Guide to ELM-alquimia-PFLOTRAN

Benjamin Sulman

_Updated 2024-05-20_

## Quick start guide:
Set up a directory to store all this stuff:

        cd $HOME
        mkdir ELM-alquimia

You can of course organize these codes however you want, as long as you make sure the paths in commands point to the right directories.

1.	Install software packages


    *	To set up environment for ELM compilation, add the following to ~/.bashrc:

            export PETSC_DIR=/ccsopen/proj/cli185/petsc-3.x-noopt/
            export PETSC_PATH=$PETSC_DIR
            export PFLOTRAN_DIR=/ccsopen/proj/cli185/b0u/pflotran
    
        Close and reopen your terminal after making these edits so they will take effect.

    * Install PFLOTRAN
        
            cd $HOME/ELM-alquimia
            git clone https://github.com/bsulman/pflotran-elm-interface.git
            pushd pflotran-elm-interface/src/pflotran
            git checkout pflotran-elm-interface
            make pflotran pflotran_rxn
            popd

    * Install python packages using the mamba python package manager

                wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
                bash Miniforge3-Linux-x86_64.sh

        You may need to close your terminal and log back in again for this to take effect before running the following commands (CIME does not seem to work for python > 3.11):

                mamba create --name myconda python=3.11 matplotlib xarray cartopy netcdf4 pandas numpy scipy ipython pygraphviz networkx cftime nc-time-axis openpyxl cffi
                mamba activate myconda
                pip install mpi4py

    * Install python scripts for handling alquimia/PFLOTRAN simulations

                git clone --recursive https://github.com/bsulman/REDOX-PFLOTRAN.git

                # Build alquimia interface
                # OPTIONAL since library is also now located at /nfs/data/ccsi/proj-shared/b0u/ELM-PFLOTRAN/alquimia
                # Note, this step will fail if environment variables PETSC_DIR and PETSC_ARCH are not correct
                # On CADES Baseline this step must be submitted using srun because it uses some MPI functions when testing the installation
                cd REDOX-PFLOTRAN/alquimia
                mkdir build
                cd build
                srun -N 1 -t 1:00:00 -A CLI185 -p batch cmake  .. \
                -DCMAKE_INSTALL_PREFIX=/ccsopen/proj/cli185/b0u/alquimia/ \
                -DCMAKE_C_COMPILER=$MPICC \
                -DCMAKE_CXX_COMPILER=$MPICXX \
                -DCMAKE_Fortran_COMPILER=$MPIF90 \
                -DCMAKE_BUILD_TYPE=Debug \
                -DXSDK_WITH_PFLOTRAN=ON \
                -DTPL_PFLOTRAN_LIBRARIES=$PFLOTRAN_DIR/libpflotranchem.a \
                -DTPL_PFLOTRAN_INCLUDE_DIRS=$PFLOTRAN_DIR \
                -DCMAKE_Fortran_FLAGS="-DPFLOTRAN_SOMDEC"

                make install


    *	Install Offline Model Testbed (OLMT) and check out alquimia branch:

            cd $HOME/ELM-alquimia
            git clone -b bsulman/coastal_main https://github.com/dmricciuto/OLMT.git

    
    *	Clone/checkout correct ELM code:

            cd $HOME/ELM-alquimia
            git -c http.sslVerify=false -c url."https://github.com/".insteadOf="git@github.com:" clone -b bsulman/lnd/tidal_multigrid --recurse-submodules https://github.com/bsulman/E3SM.git



2.	Generate PFLOTRAN input decks:

        cd $HOME/ELM-alquimia/REDOX-PFLOTRAN
        mkdir -p ELM_decks
        python network_for_ELM.py

    This command generates four PFLOTRAN input decks:
    * `CTC_alquimia_forELM.in`: Approximate recreation of normal ELM soil organic matter and inorganic N pools in PFLOTRAN without additional chemistry
    * `CTC_alquimia_forELM_adspinup.in`: Normal ELM pools with rate constants set for accelerated decomposition spinup
    * `CTC_alquimia_forELM_O2consuming.in`: ELM pools with oxygen consumption, dissolved CO2, DOM, and Fe cycling
    * `CTC_alquimia_forELM_O2consuming_adspinup.in`: ELM pools with oxygen consumption, dissolved CO2, DOM, and Fe cycling with modified rate constants for accelerated decomposition spinup

    The input deck generation uses a set of python scripts for inserting reactions and chemical species into a PFLOTRAN input deck template. The template it uses in this case is `SOMdecomp_template.txt`. The pools and reactions are all specified in the `network_for_ELM.py` script and can be changed by editing that script.

    Note: The C:N ratio of DOM is currently hard-coded into the ELM code, and reactions involving DOM are balanced in the input deck so changing either of those may cause N conservation errors and crash the model.

3.	Run simulation with OLMT. Several example simulations are in marsh_sim_cmds.txt

        cd $HOME/ELM-alquimia/OLMT
        mkdir -p ~/cases
        
        site=beo
        metdir=/gpfs/wolf2/cades/cli185/proj-shared/pt-e3sm-inputdata/atm/datm7/GSWP3_daymet/cpl_bypass_$site
        varlist="TOTVEGC,TOTSOMC,TOTLITC,SOIL1C_vr,SOIL2C_vr,SOIL3C_vr,SOIL4C_vr,LITR1C_vr,LITR2C_vr,LITR3C_vr,LEAFC,\
        soil_O2,HR,GPP,NEE,NPP,SMINN,SMINN_TO_PLANT,DIC_vr,SIC_vr,H2OSOI,H2OSFC,SOILLIQ,SOILICE,ZWT,QFLX_ADV,\
        QFLX_LAT_AQU,QFLX_EVAP_TOT,QVEGT,watsat,chem_dt,soil_salinity,SALINITY,soil_pH,DOC_vr,DIC_vr,DOC_RUNOFF,DIC_RUNOFF,SMIN_NO3_RUNOFF,\
        soil_sulfate,soil_sulfide,CH4_vr,CH4FLUX_ALQUIMIA,QDRAI,QDRAI_VR,H2OSFC_TIDE,TSOI,soil_Fe2,soil_FeOxide,soil_FeS"
        python site_fullrun.py --site AK-BEO --sitegroup NGEEArctic --caseidprefix Alaska_alquimia_7cell  \
                       --nyears_ad_spinup 100 --nyears_final_spinup 100 --tstep 1 --nyears_transient 151 \
                       --cpl_bypass --machine cades-baseline --no_dynroot --spinup_vars --gswp3 --daymet4 --nofire --nopftdyn --nopointdata \
                       --model_root $HOME/ELM-alquimia/E3SM --ccsm_input /gpfs/wolf2/cades/cli185/proj-shared/pt-e3sm-inputdata \
                       --metdir $metdir --walltime 10 \
                       --domainfile /gpfs/wolf2/cades/cli185/proj-shared/pt-e3sm-inputdata/share/domains/domain.clm/domain.lnd.7x1pt_beo-NGEE_areaC_navy.nc \
                       --surffile /gpfs/wolf2/cades/cli185/proj-shared/pt-e3sm-inputdata/lnd/clm2/surfdata_map/surfdata_7x1pt_beo-NGEE_areaC_simyr1850_c20230320.nc \
                       --caseroot ~/cases --runroot /gpfs/wolf2/cades/cli185/scratch/$USER/  --mpilib openmpi --pio_version 2 \
                       --hist_nhtfrq_trans -1 --hist_mfilt_trans 8760 --hist_mfilt_spinup 0 --hist_nhtfrq_spinup 12 --cn_only \
                       --trans_varlist $varlist \
                       --alquimia $HOME/ELM-alquimia/REDOX-PFLOTRAN/ELM_decks/CTC_alquimia_forELM_O2consuming.in \
                       --alquimia_ad $HOME/ELM-alquimia/REDOX-PFLOTRAN/ELM_decks/CTC_alquimia_forELM_O2consuming_adspinup.in \
                       --marsh --tide_forcing_file /nfs/data/ccsi/proj-shared/b0u/NGEE_ELM/BEO_hydro_BC_multicell.nc \
                       --parm_file /nfs/data/ccsi/proj-shared/b0u/NGEE_ELM/parms_BEO.txt

    This should compile ELM and run a simulation with coupler bypass and alquimia turned on, using the expanded reaction network with oxygen, DOM, and iron, going through accelerated spinup, normal spinup, and historical simulations. It will run significantly slower than a normal ELM simulation.
    
    Notes:
    * `--cn_only`: PFLOTRAN setup does not yet support P cycling so we run with only C and N
    * `--alquimia`: This flag turns on alquimia compilation in OLMT and is followed by the path to the input deck that PFLOTRAN should use. In this case, the one for the expanded decomposition network that we just generated. You can run the other reaction network produced by `network_for_ELM.py` by switching this flag to the other input deck that was generated.
    * `--alquimia_ad`: The path after this specifies the input deck to use for accelerated decomposition spinup. If not specified, it uses the same deck for both. This is necessary because the ELM-alquimia-PFLOTRAN connection uses the decomposition rate constants from the input deck and cannot automatically change them for different spinup stages.
    * `--trans_varlist`: Turns on specific outputs for the transient (historical) simulation. New variables specific to alquimia include DOC, dissolved Fe(II), dissolved O2, pH, and the actual time stepping length that alquimia uses (`chem_dt`) after using variable time stepping to ensure a valid chemistry solution. Shorter `chem_dt` due to chemistry nonconvergence is the main reason the simulation might run significantly more slowly.
    * `--pio_version 2`: Parallel input-output version. The branch of this code off E3SM v1 only works with PIO version 1 on CADES, and the current branch off E3SM v2 used in these directions only works with PIO version 2.

4. Once simulation is finished, plot ten years of output (in this case, 1880-1889). Depends on an anaconda environment "myconda3" having been set up following the REDOX-PFLOTRAN installation instructions:

        cd $HOME/ELM-alquimia/REDOX-PFLOTRAN
        conda activate myconda
        python plot_ELM_alquimia_result.py /lustre/or-scratch/cades-ccsi/$USER/Alaska_alquimia_7cell_AK-BEO_ICB20TRCNRDCTCBC/run/Alaska_alquimia_7cell_AK-BEO_ICB20TRCNRDCTCBC.h0.188?-02-01-00000.nc

    This script plots several figures on the screen so you will need to be logged into CADES with X-11 forwarding turned on (`ssh -X` ...).
    Figures:
    * **Carbon time series**: Upper panel shows total vegetation, litter, and soil organic matter pools. Lower panel shows soil surface CO<sub>2</sub> flux which is calculated from the actual equilibration of surface soil CO<sub>2</sub> concentration with the atmosphere boundary condition.
    * **Time step**: This figure shows the actual time step lengths (in seconds) that the chemistry solver used. The chemistry solver starts at the ELM time step (3600 s) and cuts the length in half if the chemistry fails to reach a valid solution or the gas transport seems too fast for the time step length. Shorter chemistry time steps make the model run slower because it is solving the chemistry more times per ELM time step.
    * **Water and oxygen**: Heat map time series (time vs. depth) and profiles of soil water content, oxygen concentration, DOC concentration, and DIC concentration. The solid, dashed, and dotted lines on the heat maps correspond with the individual profiles on the profile plots.
    * **Redox**: Heat map time series and profiles for dissolved Fe(II), Fe oxide minerals, sulfate, sulfide, and pH.
    * **Carbon**: Heat map time series and profiles for total soil carbon, DOC, and methane with depth.
    * **Hydro**: Surface water depth, porewater salinity, drainage, and lateral flow



from cffi import FFI
ffi_builder=FFI()
import os

alquimia_dir='alquimia/build/alquimia'
alquimia_include='alquimia/build/include'
pflotran_dir=os.environ['PFLOTRAN_DIR']
# pflotran_dir='/Users/b0u/Documents/Models/PFLOTRAN/pflotran-interface/src/pflotran'
# Should set it up to check that these paths actually work
petsc_paths=[]
if 'PETSC_DIR' in os.environ:
  petsc_paths.append(os.environ['PETSC_DIR'])
if 'PETSC_DIR' in os.environ and 'PETSC_ARCH' in os.environ:
  petsc_paths.append(os.path.join(os.environ['PETSC_DIR'],os.environ['PETSC_ARCH']))
if 'PETSC_PATH' in os.environ:
  petsc_paths.append(os.environ['PETSC_PATH'])

for pth in petsc_paths:
  if os.path.exists(os.path.join(pth,'lib')):
    petscpath=pth
    print(f'Found petsc library in: {petscpath}')
    break
else:
  raise RuntimeError('Failed to find PETSC in ',petsc_paths)

petsc_lib=os.path.join(petscpath,'lib')
petsc_archinclude=os.path.join(petscpath,'include')
mpi_include=os.path.join(os.environ['OPENMPI_DIR'],'include')
ffi_builder.set_source('_alquimia',
    r"""
    #include "petsc.h"
    #include "alquimia/alquimia_interface.h"
    #include "alquimia/alquimia_memory.h"
    #include "alquimia/alquimia_util.h"
    #include "alquimia/alquimia_containers.h"
    
    """,libraries=['alquimia','pflotranchem','petsc'],
        library_dirs=[alquimia_dir,pflotran_dir,petsc_lib],
        include_dirs=[os.path.join(os.environ['PETSC_DIR'],os.environ['PETSC_ARCH'],'include'),os.path.join(os.environ['PETSC_DIR'],'include'),petsc_archinclude,mpi_include,alquimia_include],
        extra_link_args=['-Wl,-rpath,alquimia/build/alquimia','-Wl,-rpath,'+petsc_lib],
        )
    
ffi_builder.cdef("""
        // Declarations shared between python and C

// alquimia_containers

  typedef struct {
    int size, capacity;
    double* data;
  } AlquimiaVectorDouble;

  typedef struct {
    int size, capacity;
    int* data;
  } AlquimiaVectorInt;

  typedef struct {
    /* NOTE: this is a vector of strings */
    int size, capacity;
    char** data;
  } AlquimiaVectorString;

  typedef struct {
    int num_primary;
    int num_sorbed;
    int num_minerals;
    int num_aqueous_complexes;
    int num_aqueous_kinetics;
    int num_surface_sites;
    int num_ion_exchange_sites;
    int num_isotherm_species;
    int num_aux_integers;
    int num_aux_doubles;
  } AlquimiaSizes;
  
  typedef struct {
    double water_density;  /* [kg/m^3] */
    double porosity;  /* [-] */
    double temperature;  /* [celsius] */
    double aqueous_pressure; /* [Pa] */
    AlquimiaVectorDouble total_mobile;  /* [molarity] */
    AlquimiaVectorDouble total_immobile;  /* [moles/m^3 bulk] */
    AlquimiaVectorDouble mineral_volume_fraction;  /* [-] */
    AlquimiaVectorDouble mineral_specific_surface_area; /* [m^2 mineral/m^3 bulk] */
    AlquimiaVectorDouble surface_site_density;  /* [moles/m^3 bulk] */
    AlquimiaVectorDouble cation_exchange_capacity;  /* [moles/m^3 bulk] */
  } AlquimiaState;
  
  typedef struct {
    double volume;  /* [m^3] */
    double saturation;  /* [m^3 liquid / m^3 pore space] */
    AlquimiaVectorDouble isotherm_kd;  /* [kg H20 / m^3 bulk] */
    AlquimiaVectorDouble freundlich_n; /* [-] */
    AlquimiaVectorDouble langmuir_b;  /* [-] */
    AlquimiaVectorDouble mineral_rate_cnst; /* [mol/m^2-sec] */    
    AlquimiaVectorDouble aqueous_kinetic_rate_cnst; /* [sec^-1] */
  } AlquimiaProperties;
  
  typedef struct {
    AlquimiaVectorInt aux_ints;  /* [-] */
    AlquimiaVectorDouble aux_doubles;  /* [-] */
  } AlquimiaAuxiliaryData;
  
  typedef struct {
    int error;
    char* message;
    bool converged;
    int num_rhs_evaluations;
    int num_jacobian_evaluations;
    int num_newton_iterations;
  } AlquimiaEngineStatus;
  
  typedef struct {
    bool thread_safe;
    bool temperature_dependent;
    bool pressure_dependent;
    bool porosity_update;
    bool operator_splitting;
    bool global_implicit;
    int index_base;
  } AlquimiaEngineFunctionality;
  
  typedef struct {
    AlquimiaVectorString primary_names;
    AlquimiaVectorInt    positivity;
    AlquimiaVectorString mineral_names;
    AlquimiaVectorString surface_site_names;
    AlquimiaVectorString ion_exchange_names;
    AlquimiaVectorString isotherm_species_names;
    AlquimiaVectorString aqueous_kinetic_names;
  } AlquimiaProblemMetaData;
  
  typedef struct {
    double pH;
    AlquimiaVectorDouble aqueous_kinetic_rate;  /* [?] */
    AlquimiaVectorDouble mineral_saturation_index;  /* [mol/sec/m^3] */
    AlquimiaVectorDouble mineral_reaction_rate;  /* [mol/sec/m^3 bulk] */
    AlquimiaVectorDouble primary_free_ion_concentration; /* [molality] */
    AlquimiaVectorDouble primary_activity_coeff; /* [-] */
    AlquimiaVectorDouble secondary_free_ion_concentration; /* [molality] */
    AlquimiaVectorDouble secondary_activity_coeff; /* [-] */
  } AlquimiaAuxiliaryOutputData;

  /* 
  ** Geochemical Conditions
  */

  typedef struct {
    char* primary_species_name;
    char* constraint_type;
    char* associated_species;
    double value;
  } AlquimiaAqueousConstraint;

  typedef struct {
    int size, capacity;
    AlquimiaAqueousConstraint* data;
  } AlquimiaAqueousConstraintVector;

  typedef struct {
    char* mineral_name;
    double volume_fraction;
    double specific_surface_area;
  } AlquimiaMineralConstraint;
  
  typedef struct {
    int size, capacity;
    AlquimiaMineralConstraint* data;
  } AlquimiaMineralConstraintVector;

  /* A geochemical condition is an array of aqueous and mineral geochemical constraints */
  typedef struct {
    char* name;
    AlquimiaAqueousConstraintVector aqueous_constraints;
    AlquimiaMineralConstraintVector mineral_constraints;
  } AlquimiaGeochemicalCondition;

  typedef struct {
    int size, capacity;
    AlquimiaGeochemicalCondition* data;
  } AlquimiaGeochemicalConditionVector;
  


// alquimia_interface.h

  /* NOTE(bja): AlquimiaData is a convenience container, and not a data
     structure that should be passed between the driver and
     engine. Conditions are not included here because they are short
     lived data structures. Need one data object per thread. */
  typedef struct {
    void* engine_state;
    AlquimiaSizes sizes;
    AlquimiaEngineFunctionality functionality;
    AlquimiaState state;
    AlquimiaProperties properties;
    AlquimiaAuxiliaryData aux_data;
    AlquimiaProblemMetaData meta_data;
    AlquimiaAuxiliaryOutputData aux_output;
  } AlquimiaData;

  /* NOTE(bja): The alquimia interface should contain nothing but
     function pointers. Inorder to thread chemistry, we need just one
     interface object. */ 
  typedef struct {
    /* read data files/structures, initialize memory, basis management
       (includes reading database, swapping basis, etc.) */
    void (*Setup)(
        const char* input_filename,
        bool hands_off,
        void* pft_engine_state,
        AlquimiaSizes* sizes,
        AlquimiaEngineFunctionality* functionality,
        AlquimiaEngineStatus* status);

    /* gracefully shutdown the engine, cleanup memory */
    void (*Shutdown)(
      void* pft_engine_state,
      AlquimiaEngineStatus* status);

    /* constrain processing for boundary/initial constraints. Called
       once for each IC/BC. */
    void (*ProcessCondition)(
        void* pft_engine_state,
        AlquimiaGeochemicalCondition* condition,
        AlquimiaProperties* props,
        AlquimiaState* state,
        AlquimiaAuxiliaryData* aux_data,
        AlquimiaEngineStatus* status);

    /* take one (or more?) reaction steps in operator split mode */
    void (*ReactionStepOperatorSplit)(
        void* pft_engine_state,
        double delta_t,
        AlquimiaProperties* props,
        AlquimiaState* state,
        AlquimiaAuxiliaryData* aux_data,
        AlquimiaEngineStatus* status);
    
    /* Access to user selected geochemical data for output, i.e. pH, 
       mineral SI, reaction rates */
    void (*GetAuxiliaryOutput)(
        void* pft_engine_state,
        AlquimiaProperties* props,
        AlquimiaState* state,
        AlquimiaAuxiliaryData* aux_data,
        AlquimiaAuxiliaryOutputData* aux_out,
        AlquimiaEngineStatus* status);
    
    void (*GetProblemMetaData)(
        void* pft_engine_state,
        AlquimiaProblemMetaData* meta_data,
        AlquimiaEngineStatus* status);
  } AlquimiaInterface;


  void CreateAlquimiaInterface(const char* const engine_name,
                               AlquimiaInterface* interface,
                               AlquimiaEngineStatus* status);

  // alquimia_memory.h
  
  /* Alquimia Vectors */
  void AllocateAlquimiaVectorDouble(const int size, AlquimiaVectorDouble* vector);
  void FreeAlquimiaVectorDouble(AlquimiaVectorDouble* vector);

  void AllocateAlquimiaVectorInt(const int size, AlquimiaVectorInt* vector);
  void FreeAlquimiaVectorInt(AlquimiaVectorInt* vector);

  void AllocateAlquimiaVectorString(const int size, AlquimiaVectorString* vector);
  void FreeAlquimiaVectorString(AlquimiaVectorString* vector);

  /* State */
  void AllocateAlquimiaState(const AlquimiaSizes* const sizes,
                             AlquimiaState* state);

  void FreeAlquimiaState(AlquimiaState* state);

  /* Auxiliary Data */ 
  void AllocateAlquimiaAuxiliaryData(const AlquimiaSizes* const sizes,
                                     AlquimiaAuxiliaryData* aux_data);
  void FreeAlquimiaAuxiliaryData(AlquimiaAuxiliaryData* aux_data);

  /* Properties */
  void AllocateAlquimiaProperties(const AlquimiaSizes* const sizes,
                                  AlquimiaProperties* props);
  void FreeAlquimiaProperties(AlquimiaProperties* props);

  /* Problem Meta Data */
  void AllocateAlquimiaProblemMetaData(const AlquimiaSizes* const sizes,
                                       AlquimiaProblemMetaData* meta_data);

  void FreeAlquimiaProblemMetaData(AlquimiaProblemMetaData* metda_data);

  /* Status */
  void AllocateAlquimiaEngineStatus(AlquimiaEngineStatus* status);

  void FreeAlquimiaEngineStatus(AlquimiaEngineStatus* status);

  /* Auxiliary Output Data */ 
  void AllocateAlquimiaAuxiliaryOutputData(const AlquimiaSizes* const sizes,
                                           AlquimiaAuxiliaryOutputData* aux_output);
  void FreeAlquimiaAuxiliaryOutputData(AlquimiaAuxiliaryOutputData* aux_output);

  /* Geochemical conditions/constraints */
  void AllocateAlquimiaGeochemicalConditionVector(const int num_conditions,
                                                  AlquimiaGeochemicalConditionVector* condition_list);
  void AllocateAlquimiaGeochemicalCondition(const int size_name,
                                            const int num_aqueous_constraints, 
                                            const int num_mineral_constraints,
                                            AlquimiaGeochemicalCondition* condition);
  void AllocateAlquimiaAqueousConstraintVector(const int num_constraints,
                                               AlquimiaAqueousConstraintVector* constraint_list);
  void AllocateAlquimiaAqueousConstraint(AlquimiaAqueousConstraint* constraint);
  void AllocateAlquimiaMineralConstraintVector(const int num_constraints,
                                               AlquimiaMineralConstraintVector* constraint_list);
  void AllocateAlquimiaMineralConstraint(AlquimiaMineralConstraint* constraint);

  void FreeAlquimiaGeochemicalConditionVector(AlquimiaGeochemicalConditionVector* condition_list);
  void FreeAlquimiaGeochemicalCondition(AlquimiaGeochemicalCondition* condition);
  void FreeAlquimiaAqueousConstraintVector(AlquimiaAqueousConstraintVector* vector);
  void FreeAlquimiaAqueousConstraint(AlquimiaAqueousConstraint* constraint);
  void FreeAlquimiaMineralConstraintVector(AlquimiaMineralConstraintVector* vector);
  void FreeAlquimiaMineralConstraint(AlquimiaMineralConstraint* constraint);


  /* Data */
  void AllocateAlquimiaData(AlquimiaData* data);
  void FreeAlquimiaData(AlquimiaData* data);
  
  // PETSC for initializing MPI
  typedef int PetscErrorCode; // from petscsystypes.h
  typedef enum { PETSC_FALSE,PETSC_TRUE } PetscBool; // from petscsystypes.h
  //PetscErrorCode PetscInitialize(int*,char***,const char[],const char[]);  // from petscsys.h
  PetscErrorCode PetscInitializeNoArguments(void);
  PetscErrorCode PetscFinalize(void);
  PetscErrorCode PetscInitialized(PetscBool*);

    """)
    
ffi_builder.compile(verbose=True)

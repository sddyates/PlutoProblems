#define  PHYSICS                 MHD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                SPHERICAL
#define  BODY_FORCE              VECTOR
#define  COOLING                 NO
#define  RECONSTRUCTION          LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 2
#define  USER_DEF_PARAMETERS     11

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  DIVB_CONTROL            DIV_CLEANING
#define  BACKGROUND_FIELD        YES
#define  RESISTIVITY             NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          YES

/* -- user-defined parameters (labels) -- */

#define  seperation              0
#define  planet_temperature      1
#define  planet_mass             2
#define  planet_radius           3
#define  planet_Bfield           4
#define  planet_surface_rho      5
#define  star_temperature        6
#define  star_mass               7
#define  star_radius             8
#define  star_Bfield             9
#define  star_surface_rho        10

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY            1.0e-15
#define  UNIT_LENGTH             6.955e+10
#define  UNIT_VELOCITY           1.0e5
#define  UNIT_MASS               (UNIT_DENSITY*pow(UNIT_LENGTH,3))
#define  UNIT_TIME               (UNIT_LENGTH/UNIT_VELOCITY)
#define  UNIT_MASS               (UNIT_DENSITY*pow(UNIT_LENGTH,3))
#define  UNIT_G                  (CONST_G*((1/pow(UNIT_LENGTH,3))*UNIT_MASS*pow(UNIT_TIME,2)))
#define  UNIT_kB                 ((CONST_kB*pow(UNIT_TIME,2))/(UNIT_DENSITY*pow(UNIT_LENGTH,5)))
#define  UNIT_B                  (sqrt(4.0*CONST_PI*UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY))
#define  tyear                   3.15569e+7
#define  tday                    8.64e+4
#define  VTK_TIME_INFO           YES
#define  VTK_VECTOR_DUMP         NO
#define  GLM_EXTENDED            NO
#define  CHOMBO_REF_VAR          RHO
#define  PPM_ORDER               4
#define  SHOW_TIME_STEPS         NO
#define  CAK                     NO
#define  CHOMBO_LOGR             YES
#define  CHOMBO_CONS_AM          NO

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING         NO
#define  WARNING_MESSAGES          NO
#define  PRINT_TO_FILE             NO
#define  INTERNAL_BOUNDARY         NO
#define  SHOCK_FLATTENING          NO
#define  CHAR_LIMITING             NO
#define  LIMITER                   VANLEER_LIM
#define  ASSIGN_VECTOR_POTENTIAL   NO

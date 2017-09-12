#define  PHYSICS                 MHD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              VECTOR
#define  COOLING                 NO
#define  RECONSTRUCTION          LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     6

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  DIVB_CONTROL            DIV_CLEANING
#define  BACKGROUND_FIELD        YES
#define  RESISTIVITY             NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  M_star                  0
#define  T_star                  1
#define  R_star                  2
#define  B_star                  3
#define  Omega                   4
#define  RHO_star                5

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY            1.0e-15
#define  UNIT_LENGTH             0.28*6.955e+10
#define  UNIT_VELOCITY           1.0e5
#define  t0                      (UNIT_LENGTH/UNIT_VELOCITY)
#define  M0                      (UNIT_DENSITY*pow(UNIT_LENGTH,3))
#define  UNIT_G                  (CONST_G*((1/pow(UNIT_LENGTH,3))*M0*pow(t0,2)))
#define  Bc                      (sqrt(4.0*CONST_PI*UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY))
#define  tyear                   3.15569e+7
#define  tday                    8.64e+4
#define  VTK_TIME_INFO           YES
#define  VTK_VECTOR_DUMP         YES
#define  GLM_EXTENDED            NO
#define  CHOMBO_REF_VAR          RHO
#define  PPM_ORDER               4
#define  SHOW_TIME_STEPS         NO
#define  CAK                     NO

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING         NO
#define  WARNING_MESSAGES          NO
#define  PRINT_TO_FILE             YES
#define  INTERNAL_BOUNDARY         YES
#define  SHOCK_FLATTENING          NO
#define  CHAR_LIMITING             NO
#define  LIMITER                   VANLEER_LIM
#define  ASSIGN_VECTOR_POTENTIAL   NO

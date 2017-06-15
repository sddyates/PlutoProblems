#define  PHYSICS                 MHD
#define  DIMENSIONS              2
#define  COMPONENTS              2
#define  GEOMETRY                POLAR
#define  BODY_FORCE              VECTOR
#define  COOLING                 POWER_LAW
#define  RECONSTRUCTION          PARABOLIC
#define  TIME_STEPPING           RK3
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     11

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  DIVB_CONTROL            CONSTRAINED_TRANSPORT
#define  BACKGROUND_FIELD        YES
#define  RESISTIVITY             NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  Cs_P                    0
#define  M_RATIO                 1
#define  L_RATIO                 2
#define  B_CGS                   3
#define  TT                      4
#define  MU                      5
#define  AA                      6
#define  Bb                      7
#define  QQ                      8
#define  aa_eff                  9
#define  BB                      10

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY            1.0e-12
#define  UNIT_LENGTH             6.955e+10*9.0
#define  UNIT_VELOCITY           1.0e+5
#define  UNIT_MASS               (UNIT_DENSITY*pow(UNIT_LENGTH,3))
#define  UNIT_TIME               (UNIT_LENGTH/UNIT_VELOCITY)
#define  UNIT_G                  (CONST_G*UNIT_DENSITY*pow(UNIT_TIME,2))
#define  UNIT_kB                 ((CONST_kB*pow(UNIT_TIME,2))/(UNIT_DENSITY*pow(UNIT_LENGTH,5)))
#define  UNIT_B                  (sqrt((4.0*CONST_PI*UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY)))
#define  UNIT_L                  (pow(UNIT_TIME,-3)*(UNIT_DENSITY*pow(UNIT_LENGTH,5)))
#define  tyear                   3.15569e+7
#define  tday                    8.64e+4
#define  L_sun                   3.846e+33
#define  VTK_VECTOR_DUMP         YES
#define  GLM_EXTENDED            YES
#define  CAK                     YES
#define  CHOMBO_LOGR             YES
#define  CHOMBO_CONS_AM          NO

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING         NO
#define  WARNING_MESSAGES          NO
#define  PRINT_TO_FILE             YES
#define  INTERNAL_BOUNDARY         YES
#define  SHOCK_FLATTENING          NO
#define  CHAR_LIMITING             NO
#define  LIMITER                   VANLEER_LIM
#define  CT_EMF_AVERAGE            UCT_HLL
#define  CT_EN_CORRECTION          NO
#define  ASSIGN_VECTOR_POTENTIAL   NO
#define  UPDATE_VECTOR_POTENTIAL   NO

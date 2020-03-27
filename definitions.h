! definitions : Directives processeurs
! en 2D : DIMENSION_GEO = 2
! en 3D : DIMENSION_GEO = 3
#define DIMENSION_GEO 2

! uniform grid = 1
! non-uniform grid = 0
#define UNIFORM_GRID 1

! continue_calculate with a result         = 1
! dont continue_calculate with a result    = 0
#define CONTINUE_CAL 0

! Cas test
! Tourbillon de Green Taylor :        T_G_T     = 1
! Cavité Entrainée :                    C_E     = 2
! Cavité Différentiellement Chauffée: C_D_C     = 3
#define TYPE_ECOULEMENT 3
#define SYMETRIE 0

! projection     incémentale = 1
! projection non-incémentale = 0 
#define INCREMENTALE 1

! Euler1 semi-implicite schema_temps = 1  ! disponible
! BDF2 Goda             schema_temps = 2  ! disponible

#define SCHEMA_TEMPS 2

! avec l équation de l enérgie : ENRG = 1 
!                      , sinon : ENRG = 0
#define ENRG 1
  
! avec la poussée d Archimède : ARCHIMEDE = 1
!                     , sinon : ARCHIMEDE = 0
#define ARCHIMEDE 1

! Newtonien                     = 0  ! disponible
! Oldroyd-B                     = 1  ! disponible
! Giesekus                      = 2  ! disponible 
!  PTT                          = 3  ! disponible  
#define VISCOELASTIC_MODELS 1

! no_exp_treatment_for_Tau      = 0
!    exp_treatment_for_Tau      = 1
#define EXP_TREATMENT 0

! traitement des termes convectifs
!if VISCOELASTIC_MODELS==0 
! Différences centrées    
!else
! WENO_3   = 10 
! WENO_5   = 11 
! WENO_7   = 12 
#define NL_TERMS 10

! central_02 = 1
! central_04 = 2
#define CENTRAL_D 2

! avec HOUC = 1 
! sans HOUC = 0
#define HOUC 1

! Quasi_Linear = 1
! Central_Difference = 0
#define QUASI_LINEAR 1

! seulement enérgie et avec u:=(u1,u2,u3) consatant
#define ONLY_ENRG 0


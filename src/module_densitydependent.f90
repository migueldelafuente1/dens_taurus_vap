!==============================================================================!
! MODULE DensityDep                                                            !
!                                                                              !
! This module contains subroutines to evaluate the spatial density to evaluate !
! density dependent interactions.                                              !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_densty_dependent                                            !
! - subroutine reduced_me_Yk                                                   !
! - function calculate_expectval_density                                       !
!                                                                              !
! Testing Methods                                                              !
! - subroutine test_printDesityWF                                              !
! - subroutine test_integrate_bulk_densities                                   !
!==============================================================================!
MODULE DensityDep

use Constants
use MathMethods
use Basis
use Wavefunctions
!use Hamiltonian
use Fields
use Lebedev


implicit none
PUBLIC

logical   :: EVAL_DENSITY_DEPENDENT = .TRUE.
logical   :: EVAL_REARRANGEMENT     = .FALSE.
logical   :: EVAL_EXPLICIT_FIELDS_DD= .FALSE.
real(r64) :: t3_DD_CONST  = 0.0!= 1350.00d+00    ! constant of the DD term [MeV]
real(r64) :: x0_DD_FACTOR = 0.0!= 1.0d+00        ! exchange factor of DD term
real(r64) :: alpha_DD     = 0.0!=  0.33333d+00   ! power of the DD term
integer, dimension(2) ::   alpha_DD_frac         ! fraction (num, den) of the power
! 0 trapezoidal, 1 Gauss-Legendre, 2 Gauss-Laguerre(r)/Legendre (t,p), 3 Laguerre-Lebedev
logical   :: export_density      = .FALSE. ! .TRUE. !

integer   :: r_dim       = 0 != 10  ! dimension of the r array
integer   :: Omega_Order = 0 != 10 ! 4

integer   :: theta_dim      ! dimension of the theta array
integer   :: phi_dim        ! dimension of the phi array
integer   :: Leb_dim
integer   :: angular_dim    ! dimension of the 1 dim array ThetaPhi

! Quadratures for integration
real(r64), dimension(:), allocatable   :: x_R
real(r64), dimension(:), allocatable   :: weight_R
real(r64), dimension(:,:), allocatable :: x_Leb
real(r64), dimension(:), allocatable   :: weight_LEB

real(r64), dimension(:), allocatable  :: r, r_export
real(r64), dimension(:), allocatable  :: theta, theta_export
real(r64), dimension(:), allocatable  :: cos_th, cos_th_export
real(r64), dimension(:), allocatable  :: phi, phi_export

complex(r64), dimension(:,:), allocatable :: rhoLR, kappaRL, kappaLR

complex(r64), dimension(:,:), allocatable :: density, density_export !(ir, iang)
complex(r64), dimension(:,:), allocatable :: dens_alpha, dens_alpm1
complex(r64), dimension(:,:,:), allocatable :: dens_pnt   !(pp, nn, pn, np, total;i_r,iang)
complex(r64), dimension(:,:), allocatable :: density_export_p, density_export_n
complex(r64), dimension(:,:), allocatable :: pairdens_export
complex(r64), dimension(:,:), allocatable :: pairdens_export_n,pairdens_export_p

!! Pre-calculated saved functions and coefficients
real(r64), dimension(:, :, :), allocatable, save :: radial_2b_sho_memo        !(ish1, ish2, ir)
real(r64), dimension(:, :, :), allocatable, save :: radial_2b_sho_export_memo !(ish1, ish2, ir)

complex(r64), dimension(:,:), allocatable,   save :: sph_harmonics_memo ! Y_(km_indx, i_ang)
complex(r64), dimension(:,:,:), allocatable, save :: sphharmDUAL_memo   ! Y*(a) Y(b)
real(r64), dimension(:,:,:), allocatable,    save :: dens_Y_KM_me       ! <a|Y_kb|b>(jm_a, jm_b, KM indx)

complex(r64), dimension(:,:,:,:), allocatable, save :: AngFunctDUAL_HF ! CGa CGb Y*(a) Y (b) [(++,+-,-+,--)], a, b, i_ang
complex(r64), dimension(:,:,:,:), allocatable, save :: AngFunctDUAL_P1 ! CGa CGb Y (a) Y (b) [(++,+-,-+,--)], a, b, i_ang
complex(r64), dimension(:,:,:,:), allocatable, save :: AngFunctDUAL_P2 ! CGa CGb Y*(a) Y*(b) [(++,+-,-+,--)], a, b, i_ang
complex(r64), dimension(:,:,:,:), allocatable, save :: BulkHF ! (tt,msms',r,ang) DEF:Sum AngFunctDUAL_HF * rho [pp,nn,pn,np, tot]  [(++,+-,-+,--)]
complex(r64), dimension(:,:,:,:), allocatable, save :: BulkP1 ! (tt,msms',r,ang) DEF:Sum AngFunctDUAL_P1(msms') - AngFunctDUAL_P1(ms'ms) * kappaLR (p, n), pn=function xM xH
complex(r64), dimension(:,:,:,:), allocatable, save :: BulkP2 ! (tt,msms',r,ang) DEF:Sum AngFunctDUAL_P2 * kappaRL
complex(r64), dimension(:,:,:,:), allocatable, save :: BulkP1_HM ! (tt,msms',r,ang) DEF:Sum AngFunctDUAL_P1(msms') - AngFunctDUAL_P1(ms'ms) * kappaLR (p, n), pn=function xM xH

complex(r64), dimension(:,:), allocatable     :: rearrangement_me  !(isp1, isp2)
complex(r64), dimension(:,:), allocatable     :: rearrang_field    !(isp1, isp2)
complex(r64), dimension(:,:,:,:), allocatable :: rea_common_RadAng !(isp1,isp2, ir,iang)
complex(r64), dimension(:,:), allocatable     :: REACommonFields   !(ir, iang))
complex(r64), dimension(:,:), allocatable     :: fixed_rearrang_field

!!! Arrays related to Differentiated density.
logical   :: EXPORT_GRAD_DD = .FALSE. ! export the Laplacian-approximation of the Rearrangement
logical :: implement_H2cpd_DD = .FALSE.

real(r64) :: last_HFB_energy

integer, dimension(:), allocatable :: HOsh_ant, HOsp_ant

!! Related to the hamiltonian and Fields
real(r64) :: REARRANGEMENT_ENERGY = 0.0d+00
integer   :: iteration = 0, global_iter_max = 0
integer   :: VSsh_dim = 0, VSsp_dim = 0, VSsp_dim2 = 0
integer   :: WBsh_dim = 0, WBsp_dim = 0
integer, dimension(:), allocatable  :: VSsh_list, WBtoHOsp_index, VSsp_VSsh,&
                                       VStoHOsp_index, VStoHOsh_index

real(r64), dimension(2) :: lambdaFer_DD ! Lag. mult. before (lambda-Z, lambda-N)


integer   :: Mphip_DD = 0, & ! number of angles in the projection for protons
             Mphin_DD = 0    !   "    "    "    "   "      "       "  neutrons
integer   :: seed_type_sym  = 0        ! (UNDEFINED)
logical   :: hasX0M1        = .FALSE.  ! |x0 - 1| > 1E-6
logical   :: evalFullSPSpace= .TRUE.   ! compute the full a,b, c,d space for explicit DD fields (cannot be set)
logical   :: exportValSpace = .FALSE.  ! the export of a reduced val.space
logical   :: evalQuasiParticleVSpace = .FALSE. ! Export for the QP sp states, not the VS procedure

integer   :: NHO_vs, NHO_co !! Major Shell number of the Valence. Sp to be exported
logical   :: NOT_DEL_FILE
logical   :: PRINT_GUTS = .TRUE.
logical   :: DOING_PROJECTION = .FALSE.
logical   :: USING_FIXED_REARRANGEMENT = .FALSE.
logical   :: EVAL_CUTOFF = .FALSE.
integer   :: CUTOFF_MODE = 0    ! 1: Energy Cutoff, 2: Kappa cutoff, 3: both
real(r64) :: CUTOFF_ENERGY_MAX = 1.0d+99
real(r64) :: CUTOFF_KAPPA      = 0.0
logical   :: CALCULATE_DD_PN_PA   = .TRUE.
logical   :: CALCULATE_DD_PN_HF   = .TRUE.

! Alternative Density dependent Profiles
integer   :: FUNCTIONAL_DENS_MODE = 1  ! 1=D1/D1S(default), 2=Nucl.Matter x0=0
real(r64) :: CONST_EDD_M2_ETA  = 0.0D+00
real(r64) :: CONST_EDD_M2_RHO0 = 0.1381553D+00  ! r0=1.2 fm
real(r64) :: CONST_EDD_M2      = 0.0D+00

logical   :: has_HEIS_MAJO_TERMS  = .FALSE.
real(r64) :: CONST_x0_EXC_HEIS = 0.0D+00
real(r64) :: CONST_x0_EXC_MAJO = 0.0D+00

logical   :: has_sevaral_DD_tems = .FALSE.
integer   :: number_DD_terms     = 1
real(r64), dimension(5,5) :: parameters_alphx0x0Hx0M  !! (alp,t3,x0, x0H,x0M)[term]

!! [END] DENSITY DEPENDENT MODIFICATIONS =====================================

CONTAINS

!-----------------------------------------------------------------------------!
! Import and define array dimensions from file DD_PARAMS.txt                  !
!-----------------------------------------------------------------------------!
subroutine import_DD_parameters

integer :: runit = 99
integer :: ios, i, seed_type_imported, aa, a, a_ant
integer :: extra_params = 0
logical :: is_exist
real(r64) :: aux_float
character(len=*), parameter :: formatST = "(1a)", &
                               formatI1 = "(1a30, 1i1)", &
                               formatI2 = "(1a30, 1i2)", &
                               formatI3 = "(1a30, 1i3)", &
                               formatF6 = "(1a30, 1f9.6)", &
                               formatEE = "(1a30, 1es12.6)", &
                               formatStrHeader = "(1a30)", &
                               formatII = "(1a30, 1i1, 99i6)", &
                               formatI3F6="(1a30, 1i3, 1es12.6)"

CHARACTER(LEN=20) :: file_input = "input_DD_PARAMS.txt"
CHARACTER(LEN=30) :: filecontents, str_
INTEGER   :: io, aux_int

inquire (file=file_input, exist=is_exist)
print *,     ""
print '(A)', "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
print '(A)', "                DENSITY DEPENDENT PARAMETERS                "
print '(A)', "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
print *,     ""

if ( is_exist ) then
  OPEN(runit, FILE=file_input, FORM="FORMATTED", STATUS="OLD", ACTION="READ")
else
  print "(a,a)", "The DD file (input_DD_PARAMS.txt) can not be found. ", &
    "The module density dependent will be skipped."
  EVAL_DENSITY_DEPENDENT = .FALSE.
  export_density = .FALSE.
  return
endif

!! Reads the input parameters
read(runit,formatST) str_
read(runit,formatI1) str_, aux_int
EVAL_DENSITY_DEPENDENT = aux_int.EQ.1
read(runit,formatI1) str_, aux_int
EVAL_REARRANGEMENT = aux_int.EQ.1
read(runit,formatI1) str_, aux_int
EVAL_EXPLICIT_FIELDS_DD = aux_int.EQ.1

read(runit,formatEE) str_, t3_DD_CONST
read(runit,formatEE) str_, x0_DD_FACTOR
read(runit,formatEE) str_, alpha_DD
read(runit,formatST) str_
read(runit,formatI3) str_, r_dim
read(runit,formatI3) str_, Omega_Order
read(runit,formatI1) str_, aux_int
export_density = aux_int.EQ.1
read(runit,formatI1) str_, aux_int
evalQuasiParticleVSpace = aux_int.GE.1
read(runit,formatI2) str_, aux_int
exportValSpace = aux_int.GE.1
!! Auxiliary additional parameters
read(runit,formatI2) str_, extra_params
if (extra_params .GT. 0) then
  print "(A,I3,A)", "  Additional ", extra_params, " Parameters-Constants."
  do i=1, extra_params
    read(runit,formatI3F6) str_, aux_int, aux_float
    call set_extra_DD_parameters(aux_int, aux_float)
  end do
end if

VSsh_dim = aux_int
if (exportValSpace) then
  if ((VSsh_dim.LE.HOsh_dim).OR.(evalQuasiParticleVSpace)) then
    print "(A,I3,2A)", "   ... Reading VS ", VSsh_dim, " sh states", &
      " (error if wrong sh dimension)"
    !print *, ""
    backspace runit
    allocate(VSsh_list(VSsh_dim))
    read(runit,formatStrHeader, advance='no') str_
    read(runit,*) VSsh_dim, (VSsh_list(i),i=1,VSsh_dim)
    call set_valence_space_to_export
  else
    !print *,"   ... Reading for FULL valence space"
    VSsh_dim  = HOsh_dim
    VSsp_dim  = HOsp_dim
    VSsp_dim2 = HOsp_dim2

    WBsh_dim  = HOsh_dim
    WBsp_dim  = HOsp_dim
    allocate(VSsh_list(VSsh_dim))
    allocate(VStoHOsp_index(VSsp_dim), VSsp_VSsh(VSsp_dim))
    allocate(WBtoHOsp_index(HOsp_dim))
    do i=1, VSsh_dim
      VSsh_list(i) = HOsh_na(i)
      VSsp_VSsh(i) = i
    end do
    do i=1, VSsp_dim
      VStoHOsp_index(i) = i
      WBtoHOsp_index(i) = i
    end do
  endif
endif

!! Allocate here the Antoine index (for the Valence Space exporting) --------
allocate(HOsh_ant(HOsh_dim))
allocate(HOsp_ant(HOsp_dim))
HOsh_ant = 0
HOsp_ant = 0
do aa = 1, HOsh_dim
  HOsh_ant(aa) = 10000*HOsh_n(aa) + 100*HOsh_l(aa) + HOsh_2j(aa)
enddo
do aa = 1, HOsp_dim
  HOsp_ant(aa) = 10000*HOsp_n(aa) + 100*HOsp_l(aa) + HOsp_2j(aa)
enddo
!!! -------------------------------------------------------------------------

read(runit,formatST) str_
rewind(runit)
close(runit)

!! set alpha more exact 1/3 double precision, and the fraction (is always dimensional)
if (abs(alpha_DD - (1.0D0/3.0D0)) < 1.0D-5) alpha_DD = 0.333333333333333
if (abs(x0_DD_FACTOR - 1.0D0) < 1.0D-5) x0_DD_FACTOR = 1.000000000000000

alpha_DD_frac = 0
if (abs(alpha_DD - (1.0/3.0)) < 1.0D-5) then
  alpha_DD_frac = (/1, 3/)
else if ((abs(alpha_DD - (1.0)) < 1.0D-5) .OR. (abs(alpha_DD) < 1.0D-4)) then
  alpha_DD_frac = (/ int(alpha_DD), 1/)
else if (abs(alpha_DD - (2.0/3.0)) < 1.0D-5) then
  alpha_DD_frac = (/2, 3/)
else if (abs(alpha_DD - (1.0/6.0)) < 1.0D-5) then
  alpha_DD_frac = (/1, 6/)
else
  print '(A,F9.5)', "[WARGNING] Awkward ALPHA constant for DD-EDF: ", alpha_DD
endif


allocate(x_R(r_dim))
allocate(weight_R(r_dim))
allocate(r(r_dim))
allocate(r_export(r_dim))
print *,           '   Density dep. Interaction values  '
print *,           '-----------------------------------------------'
print '(A,L10)',   'eval_density_dependent =', EVAL_DENSITY_DEPENDENT
print '(A,L10)',   'eval_rearrangement     =', EVAL_REARRANGEMENT
print '(A,L10)',   'eval_explicit_fieldsDD =', EVAL_EXPLICIT_FIELDS_DD
print '(A,F10.4)', 't3_DD_CONST (MeV fm-4) =', t3_DD_CONST
print '(A,F10.6)', 'x0_DD_FACTOR           =', x0_DD_FACTOR
print '(A,F10.6)', 'alpha_DD               =', alpha_DD
if (FUNCTIONAL_DENS_MODE .NE. 1) then
  print '(A,I10)',   '  MODE-density_dependent =', FUNCTIONAL_DENS_MODE
  print '(A,F10.4)', '  * Eta - factor (M2)  =', CONST_EDD_M2_ETA
  print '(A,F10.4)', '  * Dens rho0    (M2)  =', CONST_EDD_M2_RHO0
end if
print *, ''
print *,           '   Integration parameters  '
print *,           '-----------------------------------------------'
print '(A,I10)',   'r_dim              =', r_dim
print '(A,I10)',   'Omega_Order        =', Omega_Order
print '(A,L10)',   'export_density     =', export_density
print '(A,2L5)',   'eval/export Val.Sp =', evalFullSPSpace, exportValSpace
print '(A,2L10)',  'export QP Val.Sp   =', evalQuasiParticleVSpace
print *,           '-----------------------------------------------'
print *, ''

if ((.NOT.exportValSpace).AND.(implement_H2cpd_DD)) then
  if (hamil_read .EQ. 0) then
    !deallocate(hamil_H2cpd_DD) ! It wont be used
    !print "(A)", "  Hamiltonian cpd deallocated  because it will not be used!"
    continue
  endif
else if ((exportValSpace).AND.(.NOT.implement_H2cpd_DD)) then
  print "(2A)", " ERROR, do not export the matrix elements with hamiltonian",&
                " of type=1 or 2, program stops"
  STOP
else if (((exportValSpace).OR.(evalQuasiParticleVSpace)) &
         .AND.(hamil_read.EQ.1)) then
  print "(2A)", " ERROR, do not export the matrix elements when importing ",&
        " the reduced matrix elements (set hamil_read=0), program stops"
  STOP
endif
if (exportValSpace)then
  print '(A,2I4)', '    ... sh states to export DIM(sh/sp):',VSsh_dim,VSsp_dim
  allocate(VStoHOsh_index(VSsh_dim))
  do i=1,VSsh_dim
    do aa=1, HOsh_dim
      if (VSsh_list(i) .EQ. (HOsh_ant(aa))) then
        VStoHOsh_index(i) = aa
        endif
    enddo
    aa = VStoHOsh_index(i)
    print '(A,I3,2I7)',  '    ', i, VSsh_list(i), HOsh_ant(VStoHOsh_index(i))
  enddo

  NHO_vs  = 0
  NHO_co  = 99
  do aa = 1, VSsp_dim
    a   = VStoHOsp_index(aa)
    a_ant  = 2*HOsp_n(a) + HOsp_l(a)
    NHO_vs = max(a_ant, NHO_vs)     ! getting the maximum N shell of VS
    NHO_co = min(a_ant, NHO_co)     ! getting the minimum N shell of VS
  enddo
  NHO_co  = max(NHO_co - 1, 0)      ! The core is 1 shell below min_N of the VS

  !!! Exclude core for Full space M.E. exportings
  !!! NHO_co = 0 -> a 4He core, if the current nucleus is not N=Z=2, omit core.
  if (NHO_co .EQ. 0) then
    if ((dabs(valence_Z-2).GT.1.0d-6) .OR. (dabs(valence_N-2).GT.1.0d-6)) then
      NHO_co = -1
    else
      !! case where we have 4He but we don't want the core (0s_ state present)
      do aa = 1, VSsh_dim
        if (VSsh_list(aa) .EQ. 1) NHO_co = -1
      enddo
    endif
  endif

  print "(A,2I6)", '    ... N-shell HO for core(max)/vs(max):',NHO_co,NHO_vs
endif
if (EVAL_EXPLICIT_FIELDS_DD) then
  print '(A,3L10)', " [Explicit DD Field Eval.] Compute Full Valence Space =",&
    evalFullSPSpace
  print '(A)', " [DEPRECATED] Eval. DD. Explicitly = False"
  EVAL_EXPLICIT_FIELDS_DD = .FALSE.
  evalFullSPSpace = .FALSE.
endif
print *, ''
if ((FUNCTIONAL_DENS_MODE .EQ. 2) .AND. (has_HEIS_MAJO_TERMS)) then
  print *, "[ERROR] Functional (Phys.Rev.C 60 064312) AND Heis/Majo terms!"
  STOP
end if

hasX0M1 = abs(x0_DD_FACTOR - 1.0d+0) > 1.0d-6
CONST_EDD_M2 = CONST_EDD_M2_ETA / (CONST_EDD_M2_RHO0 ** alpha_DD)

print "(A)", " * Density dependent parameters imported."

print "(A,2L3)", " * [OPTIONs] Calculate DD-pn parts (HF/PA) :", &
                 CALCULATE_DD_PN_HF, CALCULATE_DD_PN_PA
print "(A,L3)",  " * [OPTIONs] Calculate DD HEISENBERG/MAJORANA :", &
                 has_HEIS_MAJO_TERMS

end subroutine import_DD_parameters


!-----------------------------------------------------------------------------!
! If there is a File called fixed_rearrangement.txt, import this matrix and   !
! set up the arguments necessary to add up as constant for the Gamma Field    !
!-----------------------------------------------------------------------------!
subroutine import_Rearrange_field_if_exist

integer   :: runit = 333, bogo_label
logical   :: is_exist
CHARACTER(LEN=25) :: file_input = "rearrangement_initial.txt"
INTEGER   :: io, aux_int
integer   :: i, j, icheck, HOsh_dim0
integer, dimension(:), allocatable :: HOsh_na0
real(r64) :: aux_real

inquire (file=file_input, exist=is_exist)
if ( is_exist ) then
  OPEN(runit, FILE=file_input, FORM="FORMATTED", STATUS="OLD", ACTION="READ")
  print "(A)", " * Initial rearrangement field PRESENT, reading from file."
else
  print "(A)", " * Initial rearrangement field NOT found."
  return
endif

read(runit,*) HOsh_dim0
allocate(HOsh_na0(HOsh_dim0))
do i = 1, HOsh_dim0
  read(runit,*) HOsh_na0(i)
enddo

!!! Stops the run if the model spaces of the wave func. and interaction differ
icheck = 0
if ( HOsh_dim0 /= HOsh_dim ) icheck = icheck + 1
do i = 1, min(HOsh_dim,HOsh_dim0)
  if ( HOsh_na0(i) /= HOsh_na(i) ) icheck = icheck + 1
enddo
if ( icheck /= 0 ) then
  print '(/,"The model space of the seed wave function is not consistent", &
        & " with the one of the interaction.")'
  print*, 'Inter:', HOsh_dim,  (HOsh_na(i), i=1,HOsh_dim)
  print*, 'State:', HOsh_dim0, (HOsh_na0(i), i=1,HOsh_dim0)
  stop
endif

allocate(fixed_rearrang_field(HOsp_dim,HOsp_dim))
fixed_rearrang_field = zzero
read(runit,*) bogo_label
do i = 1, HOsp_dim
  do j = 1, HOsp_dim
    read(runit,*) aux_real
    fixed_rearrang_field(j,i) = dcmplx(aux_real, 0.0d0)
!    if (dabs(aux_real).GT. 1.0d-07) then
!      print "(F21.15)", dreal(fixed_rearrang_field(j,i))
!    else
!      print "(D23.15)", dreal(fixed_rearrang_field(j,i))
!    endif
  enddo
enddo
USING_FIXED_REARRANGEMENT = .TRUE.
!! TODO: REARRANGE THE rearrange field if the order of the shells doesn't match

close (runit, status='keep')
deallocate(HOsh_na0)

end subroutine import_Rearrange_field_if_exist


!-----------------------------------------------------------------------------!
! Subroutine to import extra constants or modes for the calculation           !
!-----------------------------------------------------------------------------!
subroutine set_extra_DD_parameters(aux_int, aux_float)
integer,   intent(in) :: aux_int
real(r64), intent(in) :: aux_float



select case (aux_int)
  !! OPTIONS TO AVOID CALCULATING PN-DD FIELDS
  case (1)
    PRINT_GUTS = .FALSE.
  case (2)
    CALCULATE_DD_PN_PA = .FALSE. !logical(abs(aux_float) .GT. 1.0)
    print "(A,L3)", " > Calculate DD-pn PAIR Field OFF:", CALCULATE_DD_PN_PA
  case (3)
    CALCULATE_DD_PN_HF = .FALSE. !logical(abs(aux_float) .GT. 1.0)
    print "(A,L3)", " > Calculate DD-pn HF   Field OFF:", CALCULATE_DD_PN_HF
  case (4)
    EXPORT_GRAD_DD = .TRUE.
    print "(A,L3)", " > Exporting Gradient Hamiltonian:", EXPORT_GRAD_DD

  !! CUTOFF OPTIONS
  case (11)
    CUTOFF_ENERGY_MAX = aux_float
    if (CUTOFF_MODE .EQ. 2) then
      CUTOFF_MODE = 3
    else
      CUTOFF_MODE = 1
    endif
    EVAL_CUTOFF = .TRUE.
    print "(A,2I3)", " > Param CUTOFF-Energy: ", CUTOFF_MODE, CUTOFF_ENERGY_MAX
  case (12)
    CUTOFF_KAPPA  = aux_float
    if (CUTOFF_MODE .EQ. 2) then
      CUTOFF_MODE = 3
    else
      CUTOFF_MODE = 2
    end if
    EVAL_CUTOFF = .TRUE.
    print "(A,2I3)", " > Param CUTOFF-Kappa: ", CUTOFF_MODE, CUTOFF_KAPPA

  !! POTENTIALS  ---------------------------------------------------------------
  ! * Potential DD With no exchange term (x0=0) for N.matter
  !   E.Garrido, P.Sarriguren, E.Moya, N.Schuck - Phys.Rev.C 60, 064312 (1999)
  case(21)
    FUNCTIONAL_DENS_MODE = 2
    CONST_EDD_M2_ETA  = aux_float
    print "(A,F12.9)", " > Potential MODE 2: ETA  =", CONST_EDD_M2_ETA
  case(22)
    FUNCTIONAL_DENS_MODE = 2
    CONST_EDD_M2_RHO0 = aux_float
    print "(A,F12.9)", " > Potential MODE 2: RHO_0=", CONST_EDD_M2_RHO0


  case(31)
    has_HEIS_MAJO_TERMS  = .TRUE.
    CONST_x0_EXC_HEIS = aux_float
    print "(A,F12.9)", " > Exchange (spin):     Heisenberg=", CONST_x0_EXC_HEIS
  case(32)
    has_HEIS_MAJO_TERMS  = .TRUE.
    CONST_x0_EXC_MAJO = aux_float
    print "(A,F12.9)", " > Exchange (s-T):        Majorana=", CONST_x0_EXC_MAJO
  case(33:42)
    print "(3A)"," [ERROR] Several DD term implementation is not ",&
      "valid, Program to do it in repository ", &
      "[https://github.com/migueldelafuente1/densN_taurus_vap] STOP."
    STOP
  case default
    print "(2A,I4,F15.6)"," [ERROR] Invalid option SetUp Extra-argument case",&
      " not valid. Got:", aux_int, aux_float
    STOP
end select

if (has_sevaral_DD_tems)then
  parameters_alphx0x0Hx0M(1,1) = t3_DD_CONST
  parameters_alphx0x0Hx0M(1,2) = x0_DD_FACTOR
  parameters_alphx0x0Hx0M(1,3) = CONST_x0_EXC_HEIS
  parameters_alphx0x0Hx0M(1,4) = CONST_x0_EXC_MAJO
end if
end subroutine set_extra_DD_parameters


!-----------------------------------------------------------------------------!
! Subroutine to identify the index of the valence space from the basis        !
!-----------------------------------------------------------------------------!
subroutine set_valence_space_to_export
integer :: i, j, k, i_ant
integer, dimension(:), allocatable :: temp_list_index

VSsp_dim  = 0
allocate(temp_list_index(HOsp_dim))
temp_list_index = 0

do i=1, HOsp_dim
  i_ant = 10000*HOsp_n(i) + 100*HOsp_l(i) + HOsp_2j(i)
  do j=1, VSsh_dim
    if (i_ant.EQ.VSsh_list(j)) then
      VSsp_dim = VSsp_dim + 1
      temp_list_index(VSsp_dim) = i
    endif
  end do
enddo

VSsp_dim2 = VSsp_dim ** 2
if (evalQuasiParticleVSpace) then
  WBsp_dim  = HOsp_dim
  WBsh_dim  = HOsh_dim
else
  WBsp_dim  = VSsp_dim
  WBsh_dim  = VSsh_dim
end if
allocate(WBtoHOsp_index(WBsp_dim))
allocate(VStoHOsp_index(VSsp_dim), VSsp_VSsh(VSsp_dim))

do i=1, VSsp_dim
  VStoHOsp_index(i) = temp_list_index(i)
  if (.NOT.evalQuasiParticleVSpace) then !! WB for VS eval are the reduced states
    WBtoHOsp_index(i) = temp_list_index(i)
  endif

  k = VStoHOsp_index(i)
  i_ant = 10000*HOsp_n(k) + 100*HOsp_l(k) + HOsp_2j(k)
  do j = 1, VSsh_dim
    if (VSsh_list(j) .NE. i_ant) cycle
    VSsp_VSsh(i) = j
    EXIT
  end do
end do

!! in case the QP evaluation, Working basis for the hamiltonian reads all states
if (evalQuasiParticleVSpace) then
  do i = 1, HOsp_dim
    WBtoHOsp_index(i) = i
  enddo
endif

deallocate(temp_list_index)

end subroutine set_valence_space_to_export


!-----------------------------------------------------------------------------!
! subroutine set_density_dependent                                            !
!                                                                             !
! Set  the grid of radial and angular points for the grid to integrate.       !
! Export the grid and test if the sum correspond to the analytical.           !
!-----------------------------------------------------------------------------!
subroutine set_integration_grid

integer :: i_r, i_th, i_phi, i
real(r64) :: sum_, x, y

call GaussLaguerre(x_R, weight_R, r_dim, 0.5d+00)
do i_r = 1, r_dim
  r(i_r)        = HO_b * sqrt(x_R(i_r) / (2.0 + alpha_DD))
  r_export(i_r) = HO_b * sqrt(x_R(i_r))
enddo

Leb_dim   = nPoints(Omega_Order) ! := N_Omega
theta_dim = Leb_dim
phi_dim   = Leb_dim
angular_dim = Leb_dim
print '(A,i5,a,i5)','Using Gauss-Laguerre/LEBEDEV Integration method. A_dim',&
  angular_dim, "  xR_dim=", angular_dim*r_dim
print '(A,F15.13,A,F17.13)', ' b and hw in the program =', HO_b, "  ", HO_hw
allocate(x_Leb(3, Leb_dim))
allocate(weight_LEB(phi_dim))

allocate(theta(theta_dim))
allocate(theta_export(theta_dim))
allocate(cos_th(theta_dim))
allocate(cos_th_export(theta_dim))
allocate(phi(phi_dim))
allocate(phi_export(phi_dim))

call LLgrid(x_Leb, weight_LEB, Omega_Order)
sum_ = zero

if (PRINT_GUTS) then
  open(620, file='LebedevPointsWeights.gut')
  write(620,fmt='(A,I3)') "i, cos(th), phi, weight_LEB. OMEGA =", Omega_Order
endif
do i = 1, Leb_dim

  theta(i) = acos(x_Leb(3, i))
  x = x_Leb(1, i)
  y = x_Leb(2, i)
  if(y >= 0.d0) then
    phi(i) = atan2(y,x)
  else
    phi(i) = atan2(y,x) + 2.d0 * pi
  endif

  theta_export(i) = theta(i)
  cos_th(i)       = x_Leb(3, i)
  cos_th_export(i)= cos_th(i)
  phi_export(i)   = phi(i)

  if (PRINT_GUTS) then
    write(620,fmt='(I5,3D22.13)') i, cos_th(i), phi(i), weight_LEB(i)
  endif

  sum_ = sum_ + weight_LEB(i)
enddo

if (PRINT_GUTS) close(620)

end subroutine set_integration_grid

!-----------------------------------------------------------------------------!
! subroutine set_densty_dependent                                             !
!                                                                             !
! Implement the density dependent variables, the grid of integration, the     !
! radial and angular coefficients. Allocate the density dimens after the grid.!
!                                                                             !
!  Input:                                                                     !
!        seedtype = seed type with its interactions (1 cannot extract symm)   !
!        itermax  = maximum number of iterations                              !
!        Mphip = number of angles for the discretization of proton   PNR      !
!        Mphin =   "    "    "     "   "        "        "  neutrons  "       !
!-----------------------------------------------------------------------------!
subroutine set_densty_dependent(seedtype, itermax, proj_Mphip, proj_Mphin)
  integer, intent(in) :: seedtype, itermax, proj_Mphip, proj_Mphin
  complex(r64) :: x,y,z
  real(r64)    :: x1, x2,x3, y1, y2,y3, z1, z2, r, a, ALP, a2
  integer  :: i, n

  seed_type_sym = seedtype
  global_iter_max = itermax
  DOING_PROJECTION = (proj_Mphip > 1).OR.(proj_Mphin > 1)
  Mphip_DD = proj_Mphip
  Mphin_DD = proj_Mphin

  call import_DD_parameters
  call import_Rearrange_field_if_exist
  if (USING_FIXED_REARRANGEMENT) then
    hamil_H1 = hamil_H1 + fixed_rearrang_field
  endif

  if (.NOT.EVAL_DENSITY_DEPENDENT) then
    print "(A)", " * DD module is TURNED OFF, skip DD array setting."
    return
  endif

  call set_integration_grid
  call set_allocate_density_arrays

  call set_B_radial_coefficients
  call set_Radial2body_basis
  call set_SphericalHarmonic_basis
  call set_Y_KM_matrixElements
  call set_rearrangement_RadAng_fucntions

  print "(A)", " * Setting up DD module [DONE]"
  print "(A,L1)", " * DOING_PROJECTION (for DD) = ", DOING_PROJECTION

end subroutine set_densty_dependent

subroutine set_allocate_density_arrays

  allocate(rhoLR  (HOsp_dim,HOsp_dim))
  allocate(kappaLR(HOsp_dim,HOsp_dim))
  allocate(kappaRL(HOsp_dim,HOsp_dim))

  allocate(density (r_dim, angular_dim))
  allocate(dens_pnt(5, r_dim, angular_dim))
  allocate(dens_alpha(r_dim, angular_dim))
  allocate(dens_alpm1(r_dim, angular_dim))

  if (export_density) then
    allocate(density_export(r_dim, angular_dim))
    allocate(density_export_p(r_dim, angular_dim))
    allocate(density_export_n(r_dim, angular_dim))
    allocate(pairdens_export(r_dim, angular_dim))
    allocate(pairdens_export_n(r_dim, angular_dim))
    allocate(pairdens_export_p(r_dim, angular_dim))
  endif

  allocate(rearrangement_me(HOsp_dim, HOsp_dim))
  !allocate(rearrangement_me(HOsp_dim/2, HOsp_dim/2))
  allocate(rearrang_field(HOsp_dim, HOsp_dim))
  allocate(rea_common_RadAng(HOsp_dim /2, HOsp_dim /2, r_dim, angular_dim))
  allocate(REACommonFields(r_dim, angular_dim))

  !! CUTOFF ELEMENTS
  lambdaFer_DD = zero


end subroutine set_allocate_density_arrays

!-----------------------------------------------------------------------------!
! subroutine set_Radial2body_basis                                            !
!                                                                             !
! Compute and save the 2body radial wfunctions to import in the computation   !
! of the delta matrix elements and the density matrix elements.               !
! Previous r and B_coeff_radial implementation are required                   !
!-----------------------------------------------------------------------------!
subroutine set_Radial2body_basis

integer   :: a_sh, b_sh, i_r, na,la,nb,lb
real(r64) :: radial, x

allocate(radial_2b_sho_memo(HOsh_dim,HOsh_dim,r_dim))
if (export_density) then
  allocate(radial_2b_sho_export_memo(HOsh_dim,HOsh_dim,r_dim))
end if

if (PRINT_GUTS) open(629, file='R_radial2B_wf.gut')
do a_sh = 1, HOsh_dim
  na = HOsh_n(a_sh)
  la = HOsh_l(a_sh)
  do b_sh = a_sh, HOsh_dim
    nb = HOsh_n(b_sh)
    lb = HOsh_l(b_sh)

    if (PRINT_GUTS) write(629, fmt="(4I3)", advance='no') na, la, nb, lb
    do i_r = 1, r_dim
!     radial = two_sho_radial_functions_bench(a_sh, b_sh, r(i_r))
      radial = two_sho_radial_functions(a_sh, b_sh, r(i_r), .TRUE.)

      radial_2b_sho_memo(a_sh, b_sh, i_r) = radial
      radial_2b_sho_memo(b_sh, a_sh, i_r) = radial ! a_sh=b_sh overwrites

      ! Radial grid for Laguerre
      !r _lag = b *(x_R / 2+alpha)**0.5

      !! assert test R_ab = R_ba
      if (dabs(two_sho_radial_functions(a_sh, b_sh, r(i_r), .FALSE.) - &
        two_sho_radial_functions(b_sh, a_sh, r(i_r), .FALSE.)) > 1.d-12) then
        print "(A,3I4,A,2D20.13)","[ASSERT ERROR] R(ab)/=R(ba) for a,b,i_r=",&
           a_sh,b_sh,i_r," Rab/Rba=", &
           two_sho_radial_functions(a_sh, b_sh, r(i_r), .FALSE.),&
           two_sho_radial_functions(b_sh, a_sh, r(i_r), .FALSE.)

      endif
      if (PRINT_GUTS) then
        write(629,fmt='(3D22.13)',advance='no') r(i_r),radial,weight_R(i_r)
      endif

      if (export_density) then
        radial = two_sho_radial_functions(a_sh,b_sh, r_export(i_r), .FALSE.)

        radial_2b_sho_export_memo(a_sh, b_sh, i_r) = radial
        radial_2b_sho_export_memo(b_sh, a_sh, i_r) = radial
      end if
    enddo
    if (PRINT_GUTS) write(629, *) ""
  enddo
enddo

if (PRINT_GUTS) close(629)

end subroutine set_Radial2body_basis


!-----------------------------------------------------------------------------!
! function angular_momentum_index                                             !
!                                                                             !
! index accessor for angular momentum states, tuple (K, M) can be stored in   !
! order from 1 to the length of the K maximum of the array.                   !
!                                                                             !
! * For Integer : K = 0,1,2, ...  ! M=(0), (-1,0,+1), ...                     !
!       index =       (K)**2        +      K + M          + 1                 !
!               (previous elements)   (range 0 to 2K + 1)   (starting index)  !
!                                                                             !
! * For Half-Integers : K = 1,3,5, ...;   M=(-1,+1), (-3,-1,+1,+3), ...       !
!       I(K)  = (K+1)/2                                                       !
!       index = [I(K-1)**2 + I(K-1) = I*(I-1)] +   (K + M)/2           + 1    !
!                    (previous elements)        (range 0 to 2K + 1) (start)   !
!                                                                             !
! * This function can be used to define the length of the K,M array, by giving!
!   the last valid element: angular_momentum_index(K_max, K_max, half_integer)!
!-----------------------------------------------------------------------------!
function angular_momentum_index(K, M, half_integer) result (index_)
    integer, intent(in) :: K, M
    logical, intent(in) :: half_integer
    integer  :: index_

    if (half_integer) then
        index_ = ((K*K - 1)/4) + ((K + M)/2) + 1
    else
        index_ = K*(K + 1) + M + 1
    end if
    return

end function angular_momentum_index


!------------------------------------------------------------------------------!
! angular indexing for the array theta and Phi in the angular functions.       !
!                                                                              !
! i = (dim_Phi - 1)*i_theta + i_phi + 1                                        !
!------------------------------------------------------------------------------!
subroutine angular_index(index_, i_theta, i_phi)
    integer, intent(in) :: index_
    integer :: i_theta, i_phi

    if (.TRUE.) then  ! integration_method == 3
      i_theta = index_
      i_phi   = index_
    else
      !! Explicit Euclides to get the final value (Efficient?)

      !! Fortran
      i_theta = ((index_ - 1) / phi_dim) + 1
      i_phi   = MOD(index_ - 1, phi_dim) + 1
    endif

end subroutine angular_index
!-----------------------------------------------------------------------------!
! subroutine set_SphericalHarmonic_basis                                      !
!                                                                             !
! Compute and save the spherical harmonics on the theta-phi grid to avoid di- !
! rect computation. elements and the density matrix elements.                 !
!-----------------------------------------------------------------------------!
subroutine set_SphericalHarmonic_basis

integer :: i, K,M, K_max, i_th, i_phi, index_am, Kdim, la,mla, lb,mlb, index_bm
complex(r64) :: sh_val, sh_val_a, sh_val_b


K_max = 2 * MAXVAL(HOsh_l)
K_max = MAX(K_max, MAXVAL(HOsh_2j))
Kdim = angular_momentum_index(K_max, K_max, .FALSE.)

allocate(sph_harmonics_memo(Kdim, angular_dim))

if (PRINT_GUTS) open(621, file='test_sphharm_memo.gut')

do K = 0, K_max
   do M = -K, K
     index_am = angular_momentum_index(K,M,.FALSE.)

     do i = 1, angular_dim
       call angular_index(i, i_th, i_phi)
       sh_val = spherharmonic(K, M, theta(i_th), phi(i_phi))
       sph_harmonics_memo(index_am, i) = sh_val

       if (PRINT_GUTS) then
write(621, fmt='(5I4,2F20.15)') K,M, i,i_th,i_phi, dreal(sh_val), dimag(sh_val)
       endif
     enddo
   enddo
enddo

if (PRINT_GUTS) close(621)

!! Part for the non-combined spherical harmonics, allocates Y*(a) Y(b)
K_max = HO_lmax
Kdim = angular_momentum_index(K_max, K_max, .FALSE.)

allocate(sphharmDUAL_memo(Kdim, Kdim, angular_dim))

if (PRINT_GUTS) open(622, file='test_sphharmDUAL_memo.gut') ! Verified with python scipy

do la = 0, HO_lmax
  do mla = -la, la
    index_am = angular_momentum_index(la, mla, .FALSE.)

    do lb = 0, HO_lmax
      do mlb = -lb, lb
        index_bm = angular_momentum_index(lb, mlb, .FALSE.)

        do i = 1, angular_dim
          call angular_index(i, i_th, i_phi)
          sh_val_a = spherharmonic(la,mla, theta(i_th), phi(i_phi))
          sh_val_b = spherharmonic(lb,mlb, theta(i_th), phi(i_phi))

          sh_val = CONJG(sh_val_a) * sh_val_b
          sphharmDUAL_memo(index_am, index_bm, i) = sh_val

          if (PRINT_GUTS) then
            write(622, '(4I5,A,I4,2F20.15,3F20.15)') la, mla, lb, mlb, &
              " i_ang=",i,theta(i_th), &
              phi(i_phi), dREAL(sh_val), dIMAG(sh_val), weight_LEB(i)
          endif
        enddo

      enddo
    enddo
  enddo
enddo

if (PRINT_GUTS) close(622)
!! test export

end subroutine set_SphericalHarmonic_basis

!-----------------------------------------------------------------------------!
! subroutine set_Y_KM_matrixElements                                          !
!                                                                             !
! Predefine the matrix elements to be used in the 1 body density matrix,      !
! derived from spherical harmonic expansion and related to the reduced m.e.   !
! <l1,1/2,j1|| Y_K ||l2,1/2,j2>. The calculation                              !
!-----------------------------------------------------------------------------!
subroutine set_Y_KM_matrixElements
integer   :: i,j,ja,la,ma,jb,lb,mb,K,M, indx_a,indx_b,indx_km, j_max, Kdim,Jdim
integer   :: spO2  , i_an, i_th, i_phi, ms
real(r64) :: aux_val, phase, cgc_a, cgc_b

spO2 = HOsp_dim / 2

j_max = MAXVAL(HOsh_2j)
Jdim = angular_momentum_index(j_max, j_max, .TRUE.)
Kdim = angular_momentum_index(j_max, j_max, .FALSE.)
!!! Note: K_max = 2 * (j_max + j_max)[half int]

allocate(dens_Y_KM_me(Jdim, Jdim, Kdim))

Kdim = angular_momentum_index(HO_lmax, HO_lmax, .FALSE.)
allocate(AngFunctDUAL_HF(4, HOsp_dim/2, HOsp_dim/2, angular_dim))
allocate(AngFunctDUAL_P1(4, HOsp_dim/2, HOsp_dim/2, angular_dim))
allocate(AngFunctDUAL_P2(4, HOsp_dim/2, HOsp_dim/2, angular_dim))

allocate(BulkHF(5,4, r_dim, angular_dim))
allocate(BulkP1(5,4, r_dim, angular_dim))
allocate(BulkP2(5,4, r_dim, angular_dim))

AngFunctDUAL_HF = zzero
AngFunctDUAL_P1 = zzero
AngFunctDUAL_P2 = zzero

if (has_HEIS_MAJO_TERMS) then
  allocate(BulkP1_HM(5,4, r_dim, angular_dim))
  BulkP1_HM = zzero
endif

if (PRINT_GUTS) then
  open(622, file='Y_ab_K_matrix_elements.gut')
  write(622,fmt='(2A)') " la ja ma i_a|  lb jb mb i_b|   K   M i_km ",&
    "=  Dens_YKM_ab_me"
  write(622,fmt='(2A)') '-------------|--------------|---------------------',&
    '------------------------'
endif

do i = 1, spO2
  ja = HOsp_2j(i)
  la = HOsp_l(i)
  ma = HOsp_2mj(i)
  indx_a = angular_momentum_index(ja, ma, .TRUE.)

  do j = 1, spO2
    jb = HOsp_2j(j)
    lb = HOsp_l(j)
    mb = HOsp_2mj(j)
    indx_b = angular_momentum_index(jb, mb, .TRUE.)
    do K = max(0, abs(ja - jb) / 2), (ja + jb) / 2
!            M = (ma - mb)/2 !! (Fenshb)
        M = (mb - ma) /2 !! (Suhonen_Vasr)
        if (abs(M) > K) cycle
        if (MOD(la + lb + K, 2) == 1) cycle

        indx_km = angular_momentum_index(K, M, .FALSE.)

        !! Fenshbach expansion of ccgg [TESTED PASS, but conjugated]
!            call ClebschGordan(2*K, ja, jb, -2*M, ma, mb, cgc_a)
!            call ClebschGordan(ja, jb, 2*K, 1, -1, 0, cgc_b)
!            aux_val = cgc_a * cgc_b * sqrt((ja + 1)/(4*pi))
!            phase   = ((-1)**((1 + ma + mb + ja)/2))
!
!            !! Talmi formulas, (use ma - mb, there is some kind of skip)
!            call ClebschGordan(2*K, ja, jb, -2*M, ma, mb, cgc_a)
!            call ClebschGordan(ja, jb, 2*K, 1, -1, 0, cgc_b)
!            phase   = (-1)**((ja + mb) + ((jb + 3)/2)) !

        !! Varsalovich - Suhonen formulas [TEST 2 CORRECT]
        call ClebschGordan(2*K, ja, jb, 2*M, ma, mb, cgc_a)
        call ClebschGordan(ja, jb, 2*K, 1, -1, 0, cgc_b)
        aux_val = cgc_a * cgc_b * sqrt((ja + 1)/(4*pi))
        phase   = ((-1)**((ja - 1)/2 - jb + 1))

        dens_Y_KM_me(indx_a, indx_b, indx_km) = aux_val * phase

        if (PRINT_GUTS) then
          write(622,fmt='(4I3,A,4I3,A,3I3,A,F20.15,F9.3)')&
            la,ja,ma,indx_a,' | ',lb,jb,mb,indx_b,' | ',K,M,indx_km,&
            ' = ', dens_Y_KM_me(indx_a, indx_b, indx_km), phase
        endif

    enddo
    ! Sect. define uncoupled angular coefficients to precalculate fields. ==
    call set_sphhDual_precalcFields(i, ja, la, ma,  j, jb, lb, mb)

  enddo
enddo

if (PRINT_GUTS) close(622)

!TEST: export the Xi functions for a random angle
if (PRINT_GUTS) then
  open(543, file='Xi0Functions_AngFixed.gut')
  open(544, file='Xi1Functions_AngFixed.gut')
  open(545, file='Xi2Functions_AngFixed.gut')
  i_an = 53
  write(543, fmt="(A,A)")"I_ANG  ,    a  la  ja mja  ,     b   lb   jb  mjb",&
    "       msms'(++)       msms'(+-)       msms'(-+)       msms'(--) "
  write(544, fmt="(A,A)")"I_ANG  ,    a  la  ja mja  ,     b   lb   jb  mjb",&
    "       msms'(++)       msms'(+-)       msms'(-+)       msms'(--) "
  write(545, fmt="(A,A)")"I_ANG  ,    a  la  ja mja  ,     b   lb   jb  mjb",&
    "       msms'(++)       msms'(+-)       msms'(-+)       msms'(--) "


!do i_an = 1, angular_dim
do i = 1, spO2
  ja = HOsp_2j(i)
  la = HOsp_l(i)
  ma = HOsp_2mj(i)
  do j = 1, spO2
    jb = HOsp_2j(j)
    lb = HOsp_l(j)
    mb = HOsp_2mj(j)
    write(543,fmt='(I5,A,4I4,A,4I5,A)',advance='no')i_an,"  , ",i,la,ja,ma, &
      "  , ", j,lb,jb,mb, "  , "
    write(544,fmt='(I5,A,4I4,A,4I5,A)',advance='no')i_an,"  , ",i,la,ja,ma, &
      "  , ", j,lb,jb,mb, "  , "
    write(545,fmt='(I5,A,4I4,A,4I5,A)',advance='no')i_an,"  , ",i,la,ja,ma, &
      "  , ", j,lb,jb,mb, "  , "
    do ms = 1, 4
      write(543,fmt='((F17.12,SP,F17.12,"j"))',advance='no') &
        dreal(AngFunctDUAL_HF(ms,i,j,i_an)), dimag(AngFunctDUAL_HF(ms,i,j,i_an))
      write(544,fmt='((F17.12,SP,F17.12,"j"))',advance='no') &
        dreal(AngFunctDUAL_P1(ms,i,j,i_an)), dimag(AngFunctDUAL_P1(ms,i,j,i_an))
      write(545,fmt='((F17.12,SP,F17.12,"j"))',advance='no') &
        dreal(AngFunctDUAL_P2(ms,i,j,i_an)), dimag(AngFunctDUAL_P2(ms,i,j,i_an))
    end do
    write(543, fmt='(A)') ''
    write(544, fmt='(A)') ''
    write(545, fmt='(A)') ''
  end do
end do
  !write(543,fmt='(I3,F19.15,A)') i_an, weight_LEB(i_an), "  ------======------"
  !enddo
  close(543)
  close(544)
  close(545)
endif

end subroutine set_Y_KM_matrixElements

!-----------------------------------------------------------------------------!
! subroutine set_sphhDual_precalcFields                                       !
!                                                                             !
! Predefine the angular functions for the preevaluated fields, setting the    !
! sphhDUAL_** functions divided in 4 parts depending on the spin components   !
! 1:(+1/2 +1/2) 2:(+1/2 -1/2) 3:(-1/2 +1/2) 4:(-1/2 -1/2)                     !
!    the ml component is not necessary since sp coeff relate to l value and   !
!    ml = mj - ms.                                                            !
!-----------------------------------------------------------------------------!
subroutine set_sphhDual_precalcFields(a,ja,la,ma, b,jb,lb,mb)

integer, intent(in) :: a,ja,la,ma, b,jb,lb,mb
integer :: i_ang,mla,mlb, ind0_am, ind0_bm, ind1_am, ind1_bm, ind2_am, ind2_bm
real(r64) :: cgc_a_msp1, cgc_a_msm1, cgc_b_msp1, cgc_b_msm1
complex(r64) :: YaYb, YCaYb, YCaYCb
integer :: ms, i, j

!!! ASSERT ! a , b must be in the pp quadrant, the rest'll be copied afterwards
if ((a > HOsp_dim/2).OR.(b > HOsp_dim/2)) then
   print *, " [ASS. ERROR] compute_bulkDens4Fields a,b (sp) in neutron space"
   STOP
endif

call ClebschGordan(2*la, 1, ja, (ma - 1),  1, ma, cgc_a_msp1)
call ClebschGordan(2*lb, 1, jb, (mb - 1),  1, mb, cgc_b_msp1)
call ClebschGordan(2*la, 1, ja, (ma + 1), -1, ma, cgc_a_msm1)
call ClebschGordan(2*lb, 1, jb, (mb + 1), -1, mb, cgc_b_msm1)

!!  ms1_ , ms2_ = +1/2  +1/2 --------------------------------------------
mla = (ma - 1) / 2
mlb = (mb - 1) / 2
if ((abs(mla) <= la).AND.(abs(mlb) <= lb)) then
    ! Xo (Y*a, Yb)
    ind0_am = angular_momentum_index(la,  mla, .FALSE.)
    ind0_bm = angular_momentum_index(lb,  mlb, .FALSE.)
    ! X1 (Ya , Yb) = (-)mla Xo(-mla, mlb)
    ind1_am = angular_momentum_index(la, -mla, .FALSE.)
    ind1_bm = angular_momentum_index(lb,  mlb, .FALSE.)
    ! X2 (Y*a,Y*b) = (-)mlb Xo(mla, -mlb)
    ind2_am = angular_momentum_index(la,  mla, .FALSE.)
    ind2_bm = angular_momentum_index(lb, -mlb, .FALSE.)

    do i_ang = 1, angular_dim
      YCaYb    = sphharmDUAL_memo(ind0_am, ind0_bm, i_ang)
      YaYb     = ((-1)**mla)*sphharmDUAL_memo(ind1_am, ind1_bm, i_ang)
      YCaYCb   = ((-1)**mlb)*sphharmDUAL_memo(ind2_am, ind2_bm, i_ang)

      AngFunctDUAL_HF(1, a, b, i_ang) = cgc_a_msp1 * cgc_b_msp1 * YCaYb
      AngFunctDUAL_P1(1, a, b, i_ang) = cgc_a_msp1 * cgc_b_msp1 * YaYb
      AngFunctDUAL_P2(1, a, b, i_ang) = cgc_a_msp1 * cgc_b_msp1 * YCaYCb
    enddo
endif

!!  ms1_ , ms2_ = +1/2  -1/2 --------------------------------------------
mla = (ma - 1) / 2
mlb = (mb + 1) / 2
if ((abs(mla) <= la).AND.(abs(mlb) <= lb)) then
    ! Xo (Y*a, Yb)
    ind0_am = angular_momentum_index(la,  mla, .FALSE.)
    ind0_bm = angular_momentum_index(lb,  mlb, .FALSE.)
    ind1_am = angular_momentum_index(la, -mla, .FALSE.)
    ind1_bm = angular_momentum_index(lb,  mlb, .FALSE.)
    ind2_am = angular_momentum_index(la,  mla, .FALSE.)
    ind2_bm = angular_momentum_index(lb, -mlb, .FALSE.)

    do i_ang = 1, angular_dim
      YCaYb    = sphharmDUAL_memo(ind0_am, ind0_bm, i_ang)
      YaYb     = ((-1)**mla)*sphharmDUAL_memo(ind1_am, ind1_bm, i_ang)
      YCaYCb   = ((-1)**mlb)*sphharmDUAL_memo(ind2_am, ind2_bm, i_ang)

      AngFunctDUAL_HF(2, a, b, i_ang) = cgc_a_msp1 * cgc_b_msm1 * YCaYb
      AngFunctDUAL_P1(2, a, b, i_ang) = cgc_a_msp1 * cgc_b_msm1 * YaYb
      AngFunctDUAL_P2(2, a, b, i_ang) = cgc_a_msp1 * cgc_b_msm1 * YCaYCb
    enddo
endif

!!  ms1_ , ms2_ = -1/2  +1/2 --------------------------------------------
mla = (ma + 1) / 2
mlb = (mb - 1) / 2
if ((abs(mla) <= la).AND.(abs(mlb) <= lb)) then
    ! Xo (Y*a, Yb)
    ind0_am = angular_momentum_index(la,  mla, .FALSE.)
    ind0_bm = angular_momentum_index(lb,  mlb, .FALSE.)
    ind1_am = angular_momentum_index(la, -mla, .FALSE.)
    ind1_bm = angular_momentum_index(lb,  mlb, .FALSE.)
    ind2_am = angular_momentum_index(la,  mla, .FALSE.)
    ind2_bm = angular_momentum_index(lb, -mlb, .FALSE.)

    do i_ang = 1, angular_dim
      YCaYb    = sphharmDUAL_memo(ind0_am, ind0_bm, i_ang)
      YaYb     = ((-1)**mla)*sphharmDUAL_memo(ind1_am, ind1_bm, i_ang)
      YCaYCb   = ((-1)**mlb)*sphharmDUAL_memo(ind2_am, ind2_bm, i_ang)

      AngFunctDUAL_HF(3, a, b, i_ang) = cgc_a_msm1 * cgc_b_msp1 * YCaYb
      AngFunctDUAL_P1(3, a, b, i_ang) = cgc_a_msm1 * cgc_b_msp1 * YaYb
      AngFunctDUAL_P2(3, a, b, i_ang) = cgc_a_msm1 * cgc_b_msp1 * YCaYCb
    enddo
endif

!!  ms1_ , ms2_ = -1/2  -1/2 --------------------------------------------
mla = (ma + 1) / 2
mlb = (mb + 1) / 2
if ((abs(mla) <= la).AND.(abs(mlb) <= lb)) then
    ! Xo (Y*a, Yb)
    ind0_am = angular_momentum_index(la,  mla, .FALSE.)
    ind0_bm = angular_momentum_index(lb,  mlb, .FALSE.)
    ind1_am = angular_momentum_index(la, -mla, .FALSE.)
    ind1_bm = angular_momentum_index(lb,  mlb, .FALSE.)
    ind2_am = angular_momentum_index(la,  mla, .FALSE.)
    ind2_bm = angular_momentum_index(lb, -mlb, .FALSE.)

    do i_ang = 1, angular_dim
      YCaYb    = sphharmDUAL_memo(ind0_am, ind0_bm, i_ang)
      YaYb     = ((-1)**mla)*sphharmDUAL_memo(ind1_am, ind1_bm, i_ang)
      YCaYCb   = ((-1)**mlb)*sphharmDUAL_memo(ind2_am, ind2_bm, i_ang)

      AngFunctDUAL_HF(4, a, b, i_ang) = cgc_a_msm1 * cgc_b_msm1 * YCaYb
      AngFunctDUAL_P1(4, a, b, i_ang) = cgc_a_msm1 * cgc_b_msm1 * YaYb
      AngFunctDUAL_P2(4, a, b, i_ang) = cgc_a_msm1 * cgc_b_msm1 * YCaYCb
    enddo
endif

end subroutine set_sphhDual_precalcFields


!-----------------------------------------------------------------------------!
! subroutine compute_bulkDens4Fields                                          !
! a, b <int> sp states from the pp space (neutron states got here splicitly)  !
! This requires the index a:1 -> N/2 and b:1 -> N/2                           !
!-----------------------------------------------------------------------------!
!  to evaluate the the value and transposed values for all spaces             !
!-----------------------------------------------------------------------------!
subroutine compute_bulkDens4_hf(a, b, a_sh, b_sh, i_r, i_a, overlap)

integer, intent(in)      :: a, b, a_sh, b_sh, i_r, i_a
complex(r64), intent(in) :: overlap

integer      :: spO2, par_ind,    ms,ms2, a_n, b_n
complex(r64) :: roP, roN, rPN, rNP, roPt, roNt, rPNt, rNPt, &
                A_part, aux, sum_, radial_ab
logical      :: aNeQb
spO2 = HOsp_dim / 2

! assertion, a, b only from pp space:
if ((a > spO2).OR.(b > spO2)) then
  print *, "[ASSERT. ERROR] compute_bulkDens4Fields a,b (sp) in neutron space"
  STOP
endif
if (b < a) then
  print *, "[ASSERT. ERROR] compute_bulkDens4Fields a <= b!, but b < a"
  STOP
end if
aNeQb = kdelta(a, b).ne.1
a_n   = a + spO2
b_n   = b + spO2

radial_ab = radial_2b_sho_memo(a_sh, b_sh, i_r) / overlap

roP  = radial_ab * rhoLR  (b  ,a)
roN  = radial_ab * rhoLR  (b_n,a_n)
! transposed, always add
roPt = radial_ab * rhoLR  (a,  b)
roNt = radial_ab * rhoLR  (a_n,b_n)

if (CALCULATE_DD_PN_HF) then
  rPN  = radial_ab * rhoLR  (b  ,a_n)
  rNP  = radial_ab * rhoLR  (b_n,a)
  ! transposed, always add
  rPNt = radial_ab * rhoLR  (a  ,b_n)
  rNPt = radial_ab * rhoLR  (a_n,b)
else
  rPN  = zzero
  rNP  = zzero
  rPNt = zzero
  rNPt = zzero
end if


!! compute the direct bulk densities
sum_ = AngFunctDUAL_HF(1,a,b,i_a) + AngFunctDUAL_HF(4,a,b,i_a)
dens_pnt(1,i_r,i_a) = dens_pnt(1,i_r,i_a) + (sum_ * roP)
dens_pnt(2,i_r,i_a) = dens_pnt(2,i_r,i_a) + (sum_ * roN)
dens_pnt(3,i_r,i_a) = dens_pnt(3,i_r,i_a) + (sum_ * rPN)
dens_pnt(4,i_r,i_a) = dens_pnt(4,i_r,i_a) + (sum_ * rNP)
!! compute the direct bulk densities
if (aNeQb) then
  sum_ = AngFunctDUAL_HF(1,b,a,i_a) + AngFunctDUAL_HF(4,b,a,i_a)
  dens_pnt(1,i_r,i_a) = dens_pnt(1,i_r,i_a) + (sum_ * roPt)
  dens_pnt(2,i_r,i_a) = dens_pnt(2,i_r,i_a) + (sum_ * roNt)
  dens_pnt(3,i_r,i_a) = dens_pnt(3,i_r,i_a) + (sum_ * rPNt)
  dens_pnt(4,i_r,i_a) = dens_pnt(4,i_r,i_a) + (sum_ * rNPt)
endif

do ms = 1, 4
  select case (ms)
    case (2, 3)
      ms2 = 5 - ms
    case default
      ms2 = ms
  end select

  A_part = AngFunctDUAL_HF(ms2,a,b,i_a)
  BulkHF(1,ms,i_r,i_a) = BulkHF(1,ms,i_r,i_a) + (A_part * roP) !pp
  BulkHF(2,ms,i_r,i_a) = BulkHF(2,ms,i_r,i_a) + (A_part * roN) !nn
  BulkHF(3,ms,i_r,i_a) = BulkHF(3,ms,i_r,i_a) + (A_part * rPN) !pn
  BulkHF(4,ms,i_r,i_a) = BulkHF(4,ms,i_r,i_a) + (A_part * rNP) !np

  if (aNeQb) then
    A_part = AngFunctDUAL_HF(ms2,b,a,i_a)
    BulkHF(1,ms,i_r,i_a) = BulkHF(1,ms,i_r,i_a) + (A_part * roPt) !pp
    BulkHF(2,ms,i_r,i_a) = BulkHF(2,ms,i_r,i_a) + (A_part * roNt) !nn
    BulkHF(3,ms,i_r,i_a) = BulkHF(3,ms,i_r,i_a) + (A_part * rPNt) !pn
    BulkHF(4,ms,i_r,i_a) = BulkHF(4,ms,i_r,i_a) + (A_part * rNPt) !np
  endif
enddo

end subroutine compute_bulkDens4_hf

subroutine compute_bulkDens4_pair(a, b, a_sh, b_sh, i_r, i_a, overlap)

integer, intent(in)      :: a, b, a_sh, b_sh, i_r, i_a
complex(r64), intent(in) :: overlap

integer      :: spO2, par_ind,    ms,ms2, a_n, b_n
complex(r64) :: kaP,kaN,kaCcP, kaCcN,kPN,kNP,kCcNP,kCcPN, &
                kaPt,kaNt,kaCcPt, kaCcNt,kPNt,kNPt, &
                kCcNPt,kCcPNt,B1_part,B2_part, aux, sum_, radial_ab
logical      :: aNeQb
spO2 = HOsp_dim / 2

! assertion, a, b only from pp space:
if ((a > spO2).OR.(b > spO2)) then
  print *, "[ASSERT. ERROR] compute_bulkDens4Fields a,b (sp) in neutron space"
  STOP
endif
if (b < a) then
  print *, "[ASSERT. ERROR] compute_bulkDens4Fields a <= b!, but b < a"
  STOP
end if
aNeQb = kdelta(a, b).ne.1
a_n   = a + spO2
b_n   = b + spO2

radial_ab = radial_2b_sho_memo(a_sh, b_sh, i_r) / overlap

kaP    = radial_ab * kappaLR(a  ,b)
kaN    = radial_ab * kappaLR(a_n,b_n)
kaCcP  = radial_ab * kappaRL(a  ,b)
kaCcN  = radial_ab * kappaRL(a_n,b_n)
! transposed, always add
kaPt   = radial_ab * kappaLR(b  ,a)
kaNt   = radial_ab * kappaLR(b_n,a_n)
kaCcPt = radial_ab * kappaRL(b  ,a)
kaCcNt = radial_ab * kappaRL(b_n,a_n)


if (CALCULATE_DD_PN_PA) then
  kPN    = radial_ab * kappaLR(a  ,b_n)
  kNP    = radial_ab * kappaLR(a_n,b)
  kCcPN  = radial_ab * kappaRL(a  ,b_n)
  kCcNP  = radial_ab * kappaRL(a_n,b)
  ! transposed, always add
  kPNt   = radial_ab * kappaLR(b  ,a_n)
  kNPt   = radial_ab * kappaLR(b_n,a)
  kCcPNt = radial_ab * kappaRL(b  ,a_n)
  kCcNPt = radial_ab * kappaRL(b_n,a)
else
  kPN    = zzero
  kNP    = zzero
  kCcPN  = zzero
  kCcNP  = zzero
  kPNt   = zzero
  kNPt   = zzero
  kCcPNt = zzero
  kCcNPt = zzero
end if

do ms = 1, 4
  select case (ms)
    case (2, 3)
      ms2 = 5 - ms
    case default
      ms2 = ms
  end select

  B1_part = (AngFunctDUAL_P1(ms,a,b,i_a) - AngFunctDUAL_P1(ms2,a,b,i_a))
  BulkP1(1,ms,i_r,i_a) = BulkP1(1,ms,i_r,i_a) + (B1_part * kaP) !pp
  BulkP1(2,ms,i_r,i_a) = BulkP1(2,ms,i_r,i_a) + (B1_part * kaN) !nn

  B1_part = AngFunctDUAL_P1(ms ,a,b,i_a) &
             + (x0_DD_FACTOR*AngFunctDUAL_P1(ms2,a,b,i_a))
  B2_part = AngFunctDUAL_P1(ms2,a,b,i_a) &
             + (x0_DD_FACTOR*AngFunctDUAL_P1(ms ,a,b,i_a)) ! reuse of the aux variable
  BulkP1(3,ms,i_r,i_a) = BulkP1(3,ms,i_r,i_a) + (B1_part*kPN - B2_part*kNP) !pn
  BulkP1(4,ms,i_r,i_a) = BulkP1(4,ms,i_r,i_a) + (B1_part*kNP - B2_part*kPN) !np
  !! BulkP1_** are  common for the paring and rearrangement fields respectively.

  B2_part = AngFunctDUAL_P2(ms,a,b,i_a)
  BulkP2(1,ms,i_r,i_a) = BulkP2(1,ms,i_r,i_a) + (B2_part * kaCcP) !pp
  BulkP2(2,ms,i_r,i_a) = BulkP2(2,ms,i_r,i_a) + (B2_part * kaCcN) !nn
  BulkP2(3,ms,i_r,i_a) = BulkP2(3,ms,i_r,i_a) + (B2_part * kCcPN) !pn
  BulkP2(4,ms,i_r,i_a) = BulkP2(4,ms,i_r,i_a) + (B2_part * kCcNP) !np

  if (aNeQb) then
    B1_part = (AngFunctDUAL_P1(ms,b,a,i_a) - AngFunctDUAL_P1(ms2,b,a,i_a))
    BulkP1(1,ms,i_r,i_a) = BulkP1(1,ms,i_r,i_a) + (B1_part * kaPt) !pp
    BulkP1(2,ms,i_r,i_a) = BulkP1(2,ms,i_r,i_a) + (B1_part * kaNt) !nn

    B1_part = AngFunctDUAL_P1(ms ,b,a,i_a) &
               + (x0_DD_FACTOR*AngFunctDUAL_P1(ms2,b,a,i_a))
    B2_part = AngFunctDUAL_P1(ms2,b,a,i_a) &
               + (x0_DD_FACTOR*AngFunctDUAL_P1(ms ,b,a,i_a)) ! reuse of the aux variable
    BulkP1(3,ms,i_r,i_a) = BulkP1(3,ms,i_r,i_a) + (B1_part*kPNt - B2_part*kNPt) !pn
    BulkP1(4,ms,i_r,i_a) = BulkP1(4,ms,i_r,i_a) + (B1_part*kNPt - B2_part*kPNt) !np
    !! BulkP1_** are  common for the paring and rearrangement fields respectively.

    B2_part = AngFunctDUAL_P2(ms,b,a,i_a)
    BulkP2(1,ms,i_r,i_a) = BulkP2(1,ms,i_r,i_a) + (B2_part * kaCcPt) !pp
    BulkP2(2,ms,i_r,i_a) = BulkP2(2,ms,i_r,i_a) + (B2_part * kaCcNt) !nn
    BulkP2(3,ms,i_r,i_a) = BulkP2(3,ms,i_r,i_a) + (B2_part * kCcPNt) !pn
    BulkP2(4,ms,i_r,i_a) = BulkP2(4,ms,i_r,i_a) + (B2_part * kCcNPt) !np
  endif

  !! EXTENSION FOR HEISENBERG - MAJORANA PAIRING FIELDS
  if (has_HEIS_MAJO_TERMS) then

  B1_part = (AngFunctDUAL_P1(ms,a,b,i_a) - AngFunctDUAL_P1(ms2,a,b,i_a))
  BulkP1_HM(1,ms,i_r,i_a) = BulkP1_HM(1,ms,i_r,i_a) + (B1_part * kaP) !pp
  BulkP1_HM(2,ms,i_r,i_a) = BulkP1_HM(2,ms,i_r,i_a) + (B1_part * kaN) !nn

  B1_part =    CONST_x0_EXC_HEIS * AngFunctDUAL_P1(ms ,a,b,i_a) &
            - (CONST_x0_EXC_MAJO * AngFunctDUAL_P1(ms2,a,b,i_a))
  B2_part =    CONST_x0_EXC_HEIS * AngFunctDUAL_P1(ms2,a,b,i_a) &
            - (CONST_x0_EXC_MAJO * AngFunctDUAL_P1(ms ,a,b,i_a))
  BulkP1_HM(3,ms,i_r,i_a)= BulkP1_HM(3,ms,i_r,i_a) + (B1_part*kNP - B2_part*kPN) !pn
  BulkP1_HM(4,ms,i_r,i_a)= BulkP1_HM(4,ms,i_r,i_a) + (B1_part*kPN - B2_part*kNP) !np

  if (aNeQb) then
  B1_part = (AngFunctDUAL_P1(ms,b,a,i_a) - AngFunctDUAL_P1(ms2,b,a,i_a))
  BulkP1_HM(1,ms,i_r,i_a) = BulkP1_HM(1,ms,i_r,i_a) + (B1_part * kaPt) !pp
  BulkP1_HM(2,ms,i_r,i_a) = BulkP1_HM(2,ms,i_r,i_a) + (B1_part * kaNt) !nn

  B1_part =    CONST_x0_EXC_HEIS * AngFunctDUAL_P1(ms ,b,a,i_a) &
            - (CONST_x0_EXC_MAJO * AngFunctDUAL_P1(ms2,b,a,i_a))
  B2_part =    CONST_x0_EXC_HEIS * AngFunctDUAL_P1(ms2,b,a,i_a) &
            - (CONST_x0_EXC_MAJO * AngFunctDUAL_P1(ms ,b,a,i_a))
  BulkP1_HM(3,ms,i_r,i_a)= BulkP1_HM(3,ms,i_r,i_a) + (B1_part*kNPt-B2_part*kPNt) !pn
  BulkP1_HM(4,ms,i_r,i_a)= BulkP1_HM(4,ms,i_r,i_a) + (B1_part*kPNt-B2_part*kNPt) !np
  endif

  endif

enddo

end subroutine compute_bulkDens4_pair

!-----------------------------------------------------------------------------!
! Auxiliary clean the Bulk fields to be pre-calculated or get the total value !
! (after the sp loop calculation in calculate_expectval_density_)             !
!-----------------------------------------------------------------------------!
subroutine resetOrSumTotalBulkDens4Fields(reset_, i_r, i_ang)
integer, intent(in) :: i_r, i_ang
logical, intent(in) :: reset_
integer ms

if (reset_) then
  BulkHF    = zzero
  BulkP1    = zzero
  BulkP2    = zzero
  if (has_HEIS_MAJO_TERMS) BulkP1_HM = zzero
else
  do ms = 1, 4
    ! BulkHF(ms,i_r,i_ang) = BulkHF_pp(ms,i_r,i_ang) + BulkHF_nn(ms,i_r,i_ang)
    BulkHF(5,ms,i_r,i_ang) = BulkHF(1,ms,i_r,i_ang) + BulkHF(2,ms,i_r,i_ang)
    BulkP1(5,ms,i_r,i_ang) = BulkP1(1,ms,i_r,i_ang) + BulkP1(2,ms,i_r,i_ang)
    BulkP2(5,ms,i_r,i_ang) = BulkP2(1,ms,i_r,i_ang) + BulkP2(2,ms,i_r,i_ang)
    if (has_HEIS_MAJO_TERMS) then
      BulkP1_HM(5,ms,i_r,i_ang) =   BulkP1_HM(1,ms,i_r,i_ang) &
                                  + BulkP1_HM(2,ms,i_r,i_ang)
    endif
  enddo
endif

end subroutine resetOrSumTotalBulkDens4Fields

!-----------------------------------------------------------------------------!
! subroutine set_rearrangement_RadAng_fucntions                               !
!                                                                             !
! DD_matrix element and the rearrangement of the DD matrix element differs    !
! only in another delta matrix element function:                              !
!         (<k1|delta()r>|k2> = R_k1k2(r) * Angular_k1k2(Omega)                !
! Calculate them previous all the computations to save time                   !
!-----------------------------------------------------------------------------!
subroutine set_rearrangement_RadAng_fucntions

integer :: a,b,aN,bN, l1,l2,j1,j2, mj1,mj2, ind_jm_1,ind_jm_2, K, M
integer :: i_r, i_ang, ind_km, spO2
real(r64)    :: radial
complex(r64) :: ang_rea

if (.NOT.EVAL_REARRANGEMENT) then
  return
endif
rea_common_RadAng = zero
spO2 = HOsp_dim / 2

do a = 1, spO2 ! avoid pn states
  l1 = HOsp_l   (a)
  j1 = HOsp_2j  (a)
  mj1 = HOsp_2mj(a)
  ind_jm_1 = angular_momentum_index(j1, mj1, .TRUE.)
  aN = a + spO2
  do b = a, spO2
    l2 = HOsp_l  (b)
    j2 = HOsp_2j (b)
    mj2 = HOsp_2mj(b)
    ind_jm_2 = angular_momentum_index(j2, mj2, .TRUE.)
    bN = b + spO2

    do i_r = 1, r_dim
      radial = radial_2b_sho_memo(HOsp_sh(a), HOsp_sh(b),i_r)

      do i_ang = 1, angular_dim
!        !! --  Process with spatial 2- body harmonics  -----------------------!
!        ang_rea = zzero
!        do K = abs(j1 - j2) / 2, (j1 + j2) / 2
!          M = (mj2 - mj1)/2
!          if ((MOD(K + l1 + l2, 2) == 1).OR.(abs(M) > K)) cycle
!
!          ind_km   = angular_momentum_index(K,  M,  .FALSE.)
!
!          ang_rea = ang_rea + (dens_Y_KM_me(ind_jm_1, ind_jm_2, ind_km) * &
!                               sph_harmonics_memo(ind_km, i_ang))
!        enddo
!        !! ------------------------------------------------------------------ !
        ang_rea = AngFunctDUAL_HF(1,a,b,i_ang) + AngFunctDUAL_HF(4,a,b,i_ang)

        rea_common_RadAng( a, b, i_r,i_ang) = radial * ang_rea
!        rea_common_RadAng(aN,bN, i_r,i_ang) = radial * ang_rea
        if (a == b) cycle
        rea_common_RadAng( b, a, i_r,i_ang) = radial * conjg(ang_rea)
!        rea_common_RadAng(bN,aN, i_r,i_ang) = radial * conjg(ang_rea)

      enddo ! end ang
    enddo ! end radial

  enddo
enddo

end subroutine set_rearrangement_RadAng_fucntions

!-----------------------------------------------------------------------------!
! subroutine set_rearrangement_RadAng_fucntions                               !
!                                                                             !
! Update the pointers for the density matrices, their values will be stored   !
! separately from the main rho/kappa matrices from WF to include cutoff       !
! conditions and easier access. By default, these matrices are just copied    !
!-----------------------------------------------------------------------------!
subroutine update_densities_DD(UL,VL,UR,VR,rho0LR,kappa0LR,kappa0RL,ndim)

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: UL, VL, UR, VR
complex(r64), dimension(ndim,ndim), intent(in) :: rho0LR, kappa0LR, kappa0RL
complex(r64), dimension(ndim,ndim) :: URc, VRc, ULc, V2
real(r64), dimension(ndim) :: eigen_H11, ener_qp
real(r64), dimension(ndim,ndim) :: D0, rhoc, hspc, A1, A2
complex(r64) :: bdbd, bdb, bbd, bb
logical,  dimension(ndim) :: excluded_qp_indx
integer   :: info_H11, ialloc = 0
real(r64), dimension(3*ndim-1) :: work
real(r64) :: ovac0
integer   :: i, j, k, l, zn_indx, nocc0,nemp0


if ((.TRUE.).OR.(.NOT. EVAL_CUTOFF).OR.(iteration .EQ. 1)) then
  !! Just update with the main
  rhoLR   = rho0LR
  kappaLR = kappa0LR
  kappaRL = kappa0RL
else

  if (PRINT_GUTS) then
    open (623, file='cutoff_elements.gut')
    open (624, file='RHOKAPPA_withCutoff.gut')
    write(623,fmt='(A)')  "k, V^2(k,k), h_eig(k), L_Ferm, e_cutoff, excluded"
    write(623,fmt='(A,I5)') "  ITER=", iteration
    write(624,fmt='(A)')  "i, j, REAL: rho/rho0, kpa/kpa_LR, kpa/kpa_RL,  IMAG"
  endif

  eigen_H11 = zero
  ener_qp   = zero
  excluded_qp_indx = .FALSE.
  rhoLR     = zzero
  kappaLR   = zzero
  kappaRL   = zzero

  URc = conjg(UR)
  ULc = conjg(VL)
  VRc = conjg(VR)
  call zgemm('t','n',ndim,ndim,ndim,zone,VL,ndim,VR,ndim,zzero,V2,ndim)

  call calculate_fields_diag(rho0LR, kappa0LR, field_gammaLR, field_hspLR, &
                             field_deltaLR, field_deltaRL, ndim)
  field_hspRR   = real(field_hspLR) + field_gammaRR_DD + field_rearrRR_DD
  !call dsyev('v','u',ndim,field_hspRR,ndim,eigen_H11,work,3*ndim-1,info_H11)

  !!!! =====================================================================
  !!! CALCULATE THE CANONICAL BASIS
  !!! OVAC= es el el solape de la norma ERROR!!!
  call construct_canonical_basis(bogo_U0,bogo_V0,bogo_zU0c,bogo_zV0c,bogo_zD0, &
                                 ovac0,nocc0,nemp0,ndim)
  D0 = real(bogo_zD0)

  call dgemm('t','n',ndim,ndim,ndim,one,D0,ndim,dens_rhoRR,ndim,zero,A1,ndim)
  call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,D0,ndim,zero,rhoc,ndim)

  call dgemm('t','n',ndim,ndim,ndim,one,D0,ndim,field_hspRR,ndim,zero,A1,ndim)
  call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,D0,ndim,zero,hspc,ndim)

  print "(A)", "  - calculate cannonical basis [DONE]"
  ! =====================================================================================

  !!! Further reduces h in case of fully empty/occupides states
!  if ( nemp0 > 0 ) then
!    allocate (hspr(nemp0,nemp0), eigenr(nemp0),workr(3*nemp0-1))
!    hspr(1:nemp0,1:nemp0) = hspc(1:nemp0,1:nemp0)
!    call dsyev('v','u',nemp0,hspr,nemp0,eigenr,workr,3*nemp0-1,info_hsp)
!    A1 = zero
!    A2 = D0
!    do i = 1, ndim
!      A1(i,i) = one
!    enddo
!    A1(1:nemp0,1:nemp0) = hspr(1:nemp0,1:nemp0)
!    call dgemm('n','n',ndim,ndim,ndim,one,A2,ndim,A1,ndim,zero,D0,ndim)
!    deallocate(hspr, eigenr, workr)
!  endif
!
!  if ( nocc0 > 0 ) then
!    allocate (hspr(nocc0,nocc0), eigenr(nocc0),workr(3*nocc0-1))
!    hspr(1:nocc0,1:nocc0) = hspc(ndim-nocc0+1:ndim,ndim-nocc0+1:ndim)
!    call dsyev('v','u',nocc0,hspr,nocc0,eigenr,workr,3*nocc0-1,info_hsp)
!    A1 = zero
!    A2 = D0
!    do i = 1, ndim
!      A1(i,i) = one
!    enddo
!    A1(ndim-nocc0+1:ndim,ndim-nocc0+1:ndim) = hspr(1:nocc0,1:nocc0)
!    call dgemm('n','n',ndim,ndim,ndim,one,A2,ndim,A1,ndim,zero,D0,ndim)
!    deallocate(hspr, eigenr, workr)
!  endif
!
!  call dgemm('t','n',ndim,ndim,ndim,one,D0,ndim,field_hspRR,ndim,zero,A1,ndim)
!  call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,D0,ndim,zero,hspc,ndim)

  ! =====================================================================================
!
!
!  !!! Ordering of energies
!  l = 0
!  eigenh_order = 0
!  eigenh_tmp = 999
!
!  do i = 1, ndim
!    if ( abs(rhoc(i,i)) > 1.d-7 ) then
!      l = l + 1
!      eigenh_tmp(i) = hspc(i,i)
!    endif
!  enddo
!
!  do i = 1, l
!    tabmin = minloc(eigenh_tmp)
!    eigenh_order(i) = tabmin(1)
!    eigenh_tmp(tabmin(1)) = 1000
!  enddo
!
!  eigenh_tmp = 999
!
!  do i = 1, ndim
!    if ( abs(rhoc(i,i)) <= 1.d-7 ) then
!      eigenh_tmp(i) = hspc(i,i)
!    endif
!  enddo
!
!  do i = l+1, ndim
!    tabmin = minloc(eigenh_tmp)
!    eigenh_order(i) = tabmin(1)
!    eigenh_tmp(tabmin(1)) = 1000
!  enddo

  !!--------------------------------------------------------------------
!  do i = 1, ndim
!    open(ute, file='canonicalbasis.dat', status='replace', action='write', &
!           form='formatted')
!    k = eigenh_order(i)
!
!!    write(ute,format1) i, xprot, xneut, xn, xl, xpar, xj, xjz, rhoc(k,k), &
!!                       hspc(k,k)
!  enddo
  !!--------------------------------------------------------------

  do i = 1, ndim
    zn_indx = 1
    if (i > ndim/2) zn_indx = 2

    ener_qp(i) = (1 - rhoc(i,i))*hspc(i, i)!eigen_H11(i) !+ lambdaFer_DD(zn_indx)
    if (ener_qp(i) .GT. CUTOFF_ENERGY_MAX) excluded_qp_indx(i) = .TRUE.

    if (PRINT_GUTS) write(623,fmt='(I5,4F15.9,L3)') &
        i, rhoc(i,i), &
        !dreal(V2(i,i)), &
        hspc(i, i), &
        !eigen_H11(i), &
        lambdaFer_DD(zn_indx), ener_qp(i), excluded_qp_indx(i)
  enddo
  print "(A)", "  - calculate cutoffs for single particles [DONE]"
  do i = 1, ndim
    do j = 1, ndim
!      do k = 1, ndim
!        if (excluded_qp_indx(k)) cycle
!        rhoLR  (i,j) = rhoLR  (i,j) + VRc(i,k) * VL (j,k)
!        kappaLR(i,j) = kappaLR(i,j) + VRc(i,k) * UL (j,k)
!        kappaRL(i,j) = kappaRL(i,j) + VL (i,k) * URc(j,k)
!      enddo
      !! Exclusion from the rho-kappa definition
      bdbd = zzero
      bdb  = zzero
      bbd  = zzero
      bb   = zzero
      do k = 1, ndim
        if (excluded_qp_indx(k)) cycle
        do l = 1, ndim
          if (excluded_qp_indx(l)) cycle
          bdbd = bdbd + URc(j,l) * VRc(i,k)
          bdb  = bdb  + URc(j,l) * VR (i,k)
          bbd  = bbd  + VR (j,l) * VRc(i,k)
          bb   = bb   + VR (j,l) * UR (i,k)
        end do
      end do
      rhoLR  (i,j) = rhoLR  (i,j) + (bdbd + bdb + bbd + bb)
      !! kappa
      bdbd = zzero
      bdb  = zzero
      bbd  = zzero
      bb   = zzero
      do k = 1, ndim
        if (excluded_qp_indx(k)) cycle
        do l = 1, ndim
          if (excluded_qp_indx(l)) cycle
          bdbd = bdbd + VRc(j,l) * VRc(i,k)
          bdb  = bdb  + VRc(j,l) * UR (i,k)
          bbd  = bbd  + UR (j,l) * VRc(i,k)
          bb   = bb   + UR (j,l) * UR (i,k)
        end do
      end do
      kappaLR(i,j) = kappaLR(i,j) + (bdbd + bdb + bbd + bb)
      bdbd = zzero
      bdb  = zzero
      bbd  = zzero
      bb   = zzero
      do k = 1, ndim
        if (excluded_qp_indx(k)) cycle
        do l = 1, ndim
          if (excluded_qp_indx(l)) cycle
          bdbd = bdbd + URc(j,l) * URc(i,k)
          bdb  = bdb  + VR (j,l) * URc(i,k)
          bbd  = bbd  + URc(j,l) * VR (i,k)
          bb   = bb   + VR (j,l) * VR (i,k)
        end do
      end do
      kappaRL(i,j) = kappaRL(i,j) + (bdbd + bdb + bbd + bb)

      if(PRINT_GUTS) then
        write(624,fmt='(2I5,6F15.9,A,6F15.9)') i, j, &
        dreal(rhoLR(i,j)),    dreal(rhoLR(i,j)),   dreal(kappaLR (i,j)), &
        dreal(kappa0LR(i,j)), dreal(kappaRL(i,j)), dreal(kappa0RL(i,j)), &
        "    ", &
        dimag(rhoLR(i,j)),    dimag(rhoLR(i,j)),   dimag(kappaLR(i,j)), &
        dimag(kappa0LR(i,j)), dimag(kappaRL(i,j)), dimag(kappa0RL(i,j))
      endif
    end do
  end do

  if (PRINT_GUTS) close(623)
  if (PRINT_GUTS) close(624)
endif

!
!
!call zgemm('n','t',ndim,ndim,ndim,zone,VRc,ndim,VL,ndim,zzero,rhoLR,ndim)
!call zgemm('n','t',ndim,ndim,ndim,zone,VRc,ndim,UL,ndim,zzero,kappaLR,ndim)
!call zgemm('n','t',ndim,ndim,ndim,zone,VL,ndim,URc,ndim,zzero,kappaRL,ndim)
end subroutine update_densities_DD

!-----------------------------------------------------------------------------!
! subroutine calculate_expectval_density                                      !
!                                                                             !
!   Choose the correct density power when it is a complex number,             !
! the criteria is to choose the root such as: (z)^alpha = (z*)^alpha          !
!                                                                             !
!-----------------------------------------------------------------------------!
subroutine choose_riemann_fold_density(i_r, i_an)

integer, intent(in) :: i_r, i_an
integer      :: i, j
real(r64)    :: dens_R, dens_A, dens_Aa, dens_Ra, x1, x2, y1, y2, th1, th2
complex(r64) :: x

if (abs(alpha_DD - 0.3333333) .GE. 1.0d-6) then
  print "(A)", " [WARNING] If using alpha_DD /= 1/3 rethink this method 1mdd5c"
endif

dens_R = dreal(density(i_r,i_an))**2 + dimag(density(i_r,i_an))**2
dens_R = dsqrt(dens_R)
!dens_A = dacos(dreal(density(i_r, i_an)) / max(dens_R, 1.0d-30))
dens_A = datan2(dreal(density(i_r,i_an)), dimag(density(i_r,i_an)))

dens_Ra = dens_R ** alpha_DD
do i = 0, alpha_DD_frac(2) - 1
  th1 = (dens_A + (2 * pi * i)) * alpha_DD
  x1  = dcos(th1)
  y1  = dsin(th1)

  do j = 0, alpha_DD_frac(2) - 1
    th2 = (dens_A + (2 * pi * j)) * alpha_DD
    x2  = dcos(th2)
    y2  = -1.0d0 * dsin(th2)

!    print "(2(A,3F12.6))", "    ** ", th1, x1, y1, " ?= ", th2, x2, y2

!    if ((abs(x1*x2) > 0) .AND. (abs(y1*y2) > 0)) then
!      !! condition z* be in the same sector (first coincidence))
!      dens_alpha(i_r,i_an) = dCMPLX(dens_Ra * x1, dens_Ra * y1)
!
!      th1 = (dens_A + 2 * pi * i) * (alpha_DD - 1.0d0)
!      x1  = dcos(th1)
!      y1  = dsin(th1)
!      dens_Ra = dens_Ra / dens_R
!
!      x = dCMPLX(dens_Ra * x1, dens_Ra * y1)
!      if (dreal(x)**2 + dimag(x)**2 .gt. 1.0D+30) then
!        x = dCMPLX(1.0D+30*x1, 1.0D+30*y1)
!      endif
!      dens_alpm1(i_r,i_an) = x
!      return
!    end if

    !! condition for z*=z mathch (impossible)
    if (abs(x1 - x2) + abs(y1 - y2) .LT. 1.0d-4) then

      dens_alpha(i_r,i_an) = dCMPLX(dens_Ra * x1, dens_Ra * y1)

      th1 = (dens_A + 2 * pi * i) * (alpha_DD - 1.0d0)
      x1  = dcos(th1)
      y1  = dsin(th1)
      dens_Ra = dens_Ra / dens_R

      x = dCMPLX(dens_Ra * x1, dens_Ra * y1)
      if (dreal(x)**2 + dimag(x)**2 .gt. 1.0D+30) then
        x = dCMPLX(1.0D+30*x1, 1.0D+30*y1)
      endif
      dens_alpm1(i_r,i_an) = x

      return
    endif

  enddo
enddo

print "(A,2I5,2E16.9)", "[ERROR] Could not find z^for density (ir,ia)=", &
                i_r, i_an, dreal(density(i_r,i_an)), dimag(density(i_r,i_an))
print "(A)", ""
! Fold the density to the 1st quadrant. (VERSION 1 - REMOVE)
!dens_R = dreal(density(i_r,i_an))**2 + dimag((density(i_r,i_an)))**2
!dens_R = dsqrt(dens_R)
!dens_A = dacos(dreal(density(i_r, i_an)) / max(dens_R, 1.0d-30))
!if (dsin(dens_A) .lt.  zero ) then
!  dens_A = dens_A + 2*(pi - dens_A)
!endif
!
!dens_Aa = dens_A *  alpha_DD
!dens_Ra = dens_R ** alpha_DD
!
!dens_alpha(i_r,i_an) = dCMPLX(dens_Ra * (dcos(dens_Aa)), &
!                              dens_Ra * (dsin(dens_Aa)) )
!dens_Aa = dens_A *  (alpha_DD - 1)
!dens_Ra = dens_R ** (alpha_DD - 1)
!x       = dCMPLX(dens_Ra * (dcos(dens_Aa)), &
!                 dens_Ra * (dsin(dens_Aa)) )
!if (dreal(x)**2 + dimag(x)**2 .gt. 1.0D+30) then
!  x = dCMPLX(1.0D+30*cos(dens_A), 1.0D+30*sin(dens_A))
!endif
!dens_alpm1(i_r,i_an) = x

end subroutine choose_riemann_fold_density

!-----------------------------------------------------------------------------!
! subroutine calculate_expectval_density                                      !
!                                                                             !
! iopt = optimal iteration (=1 when gradient has converged)                   !
!                                                                             !
!-----------------------------------------------------------------------------!
subroutine calculate_expectval_density(overlap, iopt)

integer, intent(in) :: iopt
complex(r64), intent(in) :: overlap

integer :: a,b, a_sh, b_sh, spO2, ITER_PRNT
integer :: i_r=1, i_an=1, msp

real(r64)    :: radial_part, dens_R, dens_A, rad4Integr, dens_Aa,dens_Ra
complex(r64) :: sum_, integral_dens, sum_test, diff, x
logical      :: PRNT_

PRNT_ = (PRINT_GUTS).OR.(.FALSE.)
spO2  = HOsp_dim / 2

density   = zzero
dens_pnt  = zzero
integral_dens = zzero

call resetOrSumTotalBulkDens4Fields(.TRUE., 0, 0)

if (PRNT_) then
  open(555, file='dens_pnt.gut')
  open(556, file='Bulk_pnt.gut')
  open(551, file='Bulk1_pnt.gut')
  open(552, file='Bulk2_pnt.gut')
  write(555,fmt="(2A,2F8.4)") "i_r, i_ang      dens_pnt(pp)  nn   pn  np",&
    "  total , B alpha::", HO_b, alpha_DD
  write(556,fmt="(2A,2F8.4)") "i_r, i_ang  msp    BulkHF_pp  nn   pn  np",&
    "  total , B alpha::", HO_b, alpha_DD
  write(551,fmt="(2A,2F8.4)") "i_r, i_ang  msp    BulkP1_pp  nn   pn  np",&
    "  total , B alpha::", HO_b, alpha_DD
  write(552,fmt="(2A,2F8.4)") "i_r, i_ang  msp    BulkP2_pp  nn   pn  np",&
    "  total , B alpha::", HO_b, alpha_DD
endif

do i_r = 1, r_dim
  !! [TEST] For the density to be integrated in the variable for the m.e.
  rad4Integr =  weight_R(i_r) * exp((r(i_r)/HO_b)**2 * (1.0+alpha_DD))

  do i_an = 1, angular_dim

    do a = 1, spO2
       a_sh = HOsp_sh(a)
       do b = a, spO2         !!!!!     BENCH REQUIRES B=1      !!!!
         b_sh = HOsp_sh(b)

         call compute_bulkDens4_hf  (a, b, a_sh, b_sh, i_r, i_an, overlap)
         call compute_bulkDens4_pair(a, b, a_sh, b_sh, i_r, i_an, overlap)
      enddo ! do b
    enddo   ! do a

    dens_pnt(5,i_r,i_an) = dens_pnt(1,i_r,i_an) + dens_pnt(2,i_r,i_an)
    density(i_r,i_an)    = dens_pnt(5,i_r,i_an)

    call resetOrSumTotalBulkDens4Fields(.FALSE., i_r, i_an)

    !! [TEST] For the density to be integrated in the variable for the m.e.
    integral_dens = integral_dens + (dreal(density(i_r, i_an) * &
                                           exp( (r(i_r)/HO_b)**2)  * &
                                     weight_LEB(i_an) * rad4Integr))

    !!! calculate the density powered to alpha_DD for the matrix elements
    if (dabs(imag(density(i_r, i_an))) > 1.0e-10) then
      !! This part must not be executed for Mean Field D1S (projections yield
      !! larger values for the imaginary part), put the limit in 1E-13 precision
      !! to prompt an abnormal numerical error (r64 must have 13 decimals)
      if (.NOT.DOING_PROJECTION) then
        print *, " !!! [WARNING] density is imaginary=",imag(density(i_r,i_an))
      endif

      call choose_riemann_fold_density(i_r, i_an)

    else
      if (FUNCTIONAL_DENS_MODE .EQ. 1) then
        dens_alpha(i_r,i_an) = dreal(density(i_r,i_an)) ** alpha_DD
        dens_alpm1(i_r,i_an) = dens_alpha(i_r,i_an) / density(i_r,i_an)
        dens_alpm1(i_r,i_an) = MIN(dreal(dens_alpm1(i_r,i_an)), 1.0D+30)
      else
        dens_alpha(i_r,i_an) = 1.0d+00 - &
                          CONST_EDD_M2 * (dreal(density(i_r,i_an)) ** alpha_DD)
        dens_alpm1(i_r,i_an) = -CONST_EDD_M2 * &
                          (dreal(density(i_r,i_an)) ** (alpha_DD - 1.0d+00))
        dens_alpm1(i_r,i_an) = MIN(dreal(dens_alpm1(i_r,i_an)), 1.0D+30)
      endif
    endif

    if (PRNT_) then
      write(555,fmt="(2I5)", advance='no') i_r, i_an
      write(555,fmt='(7(F22.15,SP,F20.15,"j"))') &
        dens_pnt(1,i_r,i_an), dens_pnt(2,i_r,i_an), dens_pnt(3,i_r,i_an), &
        dens_pnt(4,i_r,i_an), dens_pnt(5,i_r,i_an), dens_alpha(i_r,i_an), &
        dens_alpm1(i_r, i_an)
      do msp =1, 4
        write(556,fmt="(3I5,A,3F9.5,A)", advance='no') i_r, i_an, msp, ',', &
          r(i_r), cos_th(i_an), phi(i_an), ','
        write(551,fmt="(3I5,A,3F9.5,A)", advance='no') i_r, i_an, msp, ',', &
          r(i_r), cos_th(i_an), phi(i_an), ','
        write(552,fmt="(3I5,A,3F9.5,A)", advance='no') i_r, i_an, msp, ',', &
          r(i_r), cos_th(i_an), phi(i_an), ','

        write(556,fmt='(5(F22.15,SP,F20.15,"j"))') &
          BulkHF(1,msp,i_r,i_an), BulkHF(2,msp,i_r,i_an), &
          BulkHF(3,msp,i_r,i_an), BulkHF(4,msp,i_r,i_an), BulkHF(5,msp,i_r,i_an)
        write(551,fmt='(5(F22.15,SP,F20.15,"j"))') &
          BulkP1(1,msp,i_r,i_an), BulkP1(2,msp,i_r,i_an), &
          BulkP1(3,msp,i_r,i_an), BulkP1(4,msp,i_r,i_an), BulkP1(5,msp,i_r,i_an)
        write(552,fmt='(5(F22.15,SP,F20.15,"j"))') &
          BulkP2(1,msp,i_r,i_an), BulkP2(2,msp,i_r,i_an), &
          BulkP2(3,msp,i_r,i_an), BulkP2(4,msp,i_r,i_an), BulkP2(5,msp,i_r,i_an)
        if (msp == 4) then
          write(556, fmt="(A)") ""
          write(551, fmt="(A)") ""
          write(552, fmt="(A)") ""
        endif
      enddo
    endif ! PRNT_ case

  enddo   !end do i_angular
enddo ! do i_r

if (PRNT_) then
  close(555)
  close(556)
  close(551)
  close(552)

  call test_integrate_bulk_densities
endif

!! [TEST] For the density to be integrated in the variable for the m.e.
ITER_PRNT = 10
if (DOING_PROJECTION) ITER_PRNT = ITER_PRNT * Mphip_DD * Mphin_DD
if ((iteration.eq.0).OR.(MOD(iteration + 1, ITER_PRNT).EQ.0)) then
  integral_dens = integral_dens * 2 * pi * (HO_b**3) / ((2.0 + alpha_DD)**1.5)
  print "(A,F13.9,A)", "      *A* ", dreal(integral_dens),"  <dens(r)> approx "
endif

end subroutine calculate_expectval_density



!------------------------------------------------------------------------------!
! subroutine calculate_common_rearrang_bulkFields                              !
!                                                                              !
! This subroutine update the common part to the new densities and bulk funct.s !
! stored in the array [REACommonFields]                                        !
! The total function is independent of sp indexes, this include only the normal!
! density-matrices that not cross the protons-neutrons.                        !
!------------------------------------------------------------------------------!
subroutine calculate_common_rearrang_bulkFields
integer      :: i_r, i_a, ms, ms2
complex(r64) :: aux_d, aux_e, aux_p, aux_pnp, aux1, aux2, aux3, aux4
real         :: X0M1

REACommonFields = zero
if (.NOT.EVAL_REARRANGEMENT) return
if (PRINT_GUTS) then
  open(665, file='BulkREA_elements.gut')
  write(665,fmt='(3A)') &
  "i_r, i_ang      B_HF_Direct_REA                               ",&
  "B_HF_Exchange_REA                          B_PairingPPNN_REA              ",&
  "              B_Pairing_PN_REA                          REA_CommonField"
endif
X0M1 = 1.0d0 - x0_DD_FACTOR

do i_r = 1, r_dim
  do i_a = 1, angular_dim

    aux_d =  dens_pnt(5,i_r,i_a)**2
    aux1  = (dens_pnt(1,i_r,i_a)**2) + (dens_pnt(2,i_r,i_a)**2)

    aux2 = zzero  ! pn np part
    if (CALCULATE_DD_PN_HF) then
      aux2  = 2.0d0 * dens_pnt(3,i_r,i_a) * dens_pnt(4,i_r,i_a)
      endif
    !! dens_pnt are complex, a test is to verify aux2 with the following to be Real
    !aux2 = 2.0d0*(dreal(dens_pnt(3,i_r,i_a))**2 - dimag(dens_pnt(3,i_r,i_a))**2)

    aux_d = aux_d - (x0_DD_FACTOR * (aux1 + aux2))

    aux_e = zzero
    aux_p = zzero
    aux_pnp = zzero
    do ms = 1, 4
      select case (ms)
        case (2, 3)
          ms2 = 5 - ms
        case default
          ms2 = ms
      end select

      !Exchange part of the HF like fields
      aux1  = BulkHF(1,ms2, i_r,i_a) * BulkHF(1,ms,i_r,i_a) !pp
      aux2  = BulkHF(2,ms2, i_r,i_a) * BulkHF(2,ms,i_r,i_a) !nn
      aux_e = aux_e + (aux1  + aux2)
        ! pn np part
      if (CALCULATE_DD_PN_HF) then
        aux1  = BulkHF(3,ms2, i_r,i_a) * BulkHF(4,ms,i_r,i_a) !pn*np
        aux2  = BulkHF(4,ms2, i_r,i_a) * BulkHF(3,ms,i_r,i_a) !np*pn
        aux_e = aux_e + (aux1 + aux2)
      end if
        !total field part
      aux1  = BulkHF(5,ms2, i_r,i_a) * BulkHF(5,ms,i_r,i_a) !tot
      aux_e = aux_e - (x0_DD_FACTOR * aux1)

      !Pairing rearrangement fields
      if (hasX0M1) then
        aux1  = BulkP2(1,ms, i_r,i_a) * BulkP1(1,ms, i_r,i_a) !pp
        aux2  = BulkP2(2,ms, i_r,i_a) * BulkP1(2,ms, i_r,i_a) !nn
        aux_p = aux_p + (aux1 + aux2)
      endif
      !pn np part (remember the 1Bpn + x0*1Bpn - 1Bnp - x0*1Bnp was done already)
      if (CALCULATE_DD_PN_PA) then
        aux1   = BulkP2(3,ms, i_r,i_a) * BulkP1(3,ms, i_r,i_a) !pn*pn
        aux2   = BulkP2(4,ms, i_r,i_a) * BulkP1(4,ms, i_r,i_a) !np*np
        aux_pnp = aux_pnp + (aux1 + aux2)
      endif

    enddo ! loop ms
    !! change 11/11/22 + sings of pairing changed to - (from -k*_ab k_cd)
    !! since (K_RL)* is the definition of the contraction, the sign is + for pairing
    aux1 = (2.0d+0*(aux_d - aux_e)) + (X0M1*aux_p) + aux_pnp

    if ((.FALSE.).AND.(dimag(aux1) > 1.0e-12)) then
      print '(A,2I4,5F20.15)',"Error!! Rearr. funct. is imaginary: ",i_r,i_a,&
        dimag(aux1), dimag(aux_d), dimag(aux_e), dimag(aux_p), dimag(aux_pnp)
    end if
    REACommonFields(i_r, i_a) = aux1   !! dreal(aux1)

    if (has_HEIS_MAJO_TERMS) call calculate_rearrang_bulkFields_HM(i_r, i_a)

    ! export Rea Fields by parts
    if (PRINT_GUTS) then
      write(665, fmt='(2I4,5(F22.15,SP,F20.15,"j"))') &
        i_r, i_a, aux_d, aux_e, aux_p, aux_pnp, REACommonFields(i_r, i_a)
    endif
  enddo
enddo
if (PRINT_GUTS) close(665)

end subroutine calculate_common_rearrang_bulkFields

!!
!! ! Extension of the rearrangement for Heisenberg-Majorana exchange operators
!!
subroutine calculate_rearrang_bulkFields_HM(i_r, i_a)
integer, intent(in) :: i_r, i_a
integer      ::  ms, ms2
complex(r64) :: aux_d, aux_e, aux_p, aux_pnp, aux1, aux2, aux3, aux4
real(r64)    :: X0MpH

X0MpH = CONST_x0_EXC_HEIS + CONST_x0_EXC_MAJO

if (.NOT. has_HEIS_MAJO_TERMS) return

aux_d =  dens_pnt(5,i_r,i_a)**2
aux1  = (dens_pnt(1,i_r,i_a)**2) + (dens_pnt(2,i_r,i_a)**2)
aux2 = zzero  ! pn np part
if (CALCULATE_DD_PN_HF) then
  aux2  = 2.0d0 * dens_pnt(3,i_r,i_a) * dens_pnt(4,i_r,i_a)
  endif
aux_d = (CONST_x0_EXC_HEIS * (aux1 + aux2)) - (CONST_x0_EXC_MAJO * aux_d)

aux_e = zzero
aux_p = zzero
aux_pnp = zzero
do ms = 1, 4
  select case (ms)
    case (2, 3)
      ms2 = 5 - ms
    case default
      ms2 = ms
  end select

  !Exchange part of the HF like fields
  aux1  = BulkHF(1,ms2, i_r,i_a) * BulkHF(1,ms,i_r,i_a) !pp
  aux2  = BulkHF(2,ms2, i_r,i_a) * BulkHF(2,ms,i_r,i_a) !nn
  aux_e = aux_e + (CONST_x0_EXC_MAJO * (aux1 + aux2))
  ! pn np part in the following <if>

  !total field part
  aux1  = BulkHF(5,ms2, i_r,i_a) * BulkHF(5,ms,i_r,i_a) !tot
  aux_e = aux_e - (CONST_x0_EXC_HEIS * aux1)

  !Pairing rearrangement fields
  aux1  = BulkP2(1,ms, i_r,i_a) * BulkP1(1,ms, i_r,i_a) !pp
  aux2  = BulkP2(2,ms, i_r,i_a) * BulkP1(2,ms, i_r,i_a) !nn
  aux_p = aux_p + (aux1 + aux2)
  !pn np part (remember the H 1Bnp - M*1Bpn - H 1Bpn  + M*1Bpn was done already)
  if (CALCULATE_DD_PN_PA) then
    aux1   = BulkP2(3,ms, i_r,i_a) * BulkP1(3,ms, i_r,i_a) !pn*pn
    aux2   = BulkP2(4,ms, i_r,i_a) * BulkP1(4,ms, i_r,i_a) !np*np
    aux_pnp = aux_pnp + (aux1 + aux2)

    aux1  = BulkHF(4,ms2, i_r,i_a) * BulkHF(3,ms,i_r,i_a) !pn*np
    aux2  = BulkHF(3,ms2, i_r,i_a) * BulkHF(4,ms,i_r,i_a) !np*pn
    aux_e = aux_e + (CONST_x0_EXC_MAJO * (aux1 + aux2))
  endif

enddo ! loop ms

aux1 = (2.0D+00*(aux_d + aux_e)) + (X0MpH*aux_p) + aux_pnp

!! Append to normal DD exchange
REACommonFields(i_r, i_a) = REACommonFields(i_r, i_a) - aux1

end subroutine calculate_rearrang_bulkFields_HM

!------------------------------------------------------------------------------!
! subroutine complete_DD_fields                                                !
!     Auxiliary subroutine to complete the fields over p/n indexes for element !
!   a, c <= HOspdim / 2 = spO2. int_hf/pa/rea are the (a,c) field final        !
!  integrals (it must contain the integral_factor and the alphaDD/4).
!------------------------------------------------------------------------------!
subroutine complete_DD_fields(int_hf, int_pa, int_rea, &
                              gammaLR, deltaLR, deltaRL, hspLR,&
                              gammaLR_DD, deltaLR_DD, deltaRL_DD, &
                              a, c, spO2, ndim)
integer, intent(in) :: ndim, spO2, a, c
complex(r64), dimension(ndim,ndim) :: gammaLR, deltaLR, deltaRL, hspLR, &
                            gammaLR_DD, deltaLR_DD, deltaRL_DD
complex(r64), dimension(4), intent(in) :: int_hf, int_pa ! all arrays are for (pp, nn, pn, np)
complex(r64), intent(in) :: int_rea
integer :: Tac, aa, cc

do Tac = 1, 4
  aa = a
  cc = c
  select case (Tac)
    case (2)        ! nn
      aa = a + spO2
      cc = c + spO2
      !rearrange values are non zero for pp and nn
      rearrang_field(aa, cc) = int_rea
      field_rearrRR_DD(aa,cc) = dreal(rearrang_field(aa,cc))
    case (3)        ! pn
      cc = c + spO2
    case (4)        ! np
      aa = a + spO2
    case default    ! pp
      rearrang_field(aa, cc) = int_rea
      field_rearrRR_DD(aa,cc) = dreal(rearrang_field(aa,cc))
  end select

  gammaLR_DD(aa, cc) = int_hf(Tac)
  deltaLR_DD(aa, cc) = int_pa(Tac)
  deltaRL_DD(aa, cc) = int_pa(Tac)

  field_gammaRR_DD(aa,cc) = dreal(gammaLR_DD(aa, cc))
  field_deltaRR_DD(aa,cc) = dreal(deltaLR_DD(aa, cc))
  !! sum to the main fields
  gammaLR(aa,cc) = gammaLR(aa,cc) + gammaLR_DD(aa,cc)
  deltaLR(aa,cc) = deltaLR(aa,cc) + deltaLR_DD(aa,cc)
  deltaRL(aa,cc) = deltaRL(aa,cc) + deltaRL_DD(aa,cc)
  hspLR  (aa,cc) = hspLR  (aa,cc) + gammaLR_DD(aa,cc)

  !! copy to the off-diagonal values and sum to the main fields
  if (a.NE.c) then ! all off-diagonal for all (diagonal is already done)
    gammaLR_DD(cc,aa) =  gammaLR_DD(aa,cc)
    deltaLR_DD(cc,aa) = -deltaLR_DD(aa,cc)
    deltaRL_DD(cc,aa) = -deltaRL_DD(aa,cc)

    field_gammaRR_DD(cc,aa) = dreal(gammaLR_DD(cc,aa))
    field_deltaRR_DD(cc,aa) = dreal(deltaLR_DD(cc,aa))

    ! adding to the program fields
    gammaLR(cc,aa) = gammaLR(cc,aa) + gammaLR_DD(cc,aa)
    deltaLR(cc,aa) = deltaLR(cc,aa) + deltaLR_DD(cc,aa)
    deltaRL(cc,aa) = deltaRL(cc,aa) + deltaRL_DD(cc,aa)
    hspLR  (cc,aa) = hspLR  (cc,aa) + gammaLR_DD(cc,aa)
  endif

  ! do rearrangement for pp, nn if proceeds
  if((EVAL_REARRANGEMENT).AND.(Tac.LT.3)) then
    hspLR(aa,cc) = hspLR(aa,cc) + rearrang_field(aa,cc)
    field_rearrRR_DD(aa,cc) = dreal(rearrang_field(aa,cc))
    ! this step is just to complete the matrix, the rearrangement is aa,cc as well
    if (a.NE.c) then
      rearrang_field(cc, aa) =  rearrang_field(aa, cc)
      hspLR(cc,aa) = hspLR(cc,aa) + rearrang_field(cc,aa)
      field_rearrRR_DD(cc,aa) = dreal(rearrang_field(cc,aa))
    endif
  endif

enddo !Tac

end subroutine complete_DD_fields

!------------------------------------------------------------------------------!
!    subroutine calculate_fields_DD_HM
! If including Heisenberg - Majorana exchange in the DD term, complete
! with this part with respect to the different angular-parts of the fields
! Notice: these fields have sign (-1), and exchange terms for HF has to be
!    inversed (* -1) since the current definition of from usual D1S interaction
!------------------------------------------------------------------------------!
subroutine calculate_fields_DD_HM(a,c, i_r,i_ang, auxHfDio, auxHfEio, aux_PEio)

integer,   intent(in) :: a, c, i_r, i_ang
complex(r64), dimension(4) :: auxHfDio, auxHfEio, aux_PEio ! all arrays are for (pp, nn, pn, np)

complex(r64), dimension(4) :: auxHfD, auxHfE, aux_PE, aux
complex(r64) :: sumD_ang
integer      :: ms, Tac
real(r64)    :: X0MpH

X0MpH = CONST_x0_EXC_HEIS + CONST_x0_EXC_MAJO

if (.NOT. (has_HEIS_MAJO_TERMS)) return

auxHfD = zzero
!! DIRECT terms for the HF field
sumD_ang  = AngFunctDUAL_HF(1,a,c,i_ang) + AngFunctDUAL_HF(4,a,c,i_ang)
auxHfD(1) = (CONST_x0_EXC_HEIS*dens_pnt(1,i_r,i_ang) - &
             CONST_x0_EXC_MAJO*dens_pnt(5,i_r,i_ang))
auxHfD(2) = (CONST_x0_EXC_HEIS*dens_pnt(2,i_r,i_ang) - &
             CONST_x0_EXC_MAJO*dens_pnt(5,i_r,i_ang))
if (CALCULATE_DD_PN_HF) then
  auxHfD(3) = -CONST_x0_EXC_HEIS * dens_pnt(3,i_r,i_ang)
  auxHfD(4) = -CONST_x0_EXC_HEIS * dens_pnt(4,i_r,i_ang)
  endif
do Tac = 1, 4
  auxHfD(Tac)   = sumD_ang * auxHfD(Tac)

  if (a.lt.10 .AND. c.le.10 .and. i_r.eq.4 .AND. i_ang.EQ.20) then
    print "(A,3I4,2F10.5)", "DIRE a,c,tac=",a,c,Tac, real(auxHfDio(Tac)),&
      real(auxHfD(Tac))
    if (Tac .EQ. 4) print *, ""
  end if

  auxHfDio(Tac) = auxHfDio(Tac) - auxHfD(Tac)
enddo

!! EXCHANGE terms for the HF fields
auxHfE = zzero
aux_PE = zzero
aux = zzero
do ms = 1, 4
  aux(ms) = CONST_x0_EXC_MAJO * BulkHF(1,ms, i_r, i_ang) !pp
  aux(ms) = aux(ms) - (CONST_x0_EXC_HEIS * BulkHF(5,ms, i_r, i_ang)) !tot
  aux(ms) = aux(ms) * AngFunctDUAL_HF(ms, a,c, i_ang)
  auxHfE(1)  = auxHfE(1) + aux(ms)

  aux(ms) = CONST_x0_EXC_MAJO * BulkHF(2,ms, i_r, i_ang) !nn
  aux(ms) = aux(ms) - (CONST_x0_EXC_HEIS * BulkHF(5,ms, i_r, i_ang))
  aux(ms) = aux(ms) * AngFunctDUAL_HF(ms, a,c, i_ang)
  auxHfE(2)  = auxHfE(2) + aux(ms)
    !pn np part
  if (CALCULATE_DD_PN_PA) then
  aux(ms) = AngFunctDUAL_HF(ms,a,c, i_ang) * BulkHF(3,ms, i_r,i_ang) !pn
  auxHfE(3)  = auxHfE(3) + (aux(ms) * CONST_x0_EXC_MAJO)
  aux(ms) = AngFunctDUAL_HF(ms,a,c, i_ang) * BulkHF(4,ms, i_r,i_ang) !np
  auxHfE(4)  = auxHfE(4) + (aux(ms) * CONST_x0_EXC_MAJO)
  endif

  !! NOTE: Angular 1, 2 functions are defined with direct form of ms,ms'
  if (hasX0M1) then
    aux(ms) = AngFunctDUAL_P2(ms,a,c,i_ang) * BulkP1_HM(1,ms,i_r,i_ang) !pp
    aux_PE(1) = aux_PE(1)  + (X0MpH * aux(ms))
    aux(ms) = AngFunctDUAL_P2(ms,a,c,i_ang) * BulkP1_HM(2,ms,i_r,i_ang) !nn
    aux_PE(2) = aux_PE(2)  + (X0MpH * aux(ms))
  endif

  !! pn np part, x0 dependence was calculated in BulkP1_**
  if (CALCULATE_DD_PN_PA) then
  aux(ms) = AngFunctDUAL_P2(ms,a,c, i_ang) * BulkP1_HM(3,ms, i_r,i_ang) !pn
  aux_PE(3)  = aux_PE(3)  + aux(ms)
  aux(ms) = AngFunctDUAL_P2(ms,a,c, i_ang) * BulkP1_HM(4,ms, i_r,i_ang) !np
  aux_PE(4)  = aux_PE(4)  + aux(ms)
  endif

enddo ! ms loop

do Tac = 1, 4
  if (a.lt.10 .AND. c.le.10 .and. i_r.eq.4 .AND. i_ang.EQ.20) then
    print "(A,3I4,4F10.5)", "EXCH-PAIR a,c,tac=",a,c,Tac, real(auxHfEio(Tac)),&
      real(auxHfE(Tac)), real(aux_PEio(Tac)), -real(aux_PE(Tac))
    if (Tac .EQ. 4) print *, ""
  end if
  auxHfEio(Tac) = auxHfEio(Tac) + auxHfE(Tac) ! + = -(from fields) * - (HF-Exch is substracted)
  aux_PEio(Tac) = aux_PEio(Tac) - aux_PE(Tac)
enddo

end subroutine calculate_fields_DD_HM
!------------------------------------------------------------------------------!
!    subroutine calculate_fields_DD                                            !
! speeds up the bench process by reading half the N/2 space, perform Trace test!
! each 10 iterations and at the first iteration.                               !
!------------------------------------------------------------------------------!
subroutine calculate_fields_DD(gammaLR,hspLR,deltaLR,deltaRL,ndim)

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim) :: gammaLR, hspLR, deltaLR, deltaRL
!! The density fields are added to the calculated with the standard hamiltonian
!! This array variables are local
complex(r64), dimension(ndim,ndim) :: gammaLR_DD, deltaLR_DD, deltaRL_DD
complex(r64), dimension(ndim,ndim) :: gam1, dltLR1, dltRL1, hsp1,A1,A2,A3

integer   :: a, b, c, d, i, j, spO2, i_r, i_ang, &
             K, M, ind_jm_a, ind_jm_c, ind_km, Tac, &
             a_sh, ja, la, ma, ta, c_sh, jc, lc, mc, tc, ms, &
             PAR, OpPAR, IP, IDN, OpIDN, R_PRINT=3, ANG_PRINT=10
real(r64) :: rad_ac, X0M1, integral_factor
complex(r64), dimension(4) :: int_hf, int_pa, &
            auxHfD, auxHfE, aux_PE, aux_pair, aux_hf ! all arrays are for (pp, nn, pn, np)

complex(r64), dimension(4,4,ndim, ndim) :: int_test_PE ! (tt)(msms')(a)(b)
complex(r64) :: sumD_ang, auxRea, int_rea, testaux,del_aux, t_rea_sum,t_gam_sum
complex(r64), dimension(4) :: aux
logical :: PRNT_, doTraceTest_

PRNT_ = (PRINT_GUTS).OR.(.FALSE.)
iteration = iteration + 1
doTraceTest_ = (iteration.eq.1).OR.(MOD(iteration, 10).EQ.0)

if(doTraceTest_) then !! copy the Fields to plot in the case of printing.
  do i = 1, ndim
    do j = 1, ndim
      gam1  (i,j) = gammaLR(i,j)
      dltLR1(i,j) = deltaLR(i,j)
      dltRL1(i,j) = deltaRL(i,j)
      hsp1  (i,j) = hspLR  (i,j)
    enddo
  enddo
endif

spO2 = HOsp_dim / 2

gammaLR_DD = zzero
deltaLR_DD = zzero
deltaRL_DD = zzero
rearrang_field = zzero

field_gammaRR_DD = zero
field_deltaRR_DD = zero
field_rearrRR_DD = zero

!! Note :: Remember that radial functions already have the factor 1/b**3
integral_factor = 0.5d0 * (HO_b**3) / ((2.0d0 + alpha_DD)**1.5d0)
integral_factor = 4.0d0 * pi * t3_DD_CONST * integral_factor
X0M1 = 1.0d0 - x0_DD_FACTOR

!! Get the rad-ang function from the current densities for rearrangement
if (EVAL_REARRANGEMENT) then
  call calculate_common_rearrang_bulkFields
endif
if (PRNT_) then
  open(555, file='BulkHF_elements_Exch.gut')
  open(556, file='BulkHF_elements_Dire.gut')
  open(557, file='BulkPA_elements_Exch.gut')
  open(558, file='BulkPA_comp_integrals.gut')

  write(555, fmt='(2A,2I4)')"[auxHf Exc] a   c  ir  ia   %% ", &
    "real(pp)  imagn(pp), nn, pn    I_r,I_ang=", R_PRINT, ANG_PRINT
  write(556, fmt='(2A,2I4)')"[auxHf Dir] a   c  ir  ia  ms   %% ", &
    "real(pp)  imagn(pp), nn, pn, Rearrange  I_r,I_ang=", R_PRINT, ANG_PRINT
  write(557, fmt='(2A,2I4)')"[aux  Pair] a   c  ir  ia  ms   %% ", &
    "real(pp)  imagn(pp), nn, pn    I_r,I_ang=", R_PRINT, ANG_PRINT
  write(558, fmt='(A)') "[pair integrals] a  c  ms  %%  I_real(pp)  nn   pn  np"
endif

if (EVAL_CUTOFF .AND. (CUTOFF_MODE.NE.1)) then
  call cutoff_by_kappa_matrix(HOsp_dim)
  call reeval_pairing_fields_after_cutoff(.TRUE., hspLR, deltaLR, deltaRL, &
                                          deltaLR_DD, deltaRL_DD, ndim)
endif

do a = 1, spO2
  !! HF field
  a_sh = HOsp_sh(a)

  do c = a, spO2
    c_sh = HOsp_sh(c)

    int_hf = zzero
    int_pa = zzero
    int_rea= zzero

    int_test_PE = zzero

    do i_r = 1, r_dim
      rad_ac = weight_R(i_r) * radial_2b_sho_memo(a_sh, c_sh, i_r)
      rad_ac = rad_ac * dexp((2.0d0+alpha_DD) * (r(i_r)/HO_b)**2)
      do i_ang = 1, angular_dim
        auxHfD = zzero
        !! DIRECT terms for the HF field
        sumD_ang  = AngFunctDUAL_HF(1,a,c,i_ang) + AngFunctDUAL_HF(4,a,c,i_ang)
        auxHfD(1) = dens_pnt(5,i_r,i_ang) - (x0_DD_FACTOR*dens_pnt(1,i_r,i_ang))
        auxHfD(2) = dens_pnt(5,i_r,i_ang) - (x0_DD_FACTOR*dens_pnt(2,i_r,i_ang))
        if (CALCULATE_DD_PN_HF) then
          auxHfD(3) = -x0_DD_FACTOR * dens_pnt(3,i_r,i_ang)
          auxHfD(4) = -x0_DD_FACTOR * dens_pnt(4,i_r,i_ang)
          endif
        do Tac = 1, 4
          auxHfD(Tac) = sumD_ang * auxHfD(Tac)
        end do

        !! EXCHANGE terms for the HF fields
        auxHfE = zzero
        aux_PE = zzero
        aux = zzero
        do ms = 1, 4
          aux(ms) = BulkHF(1,ms, i_r, i_ang) !pp
          aux(ms) = aux(ms) - (x0_DD_FACTOR * BulkHF(5,ms, i_r, i_ang)) !tot
          aux(ms) = aux(ms) * AngFunctDUAL_HF(ms, a,c, i_ang)
          auxHfE(1)  = auxHfE(1) + aux(ms)

          aux(ms) = BulkHF(2,ms, i_r, i_ang) !nn
          aux(ms) = aux(ms) - (x0_DD_FACTOR * BulkHF(5,ms, i_r, i_ang))
          aux(ms) = aux(ms) * AngFunctDUAL_HF(ms, a,c, i_ang)
          auxHfE(2)  = auxHfE(2) + aux(ms)
            !pn np part
          if (CALCULATE_DD_PN_PA) then
          aux(ms) = AngFunctDUAL_HF(ms,a,c, i_ang) * BulkHF(3,ms, i_r,i_ang) !pn
          auxHfE(3)  = auxHfE(3) + aux(ms)
          aux(ms) = AngFunctDUAL_HF(ms,a,c, i_ang) * BulkHF(4,ms, i_r,i_ang) !np
          auxHfE(4)  = auxHfE(4) + aux(ms)
          endif

          !! NOTE: Angular 1, 2 functions are defined with direct form of ms,ms'
          if (hasX0M1) then
            aux(ms) = AngFunctDUAL_P2(ms,a,c,i_ang) * BulkP1(1,ms,i_r,i_ang) !pp
            aux_PE(1) = aux_PE(1)  + (X0M1*aux(ms))
            aux(ms) = AngFunctDUAL_P2(ms,a,c,i_ang) * BulkP1(2,ms,i_r,i_ang) !nn
            aux_PE(2) = aux_PE(2)  + (X0M1*aux(ms))
          endif

          !! pn np part, x0 dependence was calculated in BulkP1_**
          if (CALCULATE_DD_PN_PA) then
          aux(ms) = AngFunctDUAL_P2(ms,a,c, i_ang) * BulkP1(3,ms, i_r,i_ang) !pn
          aux_PE(3)  = aux_PE(3)  + aux(ms)
          aux(ms) = AngFunctDUAL_P2(ms,a,c, i_ang) * BulkP1(4,ms, i_r,i_ang) !np
          aux_PE(4)  = aux_PE(4)  + aux(ms)
          endif
          !! TEST -----------------------------------------------------------
          if (PRNT_)then
            !! Integrals for the Pairing Independently
            do Tac = 1, 4
            aux(ms) = AngFunctDUAL_P2(ms,a,c,i_ang) * BulkP1(Tac,ms,i_r,i_ang)
            if ((Tac .EQ. 1).OR.(Tac .EQ. 2)) aux(ms) = aux(ms) * X0M1
            int_test_PE(Tac,ms,a,c) = int_test_PE(Tac,ms,a,c) + &
              (weight_LEB(i_ang) * rad_ac * dens_alpha(i_r,i_ang) * aux(ms))
            end do

            if ((i_ang == ANG_PRINT)) then ! (i_r == R_PRINT).AND.
              write(555, fmt='(5I4,A,4(F20.15,SP,F20.15,"j"))') &
                a,c,i_r,i_ang,ms,"%%", auxHfE(1),auxHfE(2),auxHfE(3),auxHfE(4)
              write(557, fmt='(5I4,A,4(F20.15,SP,F20.15,"j"))') &
                a,c,i_r,i_ang,ms,"%%", aux_PE(1),aux_PE(2),aux_PE(3),aux_PE(4)
            endif
          endif
          !! TEST -----------------------------------------------------------

        enddo ! ms loop

        if (has_HEIS_MAJO_TERMS) then
          if (a.eq.1 .AND. c.eq.1 .and. i_r.eq.1 .AND. i_ang.EQ.1) &
            print *, "evaluating calcu HM fields"
          call calculate_fields_DD_HM(a,c, i_r,i_ang, auxHfD, auxHfE, aux_PE)
        endif

        !! EXCHANGE Sum terms and add to the global (r,ang) value to integrate
        do Tac =  1, 4
          aux_hf(Tac)   = weight_LEB(i_ang) * rad_ac * dens_alpha(i_r,i_ang)
          aux_hf(Tac)   = (auxHfD(Tac) - auxHfE(Tac)) * aux_hf(Tac)
          int_hf(Tac)   = int_hf(Tac) + aux_hf(Tac)

          aux_pair(Tac) = weight_LEB(i_ang) * rad_ac * dens_alpha(i_r,i_ang)
          aux_pair(Tac) = aux_PE(Tac) * aux_pair(Tac)
          int_pa(Tac)   = int_pa(Tac) + aux_pair(Tac)
        enddo

        auxRea = zzero
        if (EVAL_REARRANGEMENT) then
          auxRea  = REACommonFields(i_r,i_ang) * dens_alpm1(i_r,i_ang)
          auxRea  = auxRea * rea_common_RadAng(a,c, i_r, i_ang)
          auxRea  = auxRea * dexp( (2.0d0+alpha_DD) * (r(i_r)/HO_b)**2)
          int_rea = int_rea + (auxRea * weight_R(i_r) * weight_LEB(i_ang))
        endif
        ! rearrange for pn and np are the same (pn/np are Zero)

        if ((PRNT_).AND.(i_ang == ANG_PRINT)) then ! .AND.(i_r == R_PRINT)
        write(556, fmt="(4I4,A,8F20.15)") a,c,i_r,i_ang, "%%", &
          dreal(auxHfD(1)),dimag(auxHfD(1)),dreal(auxHfD(2)),dimag(auxHfD(2)),&
          dreal(auxHfD(3)),dimag(auxHfD(3)), dreal(auxRea),dimag(auxRea)
          endif

      enddo ! loop ang
    enddo !loop r

    do Tac = 1, 4
      int_hf(Tac) = int_hf(Tac) * integral_factor
      int_pa(Tac) = int_pa(Tac) * integral_factor
    end do
    int_rea = int_rea * 0.25d+0 * integral_factor * alpha_DD

    call complete_DD_fields(int_hf, int_pa, int_rea, gammaLR, deltaLR,deltaRL,&
                            hspLR, gammaLR_DD, deltaLR_DD, deltaRL_DD, &
                            a, c, spO2, ndim)

    if (PRNT_) then
      do ms = 1, 4
        int_test_PE(Tac,ms,a,c) = int_test_PE(Tac,ms,a,c) * integral_factor
        write(558, fmt="(3I4,A)", advance='no') a,c,ms, " %%"
        do Tac = 1, 4
          write(558,fmt="(F20.15)",advance='no') dreal(int_test_PE(Tac,ms,a,c))
        enddo
        write(558, fmt='(A)') ""
      enddo

      if (a.NE.c) then
        do ms = 1, 4
        write(558, fmt="(3I4,A)", advance='no') c,a,ms, " %%"
        do Tac = 1, 4
          write(558,fmt="(F20.15)",advance='no') -dreal(int_test_PE(Tac,ms,a,c))
        enddo
        write(558, fmt='(A)') ""
        enddo
      end if
    end if

    if (.not.DOING_PROJECTION) then
    if ((dabs(imag(int_hf(1)))> 1.0d-9).OR.(dabs(imag(int_hf(2)))> 1.0d-9).OR.&
        (dabs(imag(int_hf(3)))> 1.0d-9).OR.(dabs(imag(int_hf(4)))> 1.0d-9))then
        print "(A,2I4,A,4F18.12)", "[WARNING] Imaginary part at Gamma DD (", &
            a,c, ") = ", &
            imag(int_hf(1)), imag(int_hf(2)), imag(int_hf(3)), imag(int_hf(4))
    endif
    if ((dabs(imag(int_pa(1)))> 1.0d-9).OR.(dabs(imag(int_pa(2)))> 1.0d-9).OR.&
        (dabs(imag(int_pa(3)))> 1.0d-9).OR.(dabs(imag(int_pa(4)))> 1.0d-9))then
        print "(A,2I4,A,4F18.12)", "[WARNING] Imaginary part at Delta DD (", &
            a,c, ") = ", &
            imag(int_pa(1)), imag(int_pa(2)), imag(int_pa(3)), imag(int_pa(4))
    endif
    if ((dabs(imag(int_rea))> 1.0d-9)) then
        print "(A,2I4,A,F18.12)", "[WARNING] Imaginary part Rearrang.F DD (", &
            a,c, ") = ", imag(int_rea)
    endif
    endif

  enddo
enddo

if (EVAL_CUTOFF .AND. (CUTOFF_MODE.NE.2)) then
  call cutoff_by_hspfield_matrix(hspLR, deltaLR, deltaRL, &
                                 deltaLR_DD, deltaRL_DD,ndim)
endif

!! save the last EDF HFB of the DD term
last_HFB_energy = zero
if (evalQuasiParticleVSpace .OR. exportValSpace) then

call zgemm('n','n',ndim,ndim,ndim,zone,hamil_H1,ndim,rhoLR,ndim,zzero, A1,ndim)
call zgemm('n','n',ndim,ndim,ndim,zone,gammaLR ,ndim,rhoLR,ndim,zzero, A2,ndim)
call zgemm('n','n',ndim,ndim,ndim,zone,deltaLR,ndim,kappaLR,ndim,zzero,A3,ndim)
do a=1, ndim
  last_HFB_energy =  last_HFB_energy + (A1(a,a) + 0.5d0 * (A2(a,a) - A3(a,a)))
end do
A1 = zzero
call zgemm('n','n',ndim,ndim,ndim,zone,rearrang_field,ndim,rhoLR,ndim,zzero,&
           A1,ndim)
REARRANGEMENT_ENERGY = 0.0d00
do a = 1, ndim
  REARRANGEMENT_ENERGY = REARRANGEMENT_ENERGY + dreal(A1(a, a))
end do
endif

if (PRNT_) then
  close(555)
  close(556)
  close(557)
  close(558)
endif

!! do the trace-test and printing of fields each 10 steps
if (PRNT_.OR.doTraceTest_) then
  t_rea_sum = zzero
  t_gam_sum = zzero

  if (PRNT_) then
    open (620, file='fields_matrix.gut')
    open (621, file='fields_matrix_imag.gut')
    write(620,fmt='(A,A,A)') "  i   j        gammaLR       gammaLR_DD    ",&
      "    DeltaLR        DeltaLR_DD       DeltaRL        DeltaRL_DD     ",&
      "ReaField_DD         hspLR       hspLR_dd     hspLR_rea"
    write(621,fmt='(A,A,A)') "  i   j        gammaLR       gammaLR_DD    ",&
      "    DeltaLR        DeltaLR_DD       DeltaRL        DeltaRL_DD     ",&
      "ReaField_DD         hspLR       hspLR_dd     hspLR_rea"
  endif

  do i = 1, HOsp_dim
    do j = 1, HOsp_dim
      if (PRNT_) then
        write(620, fmt='(2I4,8F16.9,2F13.6)') i, j, dreal(gam1(i,j)),&
          dreal(gammaLR_DD(i,j)), dreal(dltLR1(i,j)), dreal(deltaLR_DD(i,j)),&
          dreal(dltRL1(i,j)),dreal(deltaRL_DD(i,j)),dreal(rearrang_field(i,j)),&
          dreal(hsp1(i,j)), dreal(hsp1(i,j)+gammaLR_DD(i,j)), dreal(hspLR(i,j))
        write(621, fmt='(2I4,8F16.9,2F13.6)') i, j, dimag(gam1(i,j)),&
          dimag(gammaLR_DD(i,j)), dimag(dltLR1(i,j)), dimag(deltaLR_DD(i,j)),&
          dimag(dltRL1(i,j)),dimag(deltaRL_DD(i,j)),dimag(rearrang_field(i,j)),&
          dimag(hsp1(i,j)), dimag(hsp1(i,j)+gammaLR_DD(i,j)), dimag(hspLR(i,j))
      endif
      !! TEST: Tr(dens * Gamma) = Tr(dens * Rearrange) * alpha
      !if ((i > spO2).AND.(j > spO2)) cycle
      if (doTraceTest_) then
        t_rea_sum = t_rea_sum + (rhoLR(i,j) * rearrang_field(j,i))
        t_gam_sum = t_gam_sum + (rhoLR(i,j) * gammaLR_DD(j,i))
        t_gam_sum = t_gam_sum - 0.5d00 * (kappaRL(i,j) * deltaLR_DD(j,i))
      endif
    enddo
  enddo

  if (doTraceTest_) then
    !!
    !t_rea_sum = (0.25d00 / alpha_DD) * t_rea_sum
    del_aux = (t_gam_sum / t_rea_sum) - (alpha_DD / 2.0d00)
    if (abs(dreal(del_aux)) > 1.0D-8) then
      print '(A,2F15.6,A,F15.6)', "[Warning] Tr(rho*Gam)/= 2a*Tr(rho*Rea) =",&
        dreal(t_gam_sum),dreal(t_rea_sum), &
        ":: rhoGam/rhoRea ::", dreal(t_gam_sum / t_rea_sum)
    endif !!! ********************************************************* DELETE
  endif

  if (PRNT_) close(620)
  if (PRNT_) close(621)
endif

end subroutine calculate_fields_DD

!-----------------------------------------------------------------------------!
! subroutine to transform the final fields into the canonical basis and evaluate
! the sp energies valid to exclude the cutoff energy.
!-----------------------------------------------------------------------------!
subroutine cutoff_by_kappa_matrix(ndim)

integer, intent(in) :: ndim
!complex(r64), dimension(ndim,ndim) :: gammaLR, hspLR
!!! The density fields are added to the calculated with the standard hamiltonian
!!! This array variables are local
!complex(r64), dimension(ndim,ndim) :: gammaLR_DD

!complex(r64), dimension(ndim,ndim) :: gammaLR0, hspLR0
!real(r64), dimension(ndim,ndim)    :: gammaRR_DD_co

complex(r64), dimension(ndim,ndim)  :: rho0LR, kappa0LR, kappa0RL
real(r64), dimension(ndim/2,ndim/2) :: D0_pp, D0_nn, D0_pn
real(r64), dimension(ndim/2,ndim/2) :: rhoc_pp, rhoc_nn, rhoc_pn
real(r64), dimension(ndim/2,ndim/2) :: kapc_pp, kapc_nn, kapc_pn
!real(r64), dimension(ndim/2,ndim/2) :: hspc_pp, hspc_nn, hspc_pn
real(r64), dimension(ndim/2,ndim/2) :: bogo_U0_2,bogo_V0_2
complex(r64), dimension(ndim/2,ndim/2) :: bogo_zU0c_2,bogo_zV0c_2,bogo_zD0_2

real(r64),    dimension(ndim,ndim)  :: D0, rhoc, kapc, A1, A2!, Gamc, Delc, hspc, &
!                                       hspRR

real(r64), dimension(ndim/2,ndim/2) :: D02, rhoc2, kapc2,  A12, A22!, &
!                                       Gamc2, Delc2, hspc2
integer, dimension(2) :: k_min, k_max
logical, dimension(2) :: max_ach
real(r64) :: ovac0, e_fermi, VAL_T, aux
integer   :: i, j, k ,l, zn_indx, nocc0,nemp0, T, it, jt, n1o2

n1o2 = ndim / 2
!do i = 1, ndim
!  do j = 1, ndim
!    hspLR0  (i, j) = hspLR  (i, j) - gammaLR_DD(i, j) - rearrang_field(i, j)
!  end do
!end do

!gammaRR_DD_co = real(gammaLR_DD)
!hspRR         = real(hspLR)

!call calculate_fields_diag(rho0LR, kappa0LR, field_gammaLR, field_hspLR, &
!                           field_deltaLR, field_deltaRL, ndim)
!field_hspRR   = real(field_hspLR) + field_gammaRR_DD + field_rearrRR_DD
!call dsyev('v','u',ndim,field_hspRR,ndim,eigen_H11,work,3*ndim-1,info_H11)
!
!!!! =====================================================================
!!! CALCULATE THE CANONICAL BASIS
!call construct_canonical_basis(bogo_U0,bogo_V0,bogo_zU0c,bogo_zV0c,bogo_zD0, &
!                               ovac0,nocc0,nemp0,ndim)
!D0 = real(bogo_zD0)
!
!call dgemm('t','n',ndim,ndim,ndim,one,D0,ndim,hspRR,ndim,zero,A1,ndim)
!call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,D0,ndim,zero,hspc,ndim)

do T = 1, 3 ! pp, nn, pn
  it = 0
  jt = 0
  if ((T .EQ. 2) .OR. (T .EQ. 3)) jt = n1o2
  if  (T .EQ. 2) it = n1o2

  do i = 1, n1o2
    do j = 1, n1o2
      bogo_U0_2(i,j) = bogo_U0(i+it,j+jt)
      bogo_V0_2(i,j) = bogo_V0(i+it,j+jt)
    enddo
  enddo
  call construct_canonical_basis(bogo_U0_2,bogo_V0_2,bogo_zU0c_2,bogo_zV0c_2,&
                                 bogo_zD0_2, ovac0,nocc0,nemp0,n1o2)
  do i = 1, n1o2
    do j = 1, n1o2
      if (T.EQ.1) D0_pp(i,j) = real(bogo_zD0_2(i,j))
      if (T.EQ.2) D0_nn(i,j) = real(bogo_zD0_2(i,j))
      if (T.EQ.3) D0_pn(i,j) = real(bogo_zD0_2(i,j))
      D02(i,j) = real(bogo_zD0_2(i,j))

      if (T.EQ.1) D0(i+it,j+jt) = D0_pp(i,j)
      if (T.EQ.2) D0(i+it,j+jt) = D0_nn(i,j)
      if (T.EQ.3) then
        D0(i+it,j+jt) = D0_pn(i,j)
        D0(i+jt,j+it) = D0_pn(i,j)
      endif

      rhoc2(i,j) = dens_rhoRR   (i+it,j+jt)
      kapc2(i,j) = dens_kappaRR (i+it,j+jt)
!      hspc2(i,j) = hspRR        (i+it,j+jt)
    end do
  end do

  call dgemm('t','n',n1o2,n1o2,n1o2,one,D02,n1o2,rhoc2,n1o2,zero,A12,n1o2)
  call dgemm('n','n',n1o2,n1o2,n1o2,one,A12,n1o2,D02,n1o2,zero,rhoc2,n1o2)

  call dgemm('t','n',n1o2,n1o2,n1o2,one,D02,n1o2,kapc2,n1o2,zero,A12,n1o2)
  call dgemm('n','n',n1o2,n1o2,n1o2,one,A12,n1o2,D02,n1o2,zero,kapc2,n1o2)

!  call dgemm('t','n',n1o2,n1o2,n1o2,one,D02,n1o2,hspc2,n1o2,zero,A12,n1o2)
!  call dgemm('n','n',n1o2,n1o2,n1o2,one,A12,n1o2,D02,n1o2,zero,hspc2,n1o2)

  do i = 1, n1o2
    do j = 1, n1o2
      rhoc(i+it,j+jt) = rhoc2(i,j)
      kapc(i+it,j+jt) = kapc2(i,j)
!      hspc(i+it,j+jt) = hspc2(i,j)
      if (T.EQ.1) D0_pp(i,j) = D02(i,j)
      if (T.EQ.2) D0_nn(i,j) = D02(i,j)
      if (T.EQ.3) D0_pn(i,j) = D02(i,j)

      if (T.EQ.3) kapc_pn(i,j) = kapc2(i,j) !back up kapc_pn
    end do
  end do
end do ! T loop
!print "(A)", "  - calculate cannonical basis [DONE]"

!! Cutoff criteria from Kappa matrix.
k_min = (/0, 0/)
k_max = (/0, 0/)

max_ach = (/.FALSE., .FALSE./)
do k = 1, n1o2, 2
  do T = 1, 2
    if (abs(kapc(k + n1o2*(T-1), k + n1o2*(T-1) + 1)) .GT. CUTOFF_KAPPA) then
      if (k_min(T) .EQ. 0) then
        k_min(T) = k
        endif
      if ((.NOT.max_ach(T)) .AND. (k > 1)) then ! k in 2 steps, (k > 2)
        max_ach(T) = abs(kapc(k + n1o2*(T-1)    , k + n1o2*(T-1) + 1)) .LT. &
                     abs(kapc(k + n1o2*(T-1) - 2, k + n1o2*(T-1) - 1))
      endif
    else
      if ((.NOT.max_ach(T)) .AND. (k_min(T) .NE. 0)) then
        !case: kappa values continue increasing until a sudden drop below CUTOFF
        max_ach(T) = .TRUE.
        endif
      if (max_ach(T) .AND. (k_max(t).EQ.0)) then
        k_max(T) = k - 1
        endif
    endif
  enddo
enddo
!! remove all states but k within kappa min and max, undo the canonical transf.
if ((maxval(k_min) .GT. minval(k_max)) .OR. (minval(k_max) == 0) &
    .OR. ((k_max(1) < n1o2 / 3).OR.(k_max(2) < n1o2 / 3))) then
  k_min = (/n1o2 - 2*int(nucleus_Z), n1o2 - 2*int(nucleus_N)/)
  k_max = (/n1o2, n1o2/)
endif
print "(A,4I5,A,2L3)", " > Cutoff-kappa (p,n)(min, max):", k_min(1), k_min(2),&
            k_max(1), k_max(2), "  // max achieved=", max_ach(1), max_ach(2)

do i = 1, n1o2
  if      (i .LT. k_min(1)) then
    do j = 1, n1o2
      kapc_pn(i, j) = zero
    enddo
  else if (i .GT. k_max(1)) then
    do j = 1, n1o2
      kapc_pn(i, j) = zero
    enddo
  else ! in k_valid range
    do j = 1, max(1, k_min(2) - 1)
      kapc_pn(i, j) = zero
    enddo
    do j = min(k_max(2) + 1, n1o2), n1o2
      kapc_pn(i, j) = zero
    enddo
  end if
enddo

!! else: cannot apply the cutoff and kapc_pn remains as the initial state

open(333, file='_cannonicalFields_rhokapa.gut')
do i = 1, ndim
  do j = 1, ndim
    if (i .GT. n1o2) then
      if (j .GT. n1o2) then
        aux =  kapc(i-n1o2,j-n1o2)
      else
        aux =  kapc_pn(i-n1o2,j)
      endif
    else
      if (j .GT. n1o2) then
        aux = -kapc_pn(j-n1o2,i)
      else
        aux =  kapc(i,j)
      endif
    endif
    write(333,fmt='(2I5,4F15.5)') i, j , D0(i,j), rhoc(i,j),kapc(i,j),aux
      !,hspc(i,j)
  end do
end do
close(333)

!! Last D02 was pn transformation
call dgemm('n','n',n1o2,n1o2,n1o2,one,D0_pn,n1o2,kapc_pn,n1o2,zero,A12,n1o2)
call dgemm('n','t',n1o2,n1o2,n1o2,one,A12,n1o2,D0_pn,n1o2,zero,kapc2,n1o2)
do i = 1, n1o2
  do j = 1, n1o2
    ! no cutoff in pp-nn, copy directly
    kappaRL(i, j) = dens_kappaRR(i, j)
    kappaLR(i, j) = dens_kappaRR(i, j)
    kappaRL(i+n1o2, j+n1o2) = dens_kappaRR(i+n1o2, j+n1o2)
    kappaLR(i+n1o2, j+n1o2) = dens_kappaRR(i+n1o2, j+n1o2)
    ! pn cutoff
    kappaRL(i+n1o2, j) = -kapc2(j, i)
    kappaLR(i+n1o2, j) = -kapc2(j, i)
    kappaRL(i, j+n1o2) = kapc2(i, j)
    kappaLR(i, j+n1o2) = kapc2(i, j)
  end do
end do

!print "(A,I5)", " -cutoff subroutine DONE, iter =", iteration
end subroutine cutoff_by_kappa_matrix

subroutine cutoff_by_hspfield_matrix(hspLR, deltaLR, deltaRL, &
                                     deltaLR_DD, deltaRL_DD, ndim)

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim)  :: hspLR, deltaLR, deltaRL, &
                                       deltaLR_DD, deltaRL_DD
!! The density fields are added to the calculated with the standard hamiltonian
!! This array variables are local

complex(r64), dimension(ndim,ndim)  :: hspLR0, deltaLR0

complex(r64), dimension(ndim,ndim)  :: rho0LR, kappa0LR, kappa0RL
real(r64), dimension(ndim/2,ndim/2) :: D0_pp, D0_nn, D0_pn
real(r64), dimension(ndim/2,ndim/2) :: rhoc_pp, rhoc_nn, rhoc_pn
real(r64), dimension(ndim/2,ndim/2) :: kapc_pp, kapc_nn, kapc_pn

real(r64), dimension(ndim/2,ndim/2) :: hspc_pp, hspc_nn, hspc_pn
real(r64), dimension(ndim/2,ndim/2) :: bogo_U0_2,bogo_V0_2
complex(r64), dimension(ndim/2,ndim/2) :: bogo_zU0c_2,bogo_zV0c_2,bogo_zD0_2

real(r64),    dimension(ndim,ndim)  :: D0, rhoc, kapc, A1, A2, hspc, hspRR
real(r64), dimension(ndim/2,ndim/2) :: D02, rhoc2, kapc2,  A12, A22, hspc2

logical,   dimension(ndim) :: excluded_qp_indx
real(r64), dimension(ndim) :: ener_qp
integer, dimension(2) :: k_min, k_max
logical, dimension(2) :: max_ach
real(r64) :: ovac0, e_fermi, VAL_T, aux
integer   :: i, j, k ,l, zn_indx, nocc0,nemp0, T, it, jt, n1o2

n1o2 = ndim / 2
do i = 1, ndim
  do j = 1, ndim
    ! - gammaLR_DD(i, j) !! no, the DD gamma field is not modified with the C.O.
    hspLR0  (i, j) = hspLR  (i, j) - rearrang_field(i, j)
    deltaLR0(i, j) = deltaLR(i, j) - deltaLR_DD(i, j)
  end do
end do

hspRR = real(hspLR)

!!! Computes the expectation values of the operators
!do m = 1, constraint_dim
!  if ( m <= (constraint_dim - constraint_pair) ) then
!    call calculate_expectval_obo(dens_rhoRR,constraint_HO(1+(m-1)*ndim2), &
!                                 expvalp,expvaln,ndim)
!    expecval(m) =  expvalp + expvaln
!  else
!    call calculate_expectval_pair(dens_kappaRR,constraint_HO(1+(m-1)*ndim2), &
!                                  expecval(m),ndim)
!  endif
!enddo

do T = 1, 3 ! pp, nn, pn
  it = 0
  jt = 0
  if ((T .EQ. 2) .OR. (T .EQ. 3)) jt = n1o2
  if  (T .EQ. 2) it = n1o2

  do i = 1, n1o2
    do j = 1, n1o2
      bogo_U0_2(i,j) = bogo_U0(i+it,j+jt)
      bogo_V0_2(i,j) = bogo_V0(i+it,j+jt)
    enddo
  enddo
  call construct_canonical_basis(bogo_U0_2,bogo_V0_2,bogo_zU0c_2,bogo_zV0c_2,&
                                 bogo_zD0_2, ovac0,nocc0,nemp0,n1o2)
  do i = 1, n1o2
    do j = 1, n1o2
      if (T.EQ.1) D0_pp(i,j) = real(bogo_zD0_2(i,j))
      if (T.EQ.2) D0_nn(i,j) = real(bogo_zD0_2(i,j))
      if (T.EQ.3) D0_pn(i,j) = real(bogo_zD0_2(i,j))
      D02(i,j) = real(bogo_zD0_2(i,j))

      if (T.EQ.1) D0(i+it,j+jt) = D0_pp(i,j)
      if (T.EQ.2) D0(i+it,j+jt) = D0_nn(i,j)
      if (T.EQ.3) then
        D0(i+it,j+jt) = D0_pn(i,j)
        D0(i+jt,j+it) = D0_pn(i,j)
      endif

      rhoc2(i,j) = dens_rhoRR   (i+it,j+jt)
      kapc2(i,j) = dens_kappaRR (i+it,j+jt)
      hspc2(i,j) = hspRR        (i+it,j+jt)
      if (T .NE. 3) then
        hspc2(i,j) = hspc2(i,j) - lambdaFer_DD(T)
      endif
    enddo
  enddo

  call dgemm('t','n',n1o2,n1o2,n1o2,one,D02,n1o2,rhoc2,n1o2,zero,A12,n1o2)
  call dgemm('n','n',n1o2,n1o2,n1o2,one,A12,n1o2,D02,n1o2,zero,rhoc2,n1o2)

  call dgemm('t','n',n1o2,n1o2,n1o2,one,D02,n1o2,kapc2,n1o2,zero,A12,n1o2)
  call dgemm('n','n',n1o2,n1o2,n1o2,one,A12,n1o2,D02,n1o2,zero,kapc2,n1o2)

  call dgemm('t','n',n1o2,n1o2,n1o2,one,D02,n1o2,hspc2,n1o2,zero,A12,n1o2)
  call dgemm('n','n',n1o2,n1o2,n1o2,one,A12,n1o2,D02,n1o2,zero,hspc2,n1o2)

  do i = 1, n1o2
    do j = 1, n1o2
      rhoc(i+it,j+jt) = rhoc2(i,j)
      kapc(i+it,j+jt) = kapc2(i,j)
      hspc(i+it,j+jt) = hspc2(i,j)
      if (T.EQ.1) D0_pp(i,j) = D02(i,j)
      if (T.EQ.2) D0_nn(i,j) = D02(i,j)
      if (T.EQ.3) D0_pn(i,j) = D02(i,j)

      if (T.EQ.3) kapc_pn(i,j) = kapc2(i,j) !back up kapc_pn
    enddo
  enddo
enddo ! T loop

!! Cutoff criteria from Kappa matrix.
k_min = (/0, 0/)
k_max = (/0, 0/)
max_ach = (/.FALSE., .FALSE./)
do i = 1, ndim
  T = 1
  if (i > n1o2) T = 2

  ener_qp(i) = hspc(i,i)
  if (abs(ener_qp(i)) .GT. CUTOFF_ENERGY_MAX) then
    excluded_qp_indx(i) = .TRUE.
    if (.NOT. max_ach(T)) then
      max_ach(T) = .TRUE.
      k_min(T)   = i - n1o2*(T-1)
    endif
  else
    if ((max_ach(T)) .AND. (k_max(T) .EQ. 0)) then
      k_max(T) = i - n1o2*(T-1)
    endif
  endif
enddo
!!! remove all states but k within kappa min and max, undo the canonical transf.
!if ((maxval(k_min) .GT. minval(k_max)) .OR. (minval(k_max) == 0) &
!    .OR. ((k_max(1) < n1o2 / 3).OR.(k_max(2) < n1o2 / 3))) then
!  k_min = (/n1o2 - 2*int(nucleus_Z), n1o2 - 2*int(nucleus_N)/)
!  k_max = (/n1o2, n1o2/)
!endif
print "(A,4I5,A,2L3)", " > Cutoff-energy (p,n)(min, max):", k_min(1),k_min(2),&
            k_max(1), k_max(2), "  // max achieved=", max_ach(1), max_ach(2)

do i = 1, n1o2
  if      (i .LT. k_min(1)) then
    do j = 1, n1o2
      kapc_pn(i, j) = zero
    enddo
  else if (i .GT. k_max(1)) then
    do j = 1, n1o2
      kapc_pn(i, j) = zero
    enddo
  else ! in k_valid range
    do j = 1, max(1, k_min(2) - 1)
      kapc_pn(i, j) = zero
    enddo
    do j = min(k_max(2) + 1, n1o2), n1o2
      kapc_pn(i, j) = zero
    enddo
  end if
enddo

!! else: cannot apply the cutoff and kapc_pn remains as the initial state
!! Export
open(333, file='_cannonicalFields_coE_rhokapa.gut')
do i = 1, ndim
  do j = 1, ndim
    if (i .GT. n1o2) then
      if (j .GT. n1o2) then
        aux =  kapc(i-n1o2,j-n1o2)
      else
        aux =  kapc_pn(i-n1o2,j)
      endif
    else
      if (j .GT. n1o2) then
        aux = -kapc_pn(j-n1o2,i)
      else
        aux =  kapc(i,j)
      endif
    endif
    write(333,fmt='(2I5,5F15.5)') i, j , D0(i,j), rhoc(i,j),kapc(i,j),aux, &
                                  hspc(i,j)
  end do
end do
close(333)

!! Transform Kappa-pn to the working basis.
call dgemm('n','n',n1o2,n1o2,n1o2,one,D0_pn,n1o2,kapc_pn,n1o2,zero,A12,n1o2)
call dgemm('n','t',n1o2,n1o2,n1o2,one,A12,n1o2,D0_pn,n1o2,zero,kapc2,n1o2)
do i = 1, n1o2
  do j = 1, n1o2
    ! no cutoff in pp-nn, copy directly
    kappaRL(i, j) = dens_kappaRR(i, j)
    kappaLR(i, j) = dens_kappaRR(i, j)
    kappaRL(i+n1o2, j+n1o2) = dens_kappaRR(i+n1o2, j+n1o2)
    kappaLR(i+n1o2, j+n1o2) = dens_kappaRR(i+n1o2, j+n1o2)
    ! pn cutoff
    kappaRL(i+n1o2, j) = -kapc2(j, i)
    kappaLR(i+n1o2, j) = -kapc2(j, i)
    kappaRL(i, j+n1o2) = kapc2(i, j)
    kappaLR(i, j+n1o2) = kapc2(i, j)
  end do
end do

!!! - Copy the non-DD modified from initial kappa
do i = 1, ndim
  do j = 1, ndim
    hspLR  (i, j) = hspLR0  (i, j)
    deltaLR(i, j) = deltaLR0(i, j)
    deltaRL(i, j) = deltaLR0(i, j)
  end do
end do
call reeval_pairing_fields_after_cutoff(.FALSE., hspLR, deltaLR, deltaRL, &
                                        deltaLR_DD, deltaRL_DD, ndim)

!print "(A,I5)", " -cutoff subroutine DONE, iter =", iteration
end subroutine cutoff_by_hspfield_matrix



subroutine reeval_pairing_fields_after_cutoff(only_bulk_fields, hspLR, &
                                              deltaLR,deltaRL, &
                                              deltaLR_DD, deltaRL_DD, ndim)

logical, intent(in) :: only_bulk_fields
integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim) :: hspLR, deltaLR, deltaRL

real(r64),    dimension(ndim,ndim) :: BU_gammaRR_DD
complex(r64), dimension(ndim,ndim) :: gammaLR, gammaLR_DD,deltaLR_DD,deltaRL_DD

integer :: a,b,a_sh,b_sh,spO2,i_r,i_an, ms, Tab
complex(r64), dimension(4) :: int_hf, int_pa ! all arrays are for (pp, nn, pn, np)
complex(r64), dimension(4) :: aux, aux_PE, aux_pair
complex(r64) :: int_rea, auxRea
real(r64)    :: rad_ab, X0M1, integral_factor

spO2   = HOsp_dim / 2

BulkP1 = zzero
BulkP2 = zzero
!!-----------------------
do i_r = 1, r_dim
  do i_an = 1, angular_dim

    do a = 1, spO2
       a_sh = HOsp_sh(a)
       do b = a, spO2         !!!!!     BENCH REQUIRES B=1      !!!!
         b_sh = HOsp_sh(b)
         call compute_bulkDens4_pair(a, b, a_sh, b_sh, i_r, i_an, zone)
      enddo ! do b
    enddo   ! do a

    do ms = 1, 4
      BulkP1(5,ms,i_r,i_an) = BulkP1(1,ms,i_r,i_an) + BulkP1(2,ms,i_r,i_an)
      BulkP2(5,ms,i_r,i_an) = BulkP2(1,ms,i_r,i_an) + BulkP2(2,ms,i_r,i_an)
    enddo
  enddo
enddo

call calculate_common_rearrang_bulkFields

if (only_bulk_fields) return  !! ========================================

!! save BU for Gamma RR (ommited)
do a = 1, ndim
  do b = 1, ndim
    BU_gammaRR_DD(a, b) = field_gammaRR_DD(a, b)
  enddo
enddo
gammaLR    = zzero
gammaLR_DD = zzero
deltaLR_DD = zzero
deltaRL_DD = zzero

!!! Evaluate the pairing fields ----------------------
!! Note :: Remember that radial functions already have the factor 1/b**3
integral_factor = 0.5d0 * (HO_b**3) / ((2.0d0 + alpha_DD)**1.5d0)
integral_factor = 4.0d0 * pi * t3_DD_CONST * integral_factor
X0M1 = 1.0d0 - x0_DD_FACTOR

do a = 1, spO2
  a_sh = HOsp_sh(a)
  do b = a, spO2
    b_sh = HOsp_sh(b)

    int_pa = zzero
    int_rea= zzero

    do i_r = 1, r_dim
      rad_ab = weight_R(i_r) * radial_2b_sho_memo(a_sh, b_sh, i_r)
      rad_ab = rad_ab * dexp((2.0d0+alpha_DD) * (r(i_r)/HO_b)**2)
      do i_an = 1, angular_dim

        !! EXCHANGE terms for the HF fields
        aux_PE = zzero
        aux = zzero
        do ms = 1, 4
          !! NOTE: Angular 1, 2 functions are defined with direct form of ms,ms'
          if (hasX0M1) then
            aux(ms) = AngFunctDUAL_P2(ms,a,b,i_an) * BulkP1(1,ms,i_r,i_an) !pp
            aux_PE(1) = aux_PE(1)  + (X0M1*aux(ms))
            aux(ms) = AngFunctDUAL_P2(ms,a,b,i_an) * BulkP1(2,ms,i_r,i_an) !nn
            aux_PE(2) = aux_PE(2)  + (X0M1*aux(ms))
          endif
          !! pn np part, x0 dependence was calculated in BulkP1_**
          if (CALCULATE_DD_PN_HF) then
          aux(ms) = AngFunctDUAL_P2(ms,a,b, i_an) * BulkP1(3,ms, i_r,i_an) !pn
          aux_PE(3)  = aux_PE(3)  + aux(ms)
          aux(ms) = AngFunctDUAL_P2(ms,a,b, i_an) * BulkP1(4,ms, i_r,i_an) !np
          aux_PE(4)  = aux_PE(4)  + aux(ms)
          endif
        enddo ! ms loop

        !! EXCHANGE Sum terms and add to the global (r,ang) value to integrate
        do Tab =  1, 4
          aux_pair(Tab) = weight_LEB(i_an) * rad_ab * dens_alpha(i_r,i_an)
          aux_pair(Tab) = aux_PE(Tab) * aux_pair(Tab)
          int_pa  (Tab) = int_pa(Tab) + aux_pair(Tab)
        enddo

        auxRea = zzero
        if (EVAL_REARRANGEMENT) then
          auxRea  = REACommonFields(i_r,i_an) * dens_alpm1(i_r,i_an)
          auxRea  = auxRea * rea_common_RadAng(a,b, i_r, i_an)
          auxRea  = auxRea * dexp( (2.0d0+alpha_DD) * (r(i_r)/HO_b)**2)
          int_rea = int_rea + (auxRea * weight_R(i_r) * weight_LEB(i_an))
        endif
        ! rearrange for pn and np are the same (pn/np are Zero)

      enddo ! loop ang
    enddo !loop r

    do Tab = 1, 4
      int_pa(Tab) = int_pa(Tab) * integral_factor
    enddo
    int_rea = int_rea * 0.25d+0 * integral_factor * alpha_DD

    int_hf = zzero
    call complete_DD_fields(int_hf, int_pa, int_rea, gammaLR, deltaLR,deltaRL,&
                            hspLR, gammaLR_DD, deltaLR_DD, deltaRL_DD, &
                            a, b, spO2, ndim)

  enddo
enddo

! recover the gamma BU.
do a = 1, ndim
  do b = 1, ndim
    field_gammaRR_DD(a, b) = BU_gammaRR_DD(a, b)
  enddo
enddo

end subroutine reeval_pairing_fields_after_cutoff

!-----------------------------------------------------------------------------!
! subroutine to print the progress of iter/iter_max as a progress bar         !
!-----------------------------------------------------------------------------!
subroutine progress_bar_iteration(iter, max_iter)
integer, intent(in) :: iter, max_iter
integer :: TOTAL_SPACE=62, BAR_SPACE_LEN, completed_int, k
real    :: completed
character(len=70) :: bar=" [???.??%] |                                       &
                         &           |"
BAR_SPACE_LEN = 62 - 12
completed_int = nint(100.0d0 * iter / max_iter)
if (MOD(completed_int, 5) .NE. 0) return

completed_int = nint(1.0d0 * BAR_SPACE_LEN * iter / max_iter)

write(unit=bar(3:8),fmt="(f6.2)") 100.0d0 * iter / max_iter
do k = 1, completed_int
  bar(12+k:12+k)='*'
end do

print "(A)", bar
!! this writes in a binary ---------------------------------------------------
!write(unit=857639, fmt="(a1,a70)",advance="no") char(50), bar
!! print the progress bar
!if (iter .NE. max_iter) then
!  flush(unit=857639)
!else
!  write(unit=857639,fmt=*)
!end if
end subroutine progress_bar_iteration


!-----------------------------------------------------------------------------!
! This subroutine evaluates the energy associated to the pseudo-rearrangement !
!-----------------------------------------------------------------------------!
subroutine calculate_energy_field_rearrangement(E_core)
real(r64), intent(out) :: E_core
integer :: a, c

!! Do the trace for the energy
E_core = 0.0
do a = 1, HOsp_dim
  do c = 1, HOsp_dim
    E_core = E_core + dreal(rearrang_field(a,c)) * dens_rhoRR(c, a)
  enddo
enddo

E_core = 0.5d0 * E_core  !! the energy should be 1/2 Tr(Gamma * rho)

end subroutine calculate_energy_field_rearrangement

!-----------------------------------------------------------------------------!
! subroutine TESTS FOR THE DENSITY, SPHERICAL HARMONICS AND FUNCTIONS         !
!                                                                             !
!                                                                             !
!                                                                             !
!                                                                             !
!-----------------------------------------------------------------------------!

subroutine test_printDesityKappaWF
integer :: i, j

open(626, file = 'DIMENS_indexes_and_rhoLRkappas.gut')
write(626,'(A)') "// SHELL INDEX (i_sh, ant_label, n, l, 2j)"
do i = 1, HOsh_dim
  write(626, fmt='(5I6)') i, HOsh_na(i), HOsh_n(i), HOsh_l(i), HOsh_2j(i)
end do
write(626,'(A)')"// SINGLE PARTICLE INDEX (i_sp, i_sh, n,l,2j,2m, 2mt, i_sp_TR)"
do i=1, HOsp_dim
  write(626, fmt='(8I4)') i, HOsp_sh(i),&
        HOsp_n(i), HOsp_l(i), HOsp_2j(i), HOsp_2mj(i), HOsp_2mt(i), HOsp_tr(i)
enddo
write(626,'(A)') "// DENSITY MATRIX RHO_LR KAPPA_LR KAPPA_RL "
do i = 1, HOsp_dim
    do j = 1, HOsp_dim
    write(626, '(2I4,6F20.15)') i, j, dreal(rhoLR(i, j)), dimag(rhoLR(i, j)), &
        dreal(kappaLR(i, j)), dimag(kappaLR(i, j)), &
        dreal(kappaRL(i, j)), dimag(kappaRL(i, j))
    enddo
enddo
close(626)

!print *, "[OK] Export Rho LR matrix in file = 'DIMENS_indexes_and_rhoLRkappas.gut'"
end subroutine test_printDesityKappaWF

!------------------------------------------------------------------------------!
!  function two_angular_momentum_index(a_sh,b_sh)                              !
!          Returns an ordering for the total angular mom. values available for !
! the space, from 1 = (1/2 1/2), 2=(1,3), 3=(3,3), 4=(3,1) ...                 !
!   UPDATE: The method is used to distinguish two shell states (m independent),!
! the previous change of variables is left to remember the design for j values !
!    MAX value come from a_sh(max) b_sh(min). Ordering example:                !
!     (a)   ...   ...  ...   ...                                               !
!      3 :   5     6   (7)   ...                                               !
!      2 :   2    (3)   8    ...                                               !
!      1 :  (1)    4    9    ...                                               !
!            1     2    3    (b)                                               !
!------------------------------------------------------------------------------!
function two_shell_states_index(a_sh,b_sh) result (K)
integer, intent(in) :: a_sh, b_sh
integer :: K
integer :: la, lb

la = a_sh !(ja + 1) / 2
lb = b_sh !(jb + 1) / 2

if (lb > la) then ! la=lb case is equivalent in both
  K = ((lb-1)*(lb-1)) + la
else
  K = ((la-1)*(la-1)) + (la - lb) + la ! for completion but we expect la <= lb
endif
return

endfunction two_shell_states_index


!------------------------------------------------------------------------------!
! Integrate the matrix elements of the density  S<a|delta(r)|b> dr3 = delta_ab !
!------------------------------------------------------------------------------!
subroutine test_integrate_bulk_densities

integer :: i_r, i_a, i_t, ms, tt
real(r64)    :: radInt, int_fact, int_const
complex(r64), dimension(4)   :: int_dens  ! [pp, nn, pn, np]
complex(r64), dimension(5,4) :: int_A, int_B1, int_B2 ! [pp (++,+-,-+,--), nn,pn,np, total]
complex(r64), dimension(5) :: sum_A, sum_B1, sum_B2
character(len=28) :: fmt1 , fmt2
character(len=4), dimension(4) :: ms_str
ms_str = [character(len=4) :: " ++ ", " +- ", " -+ ", " -- "]
fmt1 = '(A,4(F17.12,SP,F16.12,"j")))'
fmt2 = '(A,5(F17.12,SP,F16.12,"j")))'

int_dens = zzero
int_A    = zzero
int_B1   = zzero
int_B2   = zzero

int_const =  2 * pi * (HO_b**3) / ((2.0 + alpha_DD)**1.5)

do i_r = 1, r_dim
  !! [TEST] For the density to be integrated in the variable for the m.e.
  radInt = int_const * weight_R(i_r) * exp((r(i_r)/HO_b)**2 * (1.0+alpha_DD))
  do i_a = 1, angular_dim
    int_fact = radInt * weight_LEB(i_a)

    do i_t = 1, 4
      int_dens(i_t) = int_dens(i_t) + (int_fact*dens_pnt(i_t, i_r, i_a))
    enddo
    !! TOTAL bulk densities
    do ms = 1, 4
      do tt = 1, 5
        int_A (tt, ms) = int_A (tt, ms) + (int_fact * BulkHF(tt,ms, i_r, i_a))
        int_B1(tt, ms) = int_B1(tt, ms) + (int_fact * BulkP1(tt,ms, i_r, i_a))
        int_B2(tt, ms) = int_B2(tt, ms) + (int_fact * BulkP2(tt,ms, i_r, i_a))
      enddo
    end do
  enddo
enddo

sum_A  = zzero
sum_B1 = zzero
sum_B2 = zzero

!!! WRITE STUFF --------------------------------------------------------------
open(321, file='test_dens_integrals.gut')
write(321, fmt='(A)') " [TEST] Integrals of Bulk densities over R3 (approx)"
write(321, fmt='(A)') "   RHO"
write(321, fmt='(12x,"real(pp)",10x,"imag(pp)",8x,"real(nn)",10x,&
  "imag(nn)",8x,"real(pn)",10x,"imag(pn)",8x,"real(np)",10x,"imag(np)")')

write(321, fmt=fmt1) "    ", int_dens(1),int_dens(2),int_dens(3),int_dens(4)
write(321, fmt='(A)') ""
write(321, fmt='(A)') "   BULK A(ms)"
write(321, fmt='(12x,"real(pp)",10x,"imag(pp)",8x,"real(nn)",10x,&
  "imag(nn)",8x,"real(pn)",10x,"imag(pn)",8x,"real(np)",10x,"imag(np)",&
  8x,"real(tt)",10x,"imag(tt)")')
do ms = 1, 4
  write(321, fmt=fmt2) ms_str(ms), &
    int_A(1,ms), int_A(2,ms), int_A(3,ms), int_A(4,ms), int_A(5,ms)
  do i_t = 1, 5
    sum_A(i_t) = sum_A(i_t) + int_A(i_t, ms)
  enddo
enddo
write(321, fmt=fmt2) " tt ", sum_A(1), sum_A(2), sum_A(3), sum_A(4), sum_A(5)
write(321, fmt='(A)') ""
write(321, fmt='(A)') "   BULK B1(ms)"
write(321, fmt='(12x,"real(pp)",10x,"imag(pp)",8x,"real(nn)",10x,&
  "imag(nn)",8x,"real(pn)",10x,"imag(pn)",8x,"real(np)",10x,"imag(np)",&
  8x,"real(tt)",10x,"imag(tt)")')
do ms = 1, 4
  write(321, fmt=fmt2) ms_str(ms), &
    int_B1(1,ms), int_B1(2,ms), int_B1(3,ms), int_B1(4,ms), int_B1(5,ms)
  do i_t = 1, 5
    sum_B1(i_t) = sum_B1(i_t) + int_B1(i_t, ms)
  enddo
enddo
write(321, fmt=fmt2) " tt ", sum_B1(1),sum_B1(2),sum_B1(3),sum_B1(4),sum_B1(5)
write(321, fmt='(A)') ""
write(321, fmt='(A)') "   BULK B2(ms)"
write(321, fmt='(12x,"real(pp)",10x,"imag(pp)",8x,"real(nn)",10x,&
  "imag(nn)",8x,"real(pn)",10x,"imag(pn)",8x,"real(np)",10x,"imag(np)",&
  8x,"real(tt)",10x,"imag(tt)")')
do ms = 1, 4
  write(321, fmt=fmt2) ms_str(ms), &
    int_B2(1,ms), int_B2(2,ms), int_B2(3,ms), int_B2(4,ms), int_B2(5,ms)
  do i_t = 1, 5
    sum_B2(i_t) = sum_B2(i_t) + int_B2(i_t, ms)
  enddo
enddo
write(321, fmt=fmt2) " tt ", sum_B2(1),sum_B2(2),sum_B2(3),sum_B2(4),sum_B2(5)

close(321) !!! ---------------------------------------------------------------

end subroutine test_integrate_bulk_densities


END MODULE DensityDep
!==============================================================================!
! End of file                                                                  !
!==============================================================================!

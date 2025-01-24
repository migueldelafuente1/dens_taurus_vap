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
use Hamiltonian
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

real(r64) :: R_MAX    = 7.5d0
integer   :: THE_grid = 1
integer   :: PHI_grid = 1

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

complex(r64), dimension(:,:,:,:), allocatable, save :: AngFunctDUAL_HF ! CGa CGb Y*(a) Y (b)
complex(r64), dimension(:,:,:,:), allocatable, save :: AngFunctDUAL_P1 ! CGa CGb Y (a) Y (b)
complex(r64), dimension(:,:,:,:), allocatable, save :: AngFunctDUAL_P2 ! CGa CGb Y*(a) Y*(b)
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
logical   :: EXPORT_PREA_DD = .TRUE.  ! export the Field-derivation of the rearrange m.e.
real(r64), dimension(:,:,:,:), allocatable  :: radial_1b_diff_memo ! (ish, i_n[-1:1], j_l[-1:1], ir)
complex(r64), dimension(:,:,:), allocatable :: partial_dens        ! (-1,0,1,2:total ,ir,iang)
real(r64), dimension(:,:), allocatable      :: hamil_GradDD_H2_byT ! 2-body grad Dens  (pppp,pnpn,pnnp,nnnn)

real(r64) :: last_HFB_energy

integer, dimension(:), allocatable :: HOsh_ant, HOsp_ant


!! Related to the hamiltonian and Fields
real(r64), dimension(:),   allocatable :: hamil_DD_H2      ! 2-body part
real(r64), dimension(:,:), allocatable :: hamil_DD_H2_byT  ! 2-body part (pppp,pnpn,pnnp,nnnn)
integer(i64) :: hamil_DD_H2dim, hamil_DD_H2dim_all         ! number of 2BME stored
integer(i16), dimension(:), allocatable :: hamil_DD_abcd   ! indices of 2BME
integer(i8) , dimension(:), allocatable :: hamil_DD_trperm ! time reversal permut.
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
if (abs(alpha_DD - (1.0/3.0)) < 1.0D-4) alpha_DD = 0.333333333333333

alpha_DD_frac = 0
if (abs(alpha_DD - (1.0/3.0)) < 1.0D-4) then
  alpha_DD_frac = (/1, 3/)
else if ((abs(alpha_DD - (1.0)) < 1.0D-4) .OR. (abs(alpha_DD) < 1.0D-4)) then
  alpha_DD_frac = (/ int(alpha_DD), 1/)
else if (abs(alpha_DD - (2.0/3.0)) < 1.0D-4) then
  alpha_DD_frac = (/2, 3/)
else if (abs(alpha_DD - (1.0/6.0)) < 1.0D-4) then
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
    deallocate(hamil_H2cpd_DD) ! It wont be used
    print "(A)", "  Hamiltonian cpd deallocated  because it will not be used!"
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
endif
print *, ''
if ((FUNCTIONAL_DENS_MODE .EQ. 2) .AND. (has_HEIS_MAJO_TERMS)) then
  print "[ERROR] Functional Form 2 (Phys.Rev.C 60 064312) and Heis/Majo terms!"
  STOP
end if

hasX0M1 = abs(x0_DD_FACTOR - 1.0d+0) > 1.0d-6
CONST_EDD_M2 = CONST_EDD_M2_ETA / (CONST_EDD_M2_RHO0 ** alpha_DD)

print "(A)", " * Density dependent parameters imported."

print "(A,2L3)", " * [OPTIONs] Calculate DD-pn parts (HF/PA) :", &
                 CALCULATE_DD_PN_HF, CALCULATE_DD_PN_PA

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
    print "(A,2I3)", " > Potential MODE 2: ETA  =", CONST_EDD_M2_ETA
  case(22)
    FUNCTIONAL_DENS_MODE = 2
    CONST_EDD_M2_RHO0 = aux_float
    print "(A,2I3)", " > Potential MODE 2: RHO_0=", CONST_EDD_M2_RHO0


  case(31)
    has_HEIS_MAJO_TERMS  = .TRUE.
    CONST_x0_EXC_HEIS = aux_float
    print "(A,F12.9)", " > Exchange (spin):       Heisenberg=", CONST_x0_EXC_HEIS
  case(32)
    has_HEIS_MAJO_TERMS  = .TRUE.
    CONST_x0_EXC_MAJO = aux_float
    print "(A,F12.9)", " > Exchange (spin-isospin): Majorana=", CONST_x0_EXC_MAJO
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
!        radial = two_sho_radial_functions_bench(a_sh, b_sh, r(i_r))
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
          write(629,fmt='(3D22.13)',advance='no') r(i_r), radial, weight_R(i_r)
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

real(r64) :: radial_part, dens_R, dens_A, rad4Integr, dens_Aa,dens_Ra
complex(r64) :: sum_, integral_dens, sum_test, diff, x
logical :: PRNT_

PRNT_ = (PRINT_GUTS).OR.(.FALSE.)
spO2 = HOsp_dim / 2

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
  print "(A,F13.9,A)", "      *A* ", dreal(integral_dens), "  <dens(r)> approx"
endif

end subroutine calculate_expectval_density


!-----------------------------------------------------------------------------!
! subroutine reconstruct_2body_DD_timerev                                     !
! implements the tReversal index whit their phase, m.e. with m_j < 0          !
!-----------------------------------------------------------------------------!
subroutine reconstruct_2body_DD_timerev

integer :: ia, ib, ic, id, ialloc=0
integer(i8)  :: perm
integer(i64) :: kk

allocate( hamil_DD_trperm(hamil_DD_H2dim), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of time reversal (DD)'

hamil_DD_trperm = 0

do kk = 1, hamil_DD_H2dim
  ia = hamil_DD_abcd(1+4*(kk-1))
  ib = hamil_DD_abcd(2+4*(kk-1))
  ic = hamil_DD_abcd(3+4*(kk-1))
  id = hamil_DD_abcd(4+4*(kk-1))

  if (evalFullSPSpace) then
    perm = 1
  else
    perm = step_reconstruct_2body_timerev(ia,ib, ic, id)
  endif

  hamil_DD_trperm(kk) = perm
enddo

end subroutine reconstruct_2body_DD_timerev


!-----------------------------------------------------------------------------!
! function step_reconstruct_2body_timerev    (private)                        !
! - This is an auxiliary function to emulate the steps to add the permutation !
! to be saved in the hamil_DD_trperm array to be used intermediately for the  !
! individual steps of the rearrangement field calculation.                    !
!-----------------------------------------------------------------------------!
function step_reconstruct_2body_timerev(ia,ib, ic, id) result (perm)
integer, intent(in) :: ia, ib, ic, id
integer             :: perm
integer             :: ta, tb, tc, td, tmp
real(r64) :: phase

perm = 0
phase = (-one)**((HOsp_2j(ia) + HOsp_2j(ib) + HOsp_2j(ic) + HOsp_2j(id))/2)

phase=phase*((-one)**(HOsp_l(ia)+HOsp_l(ib)+HOsp_l(ic)+HOsp_l(id)))
phase=phase*((-one)**((HOsp_2mj(ia)+HOsp_2mj(ib)-HOsp_2mj(ic)-HOsp_2mj(id))/2))

ta = HOsp_tr(ia)
tb = HOsp_tr(ib)
tc = HOsp_tr(ic)
td = HOsp_tr(id)

if ( ta > tb ) then
    tmp = ta
    ta = tb
    tb = tmp
    phase = -phase
    perm = perm + int(1,i8)
endif

if ( tc > td ) then
    tmp = tc
    tc = td
    td = tmp
    phase = -phase
    perm = perm + int(2,i8)
endif

if ( (ta > tc) .or. ((ta == tc) .and. (tb > td)) ) then
    perm = perm + int(4,i8)
endif

if ( phase < 0.d0 ) perm = perm - int(8,i8)

return
end function step_reconstruct_2body_timerev

!-----------------------------------------------------------------------------!
! function matrix_element_v_DD                                                !
!                                                                             !
! Computes density dependent two body matrix elements over the density average!
!    all_isos (logical) Compute 3 combinations p/n instead of the current     !
!                       ta,tb,tc,td of the sp-state.                          !
!                       v_dd_val_Real !! pppp(1), pnpn(2), pnnp(3), nnnn(4)   !
!-----------------------------------------------------------------------------!
function matrix_element_v_DD(a,b, c,d, ALL_ISOS) result (v_dd_val_Real)

integer(i32), intent(in) :: a,b,c,d
logical, intent(in) :: ALL_ISOS     !! compute 3 combinations p/n instead of the
real(r64), dimension(4) :: v_dd_val_Real !! pppp(1), pnpn(2), pnnp(3), nnnn(4)

complex(r64) :: aux, radial, aux_dir, aux_exch, ang_rea, rea,aux_rea
complex(r64), dimension(4) :: v_dd_value
real(r64)    :: delta_dir, delta_exch, angular, integral_factor,TOP, LOW
integer(i32) :: la, ja, ma, a_sh, ta, lb, jb, mb, b_sh, tb, bmax, &
                lc, jc, mc, c_sh, tc, ld, jd, md, d_sh, td,&
                i_r, i_th, i_phi, i_ang, &
                K,M, K2,M2, ind_km, ind_km_q, &
                kk1,kk2, kk1N,kk2N, l1,l2,j1,j2, mj1,mj2, ind_jm_1, ind_jm_2, &
                ind_jm_a, ind_jm_b, ind_jm_c, ind_jm_d, delta_ac_bd, delta_ad_bc
integer :: print_element = 0, HOspO2, skpd
real(r64) :: v_re_1, v_re_2

HOspO2 = HOsp_dim/2

v_dd_value = zzero
v_dd_val_Real = zero

ja = HOsp_2j(a)
ma = HOsp_2mj(a)
la = HOsp_l(a)
a_sh = HOsp_sh(a)
ind_jm_a = angular_momentum_index(ja, ma, .TRUE.)
jb = HOsp_2j(b)
mb = HOsp_2mj(b)
lb = HOsp_l(b)
b_sh = HOsp_sh(b)
ind_jm_b = angular_momentum_index(jb, mb, .TRUE.)
jc = HOsp_2j(c)
mc = HOsp_2mj(c)
lc = HOsp_l(c)
c_sh = HOsp_sh(c)
ind_jm_c = angular_momentum_index(jc, mc, .TRUE.)
jd = HOsp_2j(d)
md = HOsp_2mj(d)
ld = HOsp_l(d)
d_sh = HOsp_sh(d)
ind_jm_d = angular_momentum_index(jd, md, .TRUE.)

delta_dir  = one
delta_exch = one
if (.NOT.ALL_ISOS) then ! compute the
  !else ! evaluate the element isospin directly
  ta = HOsp_2mt(a)
  tb = HOsp_2mt(b)
  tc = HOsp_2mt(c)
  td = HOsp_2mt(d)

  delta_ac_bd = abs((ta + tc) * (tb + td) / 4)  ! a+c = (2, 0)
  delta_ad_bc = abs((ta + td) * (tb + tc) / 4)

  delta_dir  = delta_ac_bd - (x0_DD_FACTOR * delta_ad_bc)
  delta_exch = delta_ad_bc - (x0_DD_FACTOR * delta_ac_bd)

  !! Tested (abs(delta_dir)+abs(delta_exch) > 2.0d-8) & pppp/nnnn non Skipped
  if (abs(delta_dir)+abs(delta_exch) < 1.0d-6) return
endif


integral_factor = t3_DD_CONST
!! NOTE :: Remember that radial functions already have the factor 1/b**3
integral_factor = integral_factor * 0.5d0 * (HO_b**3)
integral_factor = integral_factor  / ((2.0d0 + alpha_DD)**1.5d0)
integral_factor = integral_factor * 4 * pi  ! add Lebedev norm factor


do i_r = 1, r_dim

  radial = weight_R(i_r) * radial_2b_sho_memo(a_sh, c_sh, i_r) &
                         * radial_2b_sho_memo(b_sh, d_sh, i_r) &
                         * exp((2.0d0+alpha_DD) * (r(i_r)/HO_b)**2)
  !! NOTE: the inclusion of the exponential part is necessary due the form of
  !! of the density and radial functions with the exp(-r/b^2) for stability
  !! requirement in larger shells.

  do i_ang = 1, angular_dim
      !! ====================================================================
      aux_dir  = zzero
      if (abs(delta_dir) > 1.0d-15) then
        do K = abs(ja - jc) / 2, (ja + jc) / 2
          M = (mc - ma)/2
          if ((MOD(K + la + lc, 2) == 1).OR.(abs(M) > K)) cycle

          do K2 = abs(jb - jd) / 2, (jb + jd) / 2
            !! NOTE:: in DD Hamiltonian loop, condition ma+mb=mc+md -> M=M2
            M2 = (md - mb)/2
            if ((MOD(K2 + lb + ld, 2) == 1).OR.(abs(M2) > K2)) cycle

            ind_km   = angular_momentum_index(K,  M,  .FALSE.)
            ind_km_q = angular_momentum_index(K2, M2, .FALSE.)

            aux_dir = aux_dir + (dens_Y_KM_me(ind_jm_a, ind_jm_c, ind_km)  * &
                                 dens_Y_KM_me(ind_jm_b, ind_jm_d, ind_km_q)* &
                                 sph_harmonics_memo(ind_km,   i_ang) * &
                                 sph_harmonics_memo(ind_km_q, i_ang))
          enddo
        enddo
      endif
      !! ====================================================================
      aux_exch = zzero
      if (abs(delta_exch) > 1.0d-15) then
        do K = abs(ja - jd) / 2, (ja + jd) / 2
          M = (md - ma)/2
          if ((MOD(K + la + ld, 2) == 1).OR.(abs(M) > K)) cycle

          do K2 = abs(jb - jc) / 2, (jb + jc) / 2
            M2 = (mc - mb)/2
            if ((MOD(K2 + lb + lc, 2) == 1).OR.(abs(M2) > K2)) cycle

            ind_km   = angular_momentum_index(K,  M,  .FALSE.)
            ind_km_q = angular_momentum_index(K2, M2, .FALSE.)

            aux_exch = aux_exch + (dens_Y_KM_me(ind_jm_a, ind_jm_d, ind_km)  *&
                                   dens_Y_KM_me(ind_jm_b, ind_jm_c, ind_km_q)*&
                                   sph_harmonics_memo(ind_km,   i_ang) *&
                                   sph_harmonics_memo(ind_km_q, i_ang))
          enddo
        enddo
      endif
      !! ====================================================================

      angular = weight_LEB(i_ang)

      if (ALL_ISOS) then
        !v_nnnn = v_pppp
        aux = radial * angular * (1-x0_DD_FACTOR) * ((aux_dir) - (aux_exch))
        v_dd_value(1) = v_dd_value(1) + (aux * dens_alpha(i_r, i_ang))
        v_dd_value(4) = v_dd_value(1)

        ! pn pn
        aux = radial * angular * (aux_dir + (x0_DD_FACTOR*aux_exch))
        v_dd_value(2) = v_dd_value(2) + (aux * dens_alpha(i_r, i_ang))
        ! pn np
        aux = radial * angular * ((x0_DD_FACTOR*aux_dir) + aux_exch)
        v_dd_value(3) = v_dd_value(3) - (aux * dens_alpha(i_r, i_ang))
      else
        ! first element is the one calculated (rest remain at zero))
        aux = radial * angular * ((delta_dir*aux_dir) - (delta_exch*aux_exch))
        v_dd_value(1) = v_dd_value(1) + (aux * dens_alpha(i_r, i_ang))
      end if

      !! Loop for the Rearrangement term
      if ((EVAL_REARRANGEMENT).AND.(EVAL_EXPLICIT_FIELDS_DD)) then
        !!!! if (dreal(aux)**2 + dimag(aux)**2 < 1.0d-15) cycle
        !! NOTE: don't put a skip for the |aux|<1e-15, afect the tolerance of
        !! the sum, it doesn't match exactly with the field calculation.
        aux_rea = alpha_DD * integral_factor * aux * dens_alpm1(i_r, i_ang)

        do kk1 = 1, HOspO2
          kk1N = kk1 + HOspO2
          do kk2 = 1, HOspO2
            kk2N = kk2 +HOspO2

            rea = aux_rea  * rea_common_RadAng(kk1,kk2, i_r, i_ang)
            rearrangement_me(kk1,kk2)   = rearrangement_me(kk1,kk2)   + rea
            rearrangement_me(kk1N,kk2N) = rearrangement_me(kk1N,kk2N) + rea
            enddo
        enddo
      endif

   enddo ! angular iter_
enddo    ! radial  iter_

if (EVAL_EXPLICIT_FIELDS_DD) then
  TOP = abs(maxval(real(rearrangement_me)))
  LOW = abs(minval(real(rearrangement_me)))
  if (TOP > 1.0D+10) then !((TOP < 1.0D+10).AND.(LOW > 1.E-10)) then
      print "(A,4I3,2F20.10)", "!! REA", a,b,c,d, &
          minval(real(rearrangement_me)), maxval(real(rearrangement_me))
  endif
endif

v_dd_val_Real(1) = real(v_dd_value(1), r64) * integral_factor
v_dd_val_Real(2) = real(v_dd_value(2), r64) * integral_factor
v_dd_val_Real(3) = real(v_dd_value(3), r64) * integral_factor
v_dd_val_Real(4) = real(v_dd_value(4), r64) * integral_factor

if (abs(imag(v_dd_value(1))) > 1.0d-15 ) then
    print "(A,F10.8,A,F18.15)", "  [FAIL] v_DD_abcd is not Real =", &
        real(v_dd_value(1)), " +j ", imag(v_dd_value(1))
endif

return

end function matrix_element_v_DD

!-----------------------------------------------------------------------------!
! subroutine calculate_densityDep_hamiltonian                                 !
!                                                                             !
! Computes density dependent two body matrix elements over the density average!
!                                                                             !
! * Update May 22, this subroutine fixes the DD matrix elements, that do not  !
! change along all the process. After analyse the evolution on nuclei along   !
! the process do not reach the cutoff, the dimension of H do not change in any!
! case (it do not appear any m.e. and if any turns to 0 it do not affect).    !
! * ASSERTION is necessary to check if H_abcd(kk) matches with current m.e.   !
!                                                                             !
! * Update Jun 14, fix the order of the abcd loop to the set_hamiltonian_2body!
!     MAJOR CHANGES ( cannot be used directly )                               !
! * Update May 10/23, modification of the hamiltonian evaluation for the      !
! valence space for exporting is used, the rearrangement part will not be     !
! evaluated if it is not calculated both (EVAL_EXPLICIT_FIELDS_DD = TRUE) and !
! (EVAL_REARRANGEMENT = TRUE)                                                 !
!-----------------------------------------------------------------------------!
subroutine calculate_densityDep_hamiltonian

integer(i16) :: ared, bred, cred, dred
integer(i32) :: ht, j, t, tmax, uth6=uth+8, uth7=uth+9, uth8=uth+10, ialloc=0,&
                a, ma, la, ta, b, mb, lb, tb, dmax, bmax,  bmin, cmin, dmin,&
                c, mc, lc, tc, d, md, ld, td, aa, bb, cc, dd, an, bn, cn, dn
integer(i64) :: kk, i, kkk
integer, parameter :: CONVERG_ITER = 10000
real(r64) :: xja, xjb, xjc, xjd, xjtot, xttot, phasab, phascd, Vtmp, &
             Vcut, Vdec, Vred
real(r64), dimension(4) :: me_Vdec, me_VGRc, me_VGRc1, me_VGRc2, me_VGRc3
real(r64), dimension(:), allocatable :: hamil_temp, hamil_temp_2
character(len=25) :: filename
logical   :: ALL_ISOS

ALL_ISOS = (.NOT.EVAL_EXPLICIT_FIELDS_DD)
!! NOTE: if not explicit evaluation of fields, the process was called to export v_DD
!! if ALL_ISOS = .TRUE., this subroutine can only be called once !!

if (iteration < CONVERG_ITER) then
  rewind(uth6)
  rewind(uth7)
  rewind(uth8)
  !!! Computes the two-body matrix elements in m-scheme
  open  (uth6, status='scratch', action='readwrite', access='stream', &
               form='unformatted')
  open  (uth7, status='scratch', action='readwrite', access='stream', &
               form='unformatted')
  open  (uth8, status='scratch', action='readwrite', access='stream', &
               form='unformatted')
endif

Vcut = 5.0d-14
if (ALL_ISOS) Vcut = 1.0d-9
kk = 0
NOT_DEL_FILE = .FALSE.

rearrang_field = zero

do aa = 1, WBsp_dim / 2 ! (prev = HOsp_dim)
  a  = WBtoHOsp_index(aa) !VStoHOsp_index(aa)
  an = a + (HOsp_dim / 2)
  la = HOsp_l(a)
  ma = HOsp_2mj(a)
  ta = HOsp_2mt(a)

  bmin = aa!+1
  if (evalFullSPSpace) bmin = 1
  do bb = bmin, WBsp_dim / 2 ! (prev = HOsp_dim)
    b  = WBtoHOsp_index(bb)
    bn = b  + (HOsp_dim / 2)
    lb = HOsp_l(b)
    mb = HOsp_2mj(b)
    tb = HOsp_2mt(b)

    if ((.NOT.evalFullSPSpace).AND.( ma + mb < 0 )) cycle

    cmin = aa
    if (evalFullSPSpace) cmin = 1
    do cc = cmin, WBsp_dim / 2 ! (prev = HOsp_dim)
      c  = WBtoHOsp_index(cc)
      cn = c + (HOsp_dim / 2)
      lc = HOsp_l(c)
      mc = HOsp_2mj(c)
      tc = HOsp_2mt(c)

      dmin = 1
      dmax = WBsp_dim / 2 ! (prev = HOsp_dim)
      if (.NOT.evalFullSPSpace) then
        dmin = cc!+1
        if ( cc == aa ) dmax = bb
      endif
      do dd = dmin, dmax
        d  = WBtoHOsp_index(dd)
        dn = d + (HOsp_dim / 2)
        ld = HOsp_l(d)
        md = HOsp_2mj(d)
        td = HOsp_2mt(d)

        if ( ta + tb /= tc + td ) cycle

        rearrangement_me = zero

        me_Vdec   = matrix_element_v_DD(a,b, c,d, ALL_ISOS)
        if      (EXPORT_GRAD_DD) then
          me_VGRc = matrix_element_v_gradientDD(a,b, c,d)
        else if (EXPORT_PREA_DD) then
          me_VGRc = matrix_element_pseudoRearrangement(a,b, c,d)

          !! test ANTISYMMETRY  ===================
!me_VGRc1 = matrix_element_pseudoRearrangement(a,b, d,c)
!me_VGRc2 = matrix_element_pseudoRearrangement(b,a, c,d)
!me_VGRc3 = matrix_element_pseudoRearrangement(b,a, d,c)
!
!tmax = 0
!do t = 1, 4
!  if (almost_equal(me_VGRc(t), zero, 1.0d-9)) cycle
!
!  if (.NOT. almost_equal(me_VGRc(t), -one * me_VGRc1(t), 1.0d-6)) then
!    print "(A,I3,2F15.8)"," AntySym-ERR 1 -(ab,dc)t:",t,me_VGRc(t),me_VGRc1(t)
!    tmax = tmax + 1
!  end if
!  if (.NOT. almost_equal(me_VGRc(t), -one * me_VGRc2(t), 1.0d-6)) then
!    print "(A,I3,2F15.8)"," AntySym-ERR 2 -(ba,cd)t:",t,me_VGRc(t),me_VGRc2(t)
!    tmax = tmax + 1
!  end if
!  if (.NOT. almost_equal(me_VGRc(t), me_VGRc3(t), 1.0d-6)) then
!    print "(A,I3,2F15.8)"," AntySym-ERR 3 +(ba,dc)t:",t,me_VGRc(t),me_VGRc3(t)
!    tmax = tmax + 1
!  end if
!end do
!
!if (tmax .GT. 0) print "(A,I3,A,4I5)", "   [ERROR]s:[",tmax,"] a,b,c,d",a,b,c,d

          !! ======================================
        endif

        !!! Select only matrix elements above a given cutoff to reduce the
        !!! CPU time and storage
        if (ALL_ISOS) then
          if ((maxval(me_Vdec).GE.Vcut).OR.(abs(minval(me_Vdec)).GE.Vcut)) then
            kk = kk + 1

            ared = int(a,i16)
            bred = int(b,i16)
            cred = int(c,i16)
            dred = int(d,i16)
            write(uth6) ared, bred, cred, dred
            write(uth7) me_Vdec(1), me_Vdec(2), me_Vdec(3), me_Vdec(4)
            if (EXPORT_GRAD_DD .OR. EXPORT_PREA_DD) then
              write(uth8) me_VGRc(1), me_VGRc(2), me_VGRc(3), me_VGRc(4)
              !! evaluate the rearrangement
              rearrang_field(a ,c ) = rearrang_field(a ,c ) + &
                  (me_VGRc(1)*dens_rhoLR(d ,b ) + me_VGRc(2)*dens_rhoLR(dn,bn))
              rearrang_field(an,cn) = rearrang_field(an,cn) + &
                  (me_VGRc(2)*dens_rhoLR(d ,b ) + me_VGRc(4)*dens_rhoLR(dn,bn))
              rearrang_field(a ,cn) = rearrang_field(a ,cn) + &
                  (me_VGRc(3)*dens_rhoLR(d ,bn))
              rearrang_field(an,c ) = rearrang_field(an,c ) + &
                  (me_VGRc(3)*dens_rhoLR(dn,b ))
            endif
          endif

        else !! normal case, the element is the 1st one (isospin_ form abcd)
          Vdec = me_Vdec(1)
          if ( abs(Vdec) > Vcut ) then
            kk = kk + 1

            if (iteration < CONVERG_ITER) then
              ared = int(a,i16)
              bred = int(b,i16)
              cred = int(c,i16)
              dred = int(d,i16)
              Vred = real(Vdec,r64)  !real(Vdec,r32)
              write(uth6) ared, bred, cred, dred
              write(uth7) Vred

            else
              !! ASSERT if the abcd match the previous case:
              if ((a /= hamil_DD_abcd(1+4*(kk-1))).OR. &
                  (b /= hamil_DD_abcd(2+4*(kk-1))).OR. &
                  (c /= hamil_DD_abcd(3+4*(kk-1))).OR. &
                  (d /= hamil_DD_abcd(4+4*(kk-1)))) then
                print *, "[ASSERTION ERROR]: the final v_abcd indexes change"
                print '(A,4I4,A,I9,A,I6,A,4I4)', "abcd",a,b,c,d," [kk=", kk, &
                  "] [iteration=",iteration,"] to",&
                  hamil_DD_abcd(1+4*(kk-1)),hamil_DD_abcd(2+4*(kk-1)), &
                  hamil_DD_abcd(3+4*(kk-1)),hamil_DD_abcd(4+4*(kk-1))

              endif
              !! Fix the new matrix element
              hamil_DD_H2(kk) = Vdec
            endif

            if ((EVAL_REARRANGEMENT).AND.(EVAL_EXPLICIT_FIELDS_DD)) then
              call calculate_rearrang_field_explicit(a, b, c, d, Vdec)
            endif
          endif
        endif ! select the process or to export the matrix elements

      enddo  !end loop d
    enddo  !end loop c
  enddo  !end loop b
  if (.NOT.EVAL_EXPLICIT_FIELDS_DD) call progress_bar_iteration(aa, WBsp_dim/2)
enddo  !end loop a

!!! At the first iteration, the values of the hamiltonian are saved via file
if (ALL_ISOS) then

  hamil_DD_H2dim     = kk
  hamil_DD_H2dim_all = kk

  allocate( hamil_DD_H2_byT    (4, hamil_DD_H2dim), &
            hamil_GradDD_H2_byT(4, hamil_DD_H2dim), &
            hamil_DD_abcd(4*hamil_DD_H2dim), hamil_temp(4*hamil_DD_H2dim),&
            hamil_temp_2 (4*hamil_DD_H2dim), stat=ialloc )
  if ( ialloc /= 0 ) stop 'Error during allocation of array of indices [DD]'
  rewind(uth6)
  rewind(uth7)
  rewind(uth8)
  read  (uth6) (hamil_DD_abcd(kk), kk=1, 4*hamil_DD_H2dim)
  read  (uth7) (hamil_temp(kk),    kk=1, 4*hamil_DD_H2dim)
  read  (uth8) (hamil_temp_2(kk),  kk=1, 4*hamil_DD_H2dim)
  close (uth6)
  close (uth7)
  close (uth8)

  do kk = 1, hamil_DD_H2dim
    hamil_DD_H2_byT(1, kk) =  hamil_temp(4*(kk-1) + 1) ! pp pp
    hamil_DD_H2_byT(2, kk) =  hamil_temp(4*(kk-1) + 2) ! pn pn
    hamil_DD_H2_byT(3, kk) =  hamil_temp(4*(kk-1) + 3) ! pn np
    hamil_DD_H2_byT(4, kk) =  hamil_temp(4*(kk-1) + 4) ! nn nn

    hamil_GradDD_H2_byT(1, kk) =  hamil_temp_2(4*(kk-1) + 1) ! pp pp
    hamil_GradDD_H2_byT(2, kk) =  hamil_temp_2(4*(kk-1) + 2) ! pn pn
    hamil_GradDD_H2_byT(3, kk) =  hamil_temp_2(4*(kk-1) + 3) ! pn np
    hamil_GradDD_H2_byT(4, kk) =  hamil_temp_2(4*(kk-1) + 4) ! nn nn
  enddo

  deallocate(hamil_temp, hamil_temp_2)
  if (.NOT.(EXPORT_GRAD_DD .OR. EXPORT_PREA_DD)) deallocate(hamil_GradDD_H2_byT)

!  call print_uncoupled_hamiltonian_DD(ALL_ISOS)
!  call print_uncoupled_hamiltonian_H2

else if (iteration < CONVERG_ITER) then !!! Normal Gradient DD dep. process ****

  if ((iteration > 1).AND.(EVAL_EXPLICIT_FIELDS_DD)) then
    deallocate(hamil_DD_H2, hamil_DD_abcd, hamil_DD_trperm)
  end if

  hamil_DD_H2dim = kk ! final value
  hamil_DD_H2dim_all = kk
  ! this index is used (only) in print_hamilt and cmpi in read reducced hamiltonian

  !!! Final allocation and reading of two-body matrix elements
  allocate ( hamil_DD_H2(hamil_DD_H2dim), hamil_DD_abcd(4*hamil_DD_H2dim), &
          stat=ialloc )
  if ( ialloc /= 0 ) stop 'Error during allocation of array of indices [DD]'
  rewind(uth6)
  rewind(uth7)

  read(uth6) (hamil_DD_abcd(kk), kk=1, 4*hamil_DD_H2dim)
  read(uth7) (hamil_DD_H2(kk), kk=1, hamil_DD_H2dim)

  close(uth6)
  close(uth7)
  !!! Determines the permutation needed to obtain the time-reversed two-body
  !!! matrix elements

  call reconstruct_2body_DD_timerev

  call print_uncoupled_hamiltonian_DD(.FALSE.)
endif

iteration = iteration + 1
!print *, "[OK] Testing DD hamiltonian"

end subroutine calculate_densityDep_hamiltonian

!------------------------------------------------------------------------------!
! subroutine print_uncoupled_hamiltonian_DD                                    !
!                                                                              !
! Auxiliary method to print an uncoupled list of matrix elements, depending on !
! the case of explicit internal evaluation as main process or the last export  !
! of the hamltonian_(ALL_ISOS)                                                 !
!------------------------------------------------------------------------------!
subroutine print_uncoupled_hamiltonian_DD(ALL_ISOS)

logical, intent(in) :: ALL_ISOS
real(r64) :: Vdec
character(len=20) :: filename
integer :: i, kk, a, b, c, d

filename = 'uncoupled_DD.2b'
print "(3A,3I12)", &
  "[  ] EXPORT Hamiltonian (uncoupled) for reduced Valence space [", filename,&
  "]   * VS array DD/all possible=", hamil_DD_H2dim, (VSsp_dim/2)**4

if (.NOT.ALL_ISOS) then
print '(A,F10.6,A,F10.6)'," *Top H2",MINVAL(hamil_DD_H2),' ',MAXVAL(hamil_DD_H2)
  if (iteration < 30) then
    if (iteration < 6) then ! first iterations have the larger differences
        write(filename, "(A,I0.4,A)") 'uncoupled_DD_', iteration, '.2b'
    else if (mod(iteration, 3).EQ.0) then
        write(filename, "(A,I0.4,A)") 'uncoupled_DD_', iteration, '.2b'
        !filename = 'uncoupled_DD_'//char(iteration)//'.2b'
        !print *, "printing:", filename, "   iteration:",":", char(iteration)
    endif
  elseif (iteration < 421) then
    if (mod(iteration, 30).EQ.0) then
        write(filename, "(A,I0.4,A)") 'uncoupled_DD_', iteration, '.2b'
        !filename = 'uncoupled_DD_'//char(iteration)//'.2b'
        !print *, "printing:", filename, "   iteration:",":", char(iteration)
    endif
  endif
endif

open (123, file=filename)
if (ALL_ISOS) then
  write(123,fmt='(A)')"//SING PART INDEX (i_sp, i_sh, n,l,2j,2mj,tr)"
  do kk=1, WBsp_dim
    i = WBtoHOsp_index(kk)
    write(123, fmt='(I4,6(A,I4))') i,',', HOsp_sh(i), ',', HOsp_n(i),&
      ',', HOsp_l(i),',', HOsp_2j(i),'/2,', HOsp_2mj(i),'/2,', HOsp_tr(i)
  enddo
  write(123, fmt='(3A,3I12)')"//(vs)a    b    c    d              pppp", &
    "              pnpn              pnnp              nnnn    ", &
    "* VS array DD/noDD DIM/ALL=", hamil_DD_H2dim, hamil_H2dim, (VSsp_dim/2)**4
else
  write(123,fmt='(A)')"//SING PART INDEX (sp_vs,i_sp, i_sh, n,l,2j,2m, 2mt,tr)"
  do i=1, HOsp_dim
    write(123, fmt='(I3,7(A,I4))') i,',', HOsp_sh(i), &
      ',', HOsp_n(i),',', HOsp_l(i),',', HOsp_2j(i),',', HOsp_2mj(i), &
      ',', HOsp_2mt(i),',', HOsp_tr(i)
  enddo
  write(123, fmt='(2A,2I8)')"//DD ME     a     b     c     d                ",&
    "h2bDD    DD/noDD DIM=", hamil_DD_H2dim, hamil_H2dim
endif

do kk = 1, hamil_DD_H2dim
    a = hamil_DD_abcd(1+4*(kk-1))
    b = hamil_DD_abcd(2+4*(kk-1))
    c = hamil_DD_abcd(3+4*(kk-1))
    d = hamil_DD_abcd(4+4*(kk-1))
    if (ALL_ISOS) then
      write(123, fmt='(I7,3I5,4F18.12)')  a, b, c, d, hamil_DD_H2_byT(1,kk), &
        hamil_DD_H2_byT(2,kk), hamil_DD_H2_byT(3,kk), hamil_DD_H2_byT(4,kk)
    else
      Vdec = hamil_DD_H2(kk)
      write(123, fmt='(4I6,F25.18)') a, b, c, d, Vdec
    endif

enddo
close(123)

print "(3A,3I12)", &
  "[OK] EXPORT Hamiltonian (uncoupled) for reduced Valence space [", filename,&
  "]   * VS array DD/all possible=", hamil_DD_H2dim, (VSsp_dim/2)**4

end subroutine print_uncoupled_hamiltonian_DD

!------------------------------------------------------------------------------!
! subroutine print_uncoupled_hamiltonian_DD                                    !
!                                                                              !
! Auxiliary method to print an uncoupled list of matrix elements from the H2   !
! read from the .2b file, permute and reconstruct the TR matrix elements.      !
!------------------------------------------------------------------------------!
subroutine print_uncoupled_hamiltonian_H2

integer   :: i, kk, i1, i2, i3, i4, it, perm, uth6=uth+8,uth7=uth+9, ialloc=0,&
             ared, bred, cred, dred, ndim, ndim2, k1, k2, POW10, spo2, point_,&
             j1, j2, j3, j4, tt, red_dim, sh_vs
integer(i64) :: indx_, ind_r
integer(i64), dimension(:), allocatable :: sort_indx, red_indx
integer,      dimension(:), allocatable :: temp_abcd, sort_pointer, sort_isos,&
                                           sort_red_pointer
real(r64) :: h2b, aux, i1_t, i2_t, i3_t, i4_t
real(r64), dimension(:),   allocatable :: temp_hamil
real(r64), dimension(:,:), allocatable :: temp_hamil_byT
integer,   dimension(:,:,:,:), allocatable :: red_abcd
integer, dimension(4) :: sh_curr
logical :: found
integer, dimension(:,:,:,:), allocatable:: registered_h2b ! test which m.e. is registered

print "(A,I3)", "[  ] EXPORT Hamiltonian (uncoupled) for current interaction.&
                 & sorting_dim (<5)=", floor(log10(HOsp_dim + 0.0d0)) + 1
if ((floor(log10(HOsp_dim + 0.0d0)) + 1) .GE. 5) then
  print "(A)", "[NO]  EXPORT sorting, method unavailable for sp_dim > 5. PASS"
  return
end if

allocate(registered_h2b(HOsp_dim,HOsp_dim,HOsp_dim,HOsp_dim))

open  (uth6, status='scratch', action='readwrite', access='stream', &
             form='unformatted')
open  (uth7, status='scratch', action='readwrite', access='stream', &
             form='unformatted')
! read and export all the possible matrix elements (non sorted)
ndim = 0
spo2 = WBsp_dim / 2

!! registered_h2b is an
registered_h2b = 0

open (3333, file="hamil_abcd.gut")
do kk = 1, hamil_H2dim
  i1 = hamil_abcd(1+4*(kk-1))
  i2 = hamil_abcd(2+4*(kk-1))
  i3 = hamil_abcd(3+4*(kk-1))
  i4 = hamil_abcd(4+4*(kk-1))

  if (.NOT.evalQuasiParticleVSpace) then
    !! Skip if the state is not in the WB
    found = .TRUE.
    sh_curr = (/HOsh_ant(i1), HOsh_ant(i2), HOsh_ant(i3),  HOsh_ant(i4)/)
    do k1 = 1, 4
      k2 = 0
      do sh_vs=1, VSsh_dim
        if (sh_curr(k1) .EQ. VSsh_list(sh_vs)) k2 = 1
      end do
      if (k2 .EQ. 0) then
        found = .FALSE.
        EXIT
      endif
    end do
    if (.NOT. found) cycle
  endif

  h2b  = hamil_H2(kk)
  perm = hamil_trperm(kk)

  write(3333,fmt="(2i7,A,4i4,A,F13.6)")kk,perm," indx:",i1,i2,i3,i4," =",h2b

  !!! Loop on time reversal
  do it = 1, 2
    if ( it == 2 ) then
      if ( HOsp_2mj(i1) + HOsp_2mj(i2) == 0 ) cycle
      call find_timerev(perm,i1,i2,i3,i4)
      h2b = sign(one,perm*one) * h2b
    endif

    ared = int(i1,i16)
    bred = int(i2,i16)
    cred = int(i3,i16)
    dred = int(i4,i16)
    write(uth6) ared, bred, cred, dred
    write(uth6) ared, bred, dred, cred
    write(uth6) bred, ared, cred, dred
    write(uth6) bred, ared, dred, cred
    write(uth7)  h2b, -h2b, -h2b,  h2b
    ndim = ndim + 4

    !! 1. Criteria from module_fields.calculate_fields (general)
    registered_h2b(i1,i2,i3,i4) = registered_h2b(i1,i2,i3,i4) + 1
    registered_h2b(i1,i2,i4,i3) = registered_h2b(i1,i2,i4,i3) + 1
    registered_h2b(i2,i1,i3,i4) = registered_h2b(i2,i1,i3,i4) + 1
    registered_h2b(i2,i1,i4,i3) = registered_h2b(i2,i1,i4,i3) + 1
    if ((kdelta(i1,i3) * kdelta(i2,i4)) .NE. 1) then
      registered_h2b(i3,i4,i1,i2) = registered_h2b(i3,i4,i1,i2) + 1
      registered_h2b(i3,i4,i2,i1) = registered_h2b(i3,i4,i2,i1) + 1
      registered_h2b(i4,i3,i1,i2) = registered_h2b(i4,i3,i1,i2) + 1
      registered_h2b(i4,i3,i2,i1) = registered_h2b(i4,i3,i2,i1) + 1
    endif

!    !! 2. Criteria from module_fields.calculate_fields_diag
!    registered_h2b(i1,i2,i3,i4) = registered_h2b(i1,i2,i3,i4) + 1
!    registered_h2b(i1,i2,i4,i3) = registered_h2b(i1,i2,i4,i3) + 1
!
!    if ((i1.EQ.i3) .AND. (i2.NE.i4)) then
!      registered_h2b(i3,i4,i1,i2) = registered_h2b(i3,i4,i1,i2) + 1
!    endif
!    if (i2.LE.i4) then
!      registered_h2b(i2,i1,i4,i3) = registered_h2b(i2,i1,i4,i3) + 1
!    endif
!    if (i2.LE.i3) then
!      registered_h2b(i2,i1,i3,i4) = registered_h2b(i2,i1,i3,i4) + 1
!    endif
!
!    if ((i1.NE.i3) .OR. (i2.NE.i4)) then
!      if (i4.LE.i2) then
!        registered_h2b(i4,i3,i2,i1) = registered_h2b(i4,i3,i2,i1) + 1
!      endif
!      if (i3.LE.i2) then
!        registered_h2b(i3,i4,i2,i1) = registered_h2b(i3,i4,i2,i1) + 1
!      endif
!    endif
!    !! -----------------------------------------------------------------------

    if ((kdelta(i1,i3) * kdelta(i2,i4)) .EQ. 1) cycle

    write(uth6) cred, dred, ared, bred
    write(uth6) dred, cred, ared, bred
    write(uth6) cred, dred, bred, ared
    write(uth6) dred, cred, bred, ared
    write(uth7)  h2b, -h2b, -h2b,  h2b
    ndim = ndim + 4


  enddo
end do
close(3333)
!! Allocate in a temporal hamiltonian
allocate( sort_indx(ndim), sort_pointer(ndim),&
          sort_isos(ndim), temp_hamil(ndim), temp_abcd(4*ndim),&
          red_indx (ndim), sort_red_pointer(ndim), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of array of indices [BB]'
rewind(uth6)
rewind(uth7)
read  (uth6) (temp_abcd (kk), kk=1, 4*ndim)
read  (uth7) (temp_hamil(kk), kk=1,   ndim)
close (uth6)
close (uth7)

!! Sort the indexes in order of the (a,b,c,d) by sorting a number on base HO dim
sort_indx = 0
sort_pointer = 0
sort_isos = 0
red_indx  = 0
red_dim   = 0
sort_red_pointer = 0
POW10 = floor(log10(HOsp_dim + 0.0d0)) + 1
open (3333, file="temp_abcd_init.gut")
do kk = 1, ndim
  i1 = temp_abcd(4*(kk-1) + 1)
  i2 = temp_abcd(4*(kk-1) + 2)
  i3 = temp_abcd(4*(kk-1) + 3)
  i4 = temp_abcd(4*(kk-1) + 4)

  j1 = i1
  j2 = i2
  j3 = i3
  j4 = i4

  if (i1 .GT. spo2) j1 = i1 - spo2
  if (i2 .GT. spo2) j2 = i2 - spo2
  if (i3 .GT. spo2) j3 = i3 - spo2
  if (i4 .GT. spo2) j4 = i4 - spo2

  indx_ = nint(i1*(10**(3*POW10))+i2*(10**(2*POW10))+i3*(10.0d0**(POW10)) + i4)
  ind_r = nint(j1*(10**(3*POW10))+j2*(10**(2*POW10))+j3*(10.0d0**(POW10)) + j4)

  if     ((i1 .GT. spo2).AND.(i2 .GT. spo2)) then
    tt = 4    ! nn_ nn_
  elseif ((i1 .LE. spo2).AND.(i2 .LE. spo2)) then
    tt = 1    ! pp_ pp_
  else
    tt = 4*HOsp_2mt(i1) + 2*(HOsp_2mt(i2) + HOsp_2mt(i3)) + HOsp_2mt(i4) - 1
    tt = 3 + (tt / 2)
    if ((tt.EQ.1) .OR. (tt.EQ.4)) then
      tt = 2  ! pn_ pn_
    else
      tt = 3  ! pn_ np_
    end if
  endif

  sort_indx(kk)     = indx_
  sort_pointer(kk)  = kk
  sort_isos(kk)     = tt

  if (red_dim .EQ.0) then ! 1st
    red_indx(1) = ind_r
    red_dim = 1
    sort_red_pointer(kk) = 1
  else
    found = .FALSE.
    do k1 = 1, red_dim
      if (red_indx(k1).EQ.ind_r) then
        found = .TRUE.
        sort_red_pointer(kk) = k1
        EXIT
      endif
    end do
    if (.NOT.found) then
      red_dim = red_dim + 1
      red_indx(red_dim) = ind_r
      sort_red_pointer(kk) = red_dim
    end if
  endif

  write(3333,fmt="(i6,2(A,4i4),A,2i5,A,2i10,A,F15.6)")kk,"  indx:",i1,i2,i3,i4,&
        "  r(",j1,j2,j3,j4,") tt,red:", tt,red_dim,"  hash(t/r):", indx_,ind_r,&
        " =",temp_hamil(kk)
enddo
close(3333)

deallocate(registered_h2b)

allocate(temp_hamil_byT(4, red_dim), &
         red_abcd(WBsp_dim , WBsp_dim , WBsp_dim , WBsp_dim))
temp_hamil_byT = zero
red_abcd       = 0
! bubble sorting
do k1 = 1, ndim
  do k2 = k1+1, ndim
    if (sort_indx(k1) .GT. sort_indx(k2)) then
      indx_ = sort_indx(k2)
      sort_indx(k2) = sort_indx(k1)
      sort_indx(k1) = indx_

      point_ = sort_pointer(k2)
      sort_pointer(k2) = sort_pointer(k1)
      sort_pointer(k1) = point_

      tt = sort_isos(k2)
      sort_isos(k2) = sort_isos(k1)
      sort_isos(k1) = tt

      point_ = sort_red_pointer(k2)
      sort_red_pointer(k2) = sort_red_pointer(k1)
      sort_red_pointer(k1) = point_
    end if
  enddo
  ! sort also the reduced one
  if (k1 .GT. red_dim) cycle
  do k2 = k1+1, red_dim
    if (red_indx(k1) .GT. red_indx(k2)) then
      ind_r = red_indx(k2)
      red_indx(k2) = red_indx(k1)
      red_indx(k1) = ind_r
    end if
  enddo
end do

!! Assign in the temp_hamil_byT the element by T
open (3333, file="temp_abcd_sorted.gut")
do k1 = 1, ndim
  kk = sort_pointer(k1)

  i1 = temp_abcd(4*(kk-1) + 1)
  i2 = temp_abcd(4*(kk-1) + 2)
  i3 = temp_abcd(4*(kk-1) + 3)
  i4 = temp_abcd(4*(kk-1) + 4)

  k2 = sort_red_pointer(k1) ! extract the index of the reduced space
  tt = sort_isos(k1)
  !NOTE: the sort_red_pointer still points in the reduced non zero list,
  ! here, we are not assigning k2 in order, the order will appear while reading
  temp_hamil_byT(tt, k2) = temp_hamil(kk)

  write(3333, fmt="(3i8,A,4i4,A,i2,A,F15.6)")k1,kk,k2," indx:", i1,i2,i3,i4, &
                                             "  (t:",tt,") =",temp_hamil(kk)

  if(i1 .GT. spo2) i1 = i1 - spo2
  if(i2 .GT. spo2) i2 = i2 - spo2
  if(i3 .GT. spo2) i3 = i3 - spo2
  if(i4 .GT. spo2) i4 = i4 - spo2

!  if ((red_abcd(1,k2) .NE. 0).AND.(red_abcd(1, k2) .NE. i1)) then
!    print "(A,2i5)", "[WARN] 1 NE prev assinged: ", red_abcd(1, k2) , i1
!  end if
!  if ((red_abcd(2,k2) .NE. 0).AND.(red_abcd(2, k2) .NE. i2)) then
!    print "(A,2i5)", "[WARN] 2 NE prev assinged: ", red_abcd(2, k2) , i2
!  end if
!  if ((red_abcd(3,k2) .NE. 0).AND.(red_abcd(3, k2) .NE. i3)) then
!    print "(A,2i5)", "[WARN] 3 NE prev assinged: ", red_abcd(3, k2) , i3
!  end if
!  if ((red_abcd(4,k2) .NE. 0).AND.(red_abcd(4, k2) .NE. i4)) then
!    print "(A,2i5)", "[WARN] 4 NE prev assinged: ", red_abcd(4, k2) , i4
!  end if
!  print "(A)", "-----------------------"

  red_abcd(i1,i2,i3,i4) = k2
end do
close(3333)

!! print the matrix elements in a file
open (336, file="uncoupled_BB.2b")
write(336, fmt="(A)") "//SING PART INDEX (sp_vs,i_sp, i_sh, n,l,2j,2m, tr)"
do kk = 1, WBsp_dim
  i = WBtoHOsp_index(kk)
  write(336, fmt='(I4,6(A,I4))') i,',', HOsp_sh(i), ',', HOsp_n(i),&
      ',', HOsp_l(i),',', HOsp_2j(i),'/2,', HOsp_2mj(i),'/2,', HOsp_tr(i)
enddo
write(336, fmt="(a)") "// Hamiltonian uncoupled for the m.e. given (all perm)"
write(336, fmt="(a)") "//  a    b    c    d              pppp              &
                      &pnpn              pnnp              nnnn"
do k1 = 1, red_dim
  ! red_index is sorted, so we extract the HOsp index from it,
  ind_r = red_indx(k1)
  j1    = int(ind_r / nint((10.0d0**(3*POW10)), i32))
  ind_r = MOD(ind_r,  nint((10.0d0**(3*POW10))))
  j2    = int(ind_r / nint((10.0d0**(2*POW10)), i32))
  ind_r = MOD(ind_r,  nint((10.0d0**(2*POW10))))
  j3    = int(ind_r / nint((10.0d0**(POW10)), i32))
  ind_r = MOD(ind_r,  nint((10.0d0**(POW10))))
  j4    = int(ind_r,  i32)

  k2 = red_abcd(j1,j2,j3,j4)
  if (k2.EQ.0) print "(A,4i3)", "[ERR] Invalid Access 0= red_abcd",j1,j2,j3,j4

  write(336, fmt='(4I5,4F18.12)') j1, j2, j3, j4, temp_hamil_byT(1,k2), &
    temp_hamil_byT(2,k2), temp_hamil_byT(3,k2), temp_hamil_byT(4,k2)
enddo

close(336)
deallocate(sort_indx, sort_pointer,sort_isos, temp_hamil, temp_hamil_byT, &
           temp_abcd, red_indx , sort_red_pointer, red_abcd)
print "(A)", "[OK] EXPORT Hamiltonian (uncoupled) for current interaction."
end subroutine print_uncoupled_hamiltonian_H2


!------------------------------------------------------------------------------!
! subroutine calculate_rearrang_field_explicit                                 !
!                                                                              !
! Evaluate the rearrangement field for non-null DD matrix elements just after  !
! these elements are evaluated, completing the Rearrangement Field one element !
! at a time.                                                                   !
!------------------------------------------------------------------------------!
subroutine calculate_rearrang_field_explicit(ia, ib, ic, id, Vdec)

integer, intent(in)   :: ia, ib, ic, id
real(r64), intent(in) :: Vdec
complex(r64) :: aux_rea, Fcomb, aux_reaN
                !f0, f1, f2, f3, & ! field_abcd, f_abdc, f_bacd, f_badc
                !f4, f5, f6, f7, & ! field_cdab, f_cdba, f_dcab, f_dcba
complex(r64), dimension(0:7) :: f
integer   :: HOspo2, it, k1, k2, a, b, c, d, k1N, k2N
integer   :: perm
real(r64) :: rea_h, f2r, sign_tr, rea_hN
logical :: exist_

HOspo2 = HOsp_dim / 2
perm = 1
if (.NOT.evalFullSPSpace) perm = step_reconstruct_2body_timerev(ia, ib, ic, id)
! copy indexes, find_timerev uses them as INOUT
a = ia
b = ib
c = ic
d = id

do it = 1, 2

  sign_tr = one
  if ( it == 2 ) then
    if (evalFullSPSpace) cycle
    if ( HOsp_2mj(a) + HOsp_2mj(b) == 0 ) cycle
!    write(122, fmt='(5I3,A)', advance='no') a, b, c, d, perm, " to "
    call find_timerev(perm, a, b, c, d)
    sign_tr = sign(one, perm*one)
!    write(122, fmt='(5I3)') a, b, c, d, perm

    sign_tr = sign_tr &
               * ((-1)**MOD(HOsp_l(a)+HOsp_l(b)+HOsp_l(c)+HOsp_l(d), 2))
    sign_tr = sign_tr * ((-1)**((HOsp_2mj(a) + HOsp_2mj(b) &
                                 - HOsp_2mj(c) - HOsp_2mj(d))/2))
  endif

  ! fields have to be calculated in the loop, find_timerev changes a, b, c and d
  !(Note)! For Seeds 2 and 3 kappaLR = kappaRL and do not have imaginary part.
  !! Even for seeds general with pn part has non imaginary part.
  f = zzero

  !!! Calculation of Gamma
f(0) = rhoLR(d,b)*rhoLR(c,a) - rhoLR(c,b)*rhoLR(d,a) + kappaLR(c,d)*kappaRL(a,b)

if (.NOT.evalFullSPSpace) then
 if ( ic /= id ) then
  f(1)=rhoLR(c,b)*rhoLR(d,a) - rhoLR(d,b)*rhoLR(c,a) + kappaLR(d,c)*kappaRL(a,b)
 endif
 if ( ib /= ia ) then
  f(2)=rhoLR(d,a)*rhoLR(c,b) - rhoLR(c,a)*rhoLR(d,b) + kappaLR(c,d)*kappaRL(b,a)
 endif
 if ( (ib /= ia) .and. (ic /= id) ) then
  f(3)=rhoLR(c,a)*rhoLR(d,b) - rhoLR(d,a)*rhoLR(c,b) + kappaLR(d,c)*kappaRL(b,a)
 endif

 if ( (ia /= ic) .or. (ib /= id) ) then
  !! Gamma field for bra-ket change <cd | ab>
  f(4)=rhoLR(b,d)*rhoLR(a,c) - rhoLR(a,d)*rhoLR(b,c) + kappaLR(a,b)*kappaRL(c,d)
 if (ia /= ib) then
  f(5)=rhoLR(a,d)*rhoLR(b,c) - rhoLR(b,d)*rhoLR(a,c) + kappaLR(b,a)*kappaRL(c,d)
 endif
 if (ic /= id) then
  f(6)=rhoLR(b,c)*rhoLR(a,d) - rhoLR(a,c)*rhoLR(b,d) + kappaLR(a,b)*kappaRL(d,c)
 endif
 if ( (ib /= ia) .and. (ic /= id)) then
  f(7)=rhoLR(a,c)*rhoLR(b,d) - rhoLR(b,c)*rhoLR(a,d) + kappaLR(b,a)*kappaRL(d,c)
 endif
 endif
endif ! CASE EvalFullSpace

Fcomb = (f(0) - f(1) - f(2) + f(3)) + (f(4) - f(5) - f(6) + f(7))


!  f0 = rhoLR(d,b)*rhoLR(c,a) - rhoLR(c,b)*rhoLR(d,a) + kappaLR(c,d)*kappaRL(a,b)
!  f2 = rhoLR(d,a)*rhoLR(c,b) - rhoLR(c,a)*rhoLR(d,b) + kappaLR(c,d)*kappaRL(b,a)
!  f1 = rhoLR(c,b)*rhoLR(d,a) - rhoLR(d,b)*rhoLR(c,a) + kappaLR(d,c)*kappaRL(a,b)
!  f3 = rhoLR(c,a)*rhoLR(d,b) - rhoLR(d,a)*rhoLR(c,b) + kappaLR(d,c)*kappaRL(b,a)
!
!  Fcomb = (f0 - f1 - f2 + f3)
!  if((kdelta(a,c)== 0).OR.(kdelta(b,d)==0)) then
!  f4 = rhoLR(b,d)*rhoLR(a,c) - rhoLR(a,d)*rhoLR(b,c) + kappaLR(a,b)*kappaRL(c,d)
!  f6 = rhoLR(b,c)*rhoLR(a,d) - rhoLR(a,c)*rhoLR(b,d) + kappaLR(a,b)*kappaRL(d,c)
!  f5 = rhoLR(a,d)*rhoLR(b,c) - rhoLR(b,d)*rhoLR(a,c) + kappaLR(b,a)*kappaRL(c,d)
!  f7 = rhoLR(a,c)*rhoLR(b,d) - rhoLR(b,c)*rhoLR(a,d) + kappaLR(b,a)*kappaRL(d,c)
!
!  Fcomb = Fcomb + (f4 - f5 - f6 + f7)
!  end if

  do k1 = 1, HOspo2
    k1N = k1 + HOspo2
    do k2 = 1, HOspo2
      k2N = k2 + HOspo2

      rea_h  = sign_tr * rearrangement_me(k1,k2)
      rea_hN = sign_tr * rearrangement_me(k1N,k2N)

      !!! Faster than using if ((a /= c).or.(b /= d))
      !!! Calculation of Gamma
      aux_rea  = rea_h  * Fcomb
      aux_reaN = rea_hN * Fcomb

      rearrang_field(k1,k2)   = rearrang_field(k1,k2)   + aux_rea
      rearrang_field(k1N,k2N) = rearrang_field(k1N,k2N) + aux_reaN

    enddo
  enddo


enddo

end subroutine calculate_rearrang_field_explicit

!------------------------------------------------------------------------------!
! subroutine calculate_fields_DD_explicitly                                    !
!                                                                              !
! Calculates the HFB fields h, Gamma an Delta which are then used to compute   !
! other quantities of interest (in particular the energy).                     !
!                                                                              !
! What is calculated:                                                          !
!   h^LR_{ac} = Gamma^LR_{ac} + t_{ac}                                         !
!   Gamma^LR_{ac} = sum_{bd} V_{abcd} rho^LR_{db}                              !
!   Delta^LR_{ab} = 1/2 sum_{cd}  V_{abcd} kappa^LR_{cd}                       !
!                 =     sum_{c<d} V_{abcd} kappa^LR_{cd}                       !
!   Delta^RL_{ab} = 1/2 sum_{cd}  V_{abcd} kappa^RL_{cd}                       !
!                 =     sum_{c<d} V_{abcd} kappa^RL_{cd}                       !
!                                                                              !
! The fields have the symmetry properties:                                     !
!   Delta^LR_{ab} = - Delta^LR_{ba}    (Skew symmetry)                         !
!   Delta^RL_{ab} = - Delta^RL_{ba}    (Skew symmetry)                         !
!                                                                              !
! The actual calculation below uses these symmetries and the fact that only    !
! a subset of the matrix elements of the interaction are stored.               !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        rhoLR,kappaLR,kappaRL = transition densities                          !
!                                                                              !
! Output: gammaLR,deltaLR,deltaLR = transition fields                          !
!------------------------------------------------------------------------------!
subroutine calculate_fields_DD_explicit(gammaLR,hspLR, deltaLR,deltaRL,ndim)

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim) :: gammaLR, hspLR, deltaLR, deltaRL

!! The density fields are added to the calculated with the standard hamiltonian
!! This array variables are local
complex(r64), dimension(ndim,ndim) :: gammaLR_DD, deltaLR_DD, deltaRL_DD
complex(r64) :: t_gam_sum, t_rea_sum, del_aux
integer :: i, j, ia, ib, ic, id, it, A_print, C_print, BI
integer :: perm
integer(i64) :: kk
real(r64) :: h2b, f2b, aux
character(len=10) :: g_str
character(len=70) :: mestr

!cmpi integer :: ierr=0
!cmpi complex(r64), dimension(ndim,ndim) :: gammaLR_red, deltaLR_red, &
!cmpi                                       deltaRL_red
!print *, "[  ] calculating Fields DD"

gammaLR_DD = zzero
deltaLR_DD = zzero
deltaRL_DD = zzero
!
!$OMP PARALLEL DO FIRSTPRIVATE(rhoLR,kappaLR,kappaRL) &
!$OMP             PRIVATE(ia,ib,ic,id,h2b,f2b,perm,it) &
!$OMP             REDUCTION(+:gammaLR,deltaLR,deltaRL)

A_print = 1 !5
C_print = 1 !8
write(g_str, fmt="(A,I2,I3,A)") "GM(", A_print, C_print, ")"
print "(A)", " BI = 0 Conserv, 1 parity break, 2 M break, 3 parity and M break"
print "(A)", "GAMMA_ac   TR BI [CASE]     hfb from  < a  b  c  d > (pn pn)"
print "(A)", "---------------------------------------------------------------"


!open(123, file='TRtransf_states_DD.gut')
do kk = 1, hamil_DD_H2dim
  ia = hamil_DD_abcd(1+4*(kk-1))
  ib = hamil_DD_abcd(2+4*(kk-1))
  ic = hamil_DD_abcd(3+4*(kk-1))
  id = hamil_DD_abcd(4+4*(kk-1))
  h2b  = hamil_DD_H2(kk)
  perm = hamil_DD_trperm(kk)

  !!! Loop on time reversal
  do it = 1, 2
    if ( it == 2 ) then
      if (evalFullSPSpace) cycle

      if ( HOsp_2mj(ia) + HOsp_2mj(ib) == 0 ) cycle
!      write(123, fmt='(5I3,A)', advance='no') ia, ib, ic, id, perm, " to "
      call find_timerev(perm,ia,ib,ic,id)
!      write(123, fmt='(5I3)') ia, ib, ic, id, perm
      h2b = sign(one,perm*one) * h2b
      h2b = h2b * ((-1)**((HOsp_l(ia)+HOsp_l(ib)+HOsp_l(ic)+HOsp_l(id))))
      h2b = h2b * ((-1)**((HOsp_2mj(ia) + HOsp_2mj(ib) &
                           - HOsp_2mj(ic) - HOsp_2mj(id)) / 2))
    endif

    BI = 0   !! 0 Conserving, 1 parity break, 2 M break, 3 parity and M break
    if (MOD(HOsp_l(ia)+HOsp_l(ib)+HOsp_l(ic)+HOsp_l(id), 2)==1) BI = BI + 1
    if (HOsp_2mj(ia)+HOsp_2mj(ib) /= HOsp_2mj(ic)+HOsp_2mj(id)) BI = BI + 2
    write(mestr, fmt="(4I3,2A,4(I5,A,I2,A),A,I3)") ia,ib,ic,id, " > =", &
      " ",HOsp_ant(ia),"(",HOsp_2mj(ia),")",HOsp_ant(ib),"(",HOsp_2mj(ib),&
      ")",HOsp_ant(ic),"(",HOsp_2mj(ic),")",HOsp_ant(id),"(",HOsp_2mj(id),")",&
      " perm= ", perm*(abs(1-it) + 10*(2-it))

    !!! Faster than using if ((a /= c).or.(b /= d))
    f2b = h2b * (1 - (kdelta(ia,ic) * kdelta(ib,id)))

    !!! Calculation of Gamma
    gammaLR_DD(ia,ic) = gammaLR_DD(ia,ic) + h2b * rhoLR(id,ib) ! ab cd
    !!! Calculation of Delta^10 and Delta^01
    deltaLR_DD(ib,ia) = deltaLR_DD(ib,ia) + h2b * kappaLR(id,ic)
    deltaRL_DD(ib,ia) = deltaRL_DD(ib,ia) + h2b * kappaRL(id,ic)

    if (evalFullSPSpace) cycle
    !!! Calculation of Gamma
    gammaLR_DD(ia,id) = gammaLR_DD(ia,id) - h2b * rhoLR(ic,ib) ! ab dc
    gammaLR_DD(ib,ic) = gammaLR_DD(ib,ic) - h2b * rhoLR(id,ia) ! ba cd
    gammaLR_DD(ib,id) = gammaLR_DD(ib,id) + h2b * rhoLR(ic,ia) ! ba dc

    gammaLR_DD(ic,ia) = gammaLR_DD(ic,ia) + f2b * rhoLR(ib,id) ! cd ab
    gammaLR_DD(ic,ib) = gammaLR_DD(ic,ib) - f2b * rhoLR(ia,id) ! cd ba
    gammaLR_DD(id,ia) = gammaLR_DD(id,ia) - f2b * rhoLR(ib,ic) ! dc ab
    gammaLR_DD(id,ib) = gammaLR_DD(id,ib) + f2b * rhoLR(ia,ic) ! dc ba

    !!! Calculation of Delta^10 and Delta^01
    deltaLR_DD(id,ic) = deltaLR_DD(id,ic) + f2b * kappaLR(ib,ia)
    deltaRL_DD(id,ic) = deltaRL_DD(id,ic) + f2b * kappaRL(ib,ia)

  enddo
enddo
!close(123)
!$OMP END PARALLEL DO
!!! Reduces the values for the processes in the same team
!cmpi if ( paral_myteamsize > 1 ) then
!cmpi   call mpi_reduce(gammaLR,gammaLR_red,ndim**2,mpi_double_complex, &
!cmpi                   mpi_sum,0,mpi_comm_team,ierr)
!cmpi   call mpi_reduce(deltaLR,deltaLR_red,ndim**2,mpi_double_complex, &
!cmpi                   mpi_sum,0,mpi_comm_team,ierr)
!cmpi   call mpi_reduce(deltaRL,deltaRL_red,ndim**2,mpi_double_complex, &
!cmpi                   mpi_sum,0,mpi_comm_team,ierr)
!cmpi   gammaLR = gammaLR_red
!cmpi   deltaLR = deltaLR_red
!cmpi   deltaRL = deltaRL_red
!cmpi endif


print "(3A,F15.12)", "------- Total ", g_str, " = ", gammaLR_DD(A_print,C_print)
!!! Skew symmetry ( This goes after the following loop block )
do j = 1, HOsp_dim
  do i = 1, j-1
    deltaLR_DD(i,j) = -1.0d0 * deltaLR_DD(j,i)
    deltaRL_DD(i,j) = -1.0d0 * deltaRL_DD(j,i)
  enddo
enddo

!!! h = Gamma + 1body (There is no Hamil_1 for DD))
open (620, file='fields_matrix_explicit.gut')
write(620,fmt='(A,A,A)') "  i   j        gammaLR       gammaLR_DD        ",&
  "DeltaLR        DeltaLR_DD       DeltaRL        DeltaRL_DD     ReaField_DD ",&
  "        hspLR       hspLR_dd     hspLR_rea"

t_rea_sum = zzero
t_gam_sum = zzero
do i = 1, HOsp_dim
  do j = 1, HOsp_dim

    write(620, fmt='(2I4,8F16.9)', advance='no') i,j, real(gammaLR(i,j)),&
        real(gammaLR_DD(i,j)), real(deltaLR(i,j)), real(deltaLR_DD(i,j)),&
        real(deltaRL(i,j)),real(deltaRL_DD(i,j)), real(rearrang_field(i,j)), &
        real(hspLR(i,j))

    !! introduced to update the fields with DD contributions
    !
    gammaLR(i,j) = gammaLR(i,j) + gammaLR_DD(i,j)
    deltaLR(i,j) = deltaLR(i,j) + deltaLR_DD(i,j)
    deltaRL(i,j) = deltaRL(i,j) + deltaRL_DD(i,j)
    !
    !! end
    if (EVAL_REARRANGEMENT) then
        hspLR(i,j) = hspLR(i,j) + gammaLR_DD(i,j)
        write(620, fmt='(A,F12.6)', advance='no') ' ',real(hspLR(i,j))
        hspLR(i,j) = hspLR(i,j) + (0.25d+00 * rearrang_field(i,j))
        write(620, fmt='(A,F12.6)') ' ',real(hspLR(i,j))
    else
        hspLR(i,j) = hspLR(i,j) + gammaLR_DD(i,j)
        write(620, fmt='(A,F12.6)') ' ',real(hspLR(i,j))
    endif
    !+ gamma only from DD, the other gamma was added in the routine before

    t_rea_sum = t_rea_sum + (rhoLR(i,j) * rearrang_field(j,i))
    t_gam_sum = t_gam_sum + (rhoLR(i,j) * gammaLR_DD(j,i))
    t_gam_sum = t_gam_sum - 0.5d00*(kappaRL(i,j) * deltaLR_DD(j,i))
  enddo
enddo
!!
!t_rea_sum = (0.25d00 / alpha_DD) * t_rea_sum
del_aux = (t_gam_sum / t_rea_sum) - (2.0d00 / alpha_DD)
if (dreal(del_aux)**2 + dimag(del_aux)**2 > 1.0D-15) then
  print '(A,2F20.15,A,F20.15)', "!ERROR Tr(rho*Gamma)!prop Tr(rho*Rea) ::", &
    dreal(t_gam_sum),dreal(t_rea_sum), &
    ":: rhoGamm/rhoRea ::", dreal(t_gam_sum / t_rea_sum)
endif


close(620)
!print *, "[OK] calculating Fields DD"

iteration = iteration + 1

end subroutine calculate_fields_DD_explicit


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

            if ((i_r == R_PRINT).AND.(i_ang == ANG_PRINT)) then
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

        if ((PRNT_).AND.(i_r == R_PRINT).AND.(i_ang == ANG_PRINT)) then
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
! Generalization of fields for evaluating several DD terms.
! (iterates and fix the constants)
!    Call from module_projections.
!-----------------------------------------------------------------------------!
subroutine compute_fields_DD_selector(gammaLR,hspLR,deltaLR,deltaRL,ndim)
  integer, intent(in) :: ndim
  complex(r64), dimension(ndim,ndim) :: gammaLR, hspLR, deltaLR, deltaRL
  !! The density fields are added to the calculated with the standard hamiltonian
  !! This array variables are local
  complex(r64), dimension(ndim,ndim) :: gammaLR_DD, deltaLR_DD, deltaRL_DD


  print *, " [WARNING] Unimplemented method, calling [calculate_fields_DD]"
  call calculate_fields_DD(gammaLR,hspLR,deltaLR,deltaRL,ndim)
end subroutine compute_fields_DD_selector

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
! subroutine set_Radial1b_derivates                                           !
!                                                                             !
! Set up an array for the derivatives in terms of l, n  for all the components!
! matrix: n-1 (l-1, l, l+1), n (l-1, l, l+1), n+1 (l-1, l, l+1)               !
!-----------------------------------------------------------------------------!
subroutine set_Radial1b_derivates
integer   :: a_sh, n, l, n2, l2, i_l, i_n, i_r
real(r64) :: radial

allocate(radial_1b_diff_memo(HOsh_dim, -1:1, -1:1, r_dim))
radial_1b_diff_memo = zero

do a_sh = 1, HOsh_dim
  n = HOsh_n(a_sh)
  l = HOsh_l(a_sh)

  do i_n = -1, 1
    n2 = n + i_n
    if (n2.LT.0) cycle
    do i_l = -1, 1
      l2 = l + i_l
      if (l2.LT.0) cycle

      do i_r = 1, r_dim
        radial = radial_function(n2, l2, r(i_r))  ! r = b * sqrt(x/(alpha+2))
                                                  ! already 1/b**3
        radial_1b_diff_memo(a_sh, i_n, i_l, i_r) = radial
      enddo

    end do
  end do
enddo
end subroutine set_Radial1b_derivates


!-----------------------------------------------------------------------------!
! subroutine calculate_density_laplacian                                      !
!                                                                             !
! Evaluate the gradients Grad(rho (t,ang)) and the scalar product             !
!-----------------------------------------------------------------------------!
subroutine calculate_density_laplacian(dens_rhoRR, dens_kappaRR, ndim)

integer, intent(in) :: ndim
real(r64), dimension(ndim,ndim), intent(in) :: dens_rhoRR, dens_kappaRR

integer   :: a, b, i_r,i_an, na,la,ja,mja, nb,lb,jb,mjb, a_sh, b_sh
integer   :: mla, mlb, ms, K1, M1, K2, M2, mu_, ADK2, indxa, t

real(r64) :: aux1, aux2, aux3, cgc1, cgc2, cgc3, g_kl, xikl, rad, dd_prod
real(r64), dimension(:), allocatable :: rad_diffs
complex(r64), dimension(:,:), allocatable :: rea_dens
complex(r64), dimension(:,:,:), allocatable :: prea_dir, prea_exc

complex(r64) :: ang
integer   :: ma,mb, c,lc,jc,mc, d,ld,jd,md, indx_a,indx_b,indx_c,indx_d, &
             indx_km1,indx_km2
!! Angular part is a Y_KM, up to l_max+1 (sph_harmonics_memo is up to 2*l_max)

allocate(rad_diffs(r_dim))
allocate(partial_dens(-1:2,r_dim,angular_dim))
partial_dens = zzero
!!
do a = 1, HOsp_dim
  a_sh = HOsp_sh(a)
  la = HOsp_l(a)
  na = HOsp_n(a)
  ja = HOsp_2j(a)
  mja = HOsp_2mj(a)

  do b = 1, HOsp_dim
    b_sh = HOsp_sh(b)
    lb = HOsp_l(b)
    nb = HOsp_n(b)
    jb = HOsp_2j(b)
    mjb = HOsp_2mj(b)
    !! evaluate the radial parts for the last step
    rad_diffs = zero
    do i_r = 1, r_dim
      if (na .GT. 0) then
        rad_diffs(i_r) = rad_diffs(i_r) + (&
                                     radial_1b_diff_memo(a_sh,-1,+1,i_r) * &
                                     radial_1b_diff_memo(b_sh, 0, 0,i_r) / &
                                     sqrt(na + 0.0d0))
      endif
      if (nb .GT. 0) then
        rad_diffs(i_r) = rad_diffs(i_r) + (&
                                     radial_1b_diff_memo(a_sh, 0, 0,i_r) * &
                                     radial_1b_diff_memo(b_sh,-1,+1,i_r) / &
                                     sqrt(nb + 0.0d0))
      endif
    enddo

    !! sumatory over the angular indexes
    do ms = -1, 1, 2
      mla = mja - ms
      mlb = mjb - ms

      if ((abs(mla / 2).GT.la) .OR. (abs(mlb / 2).GT.lb)) cycle

      call ClebschGordan(2*la, 1,ja, mla,ms,mja, cgc1)
      call ClebschGordan(2*lb, 1,jb, mlb,ms,mjb, cgc2)

      aux1 = ((-1)**(mla/2)) * dsqrt(((2*la) + 1)*((2*lb) + 1) / (4*pi))
      aux1 = aux1 * cgc1 * cgc2 * dens_rhoRR(a, b)

!      print "(A,3I4,F11.6)", "  1* ms,mla,b:",ms, mla,mlb, aux1

      M1 = (mjb - mja) / 2
      do K1 = abs(ja - jb)/2, (ja + jb)/2
        !! the steps of K1 have to be even
        if (abs(M1).GT.K1) cycle
        if (MOD(la+lb+K1, 2) == 1) cycle

        call ClebschGordan(2*la,2*lb,2*K1, 0,0,0, cgc1)
        call ClebschGordan(2*la,2*lb,2*K1, -mla,mlb,2*M1, cgc2)

        aux2 = cgc1 * cgc2 / dsqrt((2*K1) + 1.0d0)

        do mu_ = -1, 1, 1
          !! Components of the gradient
          M2 = (mjb - mja + 2*mu_) / 2
          do ADK2 = -1,  1, 2
            K2 = K1 + ADK2
!            print "(A,6I4,F11.6)", "  2* KM(1,2),mu,ad:",K1,M1, K2,M2, mu_,&
!                                  ADK2, aux2
            if (K2.LT.0) cycle
            if (abs(M2).GT.K2) cycle

            indxa = angular_momentum_index(K2,M2,.FALSE.)
            call ClebschGordan(2*K1,2,2*K2, 2*M1,2*mu_,2*M2, cgc3)

            !! g(K1,K2) coeff
            if (ADK2 .EQ. +1) then
              g_kl = dsqrt((K1 + 1.0d0) / ((2*K1) + 3.0d0))
              xikl = - K1
            else !! (ADK2 .EQ. -1)
              g_kl = dsqrt( K1 / ((2*K1) - 1.0d0))
              xikl = K1 + 1.0d0
            end if

            aux3 = aux1 * aux2 * g_kl * cgc3
            if (dabs(aux3) .LT. 1.0d-9) cycle

!            print "(A,I4,F11.6)", "ACCEPT  2.2 ", indxa,aux3
            do i_an = 1, angular_dim
              do i_r = 1, r_dim
                !! radial  2b functions and diff parts precalculated.
                rad = ((xikl+ (la+lb)) * HO_b / r(i_r)) - (2.0d0 * r(i_r)/HO_b)
                rad = rad * radial_2b_sho_memo(a_sh, b_sh, i_r) / HO_b
                rad = rad + (rad_diffs(i_r) / HO_b)

                partial_dens(mu_,i_r,i_an) = partial_dens(mu_,i_r,i_an) + &
                                  aux3 * rad * sph_harmonics_memo(indxa, i_an)
              enddo
            enddo !loop angular-radial

          enddo
        end do !loop mu components
      enddo

    enddo !ms

    !!
  enddo
enddo


!! test to evaluate the rearrangement density
allocate(rea_dens(r_dim, angular_dim), &
         prea_dir(2,r_dim, angular_dim), prea_exc(2,r_dim, angular_dim))
rea_dens = zzero
prea_dir = zzero
prea_exc = zzero
do a = 1, HOsp_dim
  ja = HOsp_2j(a)
  la = HOsp_l(a)
  ma = HOsp_2mj(a)
  indx_a = angular_momentum_index(ja, ma, .TRUE.)
  do b = 1, HOsp_dim
    jb = HOsp_2j(b)
    lb = HOsp_l(b)
    mb = HOsp_2mj(b)
    indx_b = angular_momentum_index(jb, mb, .TRUE.)
    do c = 1, HOsp_dim
      jc = HOsp_2j(c)
      lc = HOsp_l(c)
      mc = HOsp_2mj(c)
      indx_c = angular_momentum_index(jc, mc, .TRUE.)
      if (dabs(dens_rhoRR(c,a)) .LT. 1.0e-6) cycle

      do d = 1, HOsp_dim
        jd = HOsp_2j(d)
        ld = HOsp_l(d)
        md = HOsp_2mj(d)
        indx_d = angular_momentum_index(jd, md, .TRUE.)
        if (dabs(dens_rhoRR(d,b)) .LT. 1.0e-6) cycle

dd_prod = (2 * dens_rhoRR(d,b) * dens_rhoRR(c,a)) &
            - (dens_kappaRR(c,d) * dens_kappaRR(a,b))
!!! -------------------------------------------------------------------------
do K1 = max(0, abs(ja - jc) / 2), (ja + jc) / 2
  M1 = (mc - ma) /2 !! (Suhonen_Vasr)
  if (abs(M1) > K1) cycle
  if (MOD(la + lc + K1, 2) == 1) cycle

  indx_km1 = angular_momentum_index(K1, M1, .FALSE.)
  cgc1 = dens_Y_KM_me(indx_a,indx_c,indx_km1)
  if (dabs(cgc1) .LT. 1.0e-6) cycle
  do K2 = max(0, abs(jb - jd) / 2), (jb + jd) / 2
    M2 = (md - mb) /2 !! (Suhonen_Vasr)
    if (abs(M2) > K2) cycle
    if (MOD(lb + ld + K2, 2) == 1) cycle

    indx_km2 = angular_momentum_index(K2, M2, .FALSE.)
    cgc2 = dens_Y_KM_me(indx_b,indx_d,indx_km2)
    if (dabs(cgc2) .LT. 1.0e-6) cycle
do i_r = 1, r_dim
  rad = radial_2b_sho_memo(HOsp_sh(a), HOsp_sh(c), i_r) * &
        radial_2b_sho_memo(HOsp_sh(b), HOsp_sh(d), i_r)
  do i_an = 1, angular_dim
    ang = cgc1 * cgc2 * sph_harmonics_memo(indx_km1,i_an) * &
                        sph_harmonics_memo(indx_km2,i_an)
    rea_dens(i_r, i_an) = rea_dens(i_r, i_an) + (rad*ang*dd_prod)
  enddo
enddo !! radial - angular loop
enddo
enddo

!!! -- EXCHANGE LOOP --------------------------------------------------------

do K1 = max(0, abs(ja - jd) / 2), (ja + jd) / 2
  M1 = (md - ma) /2 !! (Suhonen_Vasr)
  if (abs(M1) > K1) cycle
  if (MOD(la + ld + K1, 2) == 1) cycle

  indx_km1 = angular_momentum_index(K1, M1, .FALSE.)
  cgc1 = dens_Y_KM_me(indx_a,indx_d,indx_km1)
  if (dabs(cgc1) .LT. 1.0e-6) cycle
  do K2 = max(0, abs(jb - jc) / 2), (jb + jc) / 2
    M2 = (mc - mb) /2 !! (Suhonen_Vasr)
    if (abs(M2) > K2) cycle
    if (MOD(lb + lc + K2, 2) == 1) cycle

    indx_km2 = angular_momentum_index(K2, M2, .FALSE.)
    cgc2 = dens_Y_KM_me(indx_b,indx_c,indx_km2)
    if (dabs(cgc2) .LT. 1.0e-6) cycle
do i_r = 1, r_dim
  rad = radial_2b_sho_memo(HOsp_sh(a), HOsp_sh(c), i_r) * &
        radial_2b_sho_memo(HOsp_sh(b), HOsp_sh(d), i_r)
  do i_an = 1, angular_dim
    ang = cgc1 * cgc2 * sph_harmonics_memo(indx_km1,i_an) * &
                        sph_harmonics_memo(indx_km2,i_an)
    rea_dens(i_r, i_an) = rea_dens(i_r, i_an) + (x0_DD_FACTOR*rad*ang*dd_prod)
  enddo
enddo !! radial - angular loop
enddo
enddo

!!! -------------------------------------------------------------------------
      end do
    end do
  end do
end do
!allocate(REACommonFields(r_dim, angular_dim))
!allocate(BulkP1(5,4,r_dim, angular_dim), BulkP2(5,4,r_dim, angular_dim))
!BulkP1 = zzero
!BulkP2 = zzero

call calculate_common_rearrang_bulkFields

!deallocate(BulkP1, BulkP2)


!! Evaluate the laplacian-derived density associated to the direct and exchange
!! functions from definition
do i_r = 1, r_dim
  do i_an = 1, angular_dim
    do mu_ = -1, 1
      partial_dens(2,i_r,i_an) = partial_dens(  2,i_r,i_an) + &
        ((-1)**mu_) * partial_dens(mu_,i_r,i_an) * partial_dens(-mu_,i_r,i_an)
    enddo
    do t=1, 2
      prea_dir(t,i_r,i_an) = prea_dir(t,i_r,i_an) + ( &
                  (dens_pnt(5,i_r,i_an) - x0_DD_FACTOR*dens_pnt(t,i_r,i_an)))
      !! Equivalent to the effect of all exchange bulk densities (does not
      !! contribute equally for each 2-body function, it depends on msms')
      do ms = 1, 4
        prea_exc(t,i_r,i_an) = prea_exc(t,i_r,i_an) - ( &
                  (BulkHF(t,ms,i_r,i_an) - x0_DD_FACTOR*BulkHF(5,ms,i_r,i_an)))
      enddo

    enddo
  enddo
enddo

!! scalar product / export for test
open (111, file='dens_differential.gut')
write(111, fmt="(A)") "  i_r i_an r(ir)    grad_den_-1     imag(grad_-1)     &
        &grad_den_0      imag(grad_0)     grad_den_+1    imag(grad_+1)    &
        &R(Laplacian)   sqrt(R(Laplac))       R(dens)   R(dens_alpha)     &
        &rea_dens *     common_rea      pseudorea_d_dir(p/n)            &
        &pseudorea_d_exch(p/n)"
do i_r = 1, r_dim
  do i_an = 1, angular_dim

    write(111,fmt='(2(I4,A),F5.2)',advance='no') i_r, ",", i_an, ",", r(i_r)
    do mu_ = -1, 1
      write(111,fmt='(A,F15.9,A,F15.9,A)',advance='no') ",", &
    dreal(partial_dens(mu_,i_r,i_an)), " ",dimag(partial_dens(mu_,i_r,i_an)),"j"
    enddo
    if (dabs(dimag(partial_dens(2,i_r,i_an))).GT.1.0e-9) then
      print "(A,2I4,D15.9,A)","[ERROR IMAG] grad diff imag > 1e-9:",&
                              i_r, i_an, partial_dens(2,i_r,i_an)
    endif

    write(111,fmt='(A,F15.9,A,F15.9)', advance='no') ",", &
      dreal(partial_dens(2,i_r,i_an)), ",",&
      dreal(partial_dens(2,i_r,i_an))**0.5d0
    write(111,fmt='(A,F15.9,A,F15.9)', advance='no') ",", &
      dreal(dens_pnt(5,i_r,i_an)), ", ", dreal(dens_alpha(i_r, i_an))
    !!! export the test for the rea_density
    write(111,fmt='(6(A,F15.9))') &
      ",", dreal(rea_dens(i_r,i_an)), " ", dreal(REACommonFields(i_r,i_an)),&
      ",", dreal(prea_dir(1,i_r, i_an)), " ", dreal(prea_dir(2,i_r, i_an)),&
      ",", dreal(prea_exc(1,i_r, i_an)), " ", dreal(prea_exc(2,i_r, i_an))
  end do
end do
close(111)
!deallocate(REACommonFields)

end subroutine calculate_density_laplacian


!-----------------------------------------------------------------------------!
! subroutine set_derivative_density_dependent                                 !
!                                                                             !
! Set up everything for the laplacian (Grad rho(r,ang))^2                     !
!-----------------------------------------------------------------------------!
subroutine set_derivative_density_dependent(dens_rhoRR, dens_kappaRR, ndim)

integer, intent(in) :: ndim
real(r64), dimension(ndim,ndim), intent(in) :: dens_rhoRR, dens_kappaRR

! 1. Set up the radial functions for cross n,l Radial functions
call set_Radial1b_derivates
print "(A)", " [DONE] Setting Radial 1b."
! 2. Calculate gradient of the last density
call calculate_density_laplacian(dens_rhoRR, dens_kappaRR, ndim)
print "(A)", " [DONE] Calculated Laplacian of the density."

end subroutine set_derivative_density_dependent

!-----------------------------------------------------------------------------!
! function matrix_element_v_gradientDD                                        !
!                                                                             !
! Computes density dependent two body matrix elements over the density average!
!    all_isos (logical) Compute 3 combinations p/n instead of the current     !
!                       ta,tb,tc,td of the sp-state.                          !
!                       v_dd_val_Real !! pppp(1), pnpn(2), pnnp(3), nnnn(4)
!             in this case, a,b,c,d are <= HOsp_dim / 2
!-----------------------------------------------------------------------------!
function matrix_element_v_gradientDD(a,b, c,d) result (v_dd_val_Real)

integer(i32), intent(in) :: a,b,c,d
real(r64), dimension(4) :: v_dd_val_Real !! pppp(1), pnpn(2), pnnp(3), nnnn(4)

integer      :: K,K2,M,M2, ind_km, ind_km_q, ind_jm_a,ind_jm_b,ind_jm_c, &
                ind_jm_d,la,lb,lc,ld, ja,jb,jc,jd, ma,mb,mc,md
complex(r64) :: aux, radial, aux_dir, aux_exch, dens_part
complex(r64), dimension(4) :: v_dd_value
real(r64)    :: angular, integral_factor
integer(i32) :: a_sh, b_sh, c_sh, d_sh, i_r, i_ang
integer      :: print_element = 0, HOspO2, skpd

HOspO2 = HOsp_dim/2

v_dd_value = zzero
v_dd_val_Real = zero

if (.NOT.EXPORT_GRAD_DD) return
if ((a.GT.HOspO2).OR.(b.GT.HOspO2).OR.(c.GT.HOspO2).OR.(d.GT.HOspO2)) then
  print "(A)", " [ERROR] (m.e. gradient), a,b,c,d > HO_sp dim /2 !!!"
  return
end if

ja = HOsp_2j(a)
ma = HOsp_2mj(a)
la = HOsp_l(a)
a_sh = HOsp_sh(a)
ind_jm_a = angular_momentum_index(ja, ma, .TRUE.)
jb = HOsp_2j(b)
mb = HOsp_2mj(b)
lb = HOsp_l(b)
b_sh = HOsp_sh(b)
ind_jm_b = angular_momentum_index(jb, mb, .TRUE.)
jc = HOsp_2j(c)
mc = HOsp_2mj(c)
lc = HOsp_l(c)
c_sh = HOsp_sh(c)
ind_jm_c = angular_momentum_index(jc, mc, .TRUE.)
jd = HOsp_2j(d)
md = HOsp_2mj(d)
ld = HOsp_l(d)
d_sh = HOsp_sh(d)
ind_jm_d = angular_momentum_index(jd, md, .TRUE.)

integral_factor = 2.0d+0 * alpha_DD * t3_DD_CONST
!! NOTE :: Remember that radial functions already have the factor 1/b**3
integral_factor = integral_factor * 0.5d0 * (HO_b**3)
integral_factor = integral_factor  / ((2.0d0 + alpha_DD)**1.5d0)
integral_factor = integral_factor * 4 * pi  ! add Lebedev norm factor

do i_r = 1, r_dim

  radial = weight_R(i_r) * radial_2b_sho_memo(a_sh, c_sh, i_r) &
                         * radial_2b_sho_memo(b_sh, d_sh, i_r) &
                         * exp((alpha_DD + 2.0d0) * (r(i_r) / HO_b)**2)
  !! NOTE: the inclusion of the exponential part included for the same reason as
  !!       in the DD and rearrangement matrix elements subroutine, (see there)
  !!                    * exp((alpha_DD + 2.0d0 + 2.0d0) * (r(i_r) / HO_b)**2)
  !! NOTE 2: Last expression is for the case of only-laplacian_ dependence,
  !!       the one with the same exponential counterpart is due to the
  !!       dimensional readjustment to imitate the rearrangement function:
  !!                  dens**(alp-1) * sqrt(Laplacian_(dens))

  do i_ang = 1, angular_dim
    !! Already deallocated
      aux_dir  = (AngFunctDUAL_HF(1,a,c,i_ang) + AngFunctDUAL_HF(4,a,c,i_ang))&
                *(AngFunctDUAL_HF(1,b,d,i_ang) + AngFunctDUAL_HF(4,b,d,i_ang))
      aux_exch = (AngFunctDUAL_HF(1,a,d,i_ang) + AngFunctDUAL_HF(4,a,d,i_ang))&
                *(AngFunctDUAL_HF(1,b,c,i_ang) + AngFunctDUAL_HF(4,b,c,i_ang))

      angular = weight_LEB(i_ang)
      dens_part = dens_alpm1(i_r,i_ang) * &
                  (dreal(partial_dens(2,i_r,i_ang))**0.5d0)

      !v_nnnn = v_pppp
      aux = radial * angular * (1-x0_DD_FACTOR) * (aux_dir - aux_exch)
      v_dd_value(1) = v_dd_value(1) + (aux * dens_part)
      v_dd_value(4) = v_dd_value(1)
      ! pn pn
      aux = radial * angular * (aux_dir + (x0_DD_FACTOR*aux_exch))
      v_dd_value(2) = v_dd_value(2) + (aux * dens_part)
      ! pn np
      aux = radial * angular * ((x0_DD_FACTOR*aux_dir) + aux_exch)
      v_dd_value(3) = v_dd_value(3) - (aux * dens_part)

  enddo  ! angular iter_
enddo    ! radial  iter_

v_dd_val_Real(1) = real(v_dd_value(1), r64) * integral_factor
v_dd_val_Real(2) = real(v_dd_value(2), r64) * integral_factor
v_dd_val_Real(3) = real(v_dd_value(3), r64) * integral_factor
v_dd_val_Real(4) = real(v_dd_value(4), r64) * integral_factor

if (abs(imag(v_dd_value(2))) > 1.0d-9 ) then
    print "(A,D15.8,A,D20.8)", "[FAIL IMAG] v_grad_DD_abcd is not Real =", &
        real(v_dd_value(2)), " +j ", imag(v_dd_value(2))
endif

!if (dabs(v_dd_val_Real(2)).GT.1.0d-6) then
!  print "(A,4I4,A,4(I6,A,I3,A),A,2D20.9)"," _Eval me(pp/pn):",a,b,c,d, " <",&
!    HOsp_ant(a), "(", HOsp_2mj(a), ")", HOsp_ant(b), "(", HOsp_2mj(b), ")", &
!    HOsp_ant(c), "(", HOsp_2mj(c), ")", HOsp_ant(d), "(", HOsp_2mj(d), ")", &
!    "> =", v_dd_val_Real(1), v_dd_val_Real(2)
!endif
return
end function matrix_element_v_gradientDD



!-----------------------------------------------------------------------------!
! This subroutine evaluates the energy associated to the pseudo-rearrangement !
!-----------------------------------------------------------------------------!
subroutine calculate_energy_field_laplacian(E_core)

real(r64), intent(out) :: E_core
complex(r64), dimension(:,:), allocatable :: psrea_field
integer :: a,c, spO2, ms, i_r, i_a, a_sh, c_sh, aa, cc, t
complex(r64), dimension(2) :: auxD, auxE
complex(r64) :: sumD_a1, sumD_a2
real(r64) :: int_const, rad_ac

spO2 = HOsp_dim / 2

!! compute the field for the pseudo-rearrangement
allocate(psrea_field(HOsp_dim, HOsp_dim))
psrea_field = zzero
int_const = 0.5d0 * (HO_b**3) / ((2.0d0 + alpha_DD)**1.5d0)
int_const = 4.0d0 * pi * int_const
!!
int_const = (2.0d0 * alpha_DD * t3_DD_CONST) * int_const
!!
do a = 1, spO2
  a_sh = HOsp_sh(a)
  do c = 1, spO2
    c_sh = HOsp_sh(c)

    do i_r = 1, r_dim
      rad_ac = weight_R(i_r) * radial_2b_sho_memo(a_sh, c_sh, i_r)
      rad_ac = rad_ac * dexp((2.0d0+alpha_DD) * (r(i_r)/HO_b)**2)
      do i_a = 1, angular_dim
        auxD = zzero
        auxE = zzero

        sumD_a1 = AngFunctDUAL_HF(1,a,c,i_a) + AngFunctDUAL_HF(4,a,c,i_a)
        do t = 1, 2
          auxD(t) = auxD(t) + sumD_a1*(dens_pnt(5, i_r, i_a) - &
                                       x0_DD_FACTOR * dens_pnt(t, i_r, i_a))
          !! exchange
          do ms = 1, 4
            sumD_a2 = AngFunctDUAL_HF(ms, a,c, i_a)
            auxE(t) = auxE(t) + sumD_a2*(BulkHF(t, ms, i_r,i_a) - &
                                         x0_DD_FACTOR * BulkHF(5, ms, i_r,i_a))
          enddo
          aa = a + ((t - 1)*spO2)
          cc = c + ((t - 1)*spO2)
          psrea_field(aa,cc) = psrea_field(aa,cc) + &
            (int_const * weight_LEB(i_a) * rad_ac * dens_alpm1(i_r,i_a) * &
            (dreal(partial_dens(2,i_r,i_a))**0.5d0)* (auxD(t) - auxE(t)))
        end do

      enddo
    enddo

  enddo
enddo

!! Do the trace for the energy
E_core = 0.0
do a = 1, HOsp_dim
  do c = 1, HOsp_dim
    E_core = E_core + dreal(psrea_field(a,c)) * dens_rhoRR(c, a)
  enddo
enddo
deallocate(psrea_field)

E_core = 0.5d0 * E_core  !! the energy should be 1/2 Tr(Gamma * rho)

end subroutine calculate_energy_field_laplacian



!-----------------------------------------------------------------------------!
! function matrix_element_pseudoRearrangement                                 !
!                                                                             !
! Computes density dependent two body matrix elements over the density average!
! The derivation of the matrix element is derived from the Fields derivation  !
!                       v_dd_val_Real !! pppp(1), pnpn(2), pnnp(3), nnnn(4)   !
!             in this case, a,b,c,d are <= HOsp_dim / 2                       !
!-----------------------------------------------------------------------------!
function matrix_element_pseudoRearrangement_v1(a,b, c,d) result (v_dd_val_Real)

integer(i32), intent(in) :: a,b,c,d
real(r64), dimension(4) :: v_dd_val_Real !! pppp(1), pnpn(2), pnnp(3), nnnn(4)

integer      :: ms, tt
complex(r64) :: aux, radial, aux_dir, aux_exch, dens_part, aux4
complex(r64), dimension(4) :: v_dd_value, term1, term2, aux1, aux2, aux3
real(r64)    :: angular, integral_factor, const_1, const_4
integer(i32) :: a_sh, b_sh, c_sh, d_sh, i_r, i_a
integer      :: HOspO2

HOspO2 = HOsp_dim/2

v_dd_value = zzero
v_dd_val_Real = zero
term1 = zzero
term2 = zzero

if (.NOT.EXPORT_PREA_DD) return
if ((a.GT.HOspO2).OR.(b.GT.HOspO2).OR.(c.GT.HOspO2).OR.(d.GT.HOspO2)) then
  print "(A)", " [ERROR] (m.e. Rea ME), a,b,c,d > HO_sp dim /2 !!!"
  return
endif

integral_factor = t3_DD_CONST
!! NOTE :: Remember that radial functions already have the factor 1/b**3
integral_factor = integral_factor * 0.5d0 * (HO_b**3)
integral_factor = integral_factor  / ((2.0d0 + alpha_DD)**1.5d0)
integral_factor = integral_factor * 4 * pi  ! add Lebedev norm factor

const_1 = 4.0d0 * alpha_DD
const_4 = alpha_DD * (alpha_DD - 1.0d0)

do i_r = 1, r_dim
  radial = weight_R(i_r) * exp((alpha_DD + 2.0d0) * (r(i_r) / HO_b)**2)
  do i_a = 1, angular_dim

    !! first term, (derivative of each rho_matrix)
    aux1 = zzero
    aux2 = zzero
    aux3 = zzero
    do tt = 1, 3 !! pnnp = 0, tt=4

      aux1(tt) = dens_pnt(5,i_r,i_a) - (x0_DD_FACTOR * dens_pnt(tt,i_r,i_a))
      aux1(tt) = aux1(tt) * rea_common_RadAng(c,d, i_r, i_a)

      do ms = 1, 4
        aux2(tt) = aux2(tt) + (((x0_DD_FACTOR * BulkHF(5, ms, i_r, i_a) * &
                                 AngFunctDUAL_HF(ms, b, d, i_a)) - &
                                (BulkHF(tt,ms,i_r,i_a) * &
                                 AngFunctDUAL_HF(ms, b, d, i_a))) * &
                               radial_2b_sho_memo(b_sh, d_sh, i_r) )
      enddo

      aux3(tt) = (aux1(tt) + aux2(tt)) * rea_common_RadAng(a,c, i_r, i_a)
      aux3(tt) = aux3(tt) * const_1 * dens_alpm1(i_r, i_a)

    enddo
    !! second term. (derivative of the inner rho^(alpha-1) )
    aux4 = rea_common_RadAng(b,d, i_r, i_a) * rea_common_RadAng(a,c, i_r, i_a)
    aux4 = aux4 * REACommonFields(i_r, i_a)
    aux4 = aux4 * dens_alpm1(i_r, i_a) / dens_pnt(5, i_r, i_a)
    aux4 = aux4 * const_4

    angular = weight_LEB(i_a)

    !v_nnnn = v_pppp
    v_dd_value(1) = v_dd_value(1) + (radial * angular * (aux3(1) + aux4))
    v_dd_value(4) = v_dd_value(4) + (radial * angular * (aux3(2) + aux4))
    ! pn pn
    v_dd_value(2) = v_dd_value(2) + (radial * angular * (aux3(3) + aux4))
    ! pn np
    v_dd_value(3) = v_dd_value(3) + 0.0d0

  enddo  ! angular iter_
enddo    ! radial  iter_

v_dd_val_Real(1) = real(v_dd_value(1), r64) * integral_factor
v_dd_val_Real(2) = real(v_dd_value(2), r64) * integral_factor
v_dd_val_Real(3) = real(v_dd_value(3), r64) * integral_factor
v_dd_val_Real(4) = real(v_dd_value(4), r64) * integral_factor

if (abs(imag(v_dd_value(2))) > 1.0d-9 ) then
    print "(A,D15.8,A,D20.8)", "[FAIL IMAG] v_prea_DD_abcd is not Real =", &
          real(v_dd_value(2)), " +j ", imag(v_dd_value(2))
endif

return
end function matrix_element_pseudoRearrangement_v1

!-----------------------------------------------------------------------------!
! function matrix_element_pseudoRearrangement                                 !
!                                                                             !
! Computes density dependent two body matrix elements over the density average!
! The derivation of the matrix element is derived from the Fields derivation  !
!                       v_dd_val_Real !! pppp(1), pnpn(2), npnp(3), nnnn(4)   !
! ## NOTE     pnnp = - pnpn(2)  and nppn = - npnp(3)
!             in this case, a,b,c,d are <= HOsp_dim / 2                       !
!-----------------------------------------------------------------------------!
function matrix_element_pseudoRearrangement_v2(a,b, c,d) result (v_dd_val_Real)

integer(i32), intent(in) :: a,b,c,d
real(r64), dimension(4) :: v_dd_val_Real !! pppp(1), pnpn(2), npnp(3), nnnn(4)

integer      :: ms, tt, t2
complex(r64) :: aux, radial, dens_part, aux5pp, aux5pn, aux5
complex(r64), dimension(4) :: v_dd_value, term1, term2, aux1, aux2, aux3, aux4
real(r64)    :: angular, integral_factor, const_1, const_5
integer(i32) :: a_sh, b_sh, c_sh, d_sh, i_r, i_a
integer      :: HOspO2

HOspO2 = HOsp_dim/2

v_dd_value = zzero
v_dd_val_Real = zero
term1 = zzero
term2 = zzero

if (.NOT.EXPORT_PREA_DD) return
if ((a.GT.HOspO2).OR.(b.GT.HOspO2).OR.(c.GT.HOspO2).OR.(d.GT.HOspO2)) then
  print "(A)", " [ERROR] (m.e. Rea ME), a,b,c,d > HO_sp dim /2 !!!"
  return
endif

integral_factor = t3_DD_CONST
!! NOTE :: Remember that radial functions already have the factor 1/b**3
integral_factor = integral_factor * 0.5d0 * (HO_b**3)
integral_factor = integral_factor  / ((2.0d0 + alpha_DD)**1.5d0)
integral_factor = integral_factor * 4 * pi  ! add Lebedev norm factor

const_1 = 4.0d0 * alpha_DD
const_5 = alpha_DD * (alpha_DD - 1.0d0)

do i_r = 1, r_dim
  radial = weight_R(i_r) * exp((alpha_DD + 2.0d0) * (r(i_r) / HO_b)**2)
  do i_a = 1, angular_dim
    angular = weight_LEB(i_a)
    !! first term, (derivative of each rho_matrix)
    aux1 = zzero
    aux2 = zzero
    aux3 = zzero
    aux4 = zzero

    !! second term, share functions for the direct term
    aux5pn = (rea_common_RadAng(a,c,i_r,i_a) * rea_common_RadAng(b,d,i_r,i_a))
    aux5pp = aux5pn - (rea_common_RadAng(a,d, i_r,i_a) * &
                       rea_common_RadAng(b,c, i_r,i_a))

    do tt = 1, 4  !! pppp  pnpn npnp nnnn
      select case(tt)
      case (1,4)
        aux1(tt) = aux5pp * dens_pnt(5, i_r, i_a)
        t2 = 1 + ((tt - 1) / 3)
        aux3(tt) = - x0_DD_FACTOR * dens_pnt(t2, i_r,i_a) * aux5pp
        aux5 = aux5pp
      case (2,3) ! pnpn, npnp
        aux1(tt) = dens_pnt(5, i_r, i_a) * aux5pn
        t2 = 4  - tt
        aux3(tt) = - x0_DD_FACTOR * dens_pnt(t2, i_r,i_a) * aux5pn
        aux5 = aux5pn
      end select

      !! exchange parts -----
      do ms = 1, 4
        select case (tt) ! ------------------
          case (1,4)
            aux = ((rea_common_RadAng(a,c, i_r,i_a) * &
                    AngFunctDUAL_HF(ms, b, d, i_a) * &
                    radial_2b_sho_memo(b_sh, d_sh, i_r)) - &
                   (rea_common_RadAng(a,d, i_r,i_a) * &
                    AngFunctDUAL_HF(ms, b, c, i_a) * &
                    radial_2b_sho_memo(b_sh, c_sh, i_r)))

            aux2(tt) = aux2(tt) + (x0_DD_FACTOR * BulkHF(5, ms, i_r, i_a) * aux)

            aux = ((rea_common_RadAng(a,c, i_r,i_a) * &
                    AngFunctDUAL_HF(ms, b, d, i_a) * &
                    radial_2b_sho_memo(b_sh, d_sh, i_r)) - &
                   (rea_common_RadAng(a,d, i_r,i_a) * &
                    AngFunctDUAL_HF(ms, b, c, i_a) * &
                    radial_2b_sho_memo(b_sh, c_sh, i_r)))
            t2 = 1 + ((tt - 1) / 3)
            aux4(tt) = aux2(tt) - (x0_DD_FACTOR * BulkHF(t2,ms, i_r, i_a) * aux)
          case (2,3)
            aux = (rea_common_RadAng(a,c, i_r,i_a) * &
                   AngFunctDUAL_HF(ms, b, d, i_a) * &
                   radial_2b_sho_memo(b_sh, d_sh, i_r))

            aux2(tt) = aux2(tt) + (x0_DD_FACTOR * BulkHF(5, ms, i_r, i_a) * aux)

            aux = (rea_common_RadAng(a,c, i_r,i_a) * &
                   AngFunctDUAL_HF(ms, b, d, i_a) * &
                   radial_2b_sho_memo(b_sh, d_sh, i_r))
            t2 = 4  - tt
            aux4(tt) = aux4(tt) - (BulkHF(t2,ms, i_r, i_a) * aux)
        end select ! ------------------
      enddo ! ms

      term1(tt) = const_1 * dens_alpm1(i_r, i_a) * &
                  (aux1(tt) + aux2(tt) + aux3(tt) + aux4(tt))
      !! second term. (derivative of the inner rho^(alpha-1) )
      term2(tt) = const_5 * aux5 * REACommonFields(i_r,i_a)&
                   * dens_alpm1(i_r,i_a) / dens_pnt(5, i_r,i_a)

      !! Final integral
      v_dd_value(tt) = v_dd_value(tt) + (radial * angular &
                                         * (term1(tt) + term2(tt)))
    enddo !tt

  enddo  ! angular iter_
enddo    ! radial  iter_

v_dd_val_Real(1) = real(v_dd_value(1), r64) * integral_factor
v_dd_val_Real(2) = real(v_dd_value(2), r64) * integral_factor
v_dd_val_Real(3) = real(v_dd_value(3), r64) * integral_factor
v_dd_val_Real(4) = real(v_dd_value(4), r64) * integral_factor

if (abs(imag(v_dd_value(2))) > 1.0d-9 ) then
    print "(A,D15.8,A,D20.8)", "[FAIL IMAG] v_prea_DD_abcd is not Real =", &
          real(v_dd_value(2)), " +j ", imag(v_dd_value(2))
endif

return
end function matrix_element_pseudoRearrangement_v2

!-----------------------------------------------------------------------------!
! function matrix_element_pseudoRearrangement                                 !
!                                                                             !
! Computes density dependent two body matrix elements over the density average!
! The derivation of the matrix element is derived from the Fields derivation  !
!                       v_dd_val_Real !! pppp(1), pnpn(2), npnp(3), nnnn(4)   !
! ## NOTE     pnnp = - pnpn(2)  and nppn = - npnp(3)
!             in this case, a,b,c,d are <= HOsp_dim / 2                       !
!-----------------------------------------------------------------------------!
function matrix_element_pseudoRearrangement(a,b, c,d) result (v_dd_val_Real)

integer(i32), intent(in) :: a,b,c,d
real(r64), dimension(4) :: v_dd_val_Real !! pppp(1), pnpn(2), pnnp(3), nnnn(4)

integer      :: ms, tt
complex(r64) :: aux, radial, aux_dir, aux_exch, dens_part, aux4
complex(r64), dimension(4) :: v_dd_value, term1, term2, aux1, aux2, aux3
real(r64)    :: angular, integral_factor, const_1, const_4
integer(i32) :: a_sh, b_sh, c_sh, d_sh, i_r, i_a
integer      :: HOspO2

HOspO2 = HOsp_dim/2

v_dd_value = zzero
v_dd_val_Real = zero
term1 = zzero
term2 = zzero

if (.NOT.EXPORT_PREA_DD) return
if ((a.GT.HOspO2).OR.(b.GT.HOspO2).OR.(c.GT.HOspO2).OR.(d.GT.HOspO2)) then
  print "(A)", " [ERROR] (m.e. Rea ME), a,b,c,d > HO_sp dim /2 !!!"
  return
endif

integral_factor = t3_DD_CONST
!! NOTE :: Remember that radial functions already have the factor 1/b**3
integral_factor = integral_factor * 0.5d0 * (HO_b**3)
integral_factor = integral_factor  / ((2.0d0 + alpha_DD)**1.5d0)
integral_factor = integral_factor * 4 * pi  ! add Lebedev norm factor

const_1 = 4.0d0 * alpha_DD
const_4 = alpha_DD !* (alpha_DD - 1.0d0)

do i_r = 1, r_dim
  radial = weight_R(i_r) * exp((alpha_DD + 2.0d0) * (r(i_r) / HO_b)**2)
  do i_a = 1, angular_dim

    !! first term, (derivative of each rho_matrix)
    aux1 = zzero
    aux2 = zzero
    aux3 = zzero

    !! second term. (derivative of the inner rho^(alpha-1) )
    aux4 = rea_common_RadAng(b,d, i_r, i_a) * rea_common_RadAng(a,c, i_r, i_a)
    aux4 = aux4 * REACommonFields(i_r, i_a)
    aux4 = aux4 * dens_alpm1(i_r, i_a) / dens_pnt(5, i_r, i_a)
    aux4 = aux4 * const_4

    angular = weight_LEB(i_a)

    !v_nnnn = v_pppp
    v_dd_value(1) = v_dd_value(1) + (radial * angular * (aux3(1) + aux4))
    v_dd_value(4) = v_dd_value(4) + (radial * angular * (aux3(2) + aux4))
    ! pn pn
    v_dd_value(2) = v_dd_value(2) + (radial * angular * (aux3(3) + aux4))
    ! pn np
    v_dd_value(3) = v_dd_value(3) + 0.0d0

  enddo  ! angular iter_
enddo    ! radial  iter_

v_dd_val_Real(1) = real(v_dd_value(1), r64) * integral_factor
v_dd_val_Real(2) = real(v_dd_value(2), r64) * integral_factor
v_dd_val_Real(3) = real(v_dd_value(3), r64) * integral_factor
v_dd_val_Real(4) = real(v_dd_value(4), r64) * integral_factor

if (abs(imag(v_dd_value(2))) > 1.0d-9 ) then
    print "(A,D15.8,A,D20.8)", "[FAIL IMAG] v_prea_DD_abcd is not Real =", &
          real(v_dd_value(2)), " +j ", imag(v_dd_value(2))
endif

return
end function matrix_element_pseudoRearrangement

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
! Print the two body matrix elements for verification                          !
!------------------------------------------------------------------------------!
subroutine are_same_states(a, b, c, d, a_eq_b, a_eq_c, c_eq_d, b_eq_d)
    integer, intent(in) :: a, b, c, d
    integer :: a_eq_b, a_eq_c, c_eq_d, b_eq_d

    a_eq_b = 0
    a_eq_c = 0
    c_eq_d = 0
    b_eq_d = 0
    if (a.eq.b) then
        a_eq_b = 1
    end if
    if (a.eq.c) then
        a_eq_c = 1
    end if
    if (c.eq.d) then
        c_eq_d = 1
    end if
    if (b.eq.d) then
        b_eq_d = 1
    end if
end subroutine are_same_states

function which_t_is_Vabcd(ta, tb, tc, td) result (t)
    !  Returns the t label position for ab_cd -> pp_pp=1, pn_pn=2, ...
    integer, intent(in) :: ta, tb, tc, td
    integer :: t
    if (ta==-1) then
        if (tb==-1) then
            t = 1
        else
            if (tc==-1) then
                t = 2
            else
                t = 3
            end if
        end if
    else
        if (tb==1) then
            t = 6
        else
            if (tc==1) then
                t = 5
            else
                t = 4
            end if
        end if
    end if
    return
end function which_t_is_Vabcd


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




subroutine test_export_pnpn_mmee_uncoupled(ndim)

integer, intent(in)     :: ndim
integer                 :: a, b, c, d, nO2, kk,i1,i2,i3,i4, &
                           ab_indx, cd_indx, it, perm, ii1,ii2,ii3,ii4, tt, sg_
real(r64) :: h2b_64
real(r32) :: h2b
real(r64), dimension(4) :: me_val
integer,   dimension(2) :: non_zero
real(r32), dimension(:,:), allocatable :: registered_h2b
logical :: REGISTER_BB
nO2 = ndim / 2

!print "(A)", "  [WARNING] Not gonna do the exporting of matrix elements"
!return
REGISTER_BB = HO_Nmax .LT. 8    !!! MZ=8 -> 88GB,  MZ=7 -> 25GB, MZ=6 -> 0.4GB

!!! Obtain the pnpn part as an (a,b) vs (c,d) matrix:
if (REGISTER_BB) then
  allocate(registered_h2b(nO2*nO2, nO2*nO2))
  registered_h2b = zero
endif

open(111, file='hamil_bb_init.gut')
do kk = 1, hamil_H2dim

  i1 = hamil_abcd(1+4*(kk-1))
  i2 = hamil_abcd(2+4*(kk-1))
  i3 = hamil_abcd(3+4*(kk-1))
  i4 = hamil_abcd(4+4*(kk-1))

!  print '(4L3,A,2L3)', (HOsp_2mt(i1).NE.-1), (HOsp_2mt(i2).NE.-1), &
!                  (HOsp_2mt(i3).NE.-1), (HOsp_2mt(i4).NE.-1), "  -> ", &
!                 ((HOsp_2mt(i1).NE.-1) .OR. (HOsp_2mt(i1).NE. 1)), &
!                 ((HOsp_2mt(i3).NE.-1) .OR. (HOsp_2mt(i4).NE. 1))
  if (((HOsp_2mt(i1).EQ.-1) .AND. (HOsp_2mt(i2).EQ.-1))) cycle
  if (((HOsp_2mt(i1).EQ. 1) .AND. (HOsp_2mt(i2).EQ. 1))) cycle
  if (((HOsp_2mt(i1).EQ. 1) .AND. (HOsp_2mt(i2).EQ.-1))) cycle
  if (((HOsp_2mt(i3).EQ. 1) .AND. (HOsp_2mt(i4).EQ.-1))) cycle
!  if (((HOsp_2mt(i3).NE.-1) .OR. (HOsp_2mt(i4).NE. 1))) cycle ! unnecessary

  h2b_64  = hamil_H2(kk)
  perm = hamil_trperm(kk)
  h2b = real(h2b_64, r32)

  !!! Loop on time reversal
  do it = 1, 2

    if ( it == 2 ) then
      if ( HOsp_2mj(i1) + HOsp_2mj(i2) == 0 ) cycle
      call find_timerev(perm,i1,i2,i3,i4)
      h2b = sign(one,perm*one) * h2b
    endif

    ii1 = i1
    ii2 = i2
    ii3 = i3
    ii4 = i4
    if (i1 > nO2) ii1 = ii1 - nO2
    if (i2 > nO2) ii2 = ii2 - nO2
    if (i3 > nO2) ii3 = ii3 - nO2
    if (i4 > nO2) ii4 = ii4 - nO2

    ab_indx = ((ii1 - 1) * nO2) + ii2
    cd_indx = ((ii3 - 1) * nO2) + ii4
    sg_ = 1
    write(111, fmt='(4I5,D15.6)') ii1, ii2, ii3, ii4, h2b

    !! 1. Criteria from module_fields.calculate_fields (general)
    if (((i1 .GT. nO2).AND.(i2 .GT. nO2)) .OR. &
        ((i1 .LE. nO2).AND.(i2 .LE. nO2))) then
      cycle    ! nn_ nn_ or ! pp_ pp_
    else
      tt = 4*HOsp_2mt(i1) + 2*(HOsp_2mt(i2) + HOsp_2mt(i3)) + HOsp_2mt(i4) - 1
      tt = 3 + (tt / 2)
      if ((tt.EQ.1) .OR. (tt.EQ.4)) then
        tt = 2  ! pn_ pn_
      else
        tt = 3  ! pn_ np_
        cd_indx = ((ii4 - 1) * nO2) + ii3
        sg_ = -1
      end if
    endif

    if (REGISTER_BB) then
      registered_h2b(ab_indx, cd_indx) = sg_ * h2b
      if ((kdelta(i1,i3) * kdelta(i2,i4)) .NE. 1) then
        registered_h2b(cd_indx, ab_indx) = sg_ * h2b
      endif
    endif

  enddo
end do !! kk loop
close(111)



open(111, file="dd_pnpn_me.gut")
open(112, file="bb_pnpn_me.gut")
!! index introduction
write(111, fmt="(A)") "%%%  MAT. ELEMS (p, n) ::  %%%%%%%%%%%%%%%%%%%%%%%%%%%%"
write(111, fmt="(I4)") ndim
write(111, fmt="(A)") "i_sp ant_inx sh  n  l  j mj"
write(112, fmt="(A)") "%%%  MAT. ELEMS (p, n) ::  %%%%%%%%%%%%%%%%%%%%%%%%%%%%"
write(112, fmt="(I4)") ndim
write(112, fmt="(A)") "i_sp ant_inx sh  n  l  j mj"
do a=1, nO2
  write(111, fmt="(I4,I8,4I3,I4)") a, HOsp_ant(a), HOsp_sh(a), HOsp_n(a), &
                                   HOsp_l(a), HOsp_2j(a), HOsp_2mj(a)
  write(112, fmt="(I4,I8,4I3,I4)") a, HOsp_ant(a), HOsp_sh(a), HOsp_n(a), &
                                   HOsp_l(a), HOsp_2j(a), HOsp_2mj(a)
end do

write(111, fmt="(A)") "%%%  MAT. ELEMS (a, b) ::  %%%%%%%%%%%%%%%%%%%%%%%%%%%%"
write(112, fmt="(A)") "%%%  MAT. ELEMS (a, b) ::  %%%%%%%%%%%%%%%%%%%%%%%%%%%%"

if (.NOT.REGISTER_BB) write(112, fmt="(A)") &
  "SKIP MZ>8 this file is not available, generate it from [hamil_bb_init.gut]"

print "(A)", " [    ] Exporting of DD non-zero PN matrix elements."

do a = 1, nO2
  do b = 1, nO2

    non_zero = 0
    !! Density-Dependent part
    if (dabs(t3_DD_CONST) .LT. 0.001) then ! skip it dd term is fixed to 0
    do c=1, no2
      do d=1, no2

        me_val = matrix_element_v_DD(a,b,c,d, .TRUE.)

        !! just check pnpn channel
        if (dabs(me_val(2)) .LT. 1.0d-09) cycle

        if (non_zero(1) .EQ. 0) then !! include the header
          write(111, fmt="(2I4,A)", advance='no') a, b, " // "
        end if
        write(111, fmt="(2I4,D15.6,A)", advance='no') c, d, me_val(2), ", "

        non_zero(1) = non_zero(1) + 1

      end do
    end do
    if (non_zero(1) .GT. 0) write(111,fmt="(A)") ""
    endif

    ab_indx = ((a - 1) * nO2) + b
    if (REGISTER_BB) then
      !! Hamiltonian part
      do c=1, no2
        do d=1, no2

          cd_indx = ((c - 1) * nO2) + d
          h2b = registered_h2b(ab_indx, cd_indx)

          !! just check pnpn channel
          if (abs(h2b) .LT. 1.0d-09) cycle

          if (non_zero(2) .EQ. 0) then !! include the header
            write(112, fmt="(2I4,A)", advance='no') a, b, " // "
          end if
          write(112, fmt="(2I4,D15.6,A)", advance='no') c, d, h2b, ", "

          non_zero(2) = non_zero(2) + 1

        end do
      end do

      if (non_zero(2) .GT. 0) write(112,fmt="(A)") ""
    endif

  end do
  print "(2(A,I6),2I9)", "   progress ... ", a,"/",b, non_zero(1),non_zero(2)
end do
close(111)
close(112)

if (REGISTER_BB) deallocate(registered_h2b)

print "(A)", " [DONE] Exporting of DD non-zero PN matrix elements."

end subroutine test_export_pnpn_mmee_uncoupled

END MODULE DensityDep
!==============================================================================!
! End of file                                                                  !
!==============================================================================!

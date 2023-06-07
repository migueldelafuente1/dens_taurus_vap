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
! - subroutine test_legendrePolynomials                                        !
! - subroutine test_sphericalHarmonics_print                                   !
! - subroutine test_sphericalHarmonics_ortogonality                            !
! - subroutine test_radial_2wf                                                 !
! - subroutine test_density_1BME_positive_def                                  !
! - subroutine test_print_density_me                                           !
!==============================================================================!
MODULE DensityDep

use Constants
use MathMethods
use Basis
use Hamiltonian
use Fields
use Lebedev

implicit none
PUBLIC

logical   :: eval_density_dependent = .TRUE.
logical   :: eval_rearrangement     = .FALSE.
logical   :: eval_explicit_fieldsDD = .FALSE.
real(r64) :: t3_DD_CONST  = 0.0!= 1350.00d+00    ! constant of the DD term [MeV]
real(r64) :: x0_DD_FACTOR = 0.0!= 1.0d+00        ! exchange factor of DD term
real(r64) :: alpha_DD     = 0.0!=  0.33333d+00   ! power of the DD term

! 0 trapezoidal, 1 Gauss-Legendre, 2 Gauss-Laguerre(r)/Legendre (t,p), 3 Laguerre-Lebedev
integer   :: integration_method  = 3
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

!! spatial steps for trapezoidal integration
real(r64) :: d_r   ! = R_MAX / (r_dim - 1)
real(r64) :: d_theta
real(r64) :: d_phi

! Quadratures for integration
real(r64), dimension(:), allocatable   :: x_R
real(r64), dimension(:), allocatable   :: weight_R
real(r64), dimension(:), allocatable   :: weight_THE
real(r64), dimension(:), allocatable   :: weight_PHI
real(r64), dimension(:,:), allocatable :: x_Leb
real(r64), dimension(:), allocatable   :: weight_LEB

real(r64), dimension(:), allocatable  :: r, r_export
real(r64), dimension(:), allocatable  :: theta, theta_export
real(r64), dimension(:), allocatable  :: cos_th, cos_th_export
real(r64), dimension(:), allocatable  :: phi, phi_export

complex(r64), dimension(:,:), allocatable :: density, density_export !(ir, iang)
complex(r64), dimension(:,:), allocatable :: dens_alpha, dens_alpm1
complex(r64), dimension(:,:,:), allocatable :: dens_pnt   !(pp, nn, pn, np, total;i_r,iang)
!complex(r64), dimension(:,:,:), allocatable :: dens_pnmix !(pn-Odd,pn-Even,np-O,np-E; i_r,iang)
complex(r64), dimension(:,:), allocatable :: density_export_p, density_export_n
complex(r64), dimension(:,:), allocatable :: pairdens_export
complex(r64), dimension(:,:), allocatable :: pairdens_export_n,pairdens_export_p

!! Pre-calculated saved functions and coefficients
real(r64), dimension(:, :, :), allocatable, save :: radial_2b_sho_memo        !(ish1, ish2, ir)
real(r64), dimension(:, :, :), allocatable, save :: radial_2b_sho_noexp_memo  !(ish1, ish2, ir)
real(r64), dimension(:, :, :), allocatable, save :: radial_2b_sho_export_memo !(ish1, ish2, ir)

complex(r64), dimension(:, :), allocatable, save    :: sph_harmonics_memo
complex(r64), dimension(:, :, :), allocatable, save :: sphharmDUAL_memo ! Y*(a) Y(b)
real(r64), dimension(:, :, :), allocatable, save    :: dens_Y_KM_me

complex(r64), dimension(:,:,:,:), allocatable, save :: AngFunctDUAL_HF ! CGa CGb Y*(a) Y (b)
complex(r64), dimension(:,:,:,:), allocatable, save :: AngFunctDUAL_P1 ! CGa CGb Y (a) Y (b)
complex(r64), dimension(:,:,:,:), allocatable, save :: AngFunctDUAL_P2 ! CGa CGb Y*(a) Y*(b)
complex(r64), dimension(:,:,:,:), allocatable, save :: BulkHF ! (tt,msms',r,ang) DEF:Sum AngFunctDUAL_HF * rho [pp,nn,pn,np, tot]  [(++,+-,-+,--)]
complex(r64), dimension(:,:,:,:), allocatable, save :: BulkP1 ! (tt,msms',r,ang) DEF:Sum AngFunctDUAL_P1(msms') - AngFunctDUAL_P1(ms'ms) * kappaLR
complex(r64), dimension(:,:,:,:), allocatable, save :: BulkP2 ! (tt,msms',r,ang) DEF:Sum AngFunctDUAL_P2 * kappaRL

complex(r64), dimension(:, :), allocatable    :: rearrangement_me  !(isp1, isp2)
complex(r64), dimension(:, :), allocatable    :: rearrang_field    !(isp1, isp2)
complex(r64), dimension(:, :, :, :), allocatable :: rea_common_RadAng !(isp1,isp2, ir,iang)
complex(r64), dimension(:, :), allocatable    :: REACommonFields   !(ir, iang))

integer, dimension(:), allocatable :: HOsh_ant, HOsp_ant

!! Related to the hamiltonian and Fields
real(r64), dimension(:), allocatable :: hamil_DD_H2        ! 2-body part
real(r64), dimension(:,:), allocatable :: hamil_DD_H2_byT  ! 2-body part
integer(i64) :: hamil_DD_H2dim, hamil_DD_H2dim_all         ! number of 2BME stored
integer(i16), dimension(:), allocatable :: hamil_DD_abcd   ! indices of 2BME
integer(i8) , dimension(:), allocatable :: hamil_DD_trperm ! time reversal permut.
integer   :: iteration = 0, global_iter_max = 0
integer   :: VSsh_dim = 0,  VSsp_dim = 0, VSsp_dim2 = 0
integer, dimension(:), allocatable     :: VSsh_list, &
                                          VStoHOsp_index, VStoHOsh_index

integer   :: seed_type_sym  = 0        ! (UNDEFINED)
logical   :: haveX0M1       = .FALSE.  ! |x0 - 1| > 1E-6
logical   :: evalFullSPSpace= .TRUE.   ! compute the full a,b, c,d space for explicit DD fields (cannot be set)
logical   :: exportVSPSpace = .FALSE.  ! the export of a reduced val.space
logical   :: evalQuasiParticleVSpace = .FALSE. ! Export for the QP sp states, not the VS procedure

integer   :: NHO_vs, NHO_co !! Major Shell number of the Valence. Sp to be exported
logical   :: NOT_DEL_FILE
logical   :: PRINT_GUTS = .FALSE.
logical   :: DOING_PROJECTION = .FALSE.

!! [END] DENSITY DEPENDENT MODIFICATIONS =====================================


CONTAINS

!-----------------------------------------------------------------------------!
! Import and define array dimensions from file DD_PARAMS.txt                  !
!-----------------------------------------------------------------------------!
subroutine import_DD_parameters

integer :: runit = 99
integer :: ios, i, seed_type_imported, aa, a, a_ant
logical :: is_exist
character(len=*), parameter :: formatST = "(1a)", &
                               formatI1 = "(1a30, 1i1)", &
                               formatI3 = "(1a30, 1i3)", &
                               formatF6 = "(1a30, 1f9.6)", &
                               formatEE = "(1a30, 1es12.6)", &
                               formatStrHeader = "(1a30)", &
                               formatII = "(1a30, 1i1, 99i6)"

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
  eval_density_dependent = .FALSE.
  export_density = .FALSE.
  return
endif

!! Reads the input parameters
read(runit,formatST) str_
read(runit,formatI1) str_, aux_int
eval_density_dependent = aux_int.EQ.1
read(runit,formatI1) str_, aux_int
eval_rearrangement = aux_int.EQ.1
read(runit,formatI1) str_, aux_int
eval_explicit_fieldsDD = aux_int.EQ.1

read(runit,formatEE) str_, t3_DD_CONST
read(runit,formatEE) str_, x0_DD_FACTOR
read(runit,formatEE) str_, alpha_DD
read(runit,formatST) str_
read(runit,formatST) str_
read(runit,formatI1) str_, integration_method
read(runit,formatI1) str_, aux_int
export_density = aux_int.EQ.1

read(runit,formatI3) str_, r_dim
read(runit,formatI3) str_, Omega_Order
read(runit,formatI3) str_, THE_grid
read(runit,formatI3) str_, PHI_grid
read(runit,formatI1) str_, aux_int
evalQuasiParticleVSpace = aux_int.GE.1
read(runit,formatI1) str_, aux_int
exportVSPSpace = aux_int.GE.1

VSsh_dim = aux_int
if (exportVSPSpace) then
  if (VSsh_dim.LE.HOsh_dim) then
    print "(A,I3,A)", "   ... Reading VS sh states", VSsh_dim, &
      " (error if wrong sh dimension)"
    print *, ""
    backspace runit
    allocate(VSsh_list(VSsh_dim))
    read(runit,formatStrHeader, advance='no') str_
    read(runit,*) VSsh_dim, (VSsh_list(i),i=1,VSsh_dim)
    call set_valence_space_to_export
  else
    print *,"   ... Reading for FULL valence space"
    VSsh_dim  = HOsh_dim
    VSsp_dim  = HOsp_dim
    VSsp_dim2 = HOsp_dim2
    allocate(VSsh_list(VSsh_dim))
    allocate(VStoHOsp_index(VSsp_dim))
    do i=1, VSsh_dim
      VSsh_list(i) = HOsh_na(i)
    end do
    do i=1, VSsp_dim
      VStoHOsp_index(i) = i
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

!! set alpha more exact 1/3 double precission
if (abs(alpha_DD - (1.0/3.0)) < 1.0D-4) alpha_DD = 0.333333333333333

d_r = R_MAX / (r_dim - 1)

allocate(x_R(r_dim))
allocate(weight_R(r_dim))
allocate(r(r_dim))
allocate(r_export(r_dim))
print *,           '   Density dep. Interaction values  '
print *,           '-----------------------------------------------'
print '(A,L10)',   'eval_density_dependent =', eval_density_dependent
print '(A,L10)',   'eval_rearrangement     =', eval_rearrangement
print '(A,L10)',   'eval_explicit_fieldsDD =', eval_explicit_fieldsDD
print '(A,F10.4)', 't3_DD_CONST (MeV fm-4) =', t3_DD_CONST
print '(A,F10.6)', 'x0_DD_FACTOR           =', x0_DD_FACTOR
print '(A,F10.6)', 'alpha_DD               =', alpha_DD
print *, ''
print *,           '   Integration parameters  '
print *,           '-----------------------------------------------'
print '(A,I10)',   'integration_method =', integration_method
print '(A,L10)',   'export_density     =', export_density
print '(A,I10)',   'r_dim              =', r_dim
print '(A,I10)',   'Omega_Order        =', Omega_Order
print '(A,I10)',   'THE_grid           =', THE_grid
print '(A,I10)',   'PHI_grid           =', PHI_grid
!print '(A,F10.6)', 'R_MAX (fm)         =', R_MAX
print '(A,2L10)',  'eval/export Val.Sp =', evalFullSPSpace, exportVSPSpace
print '(A,2L10)',  'export QP Val.Sp   =', evalQuasiParticleVSpace

if (.NOT.exportVSPSpace) then
  deallocate(hamil_H2cpd_DD) ! It wont be used
endif
if (exportVSPSpace)then
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
  print "(A,2I6)", '    ... N-shell HO for core(max)/vs(max):',NHO_co,NHO_vs
endif
if (eval_explicit_fieldsDD) then
  print '(A,3L10)', " [Explicit DD Field Eval.] Compute Full Valence Space =",&
    evalFullSPSpace
endif
print *, ''

haveX0M1 = abs(x0_DD_FACTOR - 1.0d+0) > 1.0d-6

end subroutine import_DD_parameters


!-----------------------------------------------------------------------------!
! Subroutine to identify the index of the valence space from the basis        !
!-----------------------------------------------------------------------------!
subroutine set_valence_space_to_export
integer :: i, j, i_ant
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

allocate(VStoHOsp_index(VSsp_dim))
do i=1, VSsp_dim
  VStoHOsp_index(i) = temp_list_index(i)
!  print "(A,2I6)", "VS index in HO basis:", i, VStoHOsp_index(i)
end do
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

if (integration_method == 3) then
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
  allocate(weight_THE(theta_dim))
  allocate(weight_PHI(phi_dim))
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

    weight_PHI(i)   = one
    weight_THE(i)   = weight_LEB(i) ! don't do the sqrt, weight_LEB might be < 0

    if (PRINT_GUTS) then
      write(620,fmt='(I5,3D22.13)') i, cos_th(i), phi(i), weight_LEB(i)
    endif

    sum_ = sum_ + weight_LEB(i)
  enddo

  if (PRINT_GUTS) close(620)

else ! Non Lebedev Quadratures
  theta_dim = THE_grid
  phi_dim   = PHI_grid
  angular_dim = theta_dim * phi_dim

  !! spatial steps for trapezoidal integration
  d_theta = pi / (theta_dim - 1)
  d_phi   = 2.0d0 * pi / (phi_dim - 1)

  allocate(weight_THE(theta_dim))
  allocate(weight_PHI(phi_dim))

  allocate(theta(theta_dim))
  allocate(theta_export(theta_dim))
  allocate(cos_th(theta_dim))
  allocate(cos_th_export(theta_dim))
  allocate(phi(phi_dim))
  allocate(phi_export(phi_dim))

!! RADIAL

if (integration_method == 0) then     ! Trapezoidal
  print '(A)', ' Using Trapezoidal Integration Method'
  do i_r = 1, r_dim
    r(i_r)   = (i_r - 1) * d_r
    r_export(i_r) = r(i_r)
    weight_R(i_r) = d_r
  enddo
elseif (integration_method == 1) then ! Legendre-radial
  print '(A)', ' Using Gauss-Legendre (3) Integration Method'
  call GaussLegendre(0.0d0, R_MAX, x_R,  weight_R  , r_dim)
  do i_r = 1, r_dim
    r(i_r)        = x_R(i_r)
    r_export(i_r) = x_R(i_r)
  enddo
elseif (integration_method == 2) then ! Laguerre
  print '(A)', ' Using Gauss-Laguerre/Legendre (2) Integration Method'
  call GaussLaguerre(x_R, weight_R, r_dim, 0.5d+00)
  do i_r = 1, r_dim
    r(i_r)        = HO_b * sqrt(x_R(i_r) / (2.0 + alpha_DD)) !
    r_export(i_r) = HO_b * sqrt(x_R(i_r)) !
  enddo
endif

!! ANGULAR
if (integration_method == 0) then   ! Trapezoidal
  do i_th = 1, theta_dim
    theta(i_th)  = (i_th - 1)  * d_theta
    cos_th(i_th) = cos(theta(i_th))
    theta_export(i_th)  = theta(i_th)
    cos_th_export(i_th) = cos_th(i_th)
    weight_THE(i_th)    = d_theta
  enddo
  do i_phi = 1, phi_dim
    phi(i_phi)   = (i_phi - 1) * d_phi
    phi_export(i_phi) = phi(i_phi)
    weight_PHI(i_phi) = d_phi
  enddo
else                                ! Legendre
  !call GaussLegendre( 0.0d0, pi,    theta,  weight_THE, theta_dim)
  call GaussLegendre(-1.0d0, 1.0d0, cos_th, weight_THE, theta_dim)
  do i_th = 1, theta_dim
    theta(i_th) = acos(cos_th(i_r))
    theta_export(i_th)  = theta(i_th)
    cos_th_export(i_th) = cos_th(i_th)
  enddo

  call GaussLegendre( 0.0d0, 2.0d0*pi, phi, weight_PHI, phi_dim)
  do i_phi = 1, phi_dim
    phi_export(i_phi) = phi(i_phi)
  enddo
endif

endif ! Non Lebedev Quadratures

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
  real(r64)    :: x1, x2, y1, y2, z1, z2

  print *, "Setting up DD module [   ]"
  seed_type_sym = seedtype
  global_iter_max = itermax
  DOING_PROJECTION = (proj_Mphip > 1).OR.(proj_Mphin > 1)

  call import_DD_parameters
  if (.NOT.eval_density_dependent) then
    print "(A)", "  DD module is turned off, skip DD array setting [OK]"
    return
  endif

  call set_integration_grid
  call set_allocate_density_arrays

  call set_B_radial_coefficients
  call set_Radial2body_basis
  call set_SphericalHarmonic_basis
  call set_Y_KM_matrixElements
  call set_rearrangement_RadAng_fucntions

  print "(A)", "  Setting up DD module [DONE]"
  print "(A,L1)", "  DOING_PROJECTION = ", DOING_PROJECTION

  !call test_print_basis_quantum_numbers
  !call tests_sixjsymbols
  !!! Tests over the implemented functions
  !call test_legendrePolynomials
  !call test_sphericalHarmonics_print
!  call test_sphericalHarmonics_ortogonality
  !call test_radial_2wf (DEPRECATED)
  !call test_load_computing_vs_access
  !call test_sphericalHarmonics_SumRule
  !call test_angular_function_recoupling
!  call test_dualAngularFunctionProperties
  print "(A)", "  Test DD module [DONE]"

end subroutine set_densty_dependent

subroutine set_allocate_density_arrays

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
allocate(radial_2b_sho_noexp_memo(HOsh_dim,HOsh_dim,r_dim))

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
        ! integration_method == 1 (r =x_leg over R)   * == 2 (r = b * sqrt(x_lag))
        radial = two_sho_radial_functions_bench(a_sh, b_sh, r(i_r))

        radial_2b_sho_memo(a_sh, b_sh, i_r) = radial
        radial_2b_sho_memo(b_sh, a_sh, i_r) = radial ! a_sh=b_sh overwrites

        ! Radial grid for
        !r _lag = b *(x_R / 2+alpha)**0.5
        radial = two_sho_radial_functions(a_sh, b_sh, r(i_r), .TRUE.)

        !! Note:: For using only the polynomials of the rad wf set to .FALSE.

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
        radial_2b_sho_noexp_memo(a_sh, b_sh, i_r) = radial
        radial_2b_sho_noexp_memo(b_sh, a_sh, i_r) = radial

        if (export_density) then
          if (integration_method > 1) then
            radial = two_sho_radial_functions(a_sh,b_sh, r_export(i_r), .FALSE.)
!            radial = two_sho_radial_functions_bench(a_sh, b_sh, r_export(i_r)) &
!                    * exp(r_export(i_r)**2)
          else
            !radial = two_sho_radial_functions(a_sh,b_sh, r_export(i_r), .TRUE.)
            radial = two_sho_radial_functions_bench(a_sh, b_sh, r_export(i_r))
          endif

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

    if (integration_method == 3) then
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

!!!  (OLD treatment that complete the Ang. func to directly use the full space !
!    ! == ( also fill the space for nn, pn, np space) ==
!    do i_an = 1, angular_dim
!      do ms = 1, 4
!        AngFunctDUAL_HF(ms,i+spO2,j     ,i_an) = AngFunctDUAL_HF(ms,i,j,i_an)
!        AngFunctDUAL_HF(ms,i     ,j+spO2,i_an) = AngFunctDUAL_HF(ms,i,j,i_an)
!        AngFunctDUAL_HF(ms,i+spO2,j+spO2,i_an) = AngFunctDUAL_HF(ms,i,j,i_an)
!
!        AngFunctDUAL_P1(ms,i+spO2,j     ,i_an) = AngFunctDUAL_P1(ms,i,j,i_an)
!        AngFunctDUAL_P1(ms,i     ,j+spO2,i_an) = AngFunctDUAL_P1(ms,i,j,i_an)
!        AngFunctDUAL_P1(ms,i+spO2,j+spO2,i_an) = AngFunctDUAL_P1(ms,i,j,i_an)
!
!        AngFunctDUAL_P2(ms,i+spO2,j     ,i_an) = AngFunctDUAL_P2(ms,i,j,i_an)
!        AngFunctDUAL_P2(ms,i     ,j+spO2,i_an) = AngFunctDUAL_P2(ms,i,j,i_an)
!        AngFunctDUAL_P2(ms,i+spO2,j+spO2,i_an) = AngFunctDUAL_P2(ms,i,j,i_an)
!      enddo
!    enddo !! ---------------------------------------------------------------- !

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
subroutine compute_bulkDens4Fields_bench(a, b, a_sh, b_sh, i_r, i_a, &
                                   rhoLR, kappaLR, kappaRL, ndim)
integer, intent(in)      :: a, b, a_sh, b_sh, i_r, i_a, ndim
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, kappaLR, kappaRL

integer      :: spO2, ms,ms2, la, lb, ja, jb, ma, mb, ind_jm_a, ind_jm_b
complex(r64) :: roP, roN, rPN, rNP, kaP, kaN, kaCcP, kaCcN, kPN, kNP, &
                kCcNP, kCcPN, A_part,B1_part,B2_part, aux
real(r64)    :: radial_ab
spO2 = HOsp_dim / 2

!!!  ASSERTION  !!! , a, b only from pp space:
if ((a > spO2).OR.(b > spO2)) then
   print *, " [ASS. ERROR] compute_bulkDens4Fields a,b (sp) in neutron space"
   STOP
endif

radial_ab = radial_2b_sho_noexp_memo(a_sh, b_sh, i_r)

roP = radial_ab * rhoLR  (b, a)
roN = radial_ab * rhoLR  (b +spO2, a +spO2)
rPN = radial_ab * rhoLR  (b      , a +spO2)
rNP = radial_ab * rhoLR  (b +spO2, a)

kaP = radial_ab * kappaLR(a, b)
kaN = radial_ab * kappaLR(a +spO2, b +spO2)
kPN = radial_ab * kappaLR(a      , b +spO2)
kNP = radial_ab * kappaLR(a +spO2, b)

kaCcP = radial_ab * kappaRL(a, b)
kaCcN = radial_ab * kappaRL(a +spO2, b +spO2)
kCcPN = radial_ab * kappaRL(a      , b +spO2)
kCcNP = radial_ab * kappaRL(a +spO2, b)

!---------------------- density matrix elements by spherical harmonic m.e  ---!
!aux = zzero
!la   = HOsp_l(a)
!ja   = HOsp_2j(a)
!ma   = HOsp_2mj(a)
!ind_jm_a = angular_momentum_index(ja, ma, .TRUE.)
!lb   = HOsp_l(b)
!jb   = HOsp_2j(b)
!mb   = HOsp_2mj(b)
!ind_jm_b = angular_momentum_index(jb, mb, .TRUE.)
!do K = abs(ja - jb) / 2, (ja + jb) / 2
!  M = (mb - ma)/2
!  if ((MOD(K + la + lb, 2) == 1).OR.(abs(M) > K)) cycle
!
!  ind_km = angular_momentum_index(K, M, .FALSE.)
!  aux   = aux + (dens_Y_KM_me(ind_jm_a, ind_jm_b, ind_km) * &
!                   sph_harmonics_memo(ind_km, i_a))
!enddo -----------------------------------------------------------------------!
aux = AngFunctDUAL_HF(1,a,b,i_a) + AngFunctDUAL_HF(4,a,b,i_a)

dens_pnt(1,i_r,i_a) = dens_pnt(1,i_r,i_a) + (aux * roP)
dens_pnt(2,i_r,i_a) = dens_pnt(2,i_r,i_a) + (aux * roN)
dens_pnt(3,i_r,i_a) = dens_pnt(3,i_r,i_a) + (aux * rPN)
dens_pnt(4,i_r,i_a) = dens_pnt(4,i_r,i_a) + (aux * rNP)

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

  B1_part = (AngFunctDUAL_P1(ms,a,b,i_a) - AngFunctDUAL_P1(ms2,a,b,i_a))
  BulkP1(1,ms,i_r,i_a) = BulkP1(1,ms,i_r,i_a) + (B1_part * kaP) !pp
  BulkP1(2,ms,i_r,i_a) = BulkP1(2,ms,i_r,i_a) + (B1_part * kaN) !nn

  B1_part = AngFunctDUAL_P1(ms ,a,b,i_a) &
             + (x0_DD_FACTOR*AngFunctDUAL_P1(ms2,a,b,i_a))
  B2_part = AngFunctDUAL_P1(ms2,a,b,i_a) &
             + (x0_DD_FACTOR*AngFunctDUAL_P1(ms ,a,b,i_a))
  BulkP1(3,ms,i_r,i_a) = BulkP1(3,ms,i_r,i_a) + (B1_part*kPN - B2_part*kNP) !pn
  BulkP1(4,ms,i_r,i_a) = BulkP1(4,ms,i_r,i_a) + (B1_part*kNP - B2_part*kPN) !np
  !! BulkP1_** are  common for the paring and rearrangement fields respectively.

  B2_part = AngFunctDUAL_P2(ms,a,b,i_a)
  BulkP2(1,ms,i_r,i_a) = BulkP2(1,ms,i_r,i_a) + (B2_part * kaCcP) !pp
  BulkP2(2,ms,i_r,i_a) = BulkP2(2,ms,i_r,i_a) + (B2_part * kaCcN) !nn
  BulkP2(3,ms,i_r,i_a) = BulkP2(3,ms,i_r,i_a) + (B2_part * kCcPN) !pn
  BulkP2(4,ms,i_r,i_a) = BulkP2(4,ms,i_r,i_a) + (B2_part * kCcNP) !np
enddo

end subroutine compute_bulkDens4Fields_bench

!------------------------------------------------------------------------------!
! subroutine to evaluate the the value and transposed values for all spaces    !
!------------------------------------------------------------------------------!
subroutine compute_bulkDens4Fields(a, b, a_sh, b_sh, i_r, i_a, &
                                   rhoLR, kappaLR, kappaRL, ndim)
integer, intent(in)      :: a, b, a_sh, b_sh, i_r, i_a, ndim
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, kappaLR, kappaRL

integer      :: spO2, par_ind,    ms,ms2, a_n, b_n
complex(r64) :: roP, roN, rPN, rNP, kaP,kaN,kaCcP, kaCcN,kPN,kNP,kCcNP,kCcPN, &
                roPt, roNt, rPNt, rNPt, kaPt,kaNt,kaCcPt, kaCcNt,kPNt,kNPt, &
                kCcNPt,kCcPNt, A_part,B1_part,B2_part, aux, sum_
real(r64)    :: radial_ab
logical      :: aNeQb
spO2 = HOsp_dim / 2

! assertion, a, b only from pp space:
if ((a > spO2).OR.(b > spO2)) then
   print *, " [ASS. ERROR] compute_bulkDens4Fields a,b (sp) in neutron space"
   STOP
endif
aNeQb = kdelta(a, b).ne.1
a_n   = a + spO2
b_n   = b + spO2

radial_ab = radial_2b_sho_noexp_memo(a_sh, b_sh, i_r)

roP = radial_ab * rhoLR  (b  ,a)
roN = radial_ab * rhoLR  (b_n,a_n)
rPN = radial_ab * rhoLR  (b  ,a_n)
rNP = radial_ab * rhoLR  (b_n,a)

kaP = radial_ab * kappaLR(a  ,b)
kaN = radial_ab * kappaLR(a_n,b_n)
kPN = radial_ab * kappaLR(a  ,b_n)
kNP = radial_ab * kappaLR(a_n,b)

kaCcP = radial_ab * kappaRL(a  ,b)
kaCcN = radial_ab * kappaRL(a_n,b_n)
kCcPN = radial_ab * kappaRL(a  ,b_n)
kCcNP = radial_ab * kappaRL(a_n,b)

!if (aNeQb) then !do it always, dont sum then
  roPt   = radial_ab * rhoLR  (a,  b)
  roNt   = radial_ab * rhoLR  (a_n,b_n)
  rPNt   = radial_ab * rhoLR  (a  ,b_n)
  rNPt   = radial_ab * rhoLR  (a_n,b)

  kaPt   = radial_ab * kappaLR(b  ,a)
  kaNt   = radial_ab * kappaLR(b_n,a_n)
  kPNt   = radial_ab * kappaLR(b  ,a_n)
  kNPt   = radial_ab * kappaLR(b_n,a)

  kaCcPt = radial_ab * kappaRL(b  ,a)
  kaCcNt = radial_ab * kappaRL(b_n,a_n)
  kCcPNt = radial_ab * kappaRL(b  ,a_n)
  kCcNPt = radial_ab * kappaRL(b_n,a)
!endif

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
    A_part = AngFunctDUAL_HF(ms2,b,a,i_a)
    BulkHF(1,ms,i_r,i_a) = BulkHF(1,ms,i_r,i_a) + (A_part * roPt) !pp
    BulkHF(2,ms,i_r,i_a) = BulkHF(2,ms,i_r,i_a) + (A_part * roNt) !nn
    BulkHF(3,ms,i_r,i_a) = BulkHF(3,ms,i_r,i_a) + (A_part * rPNt) !pn
    BulkHF(4,ms,i_r,i_a) = BulkHF(4,ms,i_r,i_a) + (A_part * rNPt) !np

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
enddo

end subroutine compute_bulkDens4Fields



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
else
  do ms = 1, 4
    ! BulkHF(ms,i_r,i_ang) = BulkHF_pp(ms,i_r,i_ang) + BulkHF_nn(ms,i_r,i_ang)
    BulkHF(5,ms,i_r,i_ang) = BulkHF(1,ms,i_r,i_ang) + BulkHF(2,ms,i_r,i_ang)
    BulkP1(5,ms,i_r,i_ang) = BulkP1(1,ms,i_r,i_ang) + BulkP1(2,ms,i_r,i_ang)
    BulkP2(5,ms,i_r,i_ang) = BulkP2(1,ms,i_r,i_ang) + BulkP1(2,ms,i_r,i_ang)
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

if (.NOT.eval_rearrangement) then
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
      if (integration_method < 2) then
        radial = radial_2b_sho_memo(HOsp_sh(a), HOsp_sh(b),i_r)
      else
        radial = radial_2b_sho_noexp_memo(HOsp_sh(a), HOsp_sh(b),i_r)
      end if

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
! subroutine calculate_expectval_density                                      !
!                                                                             !
! iopt = optimal iteration (=1 when gradient has converged)                   !
!
!-----------------------------------------------------------------------------!
subroutine calculate_expectval_density(rhoLR, kappaLR, kappaRL, &
                                       ndim, iopt)
integer, intent(in) :: ndim, iopt !
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, kappaRL, kappaLR

integer :: a,b, a_sh,la,ja,ma,mta, b_sh,lb,jb,mb,mtb, K, spO2, a_n, b_n
integer :: ind_jm_a, ind_jm_b, ind_km
integer :: M
integer :: i_r=1, i_an=1, msp

real(r64) :: radial_part, dens_R, dens_A, rad4Integr, cga, cgb
complex(r64) :: sum_, integral_dens, sum_test, diff
integer  :: mlb, mla, ms, ind_la, ind_lb
logical :: PRNT_
!! TEST points of the grid (with non integrable variable approximation)
complex(r64), allocatable, dimension(:,:) :: test_dens
complex(r64) :: test_2sph
allocate(test_dens(ndim, ndim))
test_dens = zzero

PRNT_ = (PRINT_GUTS).OR.(.FALSE.)
spO2 = ndim / 2

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

         call compute_bulkDens4Fields(a, b, a_sh, b_sh, i_r, i_an, &
                                      rhoLR, kappaLR, kappaRL, ndim)
!         call compute_bulkDens4Fields_bench(a, b, a_sh, b_sh, i_r, i_an, &
!                                    rhoLR, kappaLR, kappaRL, ndim) ! BENCH REQUIRES B starting at 1
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
      !!  to prompt an abnormal numerical error (r64 must have 13 decimals)
      if (.NOT.DOING_PROJECTION) then
        print *, " !!! [WARNING] density is imaginary=",imag(density(i_r,i_an))
      endif
      ! Fold the density to the 1st quadrant.
      dens_R = dreal(density(i_r,i_an))**2 + dimag((density(i_r,i_an)))**2
      dens_A = alpha_DD * dacos(dreal(density(i_r, i_an)) / dens_R)
      dens_R = dens_R ** alpha_DD

      dens_alpha(i_r,i_an) = dCMPLX(dens_R * dabs(dcos(dens_A)), &
                                    dens_R * dabs(dsin(dens_A)) )
    else
      dens_alpha(i_r,i_an) = dreal(density(i_r,i_an)) ** alpha_DD
      dens_alpm1(i_r,i_an) = dens_alpha(i_r,i_an) / density(i_r,i_an)
      dens_alpm1(i_r,i_an) = MIN(dreal(dens_alpm1(i_r,i_an)), 1.0D+30)
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
if ((iteration.eq.0).OR.(MOD(iteration + 1, 10).EQ.0)) then
  integral_dens = integral_dens * 2 * pi * (HO_b**3) / ((2.0 + alpha_DD)**1.5)
  print "(A,F13.9,A)", "      *A* ", dreal(integral_dens), "  <dens(r)> approx"
endif

!if (export_density.AND.((iteration == 1).OR.(iopt == 1))) then
!  call export_expectval_density(dens_rhoLR, ndim)
!!  call test_integrate_density_me
!end if

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
!! Note :: Remember that radial functions already have the factor 1/b**3
if (integration_method > 1) then
  integral_factor = integral_factor * 0.5d0 * (HO_b**3)
  integral_factor = integral_factor  / ((2.0d0 + alpha_DD)**1.5d0)
  if (integration_method == 3) then ! add Lebedev norm factor
    integral_factor = integral_factor * 4 * pi
  endif
elseif (integration_method == 0) then
  print *, " CANNOT EVALUATE matrix_element_v_DD with TRAPEZOIDAL integration"
  STOP
endif

do i_r = 1, r_dim
  if (integration_method > 1) then
    radial = weight_R(i_r) * radial_2b_sho_noexp_memo(a_sh, c_sh, i_r) &
                           * radial_2b_sho_noexp_memo(b_sh, d_sh, i_r)
  elseif (integration_method == 1) then
    radial = ((r(i_r))**2) * weight_R(i_r) *  &
                radial_2b_sho_memo(a_sh, c_sh, i_r) * &
                radial_2b_sho_memo(b_sh, d_sh, i_r)
  endif

  do i_ang = 1, angular_dim
      call angular_index(i_ang, i_th, i_phi)
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

      angular    = weight_THE(i_th) * weight_PHI(i_phi)

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
      if ((eval_rearrangement).AND.(eval_explicit_fieldsDD)) then
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

TOP = abs(maxval(real(rearrangement_me)))
LOW = abs(minval(real(rearrangement_me)))
if (TOP > 1.0D+10) then !((TOP < 1.0D+10).AND.(LOW > 1.E-10)) then
    print "(A,4I3,2F20.10)", "!! REA", a,b,c,d, &
        minval(real(rearrangement_me)), maxval(real(rearrangement_me))
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
! evaluated if it is not calculated both (eval_explicit_fieldsDD = TRUE) and  !
! (eval_rearrangement = TRUE)                                                 !
!-----------------------------------------------------------------------------!
subroutine calculate_densityDep_hamiltonian(dens_rhoLR, dens_kappaLR, &
                                            dens_kappaRL, ndim)
integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: dens_rhoLR
complex(r64), dimension(ndim,ndim), intent(in) :: dens_kappaLR, dens_kappaRL

integer(i16) :: ared, bred, cred, dred
integer(i32) :: ht, j, t, tmax, uth6=uth+8, uth7=uth+9, fac_ht, ialloc=0, &
                a, ma, la, ta, b, mb, lb, tb, dmax, bmax,  bmin, cmin, dmin,&
                c, mc, lc, tc, d, md, ld, td, aa, bb, cc, dd
integer(i64) :: kk, i, kkk, hdim
integer, parameter :: CONVERG_ITER = 10000
real(r64) :: xja, xjb, xjc, xjd, xjtot, xttot, phasab, phascd, Vtmp, &
             Vcut, Vdec, Vred
real(r64), dimension(4) :: me_Vdec
real(r64), dimension(:), allocatable :: hamil_temp
character(len=25) :: filename
logical   :: ALL_ISOS

ALL_ISOS = (.NOT.eval_explicit_fieldsDD)
!! NOTE: if not explicit evaluation of fields, the process was called to export v_DD
!! if ALL_ISOS = .TRUE., this subroutine can only be called once !!

if (iteration < CONVERG_ITER) then
    rewind(uth6)
    rewind(uth7)
    !!! Computes the two-body matrix elements in m-scheme
    open  (uth6, status='scratch', action='readwrite', access='stream', &
                 form='unformatted')
    open  (uth7, status='scratch', action='readwrite', access='stream', &
                 form='unformatted')
endif

hdim = VSsp_dim / 2
Vcut = 5.0d-14
if (ALL_ISOS) Vcut = 1.0d-10
kk = 0
NOT_DEL_FILE = .FALSE.

rearrang_field = zero

do aa = 1, VSsp_dim / 2 ! (prev = HOsp_dim)
  a  = VStoHOsp_index(aa)
  la = HOsp_l(a)
  ma = HOsp_2mj(a)
  ta = HOsp_2mt(a)

  bmin = aa!+1
  if (evalFullSPSpace) bmin = 1
  do bb = bmin, VSsp_dim / 2 ! (prev = HOsp_dim)
    b  = VStoHOsp_index(bb)
    lb = HOsp_l(b)
    mb = HOsp_2mj(b)
    tb = HOsp_2mt(b)

    if ((.NOT.evalFullSPSpace).AND.( ma + mb < 0 )) cycle

    cmin = aa
    if (evalFullSPSpace) cmin = 1
    do cc = cmin, VSsp_dim / 2 ! (prev = HOsp_dim)
      c  = VStoHOsp_index(cc)
      lc = HOsp_l(c)
      mc = HOsp_2mj(c)
      tc = HOsp_2mt(c)

      dmin = 1
      dmax = VSsp_dim / 2 ! (prev = HOsp_dim)
      if (.NOT.evalFullSPSpace) then
        dmin = cc!+1
        if ( cc == aa ) dmax = bb
      endif
      do dd = dmin, dmax
        d  = VStoHOsp_index(dd)
        ld = HOsp_l(d)
        md = HOsp_2mj(d)
        td = HOsp_2mt(d)

        if ( ta + tb /= tc + td ) cycle

        rearrangement_me = zero

        me_Vdec = matrix_element_v_DD(a,b, c,d, ALL_ISOS)

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

            if ((eval_rearrangement).AND.(eval_explicit_fieldsDD)) then
              call calculate_rearrang_field_explicit(a, b, c, d, Vdec,&
                                                     dens_rhoLR, dens_kappaLR,&
                                                     dens_kappaRL, ndim)
            endif
          endif
        endif ! select the process or to export the matrix elements

      enddo  !end loop b
    enddo  !end loop d
  enddo  !end loop c
enddo  !end loop a

!!! At the first iteration, the values of the hamiltonian are saved via file
if (ALL_ISOS) then

  hamil_DD_H2dim     = kk
  hamil_DD_H2dim_all = kk

  allocate( hamil_DD_H2_byT(4, hamil_DD_H2dim), &
            hamil_DD_abcd(4*hamil_DD_H2dim), hamil_temp(4*hamil_DD_H2dim),&
            stat=ialloc )
  if ( ialloc /= 0 ) stop 'Error during allocation of array of indices [DD]'
  rewind(uth6)
  rewind(uth7)

  read(uth6) (hamil_DD_abcd(kk), kk=1, 4*hamil_DD_H2dim)
  read(uth7) (hamil_temp(kk),    kk=1, 4*hamil_DD_H2dim)

  close(uth6)
  close(uth7)

  do kk = 1, hamil_DD_H2dim
    hamil_DD_H2_byT(1, kk) =  hamil_temp(4*(kk-1) + 1)
    hamil_DD_H2_byT(2, kk) =  hamil_temp(4*(kk-1) + 2)
    hamil_DD_H2_byT(3, kk) =  hamil_temp(4*(kk-1) + 3)
    hamil_DD_H2_byT(4, kk) =  hamil_temp(4*(kk-1) + 4)
  enddo
  deallocate(hamil_temp)

  call print_uncoupled_hamiltonian_DD(ALL_ISOS)

else if (iteration < CONVERG_ITER) then !!! Normal Gradient DD dep. process ****

  if ((iteration > 1).AND.(eval_explicit_fieldsDD)) then
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
write(123, fmt='(A)') "//SING PART INDEX (sp_vs,i_sp, i_sh, n,l,2j,2m, 2mt,tr)"
if (ALL_ISOS) then
  do kk=1, VSsp_dim
    i = VStoHOsp_index(kk)
    write(123, fmt='(I4,8(A,I4))') kk,',',i,',', HOsp_sh(i), &
      ',', HOsp_n(i),',', HOsp_l(i),',', HOsp_2j(i),',', HOsp_2mj(i), &
      ',', HOsp_2mt(i),',', HOsp_tr(i)
  enddo
  write(123, fmt='(3A,3I12)')"//(vs)a    b    c    d              pppp", &
    "              pnpn              pnnp              nnnn    ", &
    "* VS array DD/noDD DIM/ALL=", hamil_DD_H2dim, hamil_H2dim, (VSsp_dim/2)**4
else
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
! subroutine calculate_rearrang_field_explicit                                 !
!                                                                              !
! Evaluate the rearrangement field for non-null DD matrix elements just after  !
! these elements are evaluated, completing the Rearrangement Field one element !
! at a time.                                                                   !
!------------------------------------------------------------------------------!
subroutine calculate_rearrang_field_explicit(ia, ib, ic, id, Vdec,&
                                             rhoLR, kappaLR, kappaRL, ndim)
integer, intent(in)   :: ia, ib, ic, id, ndim
real(r64), intent(in) :: Vdec
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR
complex(r64), dimension(ndim,ndim), intent(in) :: kappaLR, kappaRL
complex(r64) :: aux_rea, Fcomb, aux_reaN
                !f0, f1, f2, f3, & ! field_abcd, f_abdc, f_bacd, f_badc
                !f4, f5, f6, f7, & ! field_cdab, f_cdba, f_dcab, f_dcba
complex(r64), dimension(0:7) :: f
integer   :: HOspo2, it, k1, k2, a, b, c, d, k1N, k2N
integer   :: perm
real(r64) :: rea_h, f2r, sign_tr, rea_hN
logical :: exist_

HOspo2 = ndim / 2
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
subroutine calculate_fields_DD_explicit(rhoLR,kappaLR,kappaRL,gammaLR,hspLR,&
                                          deltaLR,deltaRL,ndim)
integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, kappaLR, kappaRL
complex(r64), dimension(ndim,ndim)             :: gammaLR, hspLR, deltaLR, &
                                                  deltaRL
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

!    !!! Calculation of Delta^10 and Delta^01
!    deltaLR_DD(ib,ia) = deltaLR_DD(ib,ia) + h2b * kappaLR(id,ic)
!    deltaRL_DD(ib,ia) = deltaRL_DD(ib,ia) + h2b * kappaRL(id,ic)
!
!if ((ia==A_print).and.(ic==C_print)) then   !!!==========================
!  aux = dreal(h2b*rhoLR(id,ib))
!  print"(A,2I3,A,F18.15,2A)",g_str,it,BI," [  0]",aux," from <",mestr
!end if    !!!============================================
!    !!! Calculation of Gamma
!    gammaLR_DD(ia,ic)   = gammaLR_DD(ia,ic) + h2b * rhoLR(id,ib)
!    if (evalFullSPSpace) cycle
!
!    if ( ic /= id ) then
!      gammaLR_DD(ia,id) = gammaLR_DD(ia,id) - h2b * rhoLR(ic,ib)
!if ((ia==A_print).and.(id==C_print)) then   !!!===============================
!  aux = dreal(-h2b*rhoLR(ic,ib))
!  print "(A,2I3,A,F18.15,2A)",g_str,it,BI," [1cd]",aux," from <",mestr
!endif   !!!============================================
!    endif
!
!    if ( ib /= ia ) then
!      gammaLR_DD(ib,ic) = gammaLR_DD(ib,ic) - h2b * rhoLR(id,ia)
!if ((ib==A_print).and.(ic==C_print)) then   !!!===============================
!  aux = dreal(-h2b*rhoLR(id,ia))
!  print "(A,2I3,A,F18.15,2A)",g_str,it,BI," [2ba]",aux," from <",mestr
!end if   !!!============================================
!    endif
!
!    if ( (ib /= ia) .and. (ic /= id) ) then
!      gammaLR_DD(ib,id) = gammaLR_DD(ib,id) + h2b * rhoLR(ic,ia)
!if ((ib==A_print).and.(id==C_print)) then   !!!==============================
!  aux = dreal(h2b*rhoLR(ic,ia))
!  print "(A,2I3,A,F18.15,2A)",g_str,it,BI," [  3]",aux," from <",mestr
!end if   !!!============================================
!    endif
!
!    if ( (ia /= ic) .or. (ib /= id) ) then
!      !! Complete the delta Field for <cd | ab>
!      deltaLR_DD(id,ic) = deltaLR_DD(id,ic) + h2b * kappaLR(ib,ia)
!      deltaRL_DD(id,ic) = deltaRL_DD(id,ic) + h2b * kappaRL(ib,ia)
!
!      !! Gamma field for bra-ket change <cd | ab>
!      gammaLR_DD(ic,ia)   = gammaLR_DD(ic,ia) + h2b * rhoLR(ib,id)
!if ((ic==A_print).and.(ia==C_print)) then   !!!===============================
!  aux = dreal(h2b*rhoLR(ib,id))
!  print "(A,2I3,A,F18.15,2A)",g_str,it,BI," [  4]",aux," from <", mestr
!end if   !!!============================================
!
!      if ( ia /= ib ) then
!        gammaLR_DD(ic,ib) = gammaLR_DD(ic,ib) - h2b * rhoLR(ia,id)
!if ((ic==A_print).and.(ib==C_print)) then   !!!===============================
!  aux = dreal(-h2b*rhoLR(ia,id))
!  print "(A,2I3,A,F18.15,2A)",g_str,it,BI," [5ba]",aux," from <",mestr
!end if   !!!============================================
!      endif
!
!      if ( ic /= id ) then
!        gammaLR_DD(id,ia) = gammaLR_DD(id,ia) - h2b * rhoLR(ib,ic)
!if ((id==A_print).and.(ia==C_print)) then   !!!===============================
!  aux = dreal(-h2b*rhoLR(ib,ic))
!  print "(A,2I3,A,F18.15,2A)",g_str,it,BI," [6dc]",aux," from <",mestr
!end if   !!!============================================
!      endif
!
!      if ( (ib /= ia) .and. (ic /= id) ) then
!        gammaLR_DD(id,ib) = gammaLR_DD(id,ib) + h2b * rhoLR(ia,ic)
!if ((id==A_print).and.(ib==C_print)) then   !!!==============================
!  aux = dreal(h2b*rhoLR(ia,ic))
!  print "(A,2I3,A,F18.15,2A)",g_str,it,BI," [  7]",aux," from <",mestr
!end if   !!!============================================
!      endif
!
!    endif

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
    if (eval_rearrangement) then
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
if (.NOT.eval_rearrangement) return
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
      ! pn np part
    aux2  = 2 * dens_pnt(3,i_r,i_a) * dens_pnt(4,i_r,i_a)
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
      aux1  = BulkHF(3,ms2, i_r,i_a) * BulkHF(4,ms,i_r,i_a) !pn*np
      aux2  = BulkHF(4,ms2, i_r,i_a) * BulkHF(3,ms,i_r,i_a) !np*pn
      aux_e = aux_e + (aux1  + aux2)
        !total field part
      aux1  = BulkHF(5,ms2, i_r,i_a) * BulkHF(5,ms,i_r,i_a) !tot
      aux_e = aux_e - (x0_DD_FACTOR * aux1)

      !Pairing rearrangement fields
      if (haveX0M1) then
        aux1  = BulkP2(1,ms, i_r,i_a) * BulkP1(1,ms, i_r,i_a) !pp
        aux2  = BulkP2(2,ms, i_r,i_a) * BulkP1(2,ms, i_r,i_a) !nn
        aux_p = aux_p + (aux1 + aux2)
      endif
      !pn np part (remember the 1Bpn + x0*1Bpn - 1Bnp - x0*1Bnp was done already)
      aux1   = BulkP2(3,ms, i_r,i_a) * BulkP1(3,ms, i_r,i_a) !pn*np
      aux2   = BulkP2(4,ms, i_r,i_a) * BulkP1(4,ms, i_r,i_a) !np*pn
      aux_pnp = aux_pnp + (aux1 + aux2)

    enddo ! loop ms
    !! change 11/11/22 + sings of pairing changed to - (from -k*_ab k_cd)
    !! since (K_RL)* is the definition of the contraction, the sign is + for pairing
    aux1 = (2.0d+0*(aux_d - aux_e)) + (X0M1*aux_p) + aux_pnp
    if ((.FALSE.).AND.(dimag(aux1) > 1.0e-12)) then
      print '(A,2I4,5F20.15)',"Error!! Rearr. funct. is imaginary: ",i_r,i_a,&
        dimag(aux1), dimag(aux_d), dimag(aux_e), dimag(aux_p), dimag(aux_pnp)
    end if
    REACommonFields(i_r, i_a) = aux1   !! dreal(aux1)

    ! export Rea Fields by parts
    if (PRINT_GUTS) then
      write(665, fmt='(2I4,5(F22.15,SP,F20.15,"j"))') &
        i_r, i_a, aux_d, aux_e, aux_p, aux_pnp, REACommonFields(i_r, i_a)
    endif
  enddo
enddo
if (PRINT_GUTS) close(665)

end subroutine calculate_common_rearrang_bulkFields

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
    case (3)        ! pn
      cc = c + spO2
    case (4)        ! np
      aa = a + spO2
    case default    ! pp
      rearrang_field(aa, cc) = int_rea
  end select

  gammaLR_DD(aa, cc) = int_hf(Tac)
  deltaLR_DD(aa, cc) = int_pa(Tac)
  deltaRL_DD(aa, cc) = int_pa(Tac)
  !! sum to the main fields
  gammaLR(aa,cc) = gammaLR(aa,cc) + gammaLR_DD(aa,cc)
  deltaLR(aa,cc) = deltaLR(aa,cc) + deltaLR_DD(aa,cc)
  deltaRL(aa,cc) = deltaRL(aa,cc) + deltaRL_DD(aa,cc)
  hspLR  (aa,cc) = hspLR  (aa,cc) + gammaLR_DD(aa,cc)

  !! copy to the off-diagonal values and sum to the main fields
  if(a.NE.c) then ! all off-diagonal for all (diagonal is already done)
    gammaLR_DD(cc,aa) =  gammaLR_DD(aa,cc)
    deltaLR_DD(cc,aa) = -deltaLR_DD(aa,cc)
    deltaRL_DD(cc,aa) = -deltaRL_DD(aa,cc)

    ! adding to the program fields
    gammaLR(cc,aa) = gammaLR(cc,aa) + gammaLR_DD(cc,aa)
    deltaLR(cc,aa) = deltaLR(cc,aa) + deltaLR_DD(cc,aa)
    deltaRL(cc,aa) = deltaRL(cc,aa) + deltaRL_DD(cc,aa)
    hspLR  (cc,aa) = hspLR  (cc,aa) + gammaLR_DD(cc,aa)
  endif

  ! do rearrangement for pp, nn if proceeds
  if((eval_rearrangement).AND.(Tac.LT.3)) then
    hspLR(aa,cc) = hspLR(aa,cc) + rearrang_field(aa,cc)
    ! this step is just to complete the matrix, the rearrangement is aa,cc as well
    if (a.NE.c) then
      rearrang_field(cc, aa) =  rearrang_field(aa, cc)
      hspLR(cc,aa) = hspLR(cc,aa) + rearrang_field(cc,aa)
    endif
  endif

enddo !Tac

end subroutine complete_DD_fields


!------------------------------------------------------------------------------!
!    subroutine calculate_fields_DD                                            !
! speeds up the bench process by reading half the N/2 space, perform Trace test!
! each 10 iterations and at the first iteration.                               !
!------------------------------------------------------------------------------!
subroutine calculate_fields_DD(rhoLR,kappaLR,kappaRL,gammaLR,hspLR,&
                               deltaLR,deltaRL,ndim)
integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, kappaLR, kappaRL
complex(r64), dimension(ndim,ndim)             :: gammaLR, hspLR, deltaLR, &
                                                  deltaRL
!! The density fields are added to the calculated with the standard hamiltonian
!! This array variables are local
complex(r64), dimension(ndim,ndim) :: gammaLR_DD, deltaLR_DD, deltaRL_DD
complex(r64), dimension(ndim,ndim) :: gam1, dltLR1, dltRL1, hsp1

integer   :: a, b, c, d, i, j, spO2, i_r, i_ang, &
             K, M, ind_jm_a, ind_jm_c, ind_km, Tac, &
             a_sh, ja, la, ma, ta, c_sh, jc, lc, mc, tc, ms, &
             PAR, OpPAR, IP, IDN, OpIDN, R_PRINT=3, ANG_PRINT=10
real(r64) :: rad_ac, X0M1, integral_factor
complex(r64), dimension(4) :: int_hf, int_pa, &
            auxHfD, auxHfE, aux_PE, aux_pair, aux_hf ! all arrays are for (pp, nn, pn, np)

complex(r64) :: sumD_ang, auxRea, int_rea, testaux,del_aux, t_rea_sum,t_gam_sum
complex(r64), dimension(4) :: aux
complex(r64), dimension(10) :: intgama, auxgama
logical :: PRNT_, doTraceTest_

PRNT_ = (PRINT_GUTS).OR.(.FALSE.)
iteration = iteration + 1
doTraceTest_ = (iteration.eq.1).OR.(MOD(iteration, 10).EQ.0)

if (integration_method.NE.3) then
   print *,"[ERROR] Fields calculation can only be performed in Levedeb-Gauss",&
        " (integration_method=3). Program STOPS."
   STOP
endif

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

!! Note :: Remember that radial functions already have the factor 1/b**3
integral_factor = 0.5d0 * (HO_b**3) / ((2.0d0 + alpha_DD)**1.5d0)
integral_factor = 4.0d0 * pi * t3_DD_CONST * integral_factor
X0M1 = 1.0d0 - x0_DD_FACTOR

!! Get the rad-ang function from the current densities for rearrangement
if (eval_rearrangement) then
  call calculate_common_rearrang_bulkFields
endif
if (PRNT_) then
  open(555, file='BulkHF_elements_Exch.gut')
  open(556, file='BulkHF_elements_Dire.gut')
  open(557, file='BulkPA_elements_Exch.gut')

  write(555, fmt='(2A,2I4)')"[auxHf Exc] a   c  ir  ia   %% ", &
    "real(pp)  imagn(pp), nn, pn    I_r,I_ang=", R_PRINT, ANG_PRINT
  write(556, fmt='(2A,2I4)')"[auxHf Dir] a   c  ir  ia  ms   %% ", &
    "real(pp)  imagn(pp), nn, pn, Rearrange  I_r,I_ang=", R_PRINT, ANG_PRINT
  write(557, fmt='(2A,2I4)')"[aux  Pair] a   c  ir  ia  ms   %% ", &
    "real(pp)  imagn(pp), nn, pn    I_r,I_ang=", R_PRINT, ANG_PRINT
endif

do a = 1, spO2
  !! HF field
  a_sh = HOsp_sh(a)

  do c = a, spO2
    c_sh = HOsp_sh(c)

    int_hf = zzero
    int_pa = zzero
    int_rea= zzero

    intgama = zzero

    do i_r = 1, r_dim
      rad_ac = weight_R(i_r) * radial_2b_sho_noexp_memo(a_sh, c_sh, i_r)
      rad_ac = rad_ac * exp((2.0d0+alpha_DD) * (r(i_r)/HO_b)**2)
      do i_ang = 1, angular_dim
        auxHfD = zzero
        !! DIRECT terms for the HF field

        sumD_ang  = AngFunctDUAL_HF(1,a,c,i_ang) + AngFunctDUAL_HF(4,a,c,i_ang)
        auxHfD(1) = dens_pnt(5,i_r,i_ang) - (x0_DD_FACTOR*dens_pnt(1,i_r,i_ang))
        auxHfD(2) = dens_pnt(5,i_r,i_ang) - (x0_DD_FACTOR*dens_pnt(2,i_r,i_ang))
        auxHfD(3) = -x0_DD_FACTOR * dens_pnt(3,i_r,i_ang)
        auxHfD(4) = -x0_DD_FACTOR * dens_pnt(4,i_r,i_ang)
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
          aux(ms) = AngFunctDUAL_HF(ms,a,c, i_ang) * BulkHF(3,ms, i_r,i_ang) !pn
          auxHfE(3)  = auxHfE(3) + aux(ms)
          aux(ms) = AngFunctDUAL_HF(ms,a,c, i_ang) * BulkHF(4,ms, i_r,i_ang) !np
          auxHfE(4)  = auxHfE(4) + aux(ms)

          !! NOTE: Angular 1, 2 functions are defined with direct form of ms,ms'
          if (haveX0M1) then
            aux(ms) = AngFunctDUAL_P2(ms,a,c,i_ang) * BulkP1(1,ms,i_r,i_ang) !pp
            aux_PE(1) = aux_PE(1)  + (X0M1*aux(ms))
            aux(ms) = AngFunctDUAL_P2(ms,a,c,i_ang) * BulkP1(2,ms,i_r,i_ang) !nn
            aux_PE(2) = aux_PE(2)  + (X0M1*aux(ms))
          endif
          !! pn np part, x0 dependence was calculated in BulkP1_**
          aux(ms) = AngFunctDUAL_P2(ms,a,c, i_ang) * BulkP1(3,ms, i_r,i_ang) !pn
          aux_PE(3)  = aux_PE(3)  + aux(ms)
          aux(ms) = AngFunctDUAL_P2(ms,a,c, i_ang) * BulkP1(4,ms, i_r,i_ang) !np
          aux_PE(4)  = aux_PE(4)  + aux(ms)

          if ((PRNT_).AND.(i_r == R_PRINT).AND.(i_ang == ANG_PRINT)) then
            write(555, fmt='(5I4,A,4(F20.15,SP,F20.15,"j"))') a,c,i_r,i_ang,ms,&
              "%%", auxHfE(1),auxHfE(2),auxHfE(3),auxHfE(4)
            write(557, fmt='(5I4,A,4(F20.15,SP,F20.15,"j"))') a,c,i_r,i_ang,ms,&
              "%%", aux_PE(1),aux_PE(2),aux_PE(3),aux_PE(4)
          endif
          !! TEST -----------------------------------------------------------

        enddo ! ms loop

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
        if (eval_rearrangement) then
          auxRea  = REACommonFields(i_r,i_ang) * dens_alpm1(i_r,i_ang)
          auxRea  = auxRea * rea_common_RadAng(a,c, i_r, i_ang)
          auxRea  = auxRea * exp( (2.0d0+alpha_DD) * (r(i_r)/HO_b)**2)
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


if (PRNT_) then
  close(555)
  close(556)
  close(557)
endif

!! do the trace-test and printing of fields each 10 steps
if (doTraceTest_) then
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
      t_rea_sum = t_rea_sum + (rhoLR(i,j) * rearrang_field(j,i))
      t_gam_sum = t_gam_sum + (rhoLR(i,j) * gammaLR_DD(j,i))
      t_gam_sum = t_gam_sum - 0.5d00 * (kappaRL(i,j) * deltaLR_DD(j,i))
    enddo
  enddo
  !!
  !t_rea_sum = (0.25d00 / alpha_DD) * t_rea_sum
  del_aux = (t_gam_sum / t_rea_sum) - (alpha_DD / 2.0d00)
  if (abs(dreal(del_aux)) > 1.0D-8) then
    print '(A,2F15.6,A,F15.6)', "[Warning] Tr(rho*Gam)/= 2a*Tr(rho*Rea) =",&
      dreal(t_gam_sum),dreal(t_rea_sum), &
      ":: rhoGam/rhoRea ::", dreal(t_gam_sum / t_rea_sum)
  endif !!! *********************************************************** DELETE

  if (PRNT_) close(620)
  if (PRNT_) close(621)
endif


end subroutine calculate_fields_DD


!------------------------------------------------------------------------------!
!    subroutine calculate_fields_DD_bench                                      !
! Legacy method to compute directly over all states a,c, expected  !
!------------------------------------------------------------------------------!
subroutine calculate_fields_DD_bench(rhoLR,kappaLR,kappaRL,gammaLR,hspLR,&
                               deltaLR,deltaRL,ndim)
integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, kappaLR, kappaRL
complex(r64), dimension(ndim,ndim)             :: gammaLR, hspLR, deltaLR, &
                                                  deltaRL
!! The density fields are added to the calculated with the standard hamiltonian
!! This array variables are local
complex(r64), dimension(ndim,ndim) :: gammaLR_DD, deltaLR_DD, deltaRL_DD

integer   :: a, b, c, d, i, j, spO2, i_r, i_ang, &
             K, M, ind_jm_a, ind_jm_c, ind_km, Tac, &
             a_sh, ja, la, ma, ta, c_sh, jc, lc, mc, tc, ms, &
             PAR, OpPAR, IP, IDN, OpIDN, R_PRINT=3, ANG_PRINT=10, a_n, c_n
real(r64) :: rad_ac, X0M1, integral_factor
complex(r64), dimension(4) :: int_hf, int_pa, &
            auxHfD, auxHfE, aux_PE, aux_pair, aux_hf ! all arrays are for (pp, nn, pn, np)

complex(r64) :: sumD_ang, auxRea, int_rea, testaux,del_aux, t_rea_sum,t_gam_sum
complex(r64), dimension(4) :: aux
complex(r64), dimension(10) :: intgama, auxgama
logical :: PRNT_

PRNT_ = (PRINT_GUTS).OR.(.TRUE.)

if (integration_method.NE.3) then
   print *,"[ERROR] Fields calculation can only be performed in Levedeb-Gauss",&
        " (integration_method=3). Program STOPS."
   STOP
endif

spO2 = HOsp_dim / 2

gammaLR_DD = zzero
deltaLR_DD = zzero
deltaRL_DD = zzero
rearrang_field = zzero

!! Note :: Remember that radial functions already have the factor 1/b**3
integral_factor = 0.5d0 * (HO_b**3) / ((2.0d0 + alpha_DD)**1.5d0)
integral_factor = 4.0d0 * pi * t3_DD_CONST * integral_factor
X0M1 = 1.0d0 - x0_DD_FACTOR

!! Get the rad-ang function from the current densities for rearrangement
if (eval_rearrangement) then
  call calculate_common_rearrang_bulkFields
endif
if (PRNT_) then
  open(555, file='BulkHF_elements_Exch.gut')
  open(556, file='BulkHF_elements_Dire.gut')
  open(557, file='BulkPA_elements_Exch.gut')

  write(555, fmt='(2A,2I4)')"[auxHf Exc] a   c  ir  ia   %% ", &
    "real(pp)  imagn(pp), nn, pn    I_r,I_ang=", R_PRINT, ANG_PRINT
  write(556, fmt='(2A,2I4)')"[auxHf Dir] a   c  ir  ia  ms   %% ", &
    "real(pp)  imagn(pp), nn, pn, Rearrange  I_r,I_ang=", R_PRINT, ANG_PRINT
  write(557, fmt='(2A,2I4)')"[aux  Pair] a   c  ir  ia  ms   %% ", &
    "real(pp)  imagn(pp), nn, pn    I_r,I_ang=", R_PRINT, ANG_PRINT
endif

do a = 1, spO2
  !! HF field
  a_sh = HOsp_sh(a)
  a_n  = a + spO2

  do c = 1, spO2
    c_sh = HOsp_sh(c)
    c_n  = c + spO2

    int_hf = zzero
    int_pa = zzero
    int_rea= zzero

    intgama = zzero

    do i_r = 1, r_dim
      rad_ac = weight_R(i_r) * radial_2b_sho_noexp_memo(a_sh, c_sh, i_r)
      do i_ang = 1, angular_dim
        auxHfD = zzero
        !! DIRECT terms for the HF field

        sumD_ang  = AngFunctDUAL_HF(1,a,c,i_ang) + AngFunctDUAL_HF(4,a,c,i_ang)
        auxHfD(1) = dens_pnt(5,i_r,i_ang) - (x0_DD_FACTOR*dens_pnt(1,i_r,i_ang))
        auxHfD(2) = dens_pnt(5,i_r,i_ang) - (x0_DD_FACTOR*dens_pnt(2,i_r,i_ang))
        auxHfD(3) = -x0_DD_FACTOR * dens_pnt(3,i_r,i_ang)
        auxHfD(4) = -x0_DD_FACTOR * dens_pnt(4,i_r,i_ang)
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
          aux(ms) = aux(ms) - (x0_DD_FACTOR * BulkHF(5,ms, i_r, i_ang)) !tot
          aux(ms) = aux(ms) * AngFunctDUAL_HF(ms, a,c, i_ang)
          auxHfE(2)  = auxHfE(2) + aux(ms)
            !pn np part
          aux(ms) = AngFunctDUAL_HF(ms,a,c, i_ang) * BulkHF(3,ms, i_r,i_ang) !pn
          auxHfE(3)  = auxHfE(3) + aux(ms)
          aux(ms) = AngFunctDUAL_HF(ms,a,c, i_ang) * BulkHF(4,ms, i_r,i_ang) !np
          auxHfE(4)  = auxHfE(4) + aux(ms)

          !! NOTE: Angular 1, 2 functions are defined with direct form of ms,ms'
          if (haveX0M1) then
            aux(ms) = AngFunctDUAL_P2(ms,a,c,i_ang) * BulkP1(1,ms,i_r,i_ang) !pp
            aux_PE(1) = aux_PE(1)  + (X0M1*aux(ms))
            aux(ms) = AngFunctDUAL_P2(ms,a,c,i_ang) * BulkP1(2,ms,i_r,i_ang) !nn
            aux_PE(2) = aux_PE(2)  + (X0M1*aux(ms))
          endif
          !! pn np part, x0 dependence was calculated in BulkP1_**
          aux(ms) = AngFunctDUAL_P2(ms,a,c, i_ang) * BulkP1(3,ms, i_r,i_ang) !pn
          aux_PE(3)  = aux_PE(3)  + aux(ms)
          aux(ms) = AngFunctDUAL_P2(ms,a,c, i_ang) * BulkP1(4,ms, i_r,i_ang) !np
          aux_PE(4)  = aux_PE(4)  + aux(ms)

          if ((PRNT_).AND.(i_r == R_PRINT).AND.(i_ang == ANG_PRINT)) then
            write(555, fmt='(5I4,A,4(F20.15,SP,F20.15,"j"))') a,c,i_r,i_ang,ms,&
              "%%", auxHfE(1),auxHfE(2),auxHfE(3),auxHfE(4)
            write(557, fmt='(5I4,A,4(F20.15,SP,F20.15,"j"))') a,c,i_r,i_ang,ms,&
              "%%", aux_PE(1),aux_PE(2),aux_PE(3),aux_PE(4)
          endif
          !! TEST -----------------------------------------------------------

        enddo ! ms loop

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
        if (eval_rearrangement) then
          auxRea  = REACommonFields(i_r,i_ang) * dens_alpm1(i_r,i_ang)
          auxRea  = auxRea * rea_common_RadAng(a,c, i_r, i_ang)
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

    gammaLR_DD(a  ,c  ) = int_hf(1) * integral_factor
    gammaLR_DD(a_n,c_n) = int_hf(2) * integral_factor
    gammaLR_DD(a  ,c_n) = int_hf(3) * integral_factor
    gammaLR_DD(a_n,c  ) = int_hf(4) * integral_factor
    deltaLR_DD(a  ,c  ) = int_pa(1) * integral_factor
    deltaLR_DD(a_n,c_n) = int_pa(2) * integral_factor
    deltaLR_DD(a  ,c_n) = int_pa(3) * integral_factor
    deltaLR_DD(a_n,c  ) = int_pa(4) * integral_factor
    deltaRL_DD(a  ,c  ) = deltaLR_DD(a  ,c  )
    deltaRL_DD(a_n,c_n) = deltaLR_DD(a_n,c_n)
    deltaRL_DD(a  ,c_n) = deltaLR_DD(a  ,c_n)
    deltaRL_DD(a_n,c  ) = deltaLR_DD(a_n,c  )

    rearrang_field(a  ,c  ) = int_rea * 0.25d+00 * integral_factor * alpha_DD
    rearrang_field(a_n,c_n) = rearrang_field(a, c)

    if ((dabs(imag(int_hf(1)))> 1.0d-9).OR.(dabs(imag(int_hf(2)))> 1.0d-9).OR.&
        (dabs(imag(int_hf(3)))> 1.0d-9).OR.(dabs(imag(int_hf(4)))> 1.0d-9))then
        print "(A,2I4,A,4F18.12)", "[WARNING] Imaginary part at Gamma DD (", &
            a,c, ") = ", &
            imag(int_hf(1)), imag(int_hf(2)), imag(int_hf(3)), imag(int_hf(4))
    endif
    if ((dabs(imag(int_pa(1)))> 1.0d-9).OR.(dabs(imag(int_pa(2)))> 1.0d-9).OR. &
        (dabs(imag(int_pa(3)))> 1.0d-9).OR.(dabs(imag(int_pa(4)))> 1.0d-9))then
        print "(A,2I4,A,4F18.12)", "[WARNING] Imaginary part at Delta DD (", &
            a,c, ") = ", &
            imag(int_pa(1)), imag(int_pa(2)), imag(int_pa(3)), imag(int_pa(4))
    endif
    if ((dabs(imag(int_rea))> 1.0d-9)) then
        print "(A,2I4,A,F18.12)", "[WARNING] Imaginary part Rearrang.F DD (", &
            a,c, ") = ", imag(int_rea)
    endif

  enddo
enddo
if (PRNT_) then
  close(555)
  close(556)
  close(557)

  open (620, file='fields_matrix.gut')
  write(620,fmt='(A,A,A)') "  i   j        gammaLR       gammaLR_DD        ",&
    "DeltaLR        DeltaLR_DD       DeltaRL        DeltaRL_DD     ",&
    "ReaField_DD         hspLR       hspLR_dd     hspLR_rea"
endif


t_rea_sum = zzero
t_gam_sum = zzero
do i = 1, HOsp_dim
  do j = 1, HOsp_dim

    if (PRNT_) then
      write(620, fmt='(2I4,8F16.9)', advance='no') i, j, dreal(gammaLR(i,j)),&
        dreal(gammaLR_DD(i,j)), dreal(deltaLR(i,j)), dreal(deltaLR_DD(i,j)),&
        dreal(deltaRL(i,j)),dreal(deltaRL_DD(i,j)),dreal(rearrang_field(i,j)),&
        dreal(hspLR(i,j))
    endif

    !! introduced to update the fields with DD contributions
    gammaLR(i,j) = gammaLR(i,j) + gammaLR_DD(i,j)
    deltaLR(i,j) = deltaLR(i,j) + deltaLR_DD(i,j)
    deltaRL(i,j) = deltaRL(i,j) + deltaRL_DD(i,j)

    hspLR(i,j) = hspLR(i,j) + gammaLR_DD(i,j)
    if(PRNT_) write(620, fmt='(A,F12.6)', advance='no') ' ',real(hspLR(i,j))
    if (eval_rearrangement) then
      hspLR(i,j) = hspLR(i,j) + rearrang_field(i,j)
      if(PRNT_) write(620, fmt='(A,F12.6)') ' ',real(hspLR(i,j))
    else
      if(PRNT_) write(620, fmt='(A)') ''
    endif
    !+ gamma only from DD, the other gamma was added in the routine before

    !! DELETE **************************
    !! TEST: Tr(dens * Gamma) = Tr(dens * Rearrange) * alpha
    !if ((i > spO2).AND.(j > spO2)) cycle
    t_rea_sum = t_rea_sum + (rhoLR(i,j) * rearrang_field(j,i))
    t_gam_sum = t_gam_sum + (rhoLR(i,j) * gammaLR_DD(j,i))
    t_gam_sum = t_gam_sum - 0.5d00 * (kappaRL(i,j) * deltaLR_DD(j,i))
  enddo

enddo
!!
!t_rea_sum = (0.25d00 / alpha_DD) * t_rea_sum
del_aux = (t_gam_sum / t_rea_sum) - (alpha_DD / 2.0d00)
if (abs(dreal(del_aux)) > 1.0D-9) then
  print '(A,2F15.6,A,F15.6)', "[Warning] Tr(rho*Gam)!prop Tr(rho*Rea) ::", &
    dreal(t_gam_sum),dreal(t_rea_sum), &
    ":: rhoGam/rhoRea ::", dreal(t_gam_sum / t_rea_sum)
endif !!! *********************************************************** DELETE

if(PRNT_) close(620)
iteration = iteration + 1

end subroutine calculate_fields_DD_bench

!-----------------------------------------------------------------------------!
! subroutine TESTS FOR THE DENSITY, SPHERICAL HARMONICS AND FUNCTIONS         !
!                                                                             !
!                                                                             !
!                                                                             !
!                                                                             !
!-----------------------------------------------------------------------------!

subroutine test_print_basis_quantum_numbers

integer :: i, j, a_sh, na, la, ja, ma

open(621, file='DIMENS_index.gut')
write(621, *)'DIMENSIONS HO_shell dimens'
write(621, fmt='(A,I5)') 'HOsh_dim = ', HOsh_dim
write(621,fmt='(A)') '---------------------------------------------------------'
write(621, fmt='(A)') "i_sh  sh_lab sh_n sh_l sh_2j"
do i = 1, HOsh_dim
    write(621, fmt='(I3,A,I6,A,I4,A,I4,A,I4,A,I3)') i,'  ',HOsh_na(i),' ', &
        HOsh_n(i),' ',HOsh_l(i),' ',HOsh_2j(i)
end do
write(621,fmt='(A)') '---------------------------------------------------------'
write(621, fmt='(A,I3,A,I3,A,I3)') 'HOsh_dim = ', HOsh_dim, '  HO_lmax = ', &
    HO_lmax, '  HO_2jmax = ', HO_2jmax
write(621,*) ''
write(621, *)'DIMENSIONS HO_single particles '
write(621, fmt='(A,I5,A,I5)')'HOsp_dim = ', HOsp_dim,'  HOsp_dim2 = ',HOsp_dim2
write(621,fmt='(A)') '---------------------------------------------------------'
write(621, fmt='(A)') "i_sp  sp_sh sp_n  sp_l  sp_2j  sp_2mj sp_2mT sp_time_rev"
do i=1, HOsp_dim
!  write(621, fmt='(I3,A,I4,A,I5,A,I5,A,I5,A,I5,A,I5,A,I3)') i,'  ',HOsp_sh(i)&
!    ,' ', HOsp_n(i),' ',HOsp_l(i),' ',HOsp_2j(i),'  ',HOsp_2mj(i)&
!    ,'   ',HOsp_2mt(i),'        ',HOsp_tr(i)
  write(621, fmt='(I3,A,I4,A,I5,A,I5,A,I5,A,I5,A,I5,A)') i,',', HOsp_sh(i)&
    ,',', HOsp_n(i),',',HOsp_l(i),',',HOsp_2j(i),',',HOsp_2mj(i)&
    ,',',HOsp_2mt(i),'),'
enddo
close(621)

end subroutine test_print_basis_quantum_numbers



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



subroutine calculate_valenceSpaceReduced(hamilJM, dim_jm, dim_sh)
integer, intent(in) :: dim_sh, dim_jm
real(r64), dimension(4,dim_jm,dim_jm,dim_sh,dim_sh), intent(in) :: hamilJM

integer(i32) :: a, b, aa, a2, b2, a_min, a_max, b_min, b_max, ialloc=0
integer      :: a_ant,b_ant, t, tt, a_sh, b_sh, a_sh_vs, la,lb,&
                J, J_min, J_max, M,&
                ja,jb, ma,mb,ta,tb, ma2,mb2,mta2,mtb2, ja_prev, jb_prev,&
                i_jm, i_sab, Na, Nb, spO2, NormAB, delta_ab, CORE_NUMBER
real(r64) :: aux_t, aux_v, E_core, cgc1, cgc2, cgc_t1, cgc_t2, h2int
real(r64), dimension(:), allocatable :: e_sp_vs,t_sp_vs, T_core, V_core
real(r64), dimension(4) :: Vdd_dec

print *, ""
print *, " [  ] calculate_valenceSpaceReduced"

spO2 = HOsp_dim / 2
allocate(T_core(3), V_core(3))
T_core  = zero
V_core  = zero
allocate(e_sp_vs(VSsh_dim), t_sp_vs(VSsh_dim)) ! assuming different sp.energies for p-n
e_sp_vs = zero
t_sp_vs = zero

a_min   = 0
a_max   = 0
ja_prev = 0
b_min   = 0
b_max   = 0
jb_prev = 0

do a = 1, spO2

  Na = 2*HOsp_n(a) + HOsp_l(a)
  ja = HOsp_2j (a)
  ma = HOsp_2mj(a)
  a_sh = HOsp_sh(a)
  ta = HOsp_2mt(a)

  if ((a_min.EQ.0).OR.(ja.NE.ja_prev)) then ! update the first mj to evaluate b
    a_min = a
    a_max = a + ja
    ja_prev = ja

    a_sh_vs = 0 ! if it's 0, then we do not add up to the VS energy
    do aa = 1, VSsh_dim ! find the index in the VS
      if (VSsh_list(aa).EQ.HOsh_ant(a_sh)) a_sh_vs = aa
    enddo
  endif

  aux_t = hamil_H1(a, a)
  if (Na .GT. NHO_vs) then  ! outer vs outer are neglected/ useless ------------
    cycle
  else if (Na .LE. NHO_co) then    !! Kinetic Energy Core -----------------
    T_core(2+ta) = T_core(2+ta) + (aux_t / sqrt(2*ja + 1.0))
  else if ((Na .LE. NHO_vs).AND.(a_sh_vs.NE.0)) then  !! Valence Space ----
    t_sp_vs(a_sh_vs) = t_sp_vs(a_sh_vs) + (aux_t / sqrt(2*ja + 1.0))
  endif   !!    --------

  !! Calculate the 2Body Interaction for the CORE and the VALENCE
  do b = a_min, spO2
    Nb = 2*HOsp_n(b) + HOsp_l(b)
    jb = HOsp_2j (b)
    mb = HOsp_2mj(b)

    delta_ab = 0
    if ((ja.EQ.jb).AND.(la.EQ.lb).AND.(HOsp_n(a).EQ.HOsp_n(b))) delta_ab = 1

    if (Nb .GT. NHO_vs) cycle ! outer vs outer are neglected/ useless
    if ((b_min.EQ.0).OR.(jb.NE.jb_prev)) then ! update the first mj to evaluate b
      b_min = b
      b_max = b + jb
      jb_prev = jb
    endif

    M = (ma + mb) / 2
    J_min = max(M, max(abs(ja - jb)/2, 0))  ! i.e. cannot have J=1, M=+-3
    J_max = (ja + jb) / 2

    do J = J_min, J_max
      call ClebschGordan(ja,jb,2*J, ma,mb,2*M, cgc1)

      do a2 = a_min, a_max ! loop for the second CG
        ma2 = HOsp_2mj(a2)
        mb2 = 2*M - HOsp_2mj(a2)
        b2  = b_min + (jb - mb2) / 2
        if ((b2 .LT. b_min).OR.(b2 .GT. b_max)) cycle ! INVALID mb2 value

        call ClebschGordan(ja,jb,2*J, ma2, mb2,2*M, cgc2)

        Vdd_dec = matrix_element_v_DD(a, b, a2, b2, .TRUE.)

        !! T = 0
        if (delta_ab.EQ.0) NormAB = one
        if (delta_ab.EQ.1) then
          NormAB = one / 2
          if (MOD(J, 2).EQ.0) NormAB = zero
        endif
        aux_v = NormAB * cgc1 * cgc2 * sqrt(2*J + 1.0)
        if (Nb .LE. NHO_co) then !! CORE PART :
          V_core(2) = V_core(2) + (aux_v * (Vdd_dec(2) - Vdd_dec(3)))
        else if (a_sh_vs.NE.0) then ! -------- !! VALENCE SPACE SP Energies :
          aux_v = aux_v * (Vdd_dec(2) - Vdd_dec(3))
          e_sp_vs(a_sh_vs) = e_sp_vs(a_sh_vs) + aux_v
        endif

        !! T = 1
        if (delta_ab.EQ.1) then
          NormAB = one / 2
          if (MOD(J, 2).EQ.1) NormAB = zero
        endif

        aux_v = NormAB * cgc1 * cgc2 * sqrt((2*J + 1.0) * 3)
        if (Nb .LE. NHO_co) then !! CORE PART :
          V_core(1) = V_core(1) + (aux_v *  Vdd_dec(1))
          V_core(2) = V_core(2) + (aux_v * (Vdd_dec(2) + Vdd_dec(3)))
          V_core(3) = V_core(3) + (aux_v *  Vdd_dec(4))
        else if (a_sh_vs.NE.0) then ! ---------- !! VALENCE SPACE SP Energies :
          e_sp_vs(a_sh_vs) = e_sp_vs(a_sh_vs) + aux_v * Vdd_dec(1)
          e_sp_vs(a_sh_vs) = e_sp_vs(a_sh_vs) + aux_v *(Vdd_dec(2) + Vdd_dec(3))
          e_sp_vs(a_sh_vs) = e_sp_vs(a_sh_vs) + aux_v * Vdd_dec(4)
        endif

      enddo ! loop the other m_j
    enddo ! loop J

  enddo
enddo

!! SUM the DENSITY INDEPENDENT HAMILTONIAN (shell indexes)
CORE_NUMBER = 0
do a_sh = 1, HOsh_dim
  Na = 2*HOsh_n(a_sh) + HOsh_l(a_sh)
  ja = HOsh_2j (a_sh)

  if     (Na .GT. NHO_vs) then
    cycle
  elseif (Na .LE. NHO_co) then
    CORE_NUMBER = CORE_NUMBER + (ja + 1)
  endif

  a_sh_vs = 0
  do aa = 1, VSsh_dim ! find the index in the VS
    if (VSsh_list(aa).EQ.HOsh_ant(a_sh)) a_sh_vs = aa
  enddo
  do b_sh = a_sh, HOsh_dim
    Nb = 2*HOsh_n(b_sh) + HOsh_l(b_sh)
    jb = HOsh_2j(b_sh)
    if (Nb .GT. NHO_vs) cycle ! outer vs outer are neglected/ useless

    delta_ab = 0
    if ((ja.EQ.jb).AND.(la.EQ.lb).AND.(HOsh_n(a_sh).EQ.HOsh_n(b_sh))) then
      delta_ab = 1
      endif

    J_min = max(abs(ja - jb)/2, 0)
    J_max = (ja + jb) / 2
    do J = J_min, J_max

      h2int =         hamil_H2cpd_DD(1, J, a_sh, b_sh, a_sh, b_sh)
      h2int = h2int + hamil_H2cpd_DD(2, J, a_sh, b_sh, a_sh, b_sh)
      h2int = h2int + hamil_H2cpd_DD(3, J, a_sh, b_sh, a_sh, b_sh)
      h2int = h2int + hamil_H2cpd_DD(4, J, a_sh, b_sh, a_sh, b_sh)
      !! T = 1,2,3,4 (pnpn)
      if (delta_ab.EQ.0) NormAB = one

      aux_v = NormAB * sqrt(2*J + 1.0)
      if (Nb .LE. NHO_co) then !! CORE PART :
        V_core(2) = V_core(2) + (aux_v * h2int)
      else if (a_sh_vs.NE.0) then  ! --------- !! VALENCE SPACE SP Energies :
        e_sp_vs(a_sh_vs) = e_sp_vs(a_sh_vs) + (aux_v * h2int)
      endif

      !! pppp, nnnn
      if (delta_ab.EQ.1) then
        NormAB = one / 2
        if (MOD(J, 2).EQ.0) NormAB = zero
      endif
      aux_v = NormAB * sqrt(2*J + 1.0)

      h2int = hamil_H2cpd_DD(0, J, a_sh, b_sh, a_sh, b_sh)
      if (Nb .LE. NHO_co) then !! CORE PART :
        V_core(1) = V_core(1) + (aux_v * h2int)
      else if (a_sh_vs.NE.0) then  ! --------- !! VALENCE SPACE SP Energies :
        e_sp_vs(a_sh_vs) = e_sp_vs(a_sh_vs) + (aux_v * h2int)
      endif

      h2int = hamil_H2cpd_DD(5, J, a_sh, b_sh, a_sh, b_sh)
      if (Nb .LE. NHO_co) then !! CORE PART :
        V_core(3) = V_core(3) + (aux_v * h2int)
      else if (a_sh_vs.NE.0) then  ! --------- !! VALENCE SPACE SP Energies :
        e_sp_vs(a_sh_vs) = e_sp_vs(a_sh_vs) + (aux_v * h2int)
      endif

    enddo ! sum J

  enddo

  if (a_sh_vs.NE.0) then
    e_sp_vs(a_sh_vs) = t_sp_vs(a_sh_vs) + (0.5d0 * e_sp_vs(a_sh_vs)/(ja+1.0d0))
  endif

enddo

E_core = zero
do tt = 1, 3
  E_core  = E_core + T_core(tt) + (1.0d0 * V_core(tt)) !! we sum all
enddo

open (297, file="D1S_vs_scalar.sho")
write(297, fmt='(2A,F9.3,A,F10.5,A,F5.3,A,2I5)') &
  'Density 2BME on explicit HFB wf from taurus, Scalar', &
  ' PARAMS:: t3=',t3_DD_CONST,' MeV  X0=', x0_DD_FACTOR, ' ALPHA=', alpha_DD, &
  '  CORE(n,p):', CORE_NUMBER, CORE_NUMBER
write(297, fmt="(2I4,F12.6)") INT(CORE_NUMBER), INT(CORE_NUMBER), E_core

do a_sh_vs = 1, VSsh_dim
  a_ant = VSsh_list(a_sh_vs)
  write(297, fmt="(I8)", advance='no') a_ant
enddo
write(297,*) ""
do a_sh_vs = 1, VSsh_dim
  write(297, fmt="(F12.6)", advance='no') e_sp_vs(a_sh_vs)
enddo
close(297)

deallocate(T_core, V_core, e_sp_vs, t_sp_vs)

print *,  " [OK] calculate_valenceSpaceReduced"
print *, ""

end subroutine calculate_valenceSpaceReduced


!------------------------------------------------------------------------------!
! subroutine  recouple_jjLSConjugatedME                                        !
!    recouple with the conjugated elements, include the hamiltonian_ state and !
!    and verify if the current KK state is null or not.                        !
!------------------------------------------------------------------------------!
subroutine recouple_jjLSConjugatedME(a,b,c,d, a_con,b_con,c_con,d_con, &
                                     Sbra,Sket,Lbra,Lket,Jbra,Jket,Mbra,Mket,&
                                     hamilJM, dim_jm, dim_sh, TENSOR_ORD, &
                                     factor, auxHamilRed, ind_k, kval_is_zero)

integer, intent(in)  :: a,b,c,d, a_con,b_con,c_con,d_con, Sbra,Sket,Lbra,Lket,&
                        Jbra,Jket, Mbra,Mket
integer, intent(in)  :: dim_sh, dim_jm, TENSOR_ORD, ind_k
real(r64), intent(in):: factor !! factor of the 6J*9J (and all the phases)
real(r64), dimension(4,dim_jm,dim_jm,dim_sh,dim_sh), intent(in) :: hamilJM
real(r64), dimension(4,0:TENSOR_ORD,dim_jm,dim_jm) :: auxHamilRed
logical, intent(out) :: kval_is_zero

integer :: i, j, aa, bb, cc, dd, ind_jm_b, ind_jm_k, ind_sab, ind_scd, &
           ind_sab_2, ind_scd_2, tt
real(r64) :: aux_val, aux1, aux2
integer  , dimension(4)   :: ind_sh_ab, ind_sh_cd
real(r64), dimension(4)   :: aux_r_ab, aux_r_cd
logical,   dimension(4)   :: j_isitsconjugate

j_isitsconjugate = (/ a.EQ.a_con, b.EQ.b_con, c.EQ.c_con, d.EQ.d_con /)
ind_jm_b  = angular_momentum_index(Jbra, Mbra, .FALSE.)
ind_jm_k  = angular_momentum_index(Jket, Mket, .FALSE.)
aux_r_ab = zero
aux_r_cd = zero
kval_is_zero = .TRUE.
!! a, a_con, .. are shell states
!! aux_r_ab = [(a, b), (a, b_con), (a_con, b), (a_con, b_con)]

!! case 1 is always required  -------------------------------------------------
call Wigner9JCoeff(2*HOsh_l(a), 1,HOsh_2j(a), 2*HOsh_l(b), 1,HOsh_2j(b),&
                   2*Lbra, 2*Sbra, 2*Jbra, aux1)
aux_r_ab (1) = sqrt((HOsh_2j(a)+1)*(HOsh_2j(b)+1)*(2*Sbra + 1)*(2*Lbra + 1.0d0))
aux_r_ab (1) = aux_r_ab(1) * aux1
ind_sh_ab(1) = two_shell_states_index(a, b)

call Wigner9JCoeff(2*HOsh_l(c), 1,HOsh_2j(c), 2*HOsh_l(d), 1,HOsh_2j(d),&
                   2*Lket, 2*Sket, 2*Jket, aux2)
aux_r_cd (1) = sqrt((HOsh_2j(c)+1)*(HOsh_2j(d)+1)*(2*Sket + 1)*(2*Lket + 1.0d0))
aux_r_cd (1) = aux_r_cd(1) * aux2
ind_sh_cd(1) = two_shell_states_index(c, d)

!! cases for complementary j, avoid calculating twice if j=j_con_ ------------
do i = 2, 4
  aa = a
  bb = b
  cc = d
  dd = d
  if ((i .EQ. 2) .OR. (i .EQ. 4)) then
    bb = b_con
    if ((b .EQ. b_con) .OR. ((i.EQ.4).AND.(a .EQ. a_con))) then
      ind_sh_ab(i) = 0
    else
      ind_sh_ab(i) = two_shell_states_index(aa, bb)
    end if

    dd = d_con
    if ((d .EQ. d_con) .OR. ((i.EQ.4).AND.(c .EQ. c_con))) then
      ind_sh_cd(i) = 0
    else
      ind_sh_cd(i) = two_shell_states_index(cc, dd)
    end if
  endif
  !!
  if ((i .EQ. 3) .OR. (i .EQ. 4)) then
    aa = a_con
    if ((a .EQ. a_con) .OR. ((i.EQ.4).AND.(a .EQ. a_con))) then
      ind_sh_ab(i) = 0
    else
      ind_sh_ab(i) = two_shell_states_index(aa, bb)
    end if

    cc = c_con
    if ((c .EQ. c_con) .OR. ((i.EQ.4).AND.(d .EQ. d_con))) then
      ind_sh_cd(i) = 0
    else
      ind_sh_cd(i) = two_shell_states_index(cc, dd)
    end if

  endif

  !! If the index are zero, then the state is repeated with a previous one
  if (ind_sh_ab(i) .NE. 0) then
    call Wigner9JCoeff(2*HOsh_l(aa), 1,HOsh_2j(aa), &
                       2*HOsh_l(bb), 1,HOsh_2j(bb), &
                       2*Lbra, 2*Sbra, 2*Jbra, aux1)
    aux_r_ab(i) = sqrt((HOsh_2j(aa) + 1)*(HOsh_2j(bb) + 1)* &
                       (2*Sbra + 1)*(2*Lbra + 1.0d0))
    aux_r_ab(i) = aux_r_ab(i) * aux1
  endif
  if (ind_sh_cd(i) .NE. 0) then
    call Wigner9JCoeff(2*HOsh_l(cc), 1,HOsh_2j(cc), &
                       2*HOsh_l(dd), 1,HOsh_2j(dd), &
                       2*Lket, 2*Sket, 2*Jket, aux2)
    aux_r_cd(i) = sqrt((HOsh_2j(cc) + 1)*(HOsh_2j(dd) + 1)* &
                       (2*Sket + 1)*(2*Lket + 1.0d0))
    aux_r_cd(i) = aux_r_cd(i) * aux2
  endif
enddo

!! sum all non
do i =  1, 4
  aa = ind_sh_ab(i)
  if (aa .EQ. 0) cycle
  do j = 1, 4
    bb = ind_sh_cd(j)
    if (bb .EQ. 0) cycle

    do tt = 1, 4
      aux_val = aux_r_ab(i)*aux_r_cd(j)*hamilJM(tt, ind_jm_b, ind_jm_k, aa, bb)
      aux_val = aux_val * factor

      auxHamilRed(tt,ind_k, ind_jm_b,ind_jm_k) = &
          auxHamilRed(tt,ind_k, ind_jm_b,ind_jm_k) + aux_val

      if (abs(aux_val) .GT. 1.0e-10) kval_is_zero = .FALSE.
    enddo
  enddo
enddo

end subroutine recouple_jjLSConjugatedME


!------------------------------------------------------------------------------!
! subroutine print_DD_matrix_elements                                          !
!  Export of the DD + Hamil Matrix elements for a Valence space after process  !
!------------------------------------------------------------------------------!
subroutine print_DD_matrix_elements
integer(i32) :: a, b, c, d
integer      :: a_ant,b_ant,c_ant,d_ant, t, tt, &
                Jbra, Jket, Jb_min, Jb_max,Jk_min, Jk_max, Mbra, Mket,&
                ja,jb,jc,jd, ma,mb,mc,md,ta,tb,tc,td, Na,Nb,Nc,Nd,&
                ind_jm_b, ind_jm_k, ind_sab, ind_scd, delta_ab, delta_cd,&
                TENSOR_ORD = 2, KK, MM, KKmin, KKmax, &
                a_eq_b, a_eq_c, c_eq_d, b_eq_d, J1,J2,J3,J4, &
                aa, bb, cc, dd, &
                Sbra,Sket,Lbra,Lket,Lb_min,Lk_min,Lb_max,Lk_max,Sk_min,Sk_max,&
                a_con,b_con,c_con,d_con, dim_jm,dim_sh
logical      :: kval_is_zero

real(r64)    :: aux_val,aux_val2, cgc1,cgc2,norm, recoupl_factor,&
                TOL=1.0e-10, aux_1, aux_2, aux_3, aux_4, phs_pnk0
integer, dimension(:), allocatable :: reciprocal_nlj_shell
real(r64), dimension(4)            :: h2b  ! [pppp, pnpn, pnnp, nnnn]
logical      :: in_v_space
logical,   dimension(:), allocatable :: all_zero            ! all_zero(K)
real(r64), dimension(:,:,:,:,:), allocatable :: hamilJM     ! H2JM(T,JMbra, JMket, jajb, jcjd))
real(r64), dimension(:,:,:,:),   allocatable :: auxHamilRed ! H2JM(T, K, JMbra, JMket)



print *, ""
print *, "* [  ] Printing 2B Matrix elements DD from WF_HFB /dim H2_DD:", &
    hamil_DD_H2dim

open(298, file="D1S_vs_red.2b")
open(299, file="D1S_vs_scalar.2b")
open(300, file="onlyDD_D1S_scalar.2b")
open(301, file="onlyDD_D1S_k1.2b")
open(302, file="onlyDD_D1S_k2.2b")
do KK = 0, TENSOR_ORD
  write(300+KK, '(A,I2,A,F10.5,A,F10.5,A,F6.4)') &
    'Density 2BME on explicit HFB wf from taurus, Tensor=',KK, &
    ' PARAMS:: t3=',t3_DD_CONST,' MeV  X0=', x0_DD_FACTOR, ' ALPHA=', alpha_DD
enddo
write(298, '(A,A,F9.3,A,F10.5,A,F5.3,A,2F5.1)') &
    'Density 2BME on explicit HFB wf from taurus, Scalar', &
    ' PARAMS:: t3=',t3_DD_CONST,' MeV  X0=',x0_DD_FACTOR,' ALPHA=',alpha_DD, &
    '  CORE(n,p):', valence_N, valence_Z
write(299, '(A,A,F9.3,A,F10.5,A,F5.3,A,2F5.1)') &
    'Density 2BME on explicit HFB wf from taurus, Scalar', &
    ' PARAMS:: t3=',t3_DD_CONST,' MeV  X0=',x0_DD_FACTOR,' ALPHA=',alpha_DD, &
    '  CORE(n,p):', valence_N, valence_Z

!! allocate the big JM, Jm', ab, cd array for the matrix elements
dim_jm = angular_momentum_index(2*HO_2jmax,2*HO_2jmax,.FALSE.)
dim_sh = two_shell_states_index(maxval(HOsp_sh),minval(HOsp_sh))

allocate(hamilJM(4, dim_jm,dim_jm, dim_sh,dim_sh))
hamilJM = zzero
allocate(auxHamilRed(4,0:TENSOR_ORD, dim_jm,dim_jm))

print *, " ----------------------------------------------- "
print "(A,2I5)"," * Max vals: 2sh, 2jmax", maxval(HOsp_sh), 2*HO_2jmax
print "(A,2I5)","Dimensions hamilJM [4] [jm,]2 [dim_sh,]2 :", dim_jm, dim_sh
print "(A,3I5)","Dimensions auxHamilRed [4] [k] [jm,]2    :", TENSOR_ORD,dim_jm
print *, " ----------------------------------------------- "
!! define the reciprocal shells (j=l+1/2 -> j'=l-1/2) -----------------------
allocate(reciprocal_nlj_shell(VSsh_dim))
!print "(A)", "[TEST] Reciprocal shells"
do aa = 1, VSsh_dim
  a = VStoHOsh_index(aa)
  reciprocal_nlj_shell(aa) = aa !! (default value)
  do bb = 1, VSsh_dim
    b = VStoHOsh_index(bb)
    if ((HOsh_n(a) .EQ. HOsh_n(b)).AND.(HOsh_l(a) .EQ. HOsh_l(b))) then
      if ((HOsh_2j(a) .NE. HOsh_2j(b)).OR.(HOsh_l(a) .EQ. 0)) then
        reciprocal_nlj_shell(aa) = bb
        exit
      else
        cycle
      end if
    end if
  enddo
  bb = reciprocal_nlj_shell(aa)
!  print "(A,2I7,A,2I7)", "  * ", aa, VSsh_list(aa), " -> ", bb, VSsh_list(bb)
enddo      !!! --------------------------------------------------------------

!! TODO: The HamilJM can be reduced to just Core + VS -Shells
!        (saving the memory and time required by the outer shells)
do KK = 1, hamil_DD_H2dim
  a = hamil_DD_abcd(1+4*(KK-1))
  b = hamil_DD_abcd(2+4*(KK-1))
  c = hamil_DD_abcd(3+4*(KK-1))
  d = hamil_DD_abcd(4+4*(KK-1))
  do tt=1,4
    h2b(tt) = hamil_DD_H2_byT(tt, KK)
  enddo

  ! jump elements under another additional tolerace
  !if dabs(h2b) < TOL cycle  ! USELESS, elements in H_DD_H2 were non null
  delta_ab = 0
  if (a_ant==b_ant) delta_ab = 1
  delta_cd = 0
  if (c_ant==d_ant) delta_cd = 1

  ja = HOsp_2j(a)
  jb = HOsp_2j(b)
  jc = HOsp_2j(c)
  jd = HOsp_2j(d)
  ind_sab  = two_shell_states_index(HOsp_sh(a), HOsp_sh(b))
  ind_scd  = two_shell_states_index(HOsp_sh(c), HOsp_sh(d))

  ma = HOsp_2mj(a)
  mb = HOsp_2mj(b)
  mc = HOsp_2mj(c)
  md = HOsp_2mj(d)

  print "(A,I6,A,4I3,A,I2,A,8(A,2I3,A))", "kk=",kk," [",a,b,c,d,"] tt=",tt, &
    " jm_abcd=","(",ja,ma,")", "(",jb,mb,")", "(",jc,mc,")", "(",jd,md,")"

  Mbra = (ma + mb) / 2
  Mket = (mc + md) / 2

  Jb_min = abs(ja - jb) / 2
  Jb_max =    (ja + jb) / 2
  Jk_min = abs(jc - jd) / 2
  Jk_max =    (jc + jd) / 2

  print "(A,2I3,A,3I3,A,3I3,A)", " *ind_j_ab, cd=", ind_sab, ind_scd, &
    " JM,J'M', range=[", Jb_min,Jb_max,Mbra, "]   [", Jk_min,Jk_max, Mket,"]"

  do Jbra = Jb_min, Jb_max
    if (abs(Mbra) > Jbra) cycle
    call ClebschGordan(ja,jb,2*Jbra, ma,mb,2*Mbra, cgc1)
    do Jket = Jk_min, Jk_max
      if (abs(Mket) > Jket) cycle
      call ClebschGordan(jc,jd,2*Jket, mc,md,2*Mket, cgc2)

      ind_jm_b = angular_momentum_index(Jbra, Mbra, .FALSE.)
      ind_jm_k = angular_momentum_index(Jket, Mket, .FALSE.)

      norm = sqrt((1.0d0 + delta_ab*((-1)**Jbra))*(1 + delta_cd*((-1)**Jket)))
      norm = norm / ((1 + delta_ab) * (1 + delta_cd))

      print "(A,2I3,A,2I3,A,2I4,A,4F15.9)","   (jajb),JM,JM'(",Jbra,Mbra, &
          ")(",Jket,Mket,") ind_jm_bra, ket=", ind_jm_b, ind_jm_k, " +=",&
          cgc1 , cgc2, h2b(2), cgc1 * cgc2 * h2b(2)
      do tt = 1, 4
        aux_val = cgc1 * cgc2 * h2b(tt)
        if ((tt .NE. 2).AND.(tt .NE. 3)) aux_val = aux_val * norm

        hamilJM(tt,ind_jm_b, ind_jm_k, ind_sab, ind_scd) = &
                hamilJM(tt, ind_jm_b, ind_jm_k, ind_sab, ind_scd) + aux_val
      enddo
    enddo
  enddo
  print  *, ''
enddo !k
print *, " *** I have read the full uncoupled hamiltonian. Now [step 2]"

! -------------------------------------------------------------------------
!! 2. Calculate SP-Energy and core energy from Global Hamil and DD hamil
call calculate_valenceSpaceReduced(hamilJM,dim_jm,dim_sh)

print *, " *** I have evaluate the valence Space. Now [step 3]"
!return
! -------------------------------------------------------------------------
!! 3 export the matrix elemnts in antoine format, also normalize by J
allocate(all_zero(0:TENSOR_ORD))
do aa = 1, VSsh_dim
  a  = VStoHOsh_index(aa)
  ja = HOsh_2j(a)
  a_ant = VSsh_list(aa)
  do bb = aa, VSsh_dim
    b  = VStoHOsh_index(bb)
    jb = HOsh_2j(b)
    b_ant = VSsh_list(bb)
    ind_sab = two_shell_states_index(a, b)

    do cc = aa, VSsh_dim
      c  = VStoHOsh_index(cc)
      jc = HOsh_2j(c)
      c_ant = VSsh_list(cc)
      do dd = cc, VSsh_dim
        d  = VStoHOsh_index(dd)
        jd = HOsh_2j(d)
        d_ant = VSsh_list(dd)
        ind_scd = two_shell_states_index(c, d)

  !! ======= Loop  for the <ab cd> states to output =======================
  ! this only see the (n,l,j) equivalence, the particle part is in the last step

  Jb_min = abs(ja - jb) / 2
  Jb_max = (ja + jb) / 2
  Jk_min = abs(jc - jd) / 2
  Jk_max = (jc + jd) / 2

  !! ======= Extract the simpler form of the scalar D1S  ==================
  auxHamilRed = zero
  kval_is_zero = .TRUE.
  do Jbra = max(Jb_min, Jk_min), min(Jb_max, Jk_max)
    do Mbra = -Jbra, Jbra
      ind_jm_b = angular_momentum_index(Jbra, Mbra, .FALSE.)
      ind_jm_k = angular_momentum_index(Jbra, 0, .FALSE.) ! for auxHamil to save
      do t = 1, 4
        aux_val = hamilJM(t, ind_jm_b, ind_jm_b, ind_sab, ind_scd)
        if (dabs(aux_val) .GT. TOL) then
          !aux_val = aux_val * sqrt(2*Jbra + 1.0d0) ! factor for the Reduced ME
          auxHamilRed(t,0,ind_jm_k,ind_jm_k) = &
              auxHamilRed(t,0,ind_jm_k,ind_jm_k) + aux_val

          kval_is_zero = .FALSE.
        endif
      end do
    end do
  enddo
  if (.NOT.kval_is_zero) then
    write(298, fmt='(A,4I8,2I3)') ' 0 5', a_ant,b_ant,c_ant,d_ant, &
                                    max(Jb_min,Jk_min), min(Jb_max,Jk_max)

    do Jbra = max(Jb_min, Jk_min), min(Jb_max, Jk_max)

      ind_jm_b = angular_momentum_index(Jbra, 0, .FALSE.)

      aux_1 = auxHamilRed(1,0,ind_jm_b,ind_jm_b)
      write(298,fmt='(F15.10)',advance='no') &
        aux_1 !+ hamil_H2cpd_DD(0, Jbra, a,b,c,d)
      aux_2 = auxHamilRed(2,0,ind_jm_b,ind_jm_b)
      aux_3 = auxHamilRed(3,0,ind_jm_b,ind_jm_b)
      write(298,fmt='(4F15.10)',advance='no') &
        aux_2 ,&! + hamil_H2cpd_DD(1, Jbra, a,b,c,d), &
        aux_3 ,&!+ hamil_H2cpd_DD(2, Jbra, a,b,c,d), &
        aux_3 ,&!+ hamil_H2cpd_DD(3, Jbra, a,b,c,d), &
        aux_2 !+ hamil_H2cpd_DD(4, Jbra, a,b,c,d)
      aux_4 = auxHamilRed(4,0,ind_jm_b,ind_jm_b)
      write(298,fmt='(F15.10)', advance='no') &
        aux_4 !+ hamil_H2cpd_DD(5, Jbra, a,b,c,d)
      write(298,*) ''
    enddo
  endif

  !! ======= Evaluate the rearrange for tensor components on the D1S
!  print *, ""
!  print "(A,4I5,2(A,2I3))", " abcd ", a_ant,b_ant,c_ant,d_ant, " lims bra:",&
!    Jb_min,Jb_max, " ket:", Jk_min,Jk_max
  auxHamilRed = zero
  do Jbra = Jb_min, Jb_max
    Mbra = 0
    ind_jm_b = angular_momentum_index(Jbra, Mbra, .FALSE.)
    J1 = angular_momentum_index(Jbra, 0, .FALSE.)

    do KK = 0, TENSOR_ORD
      do Jket = Jk_min, Jk_max
        Mket = 0
        ind_jm_k = angular_momentum_index(Jket, Mket, .FALSE.)
        J2 = angular_momentum_index(Jket, 0, .FALSE.)

        do Sbra = 0, 1
          Lb_min = abs(Jbra - Sbra)
          Lb_max =     Jbra + Sbra

          do Lbra = Lb_min, Lb_max
            Sk_min = max(abs(Sbra -   KK), 0 )
            Sk_max = min(    Sbra +   KK , 1 )

            do Sket = Sk_min, Sk_max
              Lk_min = max(abs(Jket - Sket), abs(Jbra - Sket), abs(Lbra - KK))
              Lk_max = min(    Jket + Sket ,     Jbra + Sket ,     Lbra - KK )

              do Lket = Lk_min, Lk_max
                !!! Calculate the 6j coeffs LS-J to k on bra and ket
                call Wigner6JCoeff(2*Lbra, 2*Sbra, 2*Jbra, &
                                   2*Sket, 2*Lket, 2*KK  , aux_1)
                call Wigner6JCoeff(2*Lbra, 2*Sbra, 2*Jket, &
                                   2*Sket, 2*Lket, 2*KK  , aux_2)
                if (abs(aux_1*aux_2) .LE. TOL) cycle

                !!! Calculate the LS.J recoupling_ for coefficients
                call Wigner9JCoeff(2*HOsh_l(a),      1,     ja, &
                                   2*HOsh_l(b),      1,     jb, &
                                        2*Lbra, 2*Sbra, 2*Jbra, aux_3)
                if (abs(aux_3) .LE. TOL) cycle
                call Wigner9JCoeff(2*HOsh_l(c),      1,     jc, &
                                   2*HOsh_l(d),      1,     jd, &
                                        2*Lket, 2*Sket, 2*Jbra, aux_4)
                if (abs(aux_4) .LE. TOL) cycle

                aux_3 = aux_3 * sqrt((ja+1)*(jb+1)*(2*Sbra + 1)*(2*Lbra + 1.0))
                aux_4 = aux_4 * sqrt((jc+1)*(jd+1)*(2*Sket + 1)*(2*Lket + 1.0))

  !!! ================================================================
  ! part to recouple with the j and conjugate j elements
  a_con = reciprocal_nlj_shell(a)
  b_con = reciprocal_nlj_shell(b)
  c_con = reciprocal_nlj_shell(c)
  d_con = reciprocal_nlj_shell(d)
  recoupl_factor = ((-1)**(Jbra+Jket))*(2*KK + 1.0d0)*(2*Jket + 1.0d0)
  recoupl_factor = recoupl_factor * aux_1 * aux_2 * aux_3 * aux_4

!  print "(A,3I4,A,I4,A,3I3,F15.9)", "   > got to the recoup!: (J,S,L)bra=", &
!        Jbra,Sbra,Lbra," KK=",KK, " (J,S,L)ket=", Jket,Sket,Lket,recoupl_factor

  call recouple_jjLSConjugatedME(a,b,c,d, a_con,b_con,c_con,d_con, &
                                 Sbra,Sket,Lbra,Lket,Jbra,Jket,Mbra,Mket,&
                                 hamilJM, dim_jm, dim_sh, TENSOR_ORD, &
                                 recoupl_factor, auxHamilRed, KK, kval_is_zero)
  if (all_zero(KK)) all_zero(KK) = kval_is_zero ! modify just if all was zero
  !!! ================================================================

              enddo ! L ket_ loop

            enddo ! S ket_ loop
          enddo ! L bra_ loop
        enddo ! Sbra_ loop
      enddo ! J ket_ loop

    enddo ! k loop
  enddo ! Jbra loop

  !! WRITE the final result of the loop for each block ============
  do KK = 0, TENSOR_ORD
    if (all_zero(KK)) cycle
    if (KK > 0) then
      !! Example of J limits implemented (being * intermediate Js))
      !BRA Js:(Jbmin) * [J1](-KK)[       ]  * * (Jbmax=J2)
      !KET Js:                     (Jkmin=J3) * * [        ](+KK)[J4] * (Jkmax)
      J3 = max(Jb_min , Jk_min)
      J2 = min(Jb_max , Jk_max)
      J1 = max(J3 - KK, min(Jb_min, Jk_min))
      J4 = min(J2 + KK, max(Jb_max, Jk_max))
      write(300+KK, fmt='(A,4I8,4I3)') &
        '0 5', a_ant, b_ant, c_ant, d_ant, J1, J2, J3, J4
    else
      write(300, fmt='(A,4I8,2I3)') ' 0 5', a_ant, b_ant, c_ant, d_ant, &
         max(Jb_min, Jk_min), min(Jb_max, Jk_max)
      write(299, fmt='(A,4I8,2I3)') ' 0 5', a_ant, b_ant, c_ant, d_ant, &
         max(Jb_min, Jk_min), min(Jb_max, Jk_max)
    endif
  enddo

  do Jbra = Jb_min, Jb_max
    do Jket = Jk_min, Jk_max
      KKmin = min(abs(Jbra - Jket), TENSOR_ORD)
      KKmax = min(Jbra + Jket, TENSOR_ORD)
      if (KKmin > TENSOR_ORD) cycle

      ind_jm_b = angular_momentum_index(Jbra, 0, .FALSE.)
      ind_jm_k = angular_momentum_index(Jket, 0, .FALSE.)

!      if ((delta_ab > TOL).OR.(delta_cd > TOL)) then
!        phs_pnk0 = (-1)**(Jbra)
!        endif

      do KK = KKmin, KKmax
        if (all_zero(KK)) cycle
!  print "(A,2I4,A,I4)", " ", Jbra, Jket, "  KK=",KK
        !! Write line
        if (KK > 0) then
          write(300+KK,fmt='(2I4)',advance='no')Jbra, Jket
        end if
        do t = 1, 4

          aux_1 = auxHamilRed(t,KK,ind_jm_b,ind_jm_k)

          if (KK == 0) then
            !! select the import hamil 2B J-scheme
            select case(t)
              case (1)
                tt = 0
              case (4)
                tt = 5
              case default  !! t=2(pnpn) tt=1(,4) & t=3(pnnp) tt=2(,3)! 1,2==4,3
                tt = t - 1
            end select

!            print "(A,4I3)", "    Nshell::", Na, Nb, Nc, Nd
            aux_4 = aux_1
            if (implement_H2cpd_DD) then
              aux_4 = aux_4 + hamil_H2cpd_DD(tt, Jbra, a,b,c,d)
              endif
!              print "(A,2F15.9)","    In:",aux_3,hamil_H2cpd_DD(tt,Jbra,a,b,c,d)
            if (t == 2) then ! print the permutations (for the pnpn)
              aux_2 = auxHamilRed(2,KK,ind_jm_b,ind_jm_k)
              aux_3 = auxHamilRed(3,KK,ind_jm_b,ind_jm_k)
              write(299,fmt='(4F15.10)',advance='no') &
                aux_2 + hamil_H2cpd_DD(1, Jbra, a,b,c,d), &
                aux_3 + hamil_H2cpd_DD(2, Jbra, a,b,c,d), &
                aux_3 + hamil_H2cpd_DD(3, Jbra, a,b,c,d), &
                aux_2 + hamil_H2cpd_DD(4, Jbra, a,b,c,d)
            else if (t .EQ. 3) then
              cycle
            else
              write(299,fmt='(F15.10)',advance='no') aux_4
            endif
          endif !! K=0

          !! Parts for the Hamil DD only -----------------------------------
          if (t .EQ. 2) then ! print the permutations (for the pnpn)
            aux_2 = auxHamilRed(2,KK,ind_jm_b,ind_jm_k)
            aux_3 = auxHamilRed(3,KK,ind_jm_b,ind_jm_k)
            write(300+KK,fmt='(4F15.10)',advance='no') aux_2,aux_3,aux_3,aux_2
          else if (t .EQ. 3) then
            cycle
          else
            write(300+KK,fmt='(F15.10)',advance='no') aux_1
          endif
          !! Parts for the Hamil DD only -----------------------------------
        enddo ! t iter

        write(300+KK,*) ''
        if (KK==0) then
          write(299,*) ''
        endif

      enddo ! K tensor Loop
    enddo
  enddo
  !! ===============================================
         enddo ! shell loops
      enddo
   enddo
enddo

deallocate(hamilJM, auxHamilRed, all_zero)!, J1,J2,J3,J4)
if (implement_H2cpd_DD) deallocate(hamil_H2cpd_DD)
do KK = 0, TENSOR_ORD
  close(300 + KK)
end do
close(299)
close(298)

print *, " * [OK] Printing 2B Matrix elements DD from WF_HFB\n"
end subroutine print_DD_matrix_elements


!------------------------------------------------------------------------------!
! subroutine print_quasipartile_DD_matrix_elements                             !
!  Export of the DD + Hamil Hamiltonian in a reduced quasiparticle basis using !
!  formulas of appendix E.2 in Ring-Schuck                                     !
!------------------------------------------------------------------------------!
subroutine print_quasipartile_DD_matrix_elements(bogo_zU0, bogo_zV0, &
                                                 dens_rhoRR, dens_kappaRR, ndim)

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: bogo_zU0,bogo_zV0
real(r64), dimension(ndim,ndim), intent(in)    :: dens_rhoRR, dens_kappaRR

real(r64), dimension(:), allocatable :: eigen_H11    ! qp    "
complex(r64), dimension(ndim,ndim) :: hspRR, gammaRR, deltaRR
real(r64), dimension(ndim,ndim) :: field_hspRR, field_deltaRR

real(r64) :: xneut, xprot, xpar, xjz, xn, xj, xl, xj2, xl2
integer :: i,j, kk, k1,k2,k3,k4, l1,l2,l3,l4
integer :: iH11, ialloc
real(r64), dimension(3*ndim-1) :: work
real(r64) :: THRESHOLD = 0.0d0


character(len=*), parameter :: format1 = "(1i4,7f9.3,1x,2f12.6)"


allocate(eigen_H11(HOsp_dim), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of gradient'
eigen_H11 = zero

call calculate_fields_diag(zone*dens_rhoRR,zone*dens_kappaRR,gammaRR,hspRR, &
                             deltaRR,ndim=ndim)
field_hspRR = real(hspRR)
field_deltaRR = real(deltaRR)

! 1 Diagonalize H11 matrix to select the first states of the quasiparticles
call calculate_H11_real(ndim)

call dsyev('v','u',ndim,field_H11,ndim,eigen_H11,work,3*ndim-1,iH11)

kk=0
do i = 1, ndim
  if (eigen_H11(i).LT.THRESHOLD) then
    kk = kk + 1


  endif
enddo

deallocate(VStoHOsp_index)
VSsp_dim = kk
allocate(VStoHOsp_index(VSsp_dim))


! 2
print "(A)", " *** Start test of H11 in DD subroutine  ********** "
print "(A)", "   #      Z        N        n        l        p &
                     &       j       jz         H11"
do i = 1, ndim
  xneut = zero
  xprot = zero
  xpar  = zero
  xjz   = zero
  xj2   = zero
  xn    = zero
  xl2   = zero
  do j = 1, ndim
    xprot = xprot + field_H11(j,i)**2 * (-HOsp_2mt(j) + 1)/2.0d0
    xneut = xneut + field_H11(j,i)**2 * ( HOsp_2mt(j) + 1)/2.0d0
    xpar  = xpar  + field_H11(j,i)**2 * (-1.d0)**HOsp_l(j)
    xn    = xn    + field_H11(j,i)**2 * HOsp_n(j)
    xjz   = xjz   + field_H11(j,i)**2 * HOsp_2mj(j)/2.0d0
    xj2   = xj2   + field_H11(j,i)**2 * (HOsp_2j(j)*(HOsp_2j(j)+2))/4.0d0
    xl2   = xl2   + field_H11(j,i)**2 * (HOsp_l(j)*(HOsp_l(j)+1))
  enddo
  xj = 0.5d0 * (-1.d0 + sqrt(1+4*abs(xj2)))
  xl = 0.5d0 * (-1.d0 + sqrt(1+4*abs(xl2)))
  print format1, i, xprot, xneut, xn, xl, xpar, xj, xjz, eigen_H11(i)
enddo
print "(A)", " *** End test of H11 in DD subroutine  ********** "


end subroutine print_quasipartile_DD_matrix_elements




subroutine test_printDesityKappaWF(rhoLR, kappaLR, kappaRL, ndim)
integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, kappaLR, kappaRL

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
do i = 1, ndim
    do j = 1, ndim
    write(626, '(2I4,6F20.15)') i, j, dreal(rhoLR(i, j)), dimag(rhoLR(i, j)), &
        dreal(kappaLR(i, j)), dimag(kappaLR(i, j)), &
        dreal(kappaRL(i, j)), dimag(kappaRL(i, j))
    enddo
enddo
close(626)

!print *, "[OK] Export Rho LR matrix in file = 'DIMENS_indexes_and_rhoLRkappas.gut'"
end subroutine test_printDesityKappaWF



!==============================================================================!
! Subroutine to export a file for plotting the density in an integrable form.  !
! Whether the integral is trapezoidal, Legendre, Legendre-radial Laguerre.     !
!==============================================================================!
subroutine export_expectval_density(dens_rhoLR,dens_kappaLR,dens_kappaRL, &
                                    bogo_zU0,bogo_zV0, ndim)

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: dens_rhoLR
complex(r64), dimension(ndim,ndim), intent(in) :: dens_kappaLR,dens_kappaRL
complex(r64), dimension(ndim,ndim), intent(in) :: bogo_zU0,bogo_zV0

integer :: i,j, a_sh,la,ja,ma,mta, b_sh,lb,jb,mb,mtb, K
integer :: ind_jm_a, ind_jm_b, ind_km
integer :: M
integer :: i_r=1, i_th=1, i_phi=1, i_ang=1

real(r64) :: radial_part, weight, x, y, z, r_case_export
complex(r64) :: sum_, integral_dens, int_dens_Z, int_dens_N, aux_fi, &
        d_exp_rtp, dExpZ,dExpN, sum_p, kExpZ, kExpN

real(r64) :: st_1,st_2,st_3, st_4, fin_1,fin_2,fin_3,fin_4, dens_exact
real(r64) :: time_ij, time_r
character(len=8) :: tab = '        '

integer      :: mla, mlb, ms ! TEST, REMOVE
real(r64)    :: cga, cgb     ! TEST, REMOVE
complex(r64) :: Ya, Yb       ! TEST, REM

!! TEST points of the grid (with non integrable variable approximation)
complex(r64), allocatable, dimension(:,:) :: test_dens
allocate(test_dens(ndim, ndim))
test_dens = zzero



if (.NOT.export_density) return

if (exportVSPSpace) then !-----------------------------------------------------

call test_printDesityKappaWF(dens_rhoLR, dens_kappaLR,dens_kappaRL, ndim)
call calculate_densityDep_hamiltonian(dens_rhoLR,dens_kappaLR,dens_kappaRL,ndim)

if (.NOT.evalQuasiParticleVSpace) then
call print_DD_matrix_elements
else
call print_quasipartile_DD_matrix_elements(bogo_zU0, bogo_zV0, &
                    real(dens_rhoRR, real64), real(dens_kappaRR, real64), ndim)
endif

endif !------------------------------------------------------------------------

density_export    = zzero
density_export_n  = zzero
density_export_p  = zzero
pairdens_export   = zzero
pairdens_export_n = zzero
pairdens_export_p = zzero

integral_dens = zzero
int_dens_Z = zzero
int_dens_N = zzero

open( 613, file='export_density_rtp.txt') !====================================
write(613, fmt='(A,3I5,F10.6,I3)') &
                "RDim,CThDim,PhiDim,b lenght,integration method_", &
                r_dim, theta_dim, phi_dim, HO_b, integration_method
write(613, fmt='(A,A)') " i_r i_t i_p    r	       cos_th         phi", &
    "            REAL(dens)     IMAG(dens)     weight_prod"
open( 615, file='export_dens_pairing_rtp.txt') !===============================
write(615, fmt='(A,3I5,F10.6,I3)') &
                "RDim,CThDim,PhiDim,b lenght,integration method_", &
                r_dim, theta_dim, phi_dim, HO_b, integration_method
write(615, fmt='(A,A,A)') " i_r i_t i_p    r	       cos_th         phi", &
    "            REAL(kappaZ)   IMAG(kappaZ)    REAL(kappaN)    IMAG(kappaN)",&
    "     weight_prod"
open( 614, file='export_density_xyz.txt') !====================================
write(614, fmt='(A,3I5,F10.6,I3)') &
                "RDim,CThDim,PhiDim,b lenght_,integration_method", &
                r_dim, theta_dim, phi_dim, HO_b, integration_method
write(614, fmt='(A,A,A)') " i_r i_t i_p     X               Y               Z",&
"               REAL(densZ)     IMAG(densZ)     REAL(densN)     IMAG(densZ)",&
"     Weight_prod"
do i_r = 1, r_dim
  print '(I3,A,I3,A,I3,A,I3,A,I3,A,I3,A,F16.8,A,F12.8,A,F12.8)', &
          i_r,'/',r_dim,'r_', i_th,'/',theta_dim,'th_', i_phi,'/',phi_dim, &
          'ph :: (r=',r(i_r),',t=', theta(i_th), ',p=', phi(i_phi)

  time_ij    = zero
  call cpu_time(st_1)
  do i_ang = 1, angular_dim
    call angular_index(i_ang, i_th, i_phi)

    call cpu_time(st_2)
    do i = 1, ndim/2
       a_sh = HOsp_sh(i)
       la = HOsp_l(i)
       ja = HOsp_2j(i)
       ma = HOsp_2mj(i)
       mta= HOsp_2mt(i)
       ind_jm_a = angular_momentum_index(ja, ma, .TRUE.)

       do j = 1, ndim/2
         b_sh = HOsp_sh(j)
         lb = HOsp_l(j)
         jb = HOsp_2j(j)
         mb = HOsp_2mj(j)
         mtb= HOsp_2mt(j)
         ind_jm_b = angular_momentum_index(jb, mb, .TRUE.)

         if (mta /= mtb) cycle !! omit since the ndim / 2 only goes over protons

         sum_ = zzero
         !! --- PRECALCULATED SPEEDING ON Tested Method 2 -----
         do K = abs(ja - jb) / 2, (ja + jb) / 2
            M = (mb - ma)/2
            if ((MOD(K + la + lb, 2) == 1).OR.(abs(M) > K)) cycle

            ind_km = angular_momentum_index(K, M, .FALSE.)

            sum_ = sum_ + (dens_Y_KM_me(ind_jm_a, ind_jm_b, ind_km) * &
                           sph_harmonics_memo(ind_km, i_ang))
        enddo
        !! -------------- TEST EXPLICIT SUM OF SP.HARM a & b ---------------
        ! USED for the PAIRING tensor density Ya is CONJG(Ya) for density
        sum_p = zzero
        do mla = -(la), la
            do mlb = -(lb), lb
               do ms = -1, 1, 2
                   if (2*mla + ms /= ma) cycle
                   if (2*mlb + ms /= mb) cycle

                   call ClebschGordan(1, 2*la, ja, ms, 2*mla, ma, cga)
                   call ClebschGordan(1, 2*lb, jb, ms, 2*mlb, mb, cgb)

                   ind_km = angular_momentum_index(la, mla, .FALSE.)
                   Ya = sph_harmonics_memo(ind_km, i_ang)
                   ind_km = angular_momentum_index(lb, mlb, .FALSE.)
                   Yb = sph_harmonics_memo(ind_km, i_ang)

                   sum_p = sum_p + (cga * cgb * Ya * Yb)
               enddo
            enddo
        enddo

        radial_part = radial_2b_sho_export_memo(a_sh, b_sh, i_r)

        sum_ = sum_ * radial_part
        !! TEST
        test_dens(i,j) = test_dens(i,j) + (sum_*weight_LEB(i_ang)*weight_R(i_r))
        !! END TEST

        density_export_n(i_r, i_ang) = density_export_n(i_r, i_ang) + &
                                        (sum_ * dens_rhoLR(j+ndim/2, i+ndim/2))
        density_export_p(i_r, i_ang) = density_export_p(i_r, i_ang) + &
                                        (sum_ * dens_rhoLR(j, i))
        sum_p = sum_p * radial_part
        pairdens_export_n(i_r, i_ang) = pairdens_export_n(i_r, i_ang) + &
                                    (sum_p * dens_kappaLR(i+ndim/2, j+ndim/2))
        pairdens_export_p(i_r, i_ang) = pairdens_export_p(i_r, i_ang) + &
                                    (sum_p * dens_kappaLR(i, j))

      enddo ! do i
    enddo   ! do j

    density_export(i_r, i_ang)  = density_export_p(i_r, i_ang) + &
                                    density_export_n(i_r, i_ang)
    call cpu_time(fin_2)
    time_ij = time_ij + (fin_2 - st_2)

    weight = weight_R(i_r) * weight_THE(i_th) * weight_PHI(i_phi)
    select case(integration_method)
      case(0)
        aux_fi = weight * ((r_export(i_r))**2) * sin(theta(i_th))
        r_case_export = r_export(i_r)
        d_exp_rtp = density_export(i_r, i_ang)
        dExpZ =  density_export_p (i_r, i_ang)
        dExpN =  density_export_p (i_r, i_ang)
        kExpZ =  pairdens_export_p(i_r, i_ang)
        kExpN =  pairdens_export_n(i_r, i_ang)
      case(1)
        aux_fi = weight * ((r_export(i_r))**2)
        r_case_export = r_export(i_r)
        d_exp_rtp = density_export(i_r, i_ang)
        dExpZ =  density_export_p (i_r, i_ang)
        dExpN =  density_export_p (i_r, i_ang)
        kExpZ =  pairdens_export_p(i_r, i_ang)
        kExpN =  pairdens_export_n(i_r, i_ang)
      case(2)
        aux_fi = weight * (HO_b**3) / 2.0
        r_case_export = r_export(i_r) / HO_b
        d_exp_rtp = density_export(i_r, i_ang) / exp((r_case_export)**2)
        dExpZ =  density_export_p (i_r, i_ang) / exp((r_case_export)**2)
        dExpN =  density_export_p (i_r, i_ang) / exp((r_case_export)**2)
        kExpZ =  pairdens_export_p(i_r, i_ang) / exp((r_case_export)**2)
        kExpN =  pairdens_export_n(i_r, i_ang) / exp((r_case_export)**2)
      case(3)
        aux_fi = weight * (HO_b**3) * 4 * pi / 2.0
        r_case_export = r_export(i_r) * HO_b
        d_exp_rtp = density_export(i_r, i_ang) / exp((r_case_export)**2)
        dExpZ =  density_export_p (i_r, i_ang) / exp((r_case_export)**2)
        dExpN =  density_export_p (i_r, i_ang) / exp((r_case_export)**2)
        kExpZ =  pairdens_export_p(i_r, i_ang) / exp((r_case_export)**2)
        kExpN =  pairdens_export_n(i_r, i_ang) / exp((r_case_export)**2)
    end select ! the resultant integral must be A

    write(613, fmt='(3I4,6D20.9)') &
        i_r,i_th,i_phi, r_export(i_r), cos_th_export(i_th),phi_export(i_phi), &
        real(d_exp_rtp), max(imag(d_exp_rtp),1.0d-98), weight
    write(615, fmt='(3I4,8D20.9)') &
        i_r,i_th,i_phi, r_export(i_r), cos_th_export(i_th),phi_export(i_phi), &
        real(kExpZ),max(imag(kExpZ),1.0d-98), &
        real(kExpN),max(imag(kExpN),1.0d-98), weight

    integral_dens = integral_dens + (density_export(i_r, i_ang) * aux_fi)
    int_dens_Z = int_dens_Z + (density_export_p(i_r, i_ang) * aux_fi)
    int_dens_N = int_dens_N + (density_export_n(i_r, i_ang) * aux_fi)

    ! Export the density in cartesians (Tomas script): X, Y, Z, rhoZ, rhoN
    x = r_case_export * sin(theta_export(i_th)) * cos(phi_export(i_phi))
    y = r_case_export * sin(theta_export(i_th)) * sin(phi_export(i_phi))
    z = r_case_export * cos(theta_export(i_th))

    write(614, fmt='(3I4,8D20.9)') i_r,i_th,i_phi, x,y,z,&
         real(dExpZ), max(imag(dExpZ),1.0d-98), &
         real(dExpN), max(imag(dExpN),1.0d-98), weight

    ! -----------------------------------------------------------------------
  enddo   !end do i_angular


  call cpu_time(fin_1)
  time_r = fin_1 - st_1

  print '("    > Time [ij block] = ",f11.8," s.")',time_ij
  print '("    > Time [r step_]  = ",f11.8," s.")',time_r
enddo ! do i_r

close(613)
close(614)
close(615)

!! TEST ELEMENTS DELTA INTgRAL
if (PRINT_GUTS) then
open(557, file='test_density_delta_me.gut')
do i=1, ndim
  do j=1, ndim

    test_dens(i,j) = test_dens(i,j) * (HO_b**3) * 4 * pi / 2.0

  write(557,fmt='(2I4,2F24.16)')i,j,dreal(test_dens(i,j)),dimag(test_dens(i,j))
  end do
end do
close(557)
endif

print '(A,F18.15,A,F18.15)', "Integral density  dr^3 =", real(integral_dens), &
                                                 ' +j ', imag(integral_dens)
print '(A,F18.15)', "Int dens Protons  dr^3 =", real(int_dens_Z)
print '(A,F18.15)', "Int dens Neutrons dr^3 =", real(int_dens_N)
print *, "Integration method = ", integration_method

!export_density = .FALSE. ! Do not export more
print *, '[OK] EXPORT expected value of density: files = '
print *, "     'density_rtp.txt', 'density_xyz.txt', 'dens_pairing_rtp.txt'"

end subroutine export_expectval_density


END MODULE DensityDep
!==============================================================================!
! End of file                                                                  !
!==============================================================================!

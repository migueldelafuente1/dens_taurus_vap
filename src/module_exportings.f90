!==============================================================================!
! MODULE DensDepResultExportings                                               !
!                                                                              !
!    This module join the subroutines to evaluate the density dependent        !
! interaction and its hamiltonian and also the diagonalizing basis for the     !
! single particle hamiltonian with the j_z observable.                         !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine  export_densityAndHamiltonian                                   !
! - subroutine  calculate_valenceSpaceReduced                                  !
! - subroutine  recouple_jjLSConjugatedME                                      !
! - subroutine  print_DD_matrix_elements                                       !
! - subroutine  print_quasipartile_DD_matrix_elements                          !
! - subroutine  export_expectval_density                                       !
!                                                                              !
!==============================================================================!
MODULE DensDepResultExportings

use Constants
use MathMethods
use Basis
use Wavefunctions
use Hamiltonian
use DensityDep
use Gradient

implicit none
PUBLIC

!! Attributes

integer, private :: info_H11       ! check if problem during diag(field_H11)

real(r64), dimension(:), allocatable, private :: eigen_hsp, & ! sp energies
                                                 eigen_H11    ! qp    "
real(r64), dimension(:,:), allocatable, private :: transf_H11

integer, dimension(:), allocatable, private :: QPtoHOsp_index
!! Methods

CONTAINS

!------------------------------------------------------------------------------!
! subroutine  export_densityAndHamiltonian                                     !
!    Select the options to export after a calculation, the spatial density and !
! kappa correlation function. Export of matrix elements from the DD term, for  !
! a reduced valence space in a shell model fashion or in the Quasiparticle     !
! spherical j_z-diagonal basis for a HFB calculation                           !
!------------------------------------------------------------------------------!
subroutine export_densityAndHamiltonian(dens_rhoRR, dens_kappaRR, ndim)

integer, intent(in) :: ndim
real(r64),    dimension(ndim,ndim), intent(in) :: dens_rhoRR, dens_kappaRR
complex(r64), dimension(ndim,ndim) :: dens_rhoRRc, dens_kappaRRc
integer  :: i, j

if (.NOT.eval_density_dependent) return
print *, ""
print "(A)", " [  SR] Exporting Hamiltonian and/or spatial density on the&
             & integration grid."

if (export_density) then
  call export_expectval_density(dens_rhoRR, dens_kappaRR, dens_kappaRR, ndim)
endif

!! deallocate HF arrays from D1S to increase memory
deallocate(sphharmDUAL_memo, AngFunctDUAL_HF, AngFunctDUAL_P1, &
           AngFunctDUAL_P2, BulkHF, BulkP1, BulkP2)

if (exportValSpace) then !-----------------------------------------------------

  !! arguments for test printDensKappa must be Complex(8)
  do i = 1, ndim
    do j = 1, ndim
      dens_rhoRRc  (i,j) = dCMPLX(dens_rhoRR  (i,j), 0.0d0)
      dens_kappaRRc(i,j) = dCMPLX(dens_kappaRR(i,j), 0.0d0)
    end do
  end do

  call test_printDesityKappaWF(dens_rhoRRc, dens_kappaRRc, dens_kappaRRc, ndim)
  print "(A)", "  [  SR] Evaluating the Hamiltonian."
  if (evalQuasiParticleVSpace) print "(A)","    For full space, Be patient ..."
  call calculate_densityDep_hamiltonian(dens_rhoRRc, &
                                        dens_kappaRRc, dens_kappaRRc, ndim)
  print "(A)", "  [DONE] Evaluating the Hamiltonian."
  deallocate(rearrangement_me, rearrang_field, &
             rea_common_RadAng, REACommonFields)
  if (.NOT.evalQuasiParticleVSpace) then
    call print_DD_matrix_elements
  else
!    if (is_good_K) then
!      print "(a)", "[WARNING] Wavefunction is not axial, cannot extract &
!                    & spherical quasiparticle matrix elements. Exit."
!      return
!    endif
    call print_quasipartile_DD_matrix_elements(dens_rhoRR,dens_kappaRR,HOsp_dim)
  endif
endif

end subroutine export_densityAndHamiltonian


!------------------------------------------------------------------------------!
! subroutine  calculate_valenceSpaceReduced                                    !
!   Evaluates the single particle energies and Core energy from a HFB energy   !
! done in a shell model fashion (summing all core levels as fully occupied).   !
!------------------------------------------------------------------------------!
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
! subroutine diagonalize_H11_and_jz                                            !
!  Obtains the simultaneous diagonalization of H11 and the Jz operator,        !
!  following the same process from the module_gradient for printing eigenbasis !
!------------------------------------------------------------------------------!
subroutine diagonalize_H11_with_jz(dens_rhoRR, dens_kappaRR, ndim)

integer, intent(in) :: ndim
!complex(r64), dimension(ndim,ndim), intent(in) :: bogo_zU0,bogo_zV0
real(r64), dimension(ndim,ndim), intent(in)    :: dens_rhoRR, dens_kappaRR

integer :: i, j, k, l, m, nocc0, nemp0, evnum, ialloc=0
real(r64) :: ovac0
integer, dimension(1) :: tabmin
integer, dimension(ndim) :: eigenh_order, evdeg
real(r64), dimension(ndim) :: eigenh_tmp
real(r64), dimension(3*ndim-1) :: work
real(r64), dimension(ndim,ndim) :: D0, rhoc, hspc, A1, A2

real(r64), dimension(:,:), allocatable :: hspr
real(r64), dimension(:), allocatable :: workr, eigenr
complex(r64), dimension(ndim,ndim) :: hspRR, gammaRR, deltaRR

allocate(eigen_hsp(HOsp_dim), eigen_H11(HOsp_dim), &
         transf_H11(HOsp_dim, HOsp_dim), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of gradient'

eigen_hsp  = zero
eigen_H11  = zero
transf_H11 = zero

!!! Computes the fields
call calculate_fields_diag(zone*dens_rhoRR,zone*dens_kappaRR,gammaRR,hspRR, &
                           deltaRR,ndim=ndim)
field_hspRR   = real(hspRR)
field_deltaRR = real(deltaRR)

!cmpi if ( paral_myrank /= 0 ) return

call calculate_H11_real(ndim)

!!! hsp in canonical basis
call construct_canonical_basis(bogo_U0,bogo_V0,bogo_zU0c,bogo_zV0c,bogo_zD0, &
                               ovac0,nocc0,nemp0,ndim)
D0 = real(bogo_zD0)

call dgemm('t','n',ndim,ndim,ndim,one,D0,ndim,dens_rhoRR,ndim,zero,A1,ndim)
call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,D0,ndim,zero,rhoc,ndim)

call dgemm('t','n',ndim,ndim,ndim,one,D0,ndim,field_H11,ndim,zero,A1,ndim)
call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,D0,ndim,zero,hspc,ndim)

!!! Further reduces h in case of fully empty/occupides states
if ( nemp0 > 0 ) then
  allocate (hspr(nemp0,nemp0), eigenr(nemp0),workr(3*nemp0-1))
  hspr(1:nemp0,1:nemp0) = hspc(1:nemp0,1:nemp0)
  call dsyev('v','u',nemp0,hspr,nemp0,eigenr,workr,3*nemp0-1,info_H11)
  A1 = zero
  A2 = D0
  do i = 1, ndim
    A1(i,i) = one
  enddo
  A1(1:nemp0,1:nemp0) = hspr(1:nemp0,1:nemp0)
  call dgemm('n','n',ndim,ndim,ndim,one,A2,ndim,A1,ndim,zero,D0,ndim)
  deallocate(hspr, eigenr, workr)
endif

if ( nocc0 > 0 ) then
  allocate (hspr(nocc0,nocc0), eigenr(nocc0),workr(3*nocc0-1))
  hspr(1:nocc0,1:nocc0) = hspc(ndim-nocc0+1:ndim,ndim-nocc0+1:ndim)
  call dsyev('v','u',nocc0,hspr,nocc0,eigenr,workr,3*nocc0-1,info_H11)
  A1 = zero
  A2 = D0
  do i = 1, ndim
    A1(i,i) = one
  enddo
  A1(ndim-nocc0+1:ndim,ndim-nocc0+1:ndim) = hspr(1:nocc0,1:nocc0)
  call dgemm('n','n',ndim,ndim,ndim,one,A2,ndim,A1,ndim,zero,D0,ndim)
  deallocate(hspr, eigenr, workr)
endif

call dgemm('t','n',ndim,ndim,ndim,one,D0,ndim,field_H11,ndim,zero,A1,ndim)
call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,D0,ndim,zero,hspc,ndim)

!!! Ordering of energies
l = 0
eigenh_order = 0
eigenh_tmp = 999

do i = 1, ndim
  if ( abs(rhoc(i,i)) > 1.d-7 ) then
    l = l + 1
    eigenh_tmp(i) = hspc(i,i)
  endif
enddo

do i = 1, l
  tabmin = minloc(eigenh_tmp)
  eigenh_order(i) = tabmin(1)
  eigenh_tmp(tabmin(1)) = 1000
enddo

eigenh_tmp = 999

do i = 1, ndim
  if ( abs(rhoc(i,i)) <= 1.d-7 ) then
    eigenh_tmp(i) = hspc(i,i)
  endif
enddo

do i = l+1, ndim
  tabmin = minloc(eigenh_tmp)
  eigenh_order(i) = tabmin(1)
  eigenh_tmp(tabmin(1)) = 1000
enddo

!!! Diagonalizes hsp
call dsyev('v','u',ndim,field_H11,ndim,eigen_H11,work,3*ndim-1,info_H11)

! In the case of axial symmetry, further diagonalizes Jz in this basis
if (is_good_K) then

  ! counts the number of eigenspaces and their degeneracy
  evnum = 1
  evdeg = 0
  evdeg(1) = 1

  do i = 2, ndim
    if ( abs(eigen_H11(i-1) - eigen_H11(i)) < 1.0d-6 ) then
      evdeg(evnum) = evdeg(evnum) + 1
    else
      evnum = evnum + 1
      evdeg(evnum) = evdeg(evnum) + 1
    endif
  enddo

  ! Jz in the matrix that diagonalizes h
  D0 = field_H11
  call dgemm('t','n',ndim,ndim,ndim,one,D0,ndim,angumome_Jz(1:ndim**2),ndim, &
             zero,A1,ndim)
  call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,D0,ndim,zero,A2,ndim)

  ! block diagonalization
  j = 0
  A1 = zero

  do i = 1, evnum
    k = evdeg(i)
    allocate( hspr(k,k), eigenr(k), workr(3*k-1) )
    hspr(:,:) = A2(1+j:j+k,1+j:j+k)
    call dsyev('v','u',k,hspr,k,eigenr,workr,3*k-1,info_H11)
    A1(1+j:j+k,1+j:j+k) = hspr(1:k,1:k)
    j = j + k
    deallocate( hspr, eigenr, workr )
  enddo

  call dgemm('n','n',ndim,ndim,ndim,one,D0,ndim,A1,ndim,zero,field_H11,ndim)
endif

call sort_quasiparticle_basis(ndim)


end subroutine

!------------------------------------------------------------------------------!
! subroutine sort_quasiparticle_basis
!
!------------------------------------------------------------------------------!
subroutine sort_quasiparticle_basis(ndim)

integer, intent(in) :: ndim
real(r64) :: xn, xl2, xl, xneut, xprot, xpar, xjz, xj2, xj, fermi_p, fermi_n
real(r64), dimension(ndim)  :: qpsp_zz, qpsp_nn, qpsp_par, qpsp_jz, qpsp_n, &
                               qpsp_j, qpsp_l
integer, dimension(3) :: possible_qp_for_hosp, possible_n_for_qp ! only to export up to 3 shells
logical, dimension(:), allocatable :: QP_index_found
integer, dimension(:), allocatable :: VStoQPsp_index, &
                                      HOshells, VSshells, sortedShells
character(len=*), parameter :: format1 = "(1i4,7f9.3,1x,2f12.6)", &
                               format2 = "(1a77,/,80('-'))"
logical :: found
integer :: i,j,k,k1,k2, kk, items_found=0, items_found_2=0, alloc_it
integer :: sp_n,sp_l,sp_2j,sp_2mj,sp_2mt, sp_sh, nmaj_sh, Nsh, VSNdim, HONdim,&
           VSlim=0, HOlim=0, index_Nsh, Nsh_i, n

!!!
integer :: METHOD_SORT = 1  ! 0 vs assumed to be the first n, 1 to sort with the HO-N shell order
!!!

allocate(QP_index_found(ndim), QPtoHOsp_index(ndim))
print "(A)", " Sorting the "
! sort the shells and check if the METHOD_SORT is 1 or 0
select CASE(METHOD_SORT)
  CASE(0)
    print "(A)", " ** METHOD_SORT=0 to select the n in the quasiparticles by &
                 & the order in the HO shells, the <n>_qp order preserves"
  CASE(1)
    print "(A)", " ** METHOD_SORT=1, the shells will be ordered assuming the &
                 & valence space states are the first to be assigned."

    allocate(HOshells(HOsh_dim), VSshells(VSsh_dim), sortedShells(HOsh_dim))
    VSlim = VSsh_dim
    HOlim = HOsh_dim
    do alloc_it = 1,2
      if (alloc_it.EQ.2) then
        deallocate(HOshells, VSshells, sortedShells)
        allocate(HOshells(items_found_2), VSshells(items_found),&
                 sortedShells(items_found_2))
        VSlim = items_found
        HOlim = items_found_2
      endif
      sortedShells = -1
      VSshells     = -1
      HOshells     = -1

      ! Write the valence space oscillator N shells
      items_found = 0
      do i = 1, VSsh_dim
        kk = VStoHOsh_index(i)
        Nsh = 2*HOsh_n(kk) + HOsh_l(kk)
        found = .FALSE.
        do j = 1, VSlim
          if (Nsh .EQ. VSshells(j)) found = .TRUE.
        enddo
        if (.NOT.found) then
          items_found = items_found + 1
          VSshells(items_found) = Nsh
        endif
      enddo

      ! Write all shells
      items_found_2=0
      do kk = 1, HOsh_dim
        Nsh = 2*HOsh_n(kk) + HOsh_l(kk)
        found = .FALSE.
        do j = 1, HOlim
          if (Nsh .EQ. HOshells(j)) found = .TRUE.
        enddo
        if (.NOT.found) then
          items_found_2 = items_found_2 + 1
          HOshells(items_found_2) = Nsh
        endif
      enddo

      !Construct the shell order to print
      do i = 1, VSlim
        sortedShells(i) = VSshells(i)
      enddo
      kk = VSlim + 1
      do i=1, HOlim
        found = .FALSE.
        do j = 1, VSlim
          if (HOshells(i) .EQ. VSshells(j)) found =.TRUE.
        enddo
        if (found) cycle
        sortedShells(kk) = HOshells(i)
        kk = kk + 1
      enddo
    enddo !alloc_iter
    !! TEST PRINT
    do i = 1, VSlim
      print "(A,2i5)", " test VSshells i,N=", i, VSshells(i)
    enddo
    print *, ""
    do i = 1, HOlim
      print "(A,3i5)", " test HOshells i,N=", i, HOshells(i), sortedShells(i)
    enddo
    print *, ""

  CASE DEFAULT
    print "(A)", " [ERROR] Invalid METHOD_SORT in sort_quasiparticle_basis &
                 & subroutine, default METHOD_SORT=0 selected:"
    print "(A)", " ** METHOD_SORT=0 to select the n in the quasiparticles by &
                 & the order in the HO shells, the <n>_qp order preserves"
    METHOD_SORT = 0
end select

!!! Writes the properties of the single-particle states in a file
fermi_p = 0.d0
fermi_n = 0.d0

if ( constraint_switch(1) == 1 ) then
  fermi_p = lagrange_lambda1(1)
endif
if ( constraint_switch(2) == 1 ) then
  fermi_n = lagrange_lambda1(1 + constraint_switch(1))
endif

!!! Basis that diagonalizes h and search the QP basis
open(ute, file='eigenbasis_jzH11.dat', status='replace', action='write', &
         form='formatted')
write(ute,"(1a,1f12.6)")   "Proton  fermi energy = ",fermi_p
write(ute,"(1a,1f12.6,/)") "Neutron fermi energy = ",fermi_n
write(ute,format2) "   #      Z        N        n        l        p &
                   &       j       jz         h  "
do i = 1, ndim
  xneut = zero
  xprot = zero
  xpar  = zero
  xjz   = zero
  xj2   = zero
  xn    = zero
  xl2   = zero
  do j = 1, ndim
    transf_H11(j,i) = field_H11(j,i)

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

  qpsp_zz(i) = xprot
  qpsp_nn(i) = xneut
  qpsp_n(i) = xn
  qpsp_l(i) = xl
  qpsp_j(i) = xj
  qpsp_jz(i)  = xjz
  qpsp_par(i) = xpar

  QP_index_found(i) = .FALSE.

  write(ute,format1) i, xprot, xneut, xn, xl, xpar, xj, xjz, eigen_H11(i)
enddo
close(ute, status='keep')

!!! Search and assign the QP basis
print "(A)", " *** Print the sp states of the VS index and WB state"
do i = 1, VSsp_dim
  kk = VStoHOsp_index(i)
  print "(A,2i3,i6,A,4i3)", " HOsp(vs):", i, kk, HOsh_ant(HOsp_sh(kk)), &
              " (nljm) :: ", HOsp_n(kk), HOsp_l(kk), HOsp_2j(kk), HOsp_2mj(kk)
enddo

!! locate the order of the shells of the VS to export.
!!    n is not conserved, the N shell of the state sh_(vs) will be the Nth.




print *, ""
print "(A)", " *** Reading the full HO space to assign the QP sp states. "
!! HO element loop
do i = 1, ndim
  possible_qp_for_hosp = -1
  possible_n_for_qp    = 99
  items_found = 0

  !! QP loop to find it
  do j = 1, ndim
    ! is proton or neutron state, is l, j, jz
    if ((abs(qpsp_zz(j) - 1.0) .LT. 1.0d-6) .AND. (sp_2mt .EQ. 1)) cycle
    if ((abs(qpsp_nn(j) - 1.0) .LT. 1.0d-6) .AND. (sp_2mt .EQ.-1)) cycle
    if ( abs(qpsp_l(j) - sp_l) .GT. 1.0d-6) cycle
    if ( abs(2*qpsp_j(j) - sp_2j) .GT. 1.0d-6) cycle
    if ( abs(2*qpsp_jz(j) - sp_2mj) .GT. 1.0d-6) cycle

    if (QP_index_found(i)) cycle

    items_found = items_found + 1
    possible_qp_for_hosp(items_found) = j
    possible_n_for_qp(items_found) = qpsp_n(j)
  enddo
!  print "(A,2i3,i6,A,4i3)", " HOsp(vs):", i, kk, HOsh_ant(HOsp_sh(kk)), &
!              " (nljm) :: ", sp_n, sp_l, sp_2j, sp_2mj, sp_2mt
  print "(A,3i7)", "> i,Hoant,items_found=", i, HOsp_ant(i), items_found

  if (items_found.EQ.0) then
    print "(A,i3)", "    [ERROR] Index not found::", i
  elseif (items_found.EQ.1) then
    QPtoHOsp_index(possible_qp_for_hosp(1)) = i
    QP_index_found(possible_qp_for_hosp(1)) = .TRUE.
  else
    do k = 1, items_found
      print "(2(A,i4))", "  * posible_qp_for i=",i," :",possible_qp_for_hosp(k)
    end do
    print *, ""
    select case(METHOD_SORT)
      case (0) !--------------------------------------------------------------!
        ! Method to sort by the <n> of qp
        Nsh_i = 2* HOsp_n(i) + HOsp_l(i)
        ! order the shells with to n (that must conserve the order of Nshell)
        do k1 = 1, items_found
          do k2 = k1, items_found
            if (possible_n_for_qp(k1) > possible_n_for_qp(k2)) then
              xn = possible_n_for_qp(k2)
              possible_n_for_qp(k2) = possible_n_for_qp(k1)
              possible_n_for_qp(k1) = xn

              kk = possible_qp_for_hosp(k2)
              possible_qp_for_hosp(k2) = possible_qp_for_hosp(k1)
              possible_qp_for_hosp(k1) = kk
            endif
          end do
        end do

        ! find the place of the current Nshell of the same parity
        ! (same length as possible_qp_for_hosp)
        index_Nsh = 1
        do Nsh = 1, HOlim
          if (MOD(HOsp_l(i) + Nsh, 2) .NE. 0) cycle

          if (Nsh_i .NE. Nsh) then
            index_Nsh = index_Nsh + 1
          else
            EXIT
          endif
        end do

        QPtoHOsp_index(possible_qp_for_hosp(index_Nsh)) = i
        QP_index_found(possible_qp_for_hosp(index_Nsh)) = .TRUE.

      case (1) !--------------------------------------------------------------!
        ! 1. Find the correct value by assuming the VS state is the first in
        ! for the different quasi particles, this is intuitive if the nuclei is
        ! in a closed shell, then the lower excitation energy will be in the ones
        ! we usually consider for the VS, the rest of states will be in the core
        ! or outer, so we don't care if they are right or wrong.
        kk = 1
        do k = 1, HOlim
          Nsh = sortedShells(k)
          if (MOD(HOsh_l(i)+Nsh, 2).NE.0) cycle

          n = (Nsh - HOsp_l(i)) / 2
          print "(A,2i3)", "  Nshell, n=", Nsh, n
          if ((n - hosp_n(i)) .NE. 0) then
            kk = kk + 1
            cycle
          endif

          QPtoHOsp_index(possible_qp_for_hosp(kk)) = i
          QP_index_found(possible_qp_for_hosp(kk)) = .TRUE.
          EXIT
        enddo
    end select
  endif

enddo

! TEST
print *, ""
print "(A)", " * Results for the QP states sorted."
print "(A)","   #      Z        N        n        l        p        j       jz&
            &         h      :: qp_assigned   2mt"
do i = 1, HOsp_dim
  xprot = qpsp_zz(i)
  xneut = qpsp_nn(i)
  xn    = qpsp_n(i)
  xl    = qpsp_l(i)
  xpar  = qpsp_par(i)
  xj    = qpsp_j(i)
  xjz   = qpsp_jz(i)

  print "(A)", "a"
  kk = QPtoHOsp_index(i)
  print "(A,i9)", "kk=", kk
  print "(A,i4)", "b", HOsp_ant(kk)
  print "(A,i4)", "c", HOsp_2mt(kk)
  fermi_p = eigen_H11(i)
  print "(A)", "d"

!  print "(i5,8f6.3,A,2i7,i3)", i,xprot,xneut,xn,xl,xpar,xj,xjz,eigen_H11(i),&
!        " qp ::", kk, HOsp_ant(kk), HOsp_2mt(kk)
  print "(i5,7f6.2,f9.3)", i, xprot,xneut,xn,xl,xpar,xj,xjz, eigen_H11(i)
  print "(A,3i7)", " qp ::", kk, HOsp_ant(kk), HOsp_2mt(kk)
enddo
print "(A)", " [DONE] Results for the QP states sorted."
!! Read the indexes of the QP just to have the Valence Space


end subroutine sort_quasiparticle_basis


!------------------------------------------------------------------------------!
! subroutine print_quasipartile_DD_matrix_elements                             !
!  Export of the DD + Hamil Hamiltonian in a reduced quasiparticle basis using !
!  formulas of appendix E.2 in Ring-Schuck, modified to apply in the h_sp      !
!  diagonal basis.                                                             !
!------------------------------------------------------------------------------!
subroutine print_quasipartile_DD_matrix_elements(dens_rhoRR, dens_kappaRR, ndim)

integer, intent(in) :: ndim
!complex(r64), dimension(ndim,ndim), intent(in) :: bogo_zU0,bogo_zV0
real(r64), dimension(ndim,ndim), intent(in)    :: dens_rhoRR, dens_kappaRR

!real(r64), dimension(:), allocatable :: eigen_H11    ! qp    "
!complex(r64), dimension(ndim,ndim) :: hspRR, gammaRR, deltaRR
real(r64), dimension(ndim,ndim) :: hspRR_eigenvect, U_trans, V_trans
real(r64), dimension(ndim,ndim) :: op_xjz, op_xj2, op_l2, op_n
real(r64), dimension(ndim) :: eigen_A

real(r64) :: xneut, xprot, xpar, xjz, xn, xj, xl, xj2, xl2
integer :: i,j,ll,kk, k1,k2,k3,k4, l1,l2,l3,l4, d1,d2,d3,d4
integer :: iH11, ialloc, info_hsp
real(r64), dimension(3*ndim-1) :: work

character(len=*), parameter :: format1 = "(1i4,7f9.3,1x,2f12.6)"

print "(A)", "  1[  ] print_quasipartile_DD_matrix_elements"

call diagonalize_H11_with_jz(dens_rhoRR, dens_kappaRR, ndim)

print "(A)", "  1[OK] H11 is in diagonal with Jz."
print "(A)", "  2[  ] Sorting the QP states for the Valence sp. identification"
call sort_quasiparticle_basis(ndim)
print "(A)", "  2[OK] Sorting the QP states for the Valence sp. identification"

end subroutine print_quasipartile_DD_matrix_elements


!==============================================================================!
! Subroutine to export a file for plotting the density in an integrable form.  !
! Whether the integral is trapezoidal, Legendre, Legendre-radial Laguerre.     !
!==============================================================================!
subroutine export_expectval_density(dens_rhoLR,dens_kappaLR,dens_kappaRL,ndim)

integer, intent(in) :: ndim
real(r64), dimension(ndim,ndim), intent(in) :: dens_rhoLR
real(r64), dimension(ndim,ndim), intent(in) :: dens_kappaLR,dens_kappaRL
!complex(r64), dimension(ndim,ndim), intent(in) :: bogo_zU0,bogo_zV0

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




END MODULE DensDepResultExportings

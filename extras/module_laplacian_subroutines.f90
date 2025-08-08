!==============================================================================!
! MODULE LaplacianDensity                                                      !
!                                                                              !
!   This module stores the functions for pseudorearrangement based on the      !
!   laplacian operator, as in Dfayance interaction, does not include its       !
!   rearrangement fields, only hamiltonian terms.                              !
!                                                                              !
!==============================================================================!

MODULE LaplacianDensity

use Constants
use MathMethods
use Basis
use DensityDep
use Lebedev

implicit none
PUBLIC

real(r64), dimension(:,:,:,:), allocatable  :: radial_1b_diff_memo ! (ish, i_n[-1:1], j_l[-1:1], ir)
complex(r64), dimension(:,:,:), allocatable :: partial_dens        ! (-1,0,1,2:total ,ir,iang)


CONTAINS


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
function matrix_element_pseudoRearrangement_v3(a,b, c,d) result (v_dd_val_Real)

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
end function matrix_element_pseudoRearrangement_v3

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

integer      :: ms, ms2, tt
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

const_1 = alpha_DD / 2.0D+00

b_sh = HOsp_sh(b)
d_sh = HOsp_sh(d)

do i_r = 1, r_dim
  radial = weight_R(i_r) * exp((alpha_DD + 2.0d0) * (r(i_r) / HO_b)**2)
  do i_a = 1, angular_dim

    !! term for the hf - direct
    aux1(1) = dens_pnt(5, i_r, i_a) - x0_DD_FACTOR * dens_pnt(1, i_r, i_a)
    aux1(4) = dens_pnt(5, i_r, i_a) - x0_DD_FACTOR * dens_pnt(2, i_r, i_a)
    aux1(2) = (aux1(1) + aux1(4)) / 2.0d+00 !! it is different for pnpn npnp
    aux1(3) = - x0_DD_FACTOR * dens_pnt(3, i_r, i_a)

    aux4 = AngFunctDUAL_HF(1, b,d, i_a) + AngFunctDUAL_HF(4, b,d, i_a)
    aux4 = aux4 * radial_2b_sho_memo(b_sh, d_sh, i_r)
    do tt = 1, 4
      aux1(tt) = aux1(tt) * aux4
    end do

    aux2 = zzero
    aux3 = zzero
    do ms = 1, 4
      select case (ms)
        case (2, 3)
          ms2 = 3*(3 - ms) + 2*(ms - 2)
        case default
          ms2 = ms
      end select
      aux4 = AngFunctDUAL_HF(ms2,b,d, i_a) * radial_2b_sho_memo(b_sh,d_sh, i_r)
      aux3(1) = BulkHF(1, ms2, i_r, i_a) - x0_DD_FACTOR * BulkHF(5, ms,i_r,i_a)
      aux3(4) = BulkHF(2, ms2, i_r, i_a) - x0_DD_FACTOR * BulkHF(5, ms,i_r,i_a)
      aux3(2) = (aux3(1) + aux3(2)) / 2.0d0
      aux3(3) = BulkHF(3, ms2,i_r,i_a)

      do tt = 1, 4
        aux2(tt) = aux2(tt) + aux4 * aux3(tt)
      enddo
    enddo

    aux3 = zzero

    !! common term.
    aux4 = rea_common_RadAng(a, c, i_r, i_a) * dens_alpm1(i_r, i_a) * const_1
    angular = weight_LEB(i_a) * aux4

    !v_nnnn = v_pppp
    v_dd_value(1) = v_dd_value(1) + (radial * angular * (aux1(1) - aux2(1)))
    v_dd_value(4) = v_dd_value(4) + (radial * angular * (aux1(4) - aux2(4)))
    ! pn pn
    v_dd_value(2) = v_dd_value(2) + (radial * angular * (aux1(2) - aux2(2)))
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







END MODULE LaplacianDensity
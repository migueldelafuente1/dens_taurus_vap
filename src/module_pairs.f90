!==============================================================================!
! MODULE Pairs                                                                 !
!                                                                              !
! This module contains the variables and routines related to pair operators.   !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_pairs                                                       !
!==============================================================================!
MODULE Pairs

use Constants
use MathMethods
use Basis

implicit none
integer :: pairs_scheme
real(r64), dimension(:), allocatable :: pairs_T00_J1p1, & !T=0 MT= 0, J=1 MJ=+1
                                        pairs_T00_J10,  & !T=0 MT= 0, J=1 MJ= 0
                                        pairs_T00_J1m1, & !T=0 MT= 0, J=1 MJ=-1
                                        pairs_T1p1_J00, & !T=1 MT=+1, J=0 MJ= 0
                                        pairs_T10_J00,  & !T=1 MT= 0, J=0 MJ= 0
                                        pairs_T1m1_J00    !T=1 MT=-1, J=0 MJ= 0

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_pairs                                                         !
!                                                                              !
! Defines the HO matrix elements of the pair operators.                        !
! The operators are vectorized (for the calculation of constraints).           !
!                                                                              !
! pairs_scheme = 0 Agnostic (no special factor or selection rule)              !
!               = 1 Seniority (only pairs among particles of the same shell)   !
!                 see Eq(19.1) of "Simple Models of Complex Nuclei" by Talmi   !
!                 (ISBN: 978-3718605507)                                       !
!               = 2 LST coupling (among particles with same n and l)           !
!                 see Hinohara.2014.PhysRevC.90.031301                         !
!                 BB: it has never been benchmarked                            !
!------------------------------------------------------------------------------!
subroutine set_pairs

integer :: ia, ib, na, la, ja, ma, ta, sha, nb, lb, jb, mb, tb, shb, incr, &
           ml, ms, ialloc=0
real(r64) :: factor, facl, facs, cg1, cg2, cg3, cg4, cg5, cg6, phs_

allocate( pairs_T00_J1p1(HOsp_dim2), pairs_T1p1_J00(HOsp_dim2), &
          pairs_T00_J10(HOsp_dim2),  pairs_T10_J00(HOsp_dim2),  &
          pairs_T00_J1m1(HOsp_dim2), pairs_T1m1_J00(HOsp_dim2), &
          stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of pair operators.'

pairs_T00_J1p1 = zero
pairs_T00_J10  = zero
pairs_T00_J1m1 = zero
pairs_T1p1_J00 = zero
pairs_T10_J00  = zero
pairs_T1m1_J00 = zero

incr = 0
do ia = 1, HOsp_dim
  na  = HOsp_n(ia)
  la  = 2*HOsp_l(ia)
  ja  = HOsp_2j(ia)
  ma  = HOsp_2mj(ia)
  ta  = HOsp_2mt(ia)
  sha = HOsp_sh(ia)
  do ib = 1, HOsp_dim
    nb  = HOsp_n(ib)
    lb  = 2*HOsp_l(ib)
    jb  = HOsp_2j(ib)
    mb  = HOsp_2mj(ib)
    tb  = HOsp_2mt(ib)
    shb = HOsp_sh(ib)
    incr = incr + 1

    !!! Conditions on quantum numbers
    if ( pairs_scheme == 1 ) then
      if ( sha /= shb ) cycle
    elseif ( pairs_scheme == 2 ) then
      if ( (na /= nb) .or. (la /= lb) ) cycle
    endif

    !!! Factors for different options
    select case (pairs_scheme)
      case (0)
        factor = one
        if ( sha == shb ) factor =  1.0d0 / sqrt(2.0d0)
      case (1)
        factor = sqrt(ja + 1.0d0) / 2.0d0
      case (2)
        factor =  1.0d0 / sqrt(2.0d0)
    end select

    !!! T=0 (isoscalar)
    call ClebschGordan(1,1,0,ta,tb,0,cg1)
    if ( pairs_scheme < 2 ) then
      call ClebschGordan(ja,jb,2,ma,mb, 2,cg2)
      call ClebschGordan(ja,jb,2,ma,mb, 0,cg3)
      call ClebschGordan(ja,jb,2,ma,mb,-2,cg4)
    else
      cg2 = zero
      cg3 = zero
      cg4 = zero
      do ml = -la, la, 2
        facl = (-1)**((la-ml)/2)
        call ClebschGordan(la,1,ja, ml,1,ma,cg5)
        call ClebschGordan(lb,1,jb,-ml,1,mb,cg6)
        cg2 = cg2 + facl * cg5 * cg6
        call ClebschGordan(la,1,ja, ml,-1,ma,cg5)
        call ClebschGordan(lb,1,jb,-ml,-1,mb,cg6)
        cg4 = cg4 + facl * cg5 * cg6
        do ms = -1, 1, 2
          facs = 1.0d0 / sqrt(2.0d0)
          call ClebschGordan(la,1,ja, ml, ms,ma,cg5)
          call ClebschGordan(lb,1,jb,-ml,-ms,mb,cg6)
          cg3 = cg3 + facl * facs * cg5 * cg6
        enddo
      enddo
    endif
    pairs_T00_J1p1(incr) = factor * cg1 * cg2
    pairs_T00_J10(incr)  = factor * cg1 * cg3
    pairs_T00_J1m1(incr) = factor * cg1 * cg4

    !!! T=1 (isovector)
    call ClebschGordan(1,1,2,ta,tb, 2,cg2)
    call ClebschGordan(1,1,2,ta,tb, 0,cg3)
    call ClebschGordan(1,1,2,ta,tb,-2,cg4)
    if ( pairs_scheme < 2 ) then
      call ClebschGordan(ja,jb,0,ma,mb,0,cg1)
    else
      cg1 = zero
      do ml = -la, la, 2
        facl = (-1)**((la-ml)/2)
        do ms = -1, 1, 2
          facs = (-1)**((1-ms)/2) / sqrt(2.0d0)
          call ClebschGordan(la,1,ja, ml, ms,ma,cg5)
          call ClebschGordan(lb,1,jb,-ml,-ms,mb,cg6)
          cg1 = cg1 + facl * facs * cg5 * cg6
        enddo
      enddo
    endif
    pairs_T1p1_J00(incr) = factor * cg1 * cg2
    pairs_T10_J00(incr)  = factor * cg1 * cg3
    pairs_T1m1_J00(incr) = factor * cg1 * cg4

  enddo
enddo

end subroutine set_pairs


!------------------------------------------------------------------------------!
! subroutine calculate_pairCoupl2B_ben                                         !
!                                                                              !
! Calculates the expected value of the 2Body pair coupled Operator (delafuen)  !
! including all the terms from the rho matrices (constant Z,N,A inflation)     !
!    (Set of a logical to exclude the rho parts or not)                        !
!                                                                              !
! seniorityScheme (pairs_scheme)                                               !
!               = 0 Agnostic (no special factor or selection rule)             !
!               = 1 Seniority (only pairs among particles of the same shell)   !
!                 see Eq(19.1) of "Simple Models of Complex Nuclei" by Talmi   !
!                 (ISBN: 978-3718605507)                                       !
!               = 2 LST coupling (cannot be used, return status 0)             !
!------------------------------------------------------------------------------!
subroutine calculate_pairCoupl2B_ben(io,rhoLR,kappaLR,kappaRL, &
                                     pair_T00_J1m1,pair_T00_J10,pair_T00_J1p1,&
                                     pair_T1m1_J00,pair_T10_J00,pair_T1p1_J00,&
                                     ndim, seniorityScheme)

integer, intent(in) :: ndim, io, seniorityScheme
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, kappaLR, kappaRL
real(r64), intent(out) :: pair_T00_J1m1,pair_T00_J10,pair_T00_J1p1, &
                          pair_T1m1_J00,pair_T10_J00,pair_T1p1_J00
integer      :: hdim, a, b, a2, b2, ta, tb, ja, jb, ma, mb, ma2, mb2, M,&
                ia, ib, a0, b0
complex(r64) :: p2B_T00_J1p1, p2B_T00_J1m1, p2B_T00_J10, &
                p2B_T1p1_J00, p2B_T1m1_J00, p2B_T10_J00
real(r64) :: N_ab_J0T1, cgj1,cgj2,cgt1,cgt2, aux, aux2, aux1B ! N_ab_J1T0 == N_ab_J0T1
logical   :: INCLUDE_RHO_TERMS = .FALSE.

if ( (io /= 0) .and. (io /= 1) ) then
  print*,'Wrong argument in calculate_pairCoupl2B_bench: io = ', io
  stop
endif
if (seniorityScheme.eq.2) return

hdim = ndim/2

p2B_T00_J1p1 = zzero
p2B_T00_J1m1 = zzero
p2B_T00_J10  = zzero
p2B_T1p1_J00 = zzero
p2B_T1m1_J00 = zzero
p2B_T10_J00  = zzero
aux1B = zero

!! Operator pair2B_T00_J1p1   pair2B_T00_J1m1   pair2B_T00_J10
do a = 1, hdim
  ja = HOsp_2j (a)
  ma = HOsp_2mj(a)
  ia = (ja-ma) / 2 ! index of a in the program order [+ja, +ja-1, ...,  -ja]
  a0 = a - ia      ! limit to read the (ja, alpha')

  if (INCLUDE_RHO_TERMS) then
    aux1B = aux1B + 0.5d+0 * real(rhoLR(a,a) + rhoLR(a+hdim,a+hdim))
  endif

  do b = 1, hdim
    jb = HOsp_2j (b)
    mb = HOsp_2mj(b)
    ib = (jb-mb) / 2 ! index of a in the program order [+ja, +ja-1, ...,  -ja]
    b0 = b - ib      ! limit to read the (jb, beta')

    if (abs(ma-mb)>2) continue !! M = 0, 1, -1 only

    N_ab_J0T1 = 1.0d0
    if (seniorityScheme.eq.1) then
      if (HOsp_sh(a).ne.HOsp_sh(b)) continue
      N_ab_J0T1 = (ja + 1.0d0) * 0.5d0  ! squared
    endif

    ! index loop for the density matrices
    do a2 = a0, a0 + ja
      ma2 = HOsp_2mj(a2)
      do b2 = b0, b0 + jb
        mb2 = HOsp_2mj(b2)
        if (abs(ma2-mb2)>2) continue !! M = 0, 1, -1 only

          !! Term pn J=0, T=1, MT= +/- 1  -------------------------------------

          call ClebschGordan(ja,jb,0 , ma, -ma,0, cgj1)
          call ClebschGordan(ja,jb,0 ,ma2,-ma2,0, cgj2)
          cgt1 = 1.0d0

          aux = kappaRL(a,b) * kappaLR(a2,b2)
          if (INCLUDE_RHO_TERMS) then
            aux = aux + rhoLR(a2,a) * rhoLR(b2,b)
            aux = aux - rhoLR(b2,a) * rhoLR(a2,b)
          endif
          p2B_T1m1_J00 = p2B_T1m1_J00 + (aux * N_ab_J0T1 * cgj1 * cgj2 * cgt1)

          aux = kappaRL(a+hdim,b+hdim) * kappaLR(a2+hdim,b2+hdim)
          if (INCLUDE_RHO_TERMS) then
            aux = aux + rhoLR(a2+hdim,a+hdim) * rhoLR(b2+hdim,b+hdim)
            aux = aux - rhoLR(b2+hdim,a+hdim) * rhoLR(a2+hdim,b+hdim)
          endif
          p2B_T1p1_J00 = p2B_T1p1_J00 + (aux * N_ab_J0T1 * cgj1 * cgj2 * cgt1)

          !! Term pn J=0, T=1, MT= 0 (same clebsh gordan for J) ---------------
          cgt1 = 0.5d0 ! C10+- = C10-+

          aux =       kappaRL(a     ,b+hdim) * kappaLR(a2     ,b2+hdim)
          aux = aux + kappaRL(a     ,b+hdim) * kappaLR(a2+hdim,b2     )
          aux = aux + kappaRL(a+hdim,b     ) * kappaLR(a2     ,b2+hdim)
          aux = aux + kappaRL(a+hdim,b     ) * kappaLR(a2+hdim,b2     )

          if (INCLUDE_RHO_TERMS) then
            aux = aux + rhoLR(a2     ,a     ) * rhoLR(b2+hdim,b+hdim)
            aux = aux + rhoLR(a2+hdim,a     ) * rhoLR(b2     ,b+hdim)
            aux = aux + rhoLR(a2     ,a+hdim) * rhoLR(b2+hdim,b     )
            aux = aux + rhoLR(a2+hdim,a+hdim) * rhoLR(b2     ,b     )

            aux = aux - rhoLR(b2+hdim,a     ) * rhoLR(a2     ,b+hdim)
            aux = aux - rhoLR(b2     ,a     ) * rhoLR(a2+hdim,b+hdim)
            aux = aux - rhoLR(b2+hdim,a+hdim) * rhoLR(a2     ,b     )
            aux = aux - rhoLR(b2     ,a+hdim) * rhoLR(a2+hdim,b     )
          endif

          p2B_T10_J00 = p2B_T10_J00 + (aux * N_ab_J0T1 * cgj1 * cgj2 * cgt1)

          !! Term pn J=1,M   T=0, MT= 0  --------------------------------------
          !! M requires to be::  ma+mb = ma2 + mb2
          M = ma + mb
          if (M.ne.(ma2 + mb2)) continue ! Not interested in a pair M,M'

          call ClebschGordan(ja,jb,2 , ma,  mb,M, cgj1)
          call ClebschGordan(ja,jb,2 ,ma2, mb2,M, cgj2)
          cgt1 = 0.5d0 ! C00+- = -1 * C00-+ (sing included in the sums)

          aux =       kappaRL(a     ,b+hdim) * kappaLR(a2     ,b2+hdim)
          aux = aux - kappaRL(a     ,b+hdim) * kappaLR(a2+hdim,b2     )
          aux = aux - kappaRL(a+hdim,b     ) * kappaLR(a2     ,b2+hdim)
          aux = aux + kappaRL(a+hdim,b     ) * kappaLR(a2+hdim,b2     )

          if (INCLUDE_RHO_TERMS) then
            aux = aux + rhoLR(a2     ,a     ) * rhoLR(b2+hdim,b+hdim)
            aux = aux - rhoLR(a2+hdim,a     ) * rhoLR(b2     ,b+hdim)
            aux = aux - rhoLR(a2     ,a+hdim) * rhoLR(b2+hdim,b     )
            aux = aux + rhoLR(a2+hdim,a+hdim) * rhoLR(b2     ,b     )

            aux = aux - rhoLR(b2+hdim,a     ) * rhoLR(a2     ,b+hdim)
            aux = aux + rhoLR(b2     ,a     ) * rhoLR(a2+hdim,b+hdim)
            aux = aux + rhoLR(b2+hdim,a+hdim) * rhoLR(a2     ,b     )
            aux = aux - rhoLR(b2     ,a+hdim) * rhoLR(a2+hdim,b     )
          endif

          if      (M.eq.2) then
            p2B_T00_J1m1 = p2B_T00_J1m1 + (aux * N_ab_J0T1 * cgj1 * cgj2 * cgt1)
          else if (M.eq.0) then
            p2B_T00_J10  = p2B_T00_J10  + (aux * N_ab_J0T1 * cgj1 * cgj2 * cgt1)
          else
            p2B_T00_J1p1 = p2B_T00_J1p1 + (aux * N_ab_J0T1 * cgj1 * cgj2 * cgt1)
          end if

      enddo
    enddo

  end do
end do

aux2 = 1.0
if (seniorityScheme.eq.1) aux2 = 0.7071067811865476

call ClebschGordan(1,1,0, 1,-1,0, aux)
p2B_T00_J1p1 = p2B_T00_J1p1 * 0.50d0 * aux2 * (aux)**2.0d0
p2B_T00_J1m1 = p2B_T00_J1m1 * 0.50d0 * aux2 * (aux)**2.0d0
p2B_T00_J10  = p2B_T00_J10  * 0.50d0 * aux2 * (aux)**2.0d0

call ClebschGordan(1,1,2, 1,1,2, aux)
p2B_T1p1_J00 = p2B_T1p1_J00 * 0.25d0 * aux2 * (aux)**2.0d0
p2B_T1m1_J00 = p2B_T1m1_J00 * 0.25d0 * aux2 * (aux)**2.0d0

call ClebschGordan(1,1,2, 1,-1,0, aux)
p2B_T10_J00  = p2B_T10_J00  * 0.50d0 * aux2 * (aux)**2.0d0

pair_T00_J1p1 = real(p2B_T00_J1p1) + aux1B
pair_T00_J1m1 = real(p2B_T00_J1m1) + aux1B
pair_T00_J10  = real(p2B_T00_J10)  + aux1B
pair_T1p1_J00 = real(p2B_T1p1_J00) + aux1B
pair_T1m1_J00 = real(p2B_T1m1_J00) + aux1B
pair_T10_J00  = real(p2B_T10_J00)  + aux1B

end subroutine calculate_pairCoupl2B_ben

!------------------------------------------------------------------------------!
! subroutine calculate_particle_number                                         !
!                                                                              !
! Computes the expectation values of the proton and neutron number operators.  !
! Extended to the pn-terms for the variance evaluation.                        !
!------------------------------------------------------------------------------!
subroutine calculate_variance_components(rhoLR, kappaLR, kappaRL, &
                                         prot2, var_pn2, neut2, ndim)

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, kappaLR, kappaRL
real(r64), intent(out) :: prot2, var_pn2, neut2
integer :: hdim, j, i
complex(r64) :: tr1, tr2, tr3, tr4, tr5, tr6, prot, neut
complex(r64), dimension(ndim/2,ndim/2) :: rhoLRp, kapLRp, kapRLp, A1, A2, &
                                          rhoLRn, kapLRn, kapRLn, A3, A4, &
                                          kapLRm, kapRLm, A5, rhoLRm, A6

hdim = ndim/2

prot  = zzero
neut  = zzero
prot2 = zero
neut2 = zero

!!! N and Z
do i = 1, hdim
  prot = prot + rhoLR(i,i)
  neut = neut + rhoLR(i+hdim,i+hdim)
enddo

!!! N^2 and Z^2
do j = 1, hdim
  do i = 1, hdim
    rhoLRp(i,j) = rhoLR(i,j)
    rhoLRn(i,j) = rhoLR(i+hdim,j+hdim)
    rhoLRm(i,j) = rhoLR(i,j+hdim)
    kapLRp(i,j) = kappaLR(i,j)
    kapLRn(i,j) = kappaLR(i+hdim,j+hdim)
    kapLRm(i,j) = kappaLR(i,j+hdim)
    kapRLp(i,j) = kappaRL(i,j)
    kapRLn(i,j) = kappaRL(i+hdim,j+hdim)
    kapRLm(i,j) = kappaRL(i,j+hdim)
  enddo
enddo

call zgemm('n','n',hdim,hdim,hdim,zone,rhoLRp,hdim,rhoLRp,hdim,zzero,A1,hdim)
call zgemm('t','n',hdim,hdim,hdim,zone,kapRLp,hdim,kapLRp,hdim,zzero,A2,hdim)
call zgemm('n','n',hdim,hdim,hdim,zone,rhoLRn,hdim,rhoLRn,hdim,zzero,A3,hdim)
call zgemm('t','n',hdim,hdim,hdim,zone,kapRLn,hdim,kapLRn,hdim,zzero,A4,hdim)
call zgemm('n','n',hdim,hdim,hdim,zone,rhoLRm,hdim,rhoLRm,hdim,zzero,A6,hdim)
call zgemm('t','n',hdim,hdim,hdim,zone,kapRLm,hdim,kapLRm,hdim,zzero,A5,hdim)

tr1 = zzero
tr2 = zzero
tr3 = zzero
tr4 = zzero
tr5 = zzero
tr6 = zzero

do i = 1, hdim
  tr1 = tr1 + A1(i,i)
  tr2 = tr2 + A2(i,i)
  tr3 = tr3 + A3(i,i)
  tr4 = tr4 + A4(i,i)
  tr5 = tr5 - A5(i,i)
  tr6 = tr6 + A6(i,i)
enddo

prot2 = dreal(prot - tr1 + tr2)
neut2 = dreal(neut - tr3 + tr4)
var_pn2 = dreal(prot2 + neut2 - 2 * tr5 - 2 *tr6)


end subroutine calculate_variance_components
!!! --------------------------------------------------------------------------
END MODULE Pairs
!==============================================================================!
! End of file                                                                  !
!==============================================================================!

!==============================================================================!
! MODULE ParticleNumber                                                        !
!                                                                              !
! This module contains the variables and routines related to particle-number   !
! operators (including their square and occupation number).                    !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_particle_number                                             !
! - subroutine calculate_particle_number                                       !
! - subroutine calculate_occupation_number                                     !
!==============================================================================!
MODULE ParticleNumber

use Basis

implicit none

real(r64), dimension(:), allocatable :: partnumb_Z, & ! Proton  number operator
                                        partnumb_N, & ! Neutron    "      "
                                        partnumb_A    ! Nucleon    "      "
CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_particle_number                                               !
!                                                                              !
! Defines the HO matrix elements of the proton, neutron and nucleon number op- !
! erators. They are simply the identity matrices in their relative subspace.   !
! The operators are vectorized (for the calculation of constraints).           !
!------------------------------------------------------------------------------!
subroutine set_particle_number

integer :: i, j, k, ialloc=0

allocate( partnumb_Z(HOsp_dim2), partnumb_N(HOsp_dim2), partnumb_A(HOsp_dim2), &
          stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of particle-number operator'

partnumb_Z = zero
partnumb_N = zero
partnumb_A = zero

k = 0
do i = 1, HOsp_dim
  do j = 1, HOsp_dim
    k = k + 1
    if ( i /= j ) cycle
    if ( i <= HOsp_dim/2 ) then
      partnumb_Z(k) = one
    else
      partnumb_N(k) = one
    endif
  enddo
enddo

partnumb_A = partnumb_Z + partnumb_N

end subroutine set_particle_number

!------------------------------------------------------------------------------!
! subroutine calculate_particle_number                                         !
!                                                                              !
! Computes the expectation values of the proton and neutron number operators.  !
! < Z > = Tr(rhoLR_p) = Tr(rhoLR) from 1            to HOsp_dim/2              !
! < N > = Tr(rhoLR_n) = Tr(rhoLR) from 1+HOsp_dim/2 to HOsp_dim                !
! < Z^2 > = Tr(rhoLR_p) + Tr(rhoLR_p)^2 - Tr(rhoLR_p^2)                        !
!           + Tr(kappaRL^T_p kappaLR_p)                                        !
! < N^2 > = Tr(rhoLR_n) + Tr(rhoLR_n)^2 - Tr(rhoLR_n^2)                        !
!           + Tr(kappaRL^T_n kappaLR_n)                                        !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        io = 0 computes < Z >, < N >, < Z^2 >, < N^2 >                        !
!           = 1 computes < Z >, < N >                                          !
!        rhoLR,kappaLR,kappaRL = transition densities                          !
!                                                                              !
! Output: prot, prot2 = < Z > and < Z^2 >                                      !
!         neut, neut2 = < N > and < N^2 >                                      !
!------------------------------------------------------------------------------!
subroutine calculate_particle_number(io,rhoLR,kappaLR,kappaRL,prot,neut, &
                                     prot2,neut2,ndim)

integer, intent(in) :: ndim, io
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, kappaLR, kappaRL
complex(r64), intent(out) :: prot, neut, prot2, neut2
integer :: hdim, j, i
complex(r64) :: tr1, tr2, tr3, tr4
complex(r64), dimension(ndim/2,ndim/2) :: rhoLRp, kapLRp, kapRLp, A1, A2, &
                                          rhoLRn, kapLRn, kapRLn, A3, A4

if ( (io /= 0) .and. (io /= 1) ) then
  print*,'Wrong argument in calculate_particle_number: io = ', io
  stop
endif

hdim = ndim/2

prot  = zzero
neut  = zzero
prot2 = zzero
neut2 = zzero

!!! N and Z
do i = 1, hdim
  prot = prot + rhoLR(i,i)
  neut = neut + rhoLR(i+hdim,i+hdim)
enddo

!!! N^2 and Z^2
if ( io == 0 ) then

  do j = 1, hdim
    do i = 1, hdim
      rhoLRp(i,j) = rhoLR(i,j)
      rhoLRn(i,j) = rhoLR(i+hdim,j+hdim)
      kapLRp(i,j) = kappaLR(i,j)
      kapLRn(i,j) = kappaLR(i+hdim,j+hdim)
      kapRLp(i,j) = kappaRL(i,j)
      kapRLn(i,j) = kappaRL(i+hdim,j+hdim)
    enddo
  enddo

  call zgemm('n','n',hdim,hdim,hdim,zone,rhoLRp,hdim,rhoLRp,hdim,zzero,A1,hdim)
  call zgemm('t','n',hdim,hdim,hdim,zone,kapRLp,hdim,kapLRp,hdim,zzero,A2,hdim)
  call zgemm('n','n',hdim,hdim,hdim,zone,rhoLRn,hdim,rhoLRn,hdim,zzero,A3,hdim)
  call zgemm('t','n',hdim,hdim,hdim,zone,kapRLn,hdim,kapLRn,hdim,zzero,A4,hdim)

  tr1 = zzero
  tr2 = zzero
  tr3 = zzero
  tr4 = zzero

  do i = 1, hdim
    tr1 = tr1 + A1(i,i)
    tr2 = tr2 + A2(i,i)
    tr3 = tr3 + A3(i,i)
    tr4 = tr4 + A4(i,i)
  enddo

  prot2 = prot + prot**2 - tr1 + tr2
  neut2 = neut + neut**2 - tr3 + tr4

endif

end subroutine calculate_particle_number

!!! delafuen -----------------------------------------------------------------
!  Test Version of the Variance mean values for the PN mixed part of the wf    !
!                                                                              !
!                                                                              !
!                                                                              !
!------------------------------------------------------------------------------!

subroutine calculate_particle_number_pn(io,rhoLR,kappaLR,kappaRL, &
!                                        pn, np, pnpn, pnnp, nppn, npnp,
                                        ndim)

integer, intent(in) :: ndim, io
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, kappaLR, kappaRL
!complex(r64), intent(out) :: pn, np, pnpn, pnnp, nppn, npnp
complex(r64) :: pn, np, pnpn, pnnp, nppn, npnp, pnpn_new, var_pn
integer :: hdim, j, i
complex(r64) :: tr1, tr2, tr3, tr4, tr5, tr6, tr7, tr8, prot, neut, &
                tr9, tr10, tr11, tr12
complex(r64), dimension(ndim/2,ndim/2) :: rhoLRp, kapLRp, kapRLp, A1, A2, &
                                          rhoLRn, kapLRn, kapRLn, A3, A4, &
                                          rhoLRpn, kapLRpn, kapRLpn, A5, A6, &
                                          rhoLRnp, kapLRnp, kapRLnp, A7, A8, &
                                          A9, A10, A11, A12


if ( (io /= 0) .and. (io /= 1) ) then
  print*,'Wrong argument in calculate_particle_number: io = ', io
  stop
endif

hdim = ndim/2

pn  = zzero
np  = zzero
pnpn = zzero
pnnp = zzero
nppn = zzero
npnp = zzero
prot = zzero
neut = zzero

!!! N and Z
do i = 1, hdim
  pn = pn + rhoLR(i,i+hdim)
  np = np + rhoLR(i+hdim,i)
  prot = prot + rhoLR(i,i)
  neut = neut + rhoLR(i+hdim,i+hdim)
enddo

!!! N^2 and Z^2
if ( io == 0 ) then

  do j = 1, hdim
    do i = 1, hdim
      rhoLRp(i,j) = rhoLR(i,j)
      rhoLRn(i,j) = rhoLR(i+hdim,j+hdim)
      kapLRp(i,j) = kappaLR(i,j)
      kapLRn(i,j) = kappaLR(i+hdim,j+hdim)
      kapRLp(i,j) = kappaRL(i,j)
      kapRLn(i,j) = kappaRL(i+hdim,j+hdim)

      rhoLRpn(i,j) = rhoLR(i,j+hdim)
      rhoLRnp(i,j) = rhoLR(i+hdim,j)
      kapLRpn(i,j) = kappaLR(i,j+hdim)
      kapLRnp(i,j) = kappaLR(i+hdim,j)
      kapRLpn(i,j) = kappaRL(i,j+hdim)
      kapRLnp(i,j) = kappaRL(i+hdim,j)
    enddo
  enddo

call zgemm('n','n',hdim,hdim,hdim,zone,rhoLRpn,hdim,rhoLRpn,hdim,zzero,A1,hdim)
call zgemm('n','n',hdim,hdim,hdim,zone,rhoLRpn,hdim,rhoLRnp,hdim,zzero,A6,hdim)
call zgemm('n','n',hdim,hdim,hdim,zone,rhoLRp ,hdim,rhoLRp ,hdim,zzero,A7,hdim)
call zgemm('n','n',hdim,hdim,hdim,zone,rhoLRn ,hdim,rhoLRn ,hdim,zzero,A8,hdim)

call zgemm('t','n',hdim,hdim,hdim,zone,kapRLp ,hdim,kapLRn ,hdim,zzero,A2,hdim)
call zgemm('t','n',hdim,hdim,hdim,zone,kapRLn ,hdim,kapLRp ,hdim,zzero,A3,hdim)
call zgemm('t','n',hdim,hdim,hdim,zone,kapRLpn,hdim,kapLRnp,hdim,zzero,A4,hdim)
call zgemm('t','n',hdim,hdim,hdim,zone,kapRLpn,hdim,kapLRpn,hdim,zzero,A5,hdim)

call zgemm('n','n',hdim,hdim,hdim,zone,rhoLRp ,hdim,rhoLRn,hdim,zzero,A9 ,hdim)
call zgemm('n','n',hdim,hdim,hdim,zone,rhoLRn ,hdim,rhoLRp,hdim,zzero,A10,hdim)
call zgemm('t','n',hdim,hdim,hdim,zone,kapRLp ,hdim,kapLRp,hdim,zzero,A11,hdim)
call zgemm('t','n',hdim,hdim,hdim,zone,kapRLn ,hdim,kapLRn,hdim,zzero,A12,hdim)

  tr1 = zzero
  tr2 = zzero
  tr3 = zzero
  tr4 = zzero
  tr5 = zzero
  tr6 = zzero
  tr7 = zzero
  tr8 = zzero
  tr9 = zzero
  tr10= zzero
  tr11= zzero
  tr12= zzero

  do i = 1, hdim
    tr1 = tr1 + A1(i,i) !pn pn
    tr6 = tr6 + A6(i,i) !np np      !! CAMBIO a pn np
    tr7 = tr7 + A7(i,i) !p  p
    tr8 = tr8 + A8(i,i) !n  n
    tr9 = tr9 + A9(i,i) !nn pp
    tr10 = tr10 + A10(i,i)!pp nn

    tr2 = tr2 + A2(i,i) !RLp^T LRn
    tr3 = tr3 + A3(i,i) !RLn^T LRp
    tr4 = tr4 + A4(i,i) !RLpn^T LRnp
    tr5 = tr5 + A5(i,i) !RLnp^T LRpn     !! CAMBIO: RLpn^T LRpn (no se usa)
    tr11= tr11 + A11(i,i) !!RLp^T LRp
    tr12= tr12 + A12(i,i) !!RLn^T LRn

  enddo

print "(A)", ""
print "(2A)", "Tr:      r pp      r nn   r pp pp   r nn nn   r pn pn   ", &
              "r pn np   r pp nn  /KK   k pp pp   k nn nn   k pn np   k pp nn"
print "(A,7F10.3,A,4F10.4)","  :", real(prot), real(neut), &
      real(tr7), real(tr8), real(tr1), real(tr6), real(tr10), "  /  ", &
      real(tr11),real(tr12),real(tr4), real(tr2)

  pnpn = tr2 + (pn * np) - tr1
  pnnp = tr4 + (neut * prot) - tr6 + prot
  nppn = tr4 + (prot * neut) - tr6 + neut
  npnp = tr3 + (np * pn) - tr1

  pnpn_new = (tr4 + tr2 + tr3 + tr4) + (4*pn*np) - (tr9 + tr1 + tr1 + tr10)
  pnpn_new = 0.25 * (neut + prot + pnpn_new)

  var_pn = pnpn_new - (np * np) !- (neut + prot - tr7 - tr8 + tr11 + tr12)

endif

print "(2A)", "NP^2 vr:       tr_pn   O2_pn_sym    VAR_O_pn     M pn pn", &
              "     M np np     M pn np     M np pn    "
print "(A,7F12.6)", "       :", real(pn), real(pnpn_new), real(var_pn),&
                    real(pnpn), real(npnp), real(pnnp), real(nppn)


end subroutine calculate_particle_number_pn


!------------------------------------------------------------------------------!
! subroutine calculate_occupation_number                                       !
!                                                                              !
! Calculates the occupation numbers of the spherical shell, i.e. the number    !
! of particcle in all the 2j+1 multiplet forming a shell.                      !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        sdim = dimension of the shells                                        !
!        rhoLR = transition density                                            !
!                                                                              !
! Output: occnum = occupation numbers                                          !
!------------------------------------------------------------------------------!
subroutine calculate_occupation_number(rhoLR,occnum,ndim,sdim)

integer, intent(in) :: ndim, sdim
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR
complex(r64), dimension(sdim,2), intent(out) :: occnum
integer :: i, j, k, jmax

k = 0
occnum = zzero

do i = 1, sdim
  jmax = HOsh_2j(i) + 1
  do j = 1, jmax
    k = k + 1
    occnum(i,1) = occnum(i,1) + rhoLR(k,k)
    occnum(i,2) = occnum(i,2) + rhoLR(k+ndim/2,k+ndim/2)
  enddo
enddo

end subroutine calculate_occupation_number

END MODULE ParticleNumber
!==============================================================================!
! End of file                                                                  !
!==============================================================================!

!==============================================================================!
! MODULE MathMethods                                                           !
!                                                                              !
! This module contains the variables and routines related to mathematical fun- !
! ctions and numerical methods.                                                !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine GaussLaguerre                                                   !
! - subroutine GaussLegendre                                                   !
! - subroutine ClebschGordan                                                   !
! - function factorial                                                         !
! - function dfactorial                                                        !
! - function kdelta                                                            !
! - function genlaguerre                                                       !
! - function assolegendre                                                      !
! - function lapack_sel                                                        !
!==============================================================================!
MODULE MathMethods

use Constants

implicit none
public

CONTAINS

!------------------------------------------------------------------------------!
! function factorial                                                           !
!                                                                              !
! Computes the factorial: n! = n * (n-1) * ... * 1                             !
!------------------------------------------------------------------------------!
recursive function factorial(n) result(facto)

integer, intent(in) :: n
real(r64) :: facto

if ( n <= 0 ) then
  facto = one
else
  facto = n * factorial(n-1)
endif

end function

!------------------------------------------------------------------------------!
! function dfactorial                                                          !
!                                                                              !
! Computes the double factorial: n!! = n * (n-2) * ... * 1                     !
!------------------------------------------------------------------------------!
recursive function dfactorial(n) result(facto)

integer, intent(in) :: n
real(r64) :: facto

if ( n <= 0 ) then
  facto = one
else
  facto = n * dfactorial(n-2)
endif

end function

!------------------------------------------------------------------------------!
! function kdelta                                                              !
!                                                                              !
! Computes the Kronecker delta: \delta_ij = 1 if i = j                         !
!                                         = 0 otherwise                        !
!------------------------------------------------------------------------------!
function kdelta(i,j) result(delta)

integer, intent(in) :: i, j
integer :: delta

if ( i == j ) then
  delta = 1
else
  delta = 0
endif

end function kdelta

!------------------------------------------------------------------------------!
! function genlaguerre                                                         !
!                                                                              !
! Calculates the value of the generalized Laguerre polynomial L^a_n(x)         !
! using the reccurence relation:                                               !
! L^a_n(x) = ( (2n-1+a-x)*L^a_n-1(x) - (n-1+a)*L^a_n-2(x) ) / n                !
!------------------------------------------------------------------------------!
function genlaguerre(a,n,x) result(l1)

integer, intent(in) :: n
real(r64), intent(in) :: x, a
integer :: i
real(r64) :: l1, l2, l3

l1 = 1.d0
l2 = 0.d0

do i = 1, n
  l3 = l2
  l2 = l1
  l1 = ( (2*i-1+a-x)*l2 - (i-1+a)*l3 ) / i
enddo

end function genlaguerre

!------------------------------------------------------------------------------!
! function assolegendre                                                        !
!                                                                              !
! Calculates the value of the associated Legendre polynomial P^m_l(x) using the!
! algorithm given in the Numerical Recipes in Fortran 90, which is based on the!
! recurrence relations:                                                        !
! P^m_m(x) = (-1)**m (2m-1)!! sqrt(1-x**2)  => build P^m_m for any m           !
! P^m_m+1(x) = x (2m+1) P^m_m(x)            => build P^m_m+1 from previous one !
! (l-m+1) P^m_l+1 = x (2l+1) P^m_l(x) - (l+m) P^m_l-1(x)                       !
!                                           => build all other possible l      !
!                                                                              !
! In addition, I have added the symmetry for negative m:                       !
! P^-m_l = (-1)**m (l-m)!/(l+m)!                                               !
!------------------------------------------------------------------------------!
function assolegendre(l,n,x) result(leg)

integer, intent(in) :: l, n
real(r64), intent(in) :: x
integer :: ll, m
real(r64) :: leg, phase, pll, pmm, pmmp1, somx2

!!! Checks the validity of the arguments
if ( (l < 0) .or. (abs(m) > l) .or. (abs(x) > 1) ) then
  print "(a)", "Wrong argument(s) in function assolegendre"
endif

!!! Transforms m in -m when m is negative
m = n
phase = one

if ( m < 0 ) then
  m = abs(n)
  phase = (-1)**m * (one*factorial(l-m)) / (one*factorial(l+m))
endif

!!! Compute P^m_m
pmm = one

if ( m > 0 ) then
  somx2 = sqrt( (one-x) * (one+x) )
  pmm = dfactorial(2*m-1) * somx2**m
  if ( mod(m,2) == 1 ) pmm = -pmm
endif

if ( l == m ) then
  leg = pmm
else
  pmmp1 = x * (2*m +1) * pmm

  !!! Compute P^m_m+1
  if ( l == m+1 ) then
    leg = pmmp1
  !!! Compute P^m_l for l > m+1
  else
    do ll = m+2, l
      pll = (x*(2*ll-1)*pmmp1 - (ll+m-1)*pmm) / (ll-m)
      pmm = pmmp1
      pmmp1 = pll
    enddo
    leg = pll
  endif
endif

leg = phase * leg

end function assolegendre

!------------------------------------------------------------------------------!
! subroutine GaussLaguerre                                                     !
!                                                                              !
! Computes the abscissas and weight for the Gauss-Laguerre quadrature.         !
! The algorithm is taken from the book: Numerical Recipe in Fortran 77         !
! (ISBN: 978-0521430647).                                                      !
!------------------------------------------------------------------------------!
subroutine GaussLaguerre(x,w,n,alf)

integer, intent(in) :: n
real(r64), intent(in) :: alf
real(r64), intent(out) :: w(n), x(n)
integer, parameter :: maxit=20
real(r64), parameter :: eps=1.0d-16
integer :: i, its, j
real(r64) :: ai, p1, p2, p3, pp, z, z1

do i = 1, n
  if ( i == 1 ) then
    z = (1.0 + alf) * (3.0 + 0.92*alf) / (1.0 + 2.4*n + 1.8*alf)
  elseif ( i == 2 ) then
    z = z + (15.0 + 6.25*alf) / (1.0 + 0.9*alf + 2.5*n)
  else
    ai = i - 2
    z = z + ( ((1.0 + 2.55*ai) / (1.9*ai)) + (1.26*ai*alf / (1.0 + 3.5*ai)) ) &
            * (z - x(i-2)) / (1.0 + 0.3*alf)
  endif

  do its = 1, maxit
    p1 = 1.d0
    p2 = 0.d0
    do j = 1, n
      p3 = p2
      p2 = p1
      p1 = ( (2*j-1+alf-z)*p2 - (j-1+alf)*p3 ) / j
    enddo
    pp = (n*p1 - (n+alf)*p2) / z
    z1 = z
    z = z1 - p1/pp
    if ( abs(z - z1) <= eps) exit
  enddo

  x(i) = z
  w(i) = -exp(log_gamma(alf+n) - log_gamma(n+0.0d0)) / (pp*n*p2)
enddo

end subroutine GaussLaguerre

!------------------------------------------------------------------------------!
! subroutine GaussLegendre                                                     !
!                                                                              !
! Computes the abscissas and weight for the Gauss-Legendre quadrature.         !
! The algorithm is taken from the book: Numerical Recipe in Fortran 77         !
! (ISBN: 978-0521430647).                                                      !
!------------------------------------------------------------------------------!

subroutine GaussLegendre(x1,x2,x,w,n)

integer, intent(in) :: n
real(r64), intent(in) :: x1, x2
real(r64), intent(out) :: x(n), w(n)
real(r64), parameter :: eps=3.d-14
integer :: i, j, m
real(r64) :: p1, p2, p3, pp, xl, xm, z, z1=1.5d0

m = (n + 1) / 2
xm = 0.5d0 * (x2 + x1)
xl = 0.5d0 * (x2 - x1)

do i = 1, m
  z = cos( pi * (i-0.25d0) / (n+0.5d0) )
  do while ( abs(z-z1) .gt. eps )
    p1 = one
    p2 = zero
    do j = 1, n
      p3 = p2
      p2 = p1
      p1 = ( (2.d0*j-1.d0)*z*p2 - (j-1.d0)*p3 ) / j
    enddo
    pp = n*(z*p1-p2) / (z*z-1.d0)
    z1 = z
    z = z1 - p1/pp
  enddo
  x(i) = xm - xl*z
  x(n+1-i) = xm + xl*z
  w(i) = 2.d0 * xl / ( (1.d0-z*z)*pp*pp )
  w(n+1-i) = w(i)
enddo

end subroutine GaussLegendre

!------------------------------------------------------------------------------!
! subroutine ClebschGordan                                                     !
!                                                                              !
! Computes the ClebschGordan for the group SU(2). The algorithm used was taken !
! from technical notes from NASA written by W. F. Ford and R. C. Bruley.       !
! Ref: NASA TN D-6173                                                          !
!                                                                              !
! (j1,m1,j2,m2|j3,m3) = c * g                                                  !
! with                                                                         !
! c = D(j1j2j3) * [(j1-m1)!(j2+m2)!(j1+m1)!(j2-m2)!(j3+m3)!(j3-m3)!]**1/2      !
! g = sqrt(2*j3+1) sum_l (-1)**l [(j1+j2-j3-l)!(j1-m1-l)!(j2+m2-l)!            !
!                                 (j3-j2+m1+l)!(j3-j1-m1+l)!l!]**-1            !
! D(j1j2j3) = [(j1+j2-j3)!(j2+j3-j1)!(j3+j1-j2)!]**1/2 / [(j1+j2+j3+1)!]**1/2  !
!                                                                              !
! Be aware that all input values of j and m are supposed to be twice their     !
! values (such that we can treat easily half-integer angular momenta).         !
!------------------------------------------------------------------------------!
subroutine ClebschGordan (j1,j2,j3,m1,m2,m3,cg)

integer, intent(in) :: j1, j2, j3, m1, m2, m3
real(r64), intent(out) :: cg
integer :: n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, k1, k2, k3, k4, k5, k6, &
           l, l1, l2
real(r64) :: p, q, h, hl, hlm1

cg = zero

!!! Computes the ingredients for the factor c
n1 = 1 + (j1 + j2 - j3)/2
n2 = 1 + (j2 + j3 - j1)/2
n3 = 1 + (j3 + j1 - j2)/2
n4 = 1 + (j1 - m1)/2
n5 = 1 + (j2 + m2)/2
n6 = 1 + (j1 + m1)/2
n7 = 1 + (j2 - m2)/2
n8 = 1 + (j3 + m3)/2
n9 = 1 + (j3 - m3)/2
n10= n1 + n2 + n3 - 1

if ( (min(n1,n2,n3,n4,n5,n6,n7,n8,n9) < 1) .or. (m1+m2 /= m3) ) return

p =  log_gamma(n1+zero) + log_gamma(n2+zero) + log_gamma(n3+zero) &
   + log_gamma(n4+zero) + log_gamma(n5+zero) + log_gamma(n6+zero) &
   + log_gamma(n7+zero) + log_gamma(n8+zero) + log_gamma(n9+zero) &
   - log_gamma(n10+zero)

!!! Computes the ingredients for the factor g
k1 = n1
k2 = n4
k3 = n5
k4 = n4 - n3
k5 = n5 - n2

l1 = max(0,k4,k5)
l2 = min(k1,k2,k3)

h  = one
hl = one

do l = l1+1, l2
  hlm1 = hl
  hl = hlm1 * (l - k1) * (l - k2) * (l - k3) / ((l - k4) * (l - k5) * l)
  h = h + hl
enddo

k1 = k1 - l1
k2 = k2 - l1
k3 = k3 - l1
k4 = l1 + 1 - k4
k5 = l1 + 1 - k5
k6 = l1 + 1

q =  log_gamma(k1+zero) + log_gamma(k2+zero) + log_gamma(k3+zero) &
   + log_gamma(k4+zero) + log_gamma(k5+zero) + log_gamma(k6+zero)

!!! Computes the final value combining the two parts.
cg = dsqrt(j3 + one) * (-1)**l1 * dexp(0.5d0*p - q) * h

end subroutine ClebschGordan

!------------------------------------------------------------------------------!
! subroutine Wigner6JCoeff                                                     !
!                                                                              !
! Computes the Wigner6JCoeff from Racah coeffcient. The algorithm used was     !
! taken from technical notes from NASA written by W. F. Ford and R. C. Bruley. !
! Ref: NASA TN D-6173                                                          !
!                                                                              !
!------------------------------------------------------------------------------!
! subroutine Wigner6JCoeff                                                     !
!                                                                              !
! Computes the Wigner6JCoeff for the group SU(2). The algorithm used was taken !
! from technical notes from NASA written by W. F. Ford and R. C. Bruley.       !
! Ref: NASA TN D-6173                                                          !
!                                                                              !
! {  a  b  c  }                                                                !
! {  d  e  f  } = (-)^[a + b + d + e] W(a,b,e,d  ; c,f)                        !
!                                                                              !
! Be aware that all input values of j and m are supposed to be twice their     !
! values (such that we can treat easily half-integer angular momenta).         !
!------------------------------------------------------------------------------!
subroutine Wigner6JCoeff (a, b, c, d, e, f, c6j)

integer, intent(in) :: a, b, c, d, e, f
real(r64), intent(out) :: c6j
real(r64) :: aux

call RacahCoeff (a,b,e,d, c,f, aux)
c6j = aux
if (MOD((a + b + d + e) / 2, 2).eq.1) then
  c6j = - c6j
endif

return

end subroutine Wigner6JCoeff

!------------------------------------------------------------------------------!
!  subroutine RacahCoeff                                                       !
!                                                                              !
!  Racah coefficient, used for the six j coefficient with the same variable    !
!  naming of the reference.                                                    !
!------------------------------------------------------------------------------!
subroutine RacahCoeff (a,b,c,d, e,f, c6j)

integer, intent(in) :: a, b, c, d, e, f
real(r64), intent(out) :: c6j
integer :: j1, j2, j3, j4, j5, j6, j7, j8, k1, k2, k3, k4, l1, l2, l3, l4, &
           m1, m2, m3, m4, l, n_min, n_max
real(r64) :: p, q, h, hl, hlm1

c6j = zero
if (MOD(a+b+e, 2).ne.0) return
if (MOD(a+c+f, 2).ne.0) return
if (MOD(d+b+f, 2).ne.0) return
if (MOD(d+c+e, 2).ne.0) return

!!! Computes the ingredients for the factor c
j1 = 1 + (a + b - e)/2
j2 = 1 + (a + c - f)/2
j3 = 1 + (b + d - f)/2
j4 = 1 + (c + d - e)/2
if ( (min(j1,j2,j3,j4) < 1) ) return
k1 = 1 + (b + e - a)/2
k2 = 1 + (c + f - a)/2
k3 = 1 + (d + f - b)/2
k4 = 1 + (d + e - c)/2
if ( (min(k1,k2,k3,k4) < 1) ) return
l1 = 1 + (e + a - b)/2
l2 = 1 + (f + a - c)/2
l3 = 1 + (f + b - d)/2
l4 = 1 + (e + c - d)/2
if ( (min(l1,l2,l3,l4) < 1) ) return

m1 = j1 + k1 + l1 - 1
m2 = j2 + k2 + l2 - 1
m3 = j3 + k3 + l3 - 1
m4 = j4 + k4 + l4 - 1

p =  log_gamma(j1+zero) + log_gamma(k1+zero) + log_gamma(l1+zero) &
   + log_gamma(j2+zero) + log_gamma(k2+zero) + log_gamma(l2+zero) &
   + log_gamma(j3+zero) + log_gamma(k3+zero) + log_gamma(l3+zero) &
   + log_gamma(j4+zero) + log_gamma(k4+zero) + log_gamma(l4+zero) &
   - log_gamma(m1+zero) - log_gamma(m2+zero) - log_gamma(m3+zero) &
   - log_gamma(m4+zero)

!!! Computes the ingredients for the factor g
j5 = j1 - l2
j6 = j1 - l3
j7 = j1 + m4 - 1

n_min = max( 0,j5,j6)
n_max = min(j1,j2,j3,j4) - 1

h  = one
hl = one

if (n_min.ne.n_max) then
  do l = n_min+1, n_max
    hlm1 = hl
    hl = (l - j1) * (l - j2) * (l - j3) * (l - j4)
    hl = hl / ((l - j5) * (l - j6) * (l - j7) * l)
    hl = hlm1 * hl
    h = h + hl
  enddo
endif

j1 = j1 - n_min
j2 = j2 - n_min
j3 = j3 - n_min
j4 = j4 - n_min
j5 = n_min + 1 - j5
j6 = n_min + 1 - j6
j8 = j7 - n_min
j7 = n_min + 1

q =  log_gamma(j1+zero) + log_gamma(j2+zero) + log_gamma(j3+zero) &
   + log_gamma(j4+zero) + log_gamma(j5+zero) + log_gamma(j6+zero) &
   + log_gamma(j7+zero) - log_gamma(j8+zero)

!!! Computes the final value combining the two parts.
c6j = ((-1)**n_min) * dexp((0.5d0*p) - q) * h

end subroutine RacahCoeff

!------------------------------------------------------------------------------!
!  subroutine Wigner9JCoeff                                                    !
!                                                                              !
!    Computing 9j coefficient from Racah symbols, it validates the lsj values  !
!                     ___                                                      !
!  { l1   s1  j1 }    \                                                        !
!  { l2   s2  j2 } =  /__    (2t+1) * W(l1,l2,J,S; L,t) * W(l1,s1,J,j2 ; j1,t) !
!  {  L    S   J }    (t)                               * W(s1,S,j2,l2 ; s2,t) !
!                                                                              !
! Be aware that all input values of j and m are supposed to be twice their     !
! values (such that we can treat easily half-integer angular momenta).         !
!------------------------------------------------------------------------------!
subroutine Wigner9JCoeff (l1,s1,j1, l2,s2,j2, L,S,J, c9j)

integer, intent(in)    :: l1,s1,j1, l2,s2,j2, L,S,J
real(r64), intent(out) :: c9j
integer   :: t, t_max, t_min
real(r64) :: aux1, aux2, aux3

c9j = zero
if ((abs(l1-s1) .GT. j1).OR.(l1+s1 .LT. j1)) return
if ((abs(l2-s2) .GT. j2).OR.(l2+s2 .LT. j2)) return
if ((abs( L- S) .GT.  J).OR.(L + S .LT.  J)) return
if ((abs(l1-l2) .GT.  L).OR.(l1+l2 .LT.  L)) return
if ((abs(s1-s2) .GT.  S).OR.(s1+s2 .LT.  S)) return
if ((abs(j1-j2) .GT.  J).OR.(j1+j2 .LT.  J)) return

t_min = max(0, abs(l2-S), abs(s1-j2), abs(l1-J))
t_max = min(l2+S, s1+j2, l1+J)
t_min = t_min / 2
t_max = t_max / 2

do t = t_min, t_max
  call RacahCoeff (l1,l2, J,S, L, 2*t, aux1)
  if (abs(aux1) .LE. 1e-9) cycle
  call RacahCoeff (l1,s1,J,j2, j1,2*t, aux2)
  if (abs(aux1) .LE. 1e-9) cycle
  call RacahCoeff (s1,S,j2,l2, s2,2*t, aux3)

  c9j = c9j  + ((2*t + 1) * aux1 * aux2 * aux3)
enddo

return

end subroutine Wigner9JCoeff

!------------------------------------------------------------------------------!
! subroutine jacobi_srt  and jacobi(non sorted)                                !
!                                                                              !
! Implemented subroutine from Fortran 77 - Numerical recipes for symmetrical   !
! eigenvalue problem, to obtain the eigenvalues and eigenvectors sorted in in- !
! increasing order.                                                            !
!------------------------------------------------------------------------------!
subroutine jacobi(A, D, V, n, ndim)

integer, intent(in) :: n, ndim
!real(r64), dimension(ndim,ndim) :: Ainp
real(r64), dimension(ndim,ndim) :: V, A
real(r64), dimension(ndim)      :: D

integer   :: nrot, NMAX=500, i, ip, iq, j
real(r64) :: c, g, h, s, sm, t, tau, theta, tresh, TOL=1.0d-9
real(r64), dimension(:), allocatable :: b, z

allocate(b(NMAX), z(NMAX))

do ip = 1,n    !! Initialize to the identity matrix.
  do iq = 1,n
    v(ip,iq) = 0.0d0
  enddo
  v(ip,ip)=1.0d0
enddo

do ip = 1,n
  b(ip) = A(ip,ip)         !! Initialize b and d to the diagonal of a.
  D(ip) = b(ip)
  z(ip) = 0.0d0
enddo

nrot = 0

do i = 1, 50
  sm = 0.0d0
  do ip = 1,n-1            !! Sum off-diagonal elements.
    do iq = ip+1,n
      sm = sm + dabs(A(ip,iq))
    enddo
  enddo
  if (dabs(sm).LT.TOL) return
  if (i.LT.4) then
    tresh = 0.2d0 * sm / (n**2)
  else
    tresh = 0.0d0
  endif

  do ip = 1, n-1
    do iq = ip+1, n
      g = 100.0d0 * dabs(A(ip,iq))
      !! After 4 sweeps, skip the rotation if the off-diagonal element is small
      if ((i.GT.4) &
          .AND. (dabs( dabs(D(ip))+g - dabs(D(ip)) ).LT.TOL) &
          .AND. (dabs( dabs(D(iq))+g - dabs(D(iq)) ).LT.TOL)) then
          A(ip, iq) = 0.0d0
      elseif (dabs(A(ip,iq)) .GT. tresh) then
        h=d(iq)-d(ip)
        if(dabs( dabs(h)+g - dabs(h)) .LT. TOL)then
          t = A(ip,iq) / h
        else
          theta = 0.5d0 * h / A(ip,iq)
          t = 1.0d0 / (dabs(theta) + dsqrt(1.0d0 + theta**2))
          if (theta .LT. 0.0d0) t = -t
        endif

        c = 1.0d0 / dsqrt(1.0d0 + t**2)
        s = t*c
        tau = s / (1.0d0 + c)
        h = t * A(ip,iq)
        z(ip) = z(ip) - h
        z(iq) = z(iq) + h
        D(ip) = D(ip) - h
        D(iq) = D(iq) + h
        A(ip,iq) = 0.0d0

        do j = 1,ip-1                 !! Case of rotations 1 <= j < p
          g = A( j,ip)
          h = A( j,iq)
          A( j,ip) = g - (s*(h + g*tau))
          A( j,iq) = h + (s*(g - h*tau))
        enddo
        do j = ip+1,iq-1              !! Case of rotations p <  j < q
          g = A(ip, j)
          h = A( j,iq)
          A(ip, j) = g - (s*(h + g*tau))
          A( j,iq) = h + (s*(g - h*tau))
        enddo
        do j = iq+1,n                 !! Case of rotations q < j <= n
          g = A(ip, j)
          h = A(iq, j)
          A(ip, j) = g - (s*(h + g*tau))
          A(iq, j) = h + (s*(g - h*tau))
        enddo

        do j = 1,n
          g = V(j,ip)
          h = V(j,iq)
          V(j,ip) = g - (s*(h + g*tau))
          V(j,iq) = h + (s*(g - h*tau))
        enddo
        nrot = nrot + 1

      endif
    enddo
  enddo

  do ip = 1,n
    b(ip) = b(ip) + z(ip)
    D(ip) = b(ip)                !! Update d with the sum of t*a_pq
    z(ip) = 0.0d0                !! and reinitialize z
  enddo

enddo

print "(A,I6)", "[WARNING] Jacobi eigen. solver had too many iters 50: nrot=",&
                nrot

end subroutine jacobi




subroutine jacobi_srt(A, D, V, n, ndim)

integer, intent(in) :: n, ndim
real(r64), dimension(ndim,ndim) :: A, V
real(r64), dimension(ndim)      :: D
integer :: i,j,k
real(r64) :: p

call jacobi(A, D, V, n, ndim)

do i = 1,n-1
  k = i
  p = D(i)
  do j = i+1,n
    if (D(j) .LT. p) then
      k = j
      p = D(j)
    endif
  enddo

  if (k .NE. i) then
    D(k) = D(i)
    D(i) = p
    do j = 1,n
      p = V(j,i)
      V(j,i) = V(j,k)
      V(j,k) = p
    enddo
  endif
enddo

return
end subroutine jacobi_srt

!------------------------------------------------------------------------------!
! almost_equal (a, b, tolerance)        |a - b| < TOL                          !
!------------------------------------------------------------------------------!
function almost_equal(a, b, tolerance) result (bool_true)

real(r64), intent(in) :: a, b, tolerance
logical :: bool_true

bool_true = .FALSE.
if (dabs(a - b) .LT. tolerance) then
  bool_true = .TRUE.
endif

return
end function almost_equal


!------------------------------------------------------------------------------!
! function lapack_sel                                                          !
!                                                                              !
! Dummy logical function "select" for LAPACK (used for Schur decomposition).   !
! See LAPACK documentation for dgees.                                          !
!------------------------------------------------------------------------------!
function lapack_sel(wr,wi) result(selec)

logical :: selec
real(r64), intent(in) :: wr, wi
real(r64) :: dummy

!!! Just to remove the "unused-dummy-argument" during compilation
dummy = wr + wi

!!! Set to false because the function is only used as dummy argument
selec=.false.

end function lapack_sel

END MODULE MathMethods
!==============================================================================!
! End of file                                                                  !
!==============================================================================!

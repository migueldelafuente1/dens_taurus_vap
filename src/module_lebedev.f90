       !>      This is part of a set of subroutines that generate
       !!      Lebedev grids [1-6] for integration on a sphere. The original
       !!      C-code [1] was kindly provided by Dr. Dmitri N. Laikov and
       !!      translated into fortran by Dr. Christoph van Wuellen.
       !!      This subroutine was translated from C to fortran77 by hand.
       !!
       !!      Users of this code are asked to include reference [1] in their
       !!      publications, and in the user- and programmers-manuals
       !!      describing their codes.
       !!
       !!      This code was distributed through CCL (http://www.ccl.net/).
       !!
       !!      [1] V.I. Lebedev, and D.N. Laikov
       !!          "A quadrature formula for the sphere of the 131st
       !!           algebraic order of accuracy"
       !!          Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
       !!
       !!      [2] V.I. Lebedev
       !!          "A quadrature formula for the sphere of 59th algebraic
       !!           order of accuracy"
       !!          Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286.
       !!
       !!      [3] V.I. Lebedev, and A.L. Skorokhodov
       !!          "Quadrature formulas of orders 41, 47, and 53 for the sphere"
       !!          Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592.
       !!
       !!      [4] V.I. Lebedev
       !!          "Spherical quadrature formulas exact to orders 25-29"
       !!          Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107.
       !!
       !!      [5] V.I. Lebedev
       !!          "Quadratures on a sphere"
       !!          Computational Mathematics and Mathematical Physics, Vol. 16,
       !!          1976, pp. 10-24.
       !!
       !!      [6] V.I. Lebedev
       !!          "Values of the nodes and weights of ninth to seventeenth
       !!           order Gauss-Markov quadrature formulae invariant under the
       !!           octahedron group with inversion"
       !!          Computational Mathematics and Mathematical Physics, Vol. 15,
       !!          1975, pp. 44-51.
       !!
       module Lebedev
         implicit none

         public :: LLgrid, nPoints, MaxGrid
         private

         integer, parameter :: dp = kind(0.0d0)

         integer, parameter :: MaxGrid = 32
         integer, parameter :: nPoints(MaxGrid) = (/6,14,26,38,50,74,86,110,&
             & 146,170,194,230,266,302,350,434,590,770,974,1202,1454,1730,2030,&
             & 2354,2702,3074, 3470,3890,4334,4802,5294,5810/) !< list of grid sizes

       contains

         !> Routine to return the Lebedev-Laikov grids
         !! \param points Array of cartesian points, must be at least size (/3,
         !! nPoints(order)/)
         !! \param weights for the points
         !! \param order of the grid, between 1 and 32, if absent the largest order
         !!  that will fit into the supplied points() array will be used, but to the
         !!   maximum of 32
         subroutine LLgrid(points,weights,order)
           real(dp), intent(out)         :: points(:,:)
           real(dp), intent(out)         :: weights(:)
           integer, intent(in), optional :: order

           real(dp), allocatable :: x(:),y(:),z(:)
           integer :: iOrder, n

           if (present(order)) then
             iOrder = order
           else
             iOrder = 1
             do while(nPoints(iOrder)<size(points,dim=2) &
                      .and.iOrder<=size(nPoints))
               iOrder = iOrder + 1
             end do
           end if
           if (iOrder < 1 .or. iOrder > size(nPoints)) then
             write(*,*)"Requested Lebedev-Laikov grid not availible"
             stop
           end if

           ALLOCATE(x(nPoints(iOrder)))
           ALLOCATE(y(nPoints(iOrder)))
           ALLOCATE(z(nPoints(iOrder)))

           select case(nPoints(iOrder))
           case(6)
             call LD0006(x,y,z,weights,n)
           case(14)
             call LD0014(x,y,z,weights,n)
           case(26)
             call LD0026(x,y,z,weights,n)
           case(38)
             call LD0038(x,y,z,weights,n)
           case(50)
             call LD0050(x,y,z,weights,n)
           case(74)
             call LD0074(x,y,z,weights,n)
           case(86)
             call LD0086(x,y,z,weights,n)
           case(110)
             call LD0110(x,y,z,weights,n)
           case(146)
             call LD0146(x,y,z,weights,n)
           case(170)
             call LD0170(x,y,z,weights,n)
           case(194)
             call LD0194(x,y,z,weights,n)
           case(230)
             call LD0230(x,y,z,weights,n)
           case(266)
             call LD0266(x,y,z,weights,n)
           case(302)
             call LD0302(x,y,z,weights,n)
           case(350)
             call LD0350(x,y,z,weights,n)
           case(434)
             call LD0434(x,y,z,weights,n)
           case(590)
             call LD0590(x,y,z,weights,n)
           case(770)
             call LD0770(x,y,z,weights,n)
           case(974)
             call LD0974(x,y,z,weights,n)
           case(1202)
             call LD1202(x,y,z,weights,n)
           case(1454)
             call LD1454(x,y,z,weights,n)
           case(1730)
             call LD1730(x,y,z,weights,n)
           case(2030)
             call LD2030(x,y,z,weights,n)
           case(2354)
             call LD2354(x,y,z,weights,n)
           case(2702)
             call LD2702(x,y,z,weights,n)
           case(3074)
             call LD3074(x,y,z,weights,n)
           case(3470)
             call LD3470(x,y,z,weights,n)
           case(3890)
             call LD3890(x,y,z,weights,n)
           case(4334)
             call LD4334(x,y,z,weights,n)
           case(4802)
             call LD4802(x,y,z,weights,n)
           case(5294)
             call LD5294(x,y,z,weights,n)
           case(5810)
             call LD5810(x,y,z,weights,n)
           end select

           if (n/=nPoints(iOrder)) then
             write(*,*)"Failure in Lebedev-Laikov grid"
             stop
           end if

           points(1,1:nPoints(iOrder)) = x(:)
           points(2,1:nPoints(iOrder)) = y(:)
           points(3,1:nPoints(iOrder)) = z(:)

           DEALLOCATE(x)
           DEALLOCATE(y)
           DEALLOCATE(z)

         end subroutine LLgrid

         subroutine gen_oh(code, num, x, y, z, w, a, b, v)
           real(dp), intent(out) :: x(*),y(*),z(*),w(*)
           real(dp), intent(inout)  :: a, b
           real(dp), intent(in)    :: v
           integer, intent(in)   :: code
           integer, intent(inout)  :: num

           real(dp) :: c

           !    Given a point on a sphere (specified by a and b), generate all
           !    the equivalent points under Oh symmetry, making grid points with
           !    weight v.
           !    The variable num is increased by the number of different points
           !    generated.
           !
           !    Depending on code, there are 6...48 different but equivalent
           !    points.
           !
           !    code=1:   (0,0,1) etc                                (  6 points)
           !    code=2:   (0,a,a) etc, a=1/sqrt(2)                   ( 12 points)
           !    code=3:   (a,a,a) etc, a=1/sqrt(3)                   (  8 points)
           !    code=4:   (a,a,b) etc, b=sqrt(1-2 a^2)               ( 24 points)
           !    code=5:   (a,b,0) etc, b=sqrt(1-a^2), a input        ( 24 points)
           !    code=6:   (a,b,c) etc, c=sqrt(1-a^2-b^2), a/b input  ( 48 points)
           !

           select case(code)
           case(1)
             a=1.0_dp
             x(1) =  a
             y(1) =  0.0_dp
             z(1) =  0.0_dp
             w(1) =  v
             x(2) = -a
             y(2) =  0.0_dp
             z(2) =  0.0_dp
             w(2) =  v
             x(3) =  0.0_dp
             y(3) =  a
             z(3) =  0.0_dp
             w(3) =  v
             x(4) =  0.0_dp
             y(4) = -a
             z(4) =  0.0_dp
             w(4) =  v
             x(5) =  0.0_dp
             y(5) =  0.0_dp
             z(5) =  a
             w(5) =  v
             x(6) =  0.0_dp
             y(6) =  0.0_dp
             z(6) = -a
             w(6) =  v
             num=num+6
           case(2)
             a=sqrt(0.5_dp)
             x( 1) =  0.0_dp
             y( 1) =  a
             z( 1) =  a
             w( 1) =  v
             x( 2) =  0.0_dp
             y( 2) = -a
             z( 2) =  a
             w( 2) =  v
             x( 3) =  0.0_dp
             y( 3) =  a
             z( 3) = -a
             w( 3) =  v
             x( 4) =  0.0_dp
             y( 4) = -a
             z( 4) = -a
             w( 4) =  v
             x( 5) =  a
             y( 5) =  0.0_dp
             z( 5) =  a
             w( 5) =  v
             x( 6) = -a
             y( 6) =  0.0_dp
             z( 6) =  a
             w( 6) =  v
             x( 7) =  a
             y( 7) =  0.0_dp
             z( 7) = -a
             w( 7) =  v
             x( 8) = -a
             y( 8) =  0.0_dp
             z( 8) = -a
             w( 8) =  v
             x( 9) =  a
             y( 9) =  a
             z( 9) =  0.0_dp
             w( 9) =  v
             x(10) = -a
             y(10) =  a
             z(10) =  0.0_dp
             w(10) =  v
             x(11) =  a
             y(11) = -a
             z(11) =  0.0_dp
             w(11) =  v
             x(12) = -a
             y(12) = -a
             z(12) =  0.0_dp
             w(12) =  v
             num=num+12
           case(3)
             a = sqrt(1.0_dp/3.0_dp)
             x(1) =  a
             y(1) =  a
             z(1) =  a
             w(1) =  v
             x(2) = -a
             y(2) =  a
             z(2) =  a
             w(2) =  v
             x(3) =  a
             y(3) = -a
             z(3) =  a
             w(3) =  v
             x(4) = -a
             y(4) = -a
             z(4) =  a
             w(4) =  v
             x(5) =  a
             y(5) =  a
             z(5) = -a
             w(5) =  v
             x(6) = -a
             y(6) =  a
             z(6) = -a
             w(6) =  v
             x(7) =  a
             y(7) = -a
             z(7) = -a
             w(7) =  v
             x(8) = -a
             y(8) = -a
             z(8) = -a
             w(8) =  v
             num=num+8
           case(4)
             b = sqrt(1.0_dp - 2.0_dp*a*a)
             x( 1) =  a
             y( 1) =  a
             z( 1) =  b
             w( 1) =  v
             x( 2) = -a
             y( 2) =  a
             z( 2) =  b
             w( 2) =  v
             x( 3) =  a
             y( 3) = -a
             z( 3) =  b
             w( 3) =  v
             x( 4) = -a
             y( 4) = -a
             z( 4) =  b
             w( 4) =  v
             x( 5) =  a
             y( 5) =  a
             z( 5) = -b
             w( 5) =  v
             x( 6) = -a
             y( 6) =  a
             z( 6) = -b
             w( 6) =  v
             x( 7) =  a
             y( 7) = -a
             z( 7) = -b
             w( 7) =  v
             x( 8) = -a
             y( 8) = -a
             z( 8) = -b
             w( 8) =  v
             x( 9) =  a
             y( 9) =  b
             z( 9) =  a
             w( 9) =  v
             x(10) = -a
             y(10) =  b
             z(10) =  a
             w(10) =  v
             x(11) =  a
             y(11) = -b
             z(11) =  a
             w(11) =  v
             x(12) = -a
             y(12) = -b
             z(12) =  a
             w(12) =  v
             x(13) =  a
             y(13) =  b
             z(13) = -a
             w(13) =  v
             x(14) = -a
             y(14) =  b
             z(14) = -a
             w(14) =  v
             x(15) =  a
             y(15) = -b
             z(15) = -a
             w(15) =  v
             x(16) = -a
             y(16) = -b
             z(16) = -a
             w(16) =  v
             x(17) =  b
             y(17) =  a
             z(17) =  a
             w(17) =  v
             x(18) = -b
             y(18) =  a
             z(18) =  a
             w(18) =  v
             x(19) =  b
             y(19) = -a
             z(19) =  a
             w(19) =  v
             x(20) = -b
             y(20) = -a
             z(20) =  a
             w(20) =  v
             x(21) =  b
             y(21) =  a
             z(21) = -a
             w(21) =  v
             x(22) = -b
             y(22) =  a
             z(22) = -a
             w(22) =  v
             x(23) =  b
             y(23) = -a
             z(23) = -a
             w(23) =  v
             x(24) = -b
             y(24) = -a
             z(24) = -a
             w(24) =  v
             num=num+24
           case(5)
             b=sqrt(1.0_dp-a*a)
             x( 1) =  a
             y( 1) =  b
             z( 1) =  0.0_dp
             w( 1) =  v
             x( 2) = -a
             y( 2) =  b
             z( 2) =  0.0_dp
             w( 2) =  v
             x( 3) =  a
             y( 3) = -b
             z( 3) =  0.0_dp
             w( 3) =  v
             x( 4) = -a
             y( 4) = -b
             z( 4) =  0.0_dp
             w( 4) =  v
             x( 5) =  b
             y( 5) =  a
             z( 5) =  0.0_dp
             w( 5) =  v
             x( 6) = -b
             y( 6) =  a
             z( 6) =  0.0_dp
             w( 6) =  v
             x( 7) =  b
             y( 7) = -a
             z( 7) =  0.0_dp
             w( 7) =  v
             x( 8) = -b
             y( 8) = -a
             z( 8) =  0.0_dp
             w( 8) =  v
             x( 9) =  a
             y( 9) =  0.0_dp
             z( 9) =  b
             w( 9) =  v
             x(10) = -a
             y(10) =  0.0_dp
             z(10) =  b
             w(10) =  v
             x(11) =  a
             y(11) =  0.0_dp
             z(11) = -b
             w(11) =  v
             x(12) = -a
             y(12) =  0.0_dp
             z(12) = -b
             w(12) =  v
             x(13) =  b
             y(13) =  0.0_dp
             z(13) =  a
             w(13) =  v
             x(14) = -b
             y(14) =  0.0_dp
             z(14) =  a
             w(14) =  v
             x(15) =  b
             y(15) =  0.0_dp
             z(15) = -a
             w(15) =  v
             x(16) = -b
             y(16) =  0.0_dp
             z(16) = -a
             w(16) =  v
             x(17) =  0.0_dp
             y(17) =  a
             z(17) =  b
             w(17) =  v
             x(18) =  0.0_dp
             y(18) = -a
             z(18) =  b
             w(18) =  v
             x(19) =  0.0_dp
             y(19) =  a
             z(19) = -b
             w(19) =  v
             x(20) =  0.0_dp
             y(20) = -a
             z(20) = -b
             w(20) =  v
             x(21) =  0.0_dp
             y(21) =  b
             z(21) =  a
             w(21) =  v
             x(22) =  0.0_dp
             y(22) = -b
             z(22) =  a
             w(22) =  v
             x(23) =  0.0_dp
             y(23) =  b
             z(23) = -a
             w(23) =  v
             x(24) =  0.0_dp
             y(24) = -b
             z(24) = -a
             w(24) =  v
             num=num+24
           case(6)
             c=sqrt(1.0_dp - a*a - b*b)
             x( 1) =  a
             y( 1) =  b
             z( 1) =  c
             w( 1) =  v
             x( 2) = -a
             y( 2) =  b
             z( 2) =  c
             w( 2) =  v
             x( 3) =  a
             y( 3) = -b
             z( 3) =  c
             w( 3) =  v
             x( 4) = -a
             y( 4) = -b
             z( 4) =  c
             w( 4) =  v
             x( 5) =  a
             y( 5) =  b
             z( 5) = -c
             w( 5) =  v
             x( 6) = -a
             y( 6) =  b
             z( 6) = -c
             w( 6) =  v
             x( 7) =  a
             y( 7) = -b
             z( 7) = -c
             w( 7) =  v
             x( 8) = -a
             y( 8) = -b
             z( 8) = -c
             w( 8) =  v
             x( 9) =  a
             y( 9) =  c
             z( 9) =  b
             w( 9) =  v
             x(10) = -a
             y(10) =  c
             z(10) =  b
             w(10) =  v
             x(11) =  a
             y(11) = -c
             z(11) =  b
             w(11) =  v
             x(12) = -a
             y(12) = -c
             z(12) =  b
             w(12) =  v
             x(13) =  a
             y(13) =  c
             z(13) = -b
             w(13) =  v
             x(14) = -a
             y(14) =  c
             z(14) = -b
             w(14) =  v
             x(15) =  a
             y(15) = -c
             z(15) = -b
             w(15) =  v
             x(16) = -a
             y(16) = -c
             z(16) = -b
             w(16) =  v
             x(17) =  b
             y(17) =  a
             z(17) =  c
             w(17) =  v
             x(18) = -b
             y(18) =  a
             z(18) =  c
             w(18) =  v
             x(19) =  b
             y(19) = -a
             z(19) =  c
             w(19) =  v
             x(20) = -b
             y(20) = -a
             z(20) =  c
             w(20) =  v
             x(21) =  b
             y(21) =  a
             z(21) = -c
             w(21) =  v
             x(22) = -b
             y(22) =  a
             z(22) = -c
             w(22) =  v
             x(23) =  b
             y(23) = -a
             z(23) = -c
             w(23) =  v
             x(24) = -b
             y(24) = -a
             z(24) = -c
             w(24) =  v
             x(25) =  b
             y(25) =  c
             z(25) =  a
             w(25) =  v
             x(26) = -b
             y(26) =  c
             z(26) =  a
             w(26) =  v
             x(27) =  b
             y(27) = -c
             z(27) =  a
             w(27) =  v
             x(28) = -b
             y(28) = -c
             z(28) =  a
             w(28) =  v
             x(29) =  b
             y(29) =  c
             z(29) = -a
             w(29) =  v
             x(30) = -b
             y(30) =  c
             z(30) = -a
             w(30) =  v
             x(31) =  b
             y(31) = -c
             z(31) = -a
             w(31) =  v
             x(32) = -b
             y(32) = -c
             z(32) = -a
             w(32) =  v
             x(33) =  c
             y(33) =  a
             z(33) =  b
             w(33) =  v
             x(34) = -c
             y(34) =  a
             z(34) =  b
             w(34) =  v
             x(35) =  c
             y(35) = -a
             z(35) =  b
             w(35) =  v
             x(36) = -c
             y(36) = -a
             z(36) =  b
             w(36) =  v
             x(37) =  c
             y(37) =  a
             z(37) = -b
             w(37) =  v
             x(38) = -c
             y(38) =  a
             z(38) = -b
             w(38) =  v
             x(39) =  c
             y(39) = -a
             z(39) = -b
             w(39) =  v
             x(40) = -c
             y(40) = -a
             z(40) = -b
             w(40) =  v
             x(41) =  c
             y(41) =  b
             z(41) =  a
             w(41) =  v
             x(42) = -c
             y(42) =  b
             z(42) =  a
             w(42) =  v
             x(43) =  c
             y(43) = -b
             z(43) =  a
             w(43) =  v
             x(44) = -c
             y(44) = -b
             z(44) =  a
             w(44) =  v
             x(45) =  c
             y(45) =  b
             z(45) = -a
             w(45) =  v
             x(46) = -c
             y(46) =  b
             z(46) = -a
             w(46) =  v
             x(47) =  c
             y(47) = -b
             z(47) = -a
             w(47) =  v
             x(48) = -c
             y(48) = -b
             z(48) = -a
             w(48) =  v
             num=num+48
           case default
             write(*,*)"Gen_Oh: Invalid Code"
             stop
           end select
         end subroutine gen_oh

         subroutine LD0006(X,Y,Z,W,N)
           real(dp), intent(out) :: X(   6)
           real(dp), intent(out) :: Y(   6)
           real(dp), intent(out) :: Z(   6)
           real(dp), intent(out) :: W(   6)
           integer, intent(out) :: N

           real(dp) :: A,B,V

           N=1
           V=0.1666666666666667_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
           end subroutine LD0006

         subroutine LD0014(X,Y,Z,W,N)
           real(dp), intent(out) :: X(  14)
           real(dp), intent(out) :: Y(  14)
           real(dp), intent(out) :: Z(  14)
           real(dp), intent(out) :: W(  14)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.6666666666666667_dp*0.1_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.7500000000000000_dp*0.1_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
           end subroutine LD0014

         subroutine LD0026(X,Y,Z,W,N)
           real(dp), intent(out) :: X(  26)
           real(dp), intent(out) :: Y(  26)
           real(dp), intent(out) :: Z(  26)
           real(dp), intent(out) :: W(  26)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.4761904761904762_dp*0.1_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.3809523809523810_dp*0.1_dp
           Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.3214285714285714_dp*0.1_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
           end subroutine LD0026

         subroutine LD0038(X,Y,Z,W,N)
           real(dp), intent(out) :: X(  38)
           real(dp), intent(out) :: Y(  38)
           real(dp), intent(out) :: Z(  38)
           real(dp), intent(out) :: W(  38)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.9523809523809524_dp*0.01_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.3214285714285714_dp*0.1_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4597008433809831_dp
           V=0.2857142857142857_dp*0.1_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
           end subroutine LD0038

         subroutine LD0050(X,Y,Z,W,N)
           real(dp), intent(out) :: X(  50)
           real(dp), intent(out) :: Y(  50)
           real(dp), intent(out) :: Z(  50)
           real(dp), intent(out) :: W(  50)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.1269841269841270_dp*0.1_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.2257495590828924_dp*0.1_dp
           Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.2109375000000000_dp*0.1_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3015113445777636_dp
           V=0.2017333553791887_dp*0.1_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
           end subroutine LD0050

         subroutine LD0074(X,Y,Z,W,N)
           real(dp), intent(out) :: X(  74)
           real(dp), intent(out) :: Y(  74)
           real(dp), intent(out) :: Z(  74)
           real(dp), intent(out) :: W(  74)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.5130671797338464_dp*0.001_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.1660406956574204_dp*0.1_dp
           Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=-0.2958603896103896_dp*0.1_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4803844614152614_dp
           V=0.2657620708215946_dp*0.1_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3207726489807764_dp
           V=0.1652217099371571_dp*0.1_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
           end subroutine LD0074

         subroutine LD0086(X,Y,Z,W,N)
           real(dp), intent(out) :: X(  86)
           real(dp), intent(out) :: Y(  86)
           real(dp), intent(out) :: Z(  86)
           real(dp), intent(out) :: W(  86)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.1154401154401154_dp*0.1_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.1194390908585628_dp*0.1_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3696028464541502_dp
           V=0.1111055571060340_dp*0.1_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6943540066026664_dp
           V=0.1187650129453714_dp*0.1_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3742430390903412_dp
           V=0.1181230374690448_dp*0.1_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
           end subroutine LD0086

         subroutine LD0110(X,Y,Z,W,N)
           real(dp), intent(out) :: X( 110)
           real(dp), intent(out) :: Y( 110)
           real(dp), intent(out) :: Z( 110)
           real(dp), intent(out) :: W( 110)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.3828270494937162_dp*0.01_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.9793737512487512_dp*0.01_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1851156353447362_dp
           V=0.8211737283191111_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6904210483822922_dp
           V=0.9942814891178103_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3956894730559419_dp
           V=0.9595471336070963_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4783690288121502_dp
           V=0.9694996361663028_dp*0.01_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
           end subroutine LD0110

         subroutine LD0146(X,Y,Z,W,N)
           real(dp), intent(out) :: X( 146)
           real(dp), intent(out) :: Y( 146)
           real(dp), intent(out) :: Z( 146)
           real(dp), intent(out) :: W( 146)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.5996313688621381_dp*0.001_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.7372999718620756_dp*0.01_dp
           Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.7210515360144488_dp*0.01_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6764410400114264_dp
           V=0.7116355493117555_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4174961227965453_dp
           V=0.6753829486314477_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1574676672039082_dp
           V=0.7574394159054034_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1403553811713183_dp
           B=0.4493328323269557_dp
           V=0.6991087353303262_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
           end subroutine LD0146

         subroutine LD0170(X,Y,Z,W,N)
           real(dp), intent(out) :: X( 170)
           real(dp), intent(out) :: Y( 170)
           real(dp), intent(out) :: Z( 170)
           real(dp), intent(out) :: W( 170)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.5544842902037365_dp*0.01_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.6071332770670752_dp*0.01_dp
           Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.6383674773515093_dp*0.01_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2551252621114134_dp
           V=0.5183387587747790_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6743601460362766_dp
           V=0.6317929009813725_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4318910696719410_dp
           V=0.6201670006589077_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2613931360335988_dp
           V=0.5477143385137348_dp*0.01_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4990453161796037_dp
           B=0.1446630744325115_dp
           V=0.5968383987681156_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
           end subroutine LD0170

         subroutine LD0194(X,Y,Z,W,N)
           real(dp), intent(out) :: X( 194)
           real(dp), intent(out) :: Y( 194)
           real(dp), intent(out) :: Z( 194)
           real(dp), intent(out) :: W( 194)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.1782340447244611_dp*0.01_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.5716905949977102_dp*0.01_dp
           Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.5573383178848738_dp*0.01_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6712973442695226_dp
           V=0.5608704082587997_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2892465627575439_dp
           V=0.5158237711805383_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4446933178717437_dp
           V=0.5518771467273614_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1299335447650067_dp
           V=0.4106777028169394_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3457702197611283_dp
           V=0.5051846064614808_dp*0.01_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1590417105383530_dp
           B=0.8360360154824589_dp
           V=0.5530248916233094_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
           end subroutine LD0194

         subroutine LD0230(X,Y,Z,W,N)
           real(dp), intent(out) :: X( 230)
           real(dp), intent(out) :: Y( 230)
           real(dp), intent(out) :: Z( 230)
           real(dp), intent(out) :: W( 230)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=-0.5522639919727325_dp*0.1_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.4450274607445226_dp*0.01_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4492044687397611_dp
           V=0.4496841067921404_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2520419490210201_dp
           V=0.5049153450478750_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6981906658447242_dp
           V=0.3976408018051883_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6587405243460960_dp
           V=0.4401400650381014_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4038544050097660_dp*0.1_dp
           V=0.1724544350544401_dp*0.1_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5823842309715585_dp
           V=0.4231083095357343_dp*0.01_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3545877390518688_dp
           V=0.5198069864064399_dp*0.01_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2272181808998187_dp
           B=0.4864661535886647_dp
           V=0.4695720972568883_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
           end subroutine LD0230

         subroutine LD0266(X,Y,Z,W,N)
           real(dp), intent(out) :: X( 266)
           real(dp), intent(out) :: Y( 266)
           real(dp), intent(out) :: Z( 266)
           real(dp), intent(out) :: W( 266)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=-0.1313769127326952_dp*0.01_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=-0.2522728704859336_dp*0.01_dp
           Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.4186853881700583_dp*0.01_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7039373391585475_dp
           V=0.5315167977810885_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1012526248572414_dp
           V=0.4047142377086219_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4647448726420539_dp
           V=0.4112482394406990_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3277420654971629_dp
           V=0.3595584899758782_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6620338663699974_dp
           V=0.4256131351428158_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.8506508083520399_dp
           V=0.4229582700647240_dp*0.01_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3233484542692899_dp
           B=0.1153112011009701_dp
           V=0.4080914225780505_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2314790158712601_dp
           B=0.5244939240922365_dp
           V=0.4071467593830964_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
           end subroutine LD0266

         subroutine LD0302(X,Y,Z,W,N)
           real(dp), intent(out) :: X( 302)
           real(dp), intent(out) :: Y( 302)
           real(dp), intent(out) :: Z( 302)
           real(dp), intent(out) :: W( 302)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.8545911725128148_dp*0.001_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.3599119285025571_dp*0.01_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3515640345570105_dp
           V=0.3449788424305883_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6566329410219612_dp
           V=0.3604822601419882_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4729054132581005_dp
           V=0.3576729661743367_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.9618308522614784_dp*0.1_dp
           V=0.2352101413689164_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2219645236294178_dp
           V=0.3108953122413675_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7011766416089545_dp
           V=0.3650045807677255_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2644152887060663_dp
           V=0.2982344963171804_dp*0.01_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5718955891878961_dp
           V=0.3600820932216460_dp*0.01_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2510034751770465_dp
           B=0.8000727494073952_dp
           V=0.3571540554273387_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1233548532583327_dp
           B=0.4127724083168531_dp
           V=0.3392312205006170_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
           end subroutine LD0302

         subroutine LD0350(X,Y,Z,W,N)
           real(dp), intent(out) :: X( 350)
           real(dp), intent(out) :: Y( 350)
           real(dp), intent(out) :: Z( 350)
           real(dp), intent(out) :: W( 350)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.3006796749453936_dp*0.01_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.3050627745650771_dp*0.01_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7068965463912316_dp
           V=0.1621104600288991_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4794682625712025_dp
           V=0.3005701484901752_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1927533154878019_dp
           V=0.2990992529653774_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6930357961327123_dp
           V=0.2982170644107595_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3608302115520091_dp
           V=0.2721564237310992_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6498486161496169_dp
           V=0.3033513795811141_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1932945013230339_dp
           V=0.3007949555218533_dp*0.01_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3800494919899303_dp
           V=0.2881964603055307_dp*0.01_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2899558825499574_dp
           B=0.7934537856582316_dp
           V=0.2958357626535696_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.9684121455103957_dp*0.1_dp
           B=0.8280801506686862_dp
           V=0.3036020026407088_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1833434647041659_dp
           B=0.9074658265305127_dp
           V=0.2832187403926303_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
           end subroutine LD0350

         subroutine LD0434(X,Y,Z,W,N)
           real(dp), intent(out) :: X( 434)
           real(dp), intent(out) :: Y( 434)
           real(dp), intent(out) :: Z( 434)
           real(dp), intent(out) :: W( 434)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.5265897968224436_dp*0.001_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.2548219972002607_dp*0.01_dp
           Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.2512317418927307_dp*0.01_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6909346307509111_dp
           V=0.2530403801186355_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1774836054609158_dp
           V=0.2014279020918528_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4914342637784746_dp
           V=0.2501725168402936_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6456664707424256_dp
           V=0.2513267174597564_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2861289010307638_dp
           V=0.2302694782227416_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7568084367178018_dp*0.1_dp
           V=0.1462495621594614_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3927259763368002_dp
           V=0.2445373437312980_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.8818132877794288_dp
           V=0.2417442375638981_dp*0.01_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.9776428111182649_dp
           V=0.1910951282179532_dp*0.01_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2054823696403044_dp
           B=0.8689460322872412_dp
           V=0.2416930044324775_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5905157048925271_dp
           B=0.7999278543857286_dp
           V=0.2512236854563495_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5550152361076807_dp
           B=0.7717462626915901_dp
           V=0.2496644054553086_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.9371809858553722_dp
           B=0.3344363145343455_dp
           V=0.2236607760437849_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
           end subroutine LD0434

         subroutine LD0590(X,Y,Z,W,N)
           real(dp), intent(out) :: X( 590)
           real(dp), intent(out) :: Y( 590)
           real(dp), intent(out) :: Z( 590)
           real(dp), intent(out) :: W( 590)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.3095121295306187_dp*0.001_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.1852379698597489_dp*0.01_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7040954938227469_dp
           V=0.1871790639277744_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6807744066455243_dp
           V=0.1858812585438317_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6372546939258752_dp
           V=0.1852028828296213_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5044419707800358_dp
           V=0.1846715956151242_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4215761784010967_dp
           V=0.1818471778162769_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3317920736472123_dp
           V=0.1749564657281154_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2384736701421887_dp
           V=0.1617210647254411_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1459036449157763_dp
           V=0.1384737234851692_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6095034115507196_dp*0.1_dp
           V=0.9764331165051050_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6116843442009876_dp
           V=0.1857161196774078_dp*0.01_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3964755348199858_dp
           V=0.1705153996395864_dp*0.01_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1724782009907724_dp
           V=0.1300321685886048_dp*0.01_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5610263808622060_dp
           B=0.3518280927733519_dp
           V=0.1842866472905286_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4742392842551980_dp
           B=0.2634716655937950_dp
           V=0.1802658934377451_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5984126497885380_dp
           B=0.1816640840360209_dp
           V=0.1849830560443660_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3791035407695563_dp
           B=0.1720795225656878_dp
           V=0.1713904507106709_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2778673190586244_dp
           B=0.8213021581932511_dp*0.1_dp
           V=0.1555213603396808_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5033564271075117_dp
           B=0.8999205842074875_dp*0.1_dp
           V=0.1802239128008525_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
         end subroutine LD0590

         subroutine LD0770(X,Y,Z,W,N)
           real(dp), intent(out) :: X( 770)
           real(dp), intent(out) :: Y( 770)
           real(dp), intent(out) :: Z( 770)
           real(dp), intent(out) :: W( 770)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.2192942088181184_dp*0.001_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.1436433617319080_dp*0.01_dp
           Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.1421940344335877_dp*0.01_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5087204410502360_dp*0.1_dp
           V=0.6798123511050502_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1228198790178831_dp
           V=0.9913184235294912_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2026890814408786_dp
           V=0.1180207833238949_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2847745156464294_dp
           V=0.1296599602080921_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3656719078978026_dp
           V=0.1365871427428316_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4428264886713469_dp
           V=0.1402988604775325_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5140619627249735_dp
           V=0.1418645563595609_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6306401219166803_dp
           V=0.1421376741851662_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6716883332022612_dp
           V=0.1423996475490962_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6979792685336881_dp
           V=0.1431554042178567_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1446865674195309_dp
           V=0.9254401499865368_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3390263475411216_dp
           V=0.1250239995053509_dp*0.01_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5335804651263506_dp
           V=0.1394365843329230_dp*0.01_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6944024393349413_dp*0.1_dp
           B=0.2355187894242326_dp
           V=0.1127089094671749_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2269004109529460_dp
           B=0.4102182474045730_dp
           V=0.1345753760910670_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.8025574607775339_dp*0.1_dp
           B=0.6214302417481605_dp
           V=0.1424957283316783_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1467999527896572_dp
           B=0.3245284345717394_dp
           V=0.1261523341237750_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1571507769824727_dp
           B=0.5224482189696630_dp
           V=0.1392547106052696_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2365702993157246_dp
           B=0.6017546634089558_dp
           V=0.1418761677877656_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7714815866765732_dp*0.1_dp
           B=0.4346575516141163_dp
           V=0.1338366684479554_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3062936666210730_dp
           B=0.4908826589037616_dp
           V=0.1393700862676131_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3822477379524787_dp
           B=0.5648768149099500_dp
           V=0.1415914757466932_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
         end subroutine LD0770

         subroutine LD0974(X,Y,Z,W,N)
           real(dp), intent(out) :: X( 974)
           real(dp), intent(out) :: Y( 974)
           real(dp), intent(out) :: Z( 974)
           real(dp), intent(out) :: W( 974)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.1438294190527431_dp*0.001_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.1125772288287004_dp*0.01_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4292963545341347_dp*0.1_dp
           V=0.4948029341949241_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1051426854086404_dp
           V=0.7357990109125470_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1750024867623087_dp
           V=0.8889132771304384_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2477653379650257_dp
           V=0.9888347838921435_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3206567123955957_dp
           V=0.1053299681709471_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3916520749849983_dp
           V=0.1092778807014578_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4590825874187624_dp
           V=0.1114389394063227_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5214563888415861_dp
           V=0.1123724788051555_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6253170244654199_dp
           V=0.1125239325243814_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6637926744523170_dp
           V=0.1126153271815905_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6910410398498301_dp
           V=0.1130286931123841_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7052907007457760_dp
           V=0.1134986534363955_dp*0.01_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1236686762657990_dp
           V=0.6823367927109931_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2940777114468387_dp
           V=0.9454158160447096_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4697753849207649_dp
           V=0.1074429975385679_dp*0.01_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6334563241139567_dp
           V=0.1129300086569132_dp*0.01_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5974048614181342_dp*0.1_dp
           B=0.2029128752777523_dp
           V=0.8436884500901954_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1375760408473636_dp
           B=0.4602621942484054_dp
           V=0.1075255720448885_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3391016526336286_dp
           B=0.5030673999662036_dp
           V=0.1108577236864462_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1271675191439820_dp
           B=0.2817606422442134_dp
           V=0.9566475323783357_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2693120740413512_dp
           B=0.4331561291720157_dp
           V=0.1080663250717391_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1419786452601918_dp
           B=0.6256167358580814_dp
           V=0.1126797131196295_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6709284600738255_dp*0.1_dp
           B=0.3798395216859157_dp
           V=0.1022568715358061_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7057738183256172_dp*0.1_dp
           B=0.5517505421423520_dp
           V=0.1108960267713108_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2783888477882155_dp
           B=0.6029619156159187_dp
           V=0.1122790653435766_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1979578938917407_dp
           B=0.3589606329589096_dp
           V=0.1032401847117460_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2087307061103274_dp
           B=0.5348666438135476_dp
           V=0.1107249382283854_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4055122137872836_dp
           B=0.5674997546074373_dp
           V=0.1121780048519972_dp*0.01_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
         end subroutine LD0974

         subroutine LD1202(X,Y,Z,W,N)
           real(dp), intent(out) :: X(1202)
           real(dp), intent(out) :: Y(1202)
           real(dp), intent(out) :: Z(1202)
           real(dp), intent(out) :: W(1202)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.1105189233267572_dp*0.001_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.9205232738090741_dp*0.001_dp
           Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.9133159786443561_dp*0.001_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3712636449657089_dp*0.1_dp
           V=0.3690421898017899_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.9140060412262223_dp*0.1_dp
           V=0.5603990928680660_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1531077852469906_dp
           V=0.6865297629282609_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2180928891660612_dp
           V=0.7720338551145630_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2839874532200175_dp
           V=0.8301545958894795_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3491177600963764_dp
           V=0.8686692550179628_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4121431461444309_dp
           V=0.8927076285846890_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4718993627149127_dp
           V=0.9060820238568219_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5273145452842337_dp
           V=0.9119777254940867_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6209475332444019_dp
           V=0.9128720138604181_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6569722711857291_dp
           V=0.9130714935691735_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6841788309070143_dp
           V=0.9152873784554116_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7012604330123631_dp
           V=0.9187436274321654_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1072382215478166_dp
           V=0.5176977312965694_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2582068959496968_dp
           V=0.7331143682101417_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4172752955306717_dp
           V=0.8463232836379928_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5700366911792503_dp
           V=0.9031122694253992_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.9827986018263947_dp
           B=0.1771774022615325_dp
           V=0.6485778453163257_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.9624249230326228_dp
           B=0.2475716463426288_dp
           V=0.7435030910982369_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.9402007994128811_dp
           B=0.3354616289066489_dp
           V=0.7998527891839054_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.9320822040143202_dp
           B=0.3173615246611977_dp
           V=0.8101731497468018_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.9043674199393299_dp
           B=0.4090268427085357_dp
           V=0.8483389574594331_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.8912407560074747_dp
           B=0.3854291150669224_dp
           V=0.8556299257311812_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.8676435628462708_dp
           B=0.4932221184851285_dp
           V=0.8803208679738260_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.8581979986041619_dp
           B=0.4785320675922435_dp
           V=0.8811048182425720_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.8396753624049856_dp
           B=0.4507422593157064_dp
           V=0.8850282341265444_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.8165288564022188_dp
           B=0.5632123020762100_dp
           V=0.9021342299040653_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.8015469370783529_dp
           B=0.5434303569693900_dp
           V=0.9010091677105086_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7773563069070351_dp
           B=0.5123518486419871_dp
           V=0.9022692938426915_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7661621213900394_dp
           B=0.6394279634749102_dp
           V=0.9158016174693465_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7553584143533510_dp
           B=0.6269805509024392_dp
           V=0.9131578003189435_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7344305757559503_dp
           B=0.6031161693096310_dp
           V=0.9107813579482705_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7043837184021765_dp
           B=0.5693702498468441_dp
           V=0.9105760258970126_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
         end subroutine LD1202

         subroutine LD1454(X,Y,Z,W,N)
           real(dp), intent(out) :: X(1454)
           real(dp), intent(out) :: Y(1454)
           real(dp), intent(out) :: Z(1454)
           real(dp), intent(out) :: W(1454)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.7777160743261247_dp*0.0001_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.7557646413004701_dp*0.001_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3229290663413854_dp*0.1_dp
           V=0.2841633806090617_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.8036733271462222_dp*0.1_dp
           V=0.4374419127053555_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1354289960531653_dp
           V=0.5417174740872172_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1938963861114426_dp
           V=0.6148000891358593_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2537343715011275_dp
           V=0.6664394485800705_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3135251434752570_dp
           V=0.7025039356923220_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3721558339375338_dp
           V=0.7268511789249627_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4286809575195696_dp
           V=0.7422637534208629_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4822510128282994_dp
           V=0.7509545035841214_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5320679333566263_dp
           V=0.7548535057718401_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6172998195394274_dp
           V=0.7554088969774001_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6510679849127481_dp
           V=0.7553147174442808_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6777315251687360_dp
           V=0.7564767653292297_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6963109410648741_dp
           V=0.7587991808518730_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7058935009831749_dp
           V=0.7608261832033027_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.9955546194091857_dp
           V=0.4021680447874916_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.9734115901794209_dp
           V=0.5804871793945964_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.9275693732388626_dp
           V=0.6792151955945159_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.8568022422795103_dp
           V=0.7336741211286294_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7623495553719372_dp
           V=0.7581866300989608_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5707522908892223_dp
           B=0.4387028039889501_dp
           V=0.7538257859800743_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5196463388403083_dp
           B=0.3858908414762617_dp
           V=0.7483517247053123_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4646337531215351_dp
           B=0.3301937372343854_dp
           V=0.7371763661112059_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4063901697557691_dp
           B=0.2725423573563777_dp
           V=0.7183448895756934_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3456329466643087_dp
           B=0.2139510237495250_dp
           V=0.6895815529822191_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2831395121050332_dp
           B=0.1555922309786647_dp
           V=0.6480105801792886_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2197682022925330_dp
           B=0.9892878979686097_dp*0.1_dp
           V=0.5897558896594636_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1564696098650355_dp
           B=0.4598642910675510_dp*0.1_dp
           V=0.5095708849247346_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6027356673721295_dp
           B=0.3376625140173426_dp
           V=0.7536906428909755_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5496032320255096_dp
           B=0.2822301309727988_dp
           V=0.7472505965575118_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4921707755234567_dp
           B=0.2248632342592540_dp
           V=0.7343017132279698_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4309422998598483_dp
           B=0.1666224723456479_dp
           V=0.7130871582177445_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3664108182313672_dp
           B=0.1086964901822169_dp
           V=0.6817022032112776_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2990189057758436_dp
           B=0.5251989784120085_dp*0.1_dp
           V=0.6380941145604121_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6268724013144998_dp
           B=0.2297523657550023_dp
           V=0.7550381377920310_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5707324144834607_dp
           B=0.1723080607093800_dp
           V=0.7478646640144802_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5096360901960365_dp
           B=0.1140238465390513_dp
           V=0.7335918720601220_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4438729938312456_dp
           B=0.5611522095882537_dp*0.1_dp
           V=0.7110120527658118_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6419978471082389_dp
           B=0.1164174423140873_dp
           V=0.7571363978689501_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5817218061802611_dp
           B=0.5797589531445219_dp*0.1_dp
           V=0.7489908329079234_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
         end subroutine LD1454

         subroutine LD1730(X,Y,Z,W,N)
           real(dp), intent(out) :: X(1730)
           real(dp), intent(out) :: Y(1730)
           real(dp), intent(out) :: Z(1730)
           real(dp), intent(out) :: W(1730)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.6309049437420976_dp*0.0001_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.6398287705571748_dp*0.001_dp
           Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.6357185073530720_dp*0.001_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2860923126194662_dp*0.1_dp
           V=0.2221207162188168_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7142556767711522_dp*0.1_dp
           V=0.3475784022286848_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1209199540995559_dp
           V=0.4350742443589804_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1738673106594379_dp
           V=0.4978569136522127_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2284645438467734_dp
           V=0.5435036221998053_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2834807671701512_dp
           V=0.5765913388219542_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3379680145467339_dp
           V=0.6001200359226003_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3911355454819537_dp
           V=0.6162178172717512_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4422860353001403_dp
           V=0.6265218152438485_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4907781568726057_dp
           V=0.6323987160974212_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5360006153211468_dp
           V=0.6350767851540569_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6142105973596603_dp
           V=0.6354362775297107_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6459300387977504_dp
           V=0.6352302462706235_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6718056125089225_dp
           V=0.6358117881417972_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6910888533186254_dp
           V=0.6373101590310117_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7030467416823252_dp
           V=0.6390428961368665_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.8354951166354646_dp*0.1_dp
           V=0.3186913449946576_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2050143009099486_dp
           V=0.4678028558591711_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3370208290706637_dp
           V=0.5538829697598626_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4689051484233963_dp
           V=0.6044475907190476_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5939400424557334_dp
           V=0.6313575103509012_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1394983311832261_dp
           B=0.4097581162050343_dp*0.1_dp
           V=0.4078626431855630_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1967999180485014_dp
           B=0.8851987391293348_dp*0.1_dp
           V=0.4759933057812725_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2546183732548967_dp
           B=0.1397680182969819_dp
           V=0.5268151186413440_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3121281074713875_dp
           B=0.1929452542226526_dp
           V=0.5643048560507316_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3685981078502492_dp
           B=0.2467898337061562_dp
           V=0.5914501076613073_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4233760321547856_dp
           B=0.3003104124785409_dp
           V=0.6104561257874195_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4758671236059246_dp
           B=0.3526684328175033_dp
           V=0.6230252860707806_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5255178579796463_dp
           B=0.4031134861145713_dp
           V=0.6305618761760796_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5718025633734589_dp
           B=0.4509426448342351_dp
           V=0.6343092767597889_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2686927772723415_dp
           B=0.4711322502423248_dp*0.1_dp
           V=0.5176268945737826_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3306006819904809_dp
           B=0.9784487303942695_dp*0.1_dp
           V=0.5564840313313692_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3904906850594983_dp
           B=0.1505395810025273_dp
           V=0.5856426671038980_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4479957951904390_dp
           B=0.2039728156296050_dp
           V=0.6066386925777091_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5027076848919780_dp
           B=0.2571529941121107_dp
           V=0.6208824962234458_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5542087392260217_dp
           B=0.3092191375815670_dp
           V=0.6296314297822907_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6020850887375187_dp
           B=0.3593807506130276_dp
           V=0.6340423756791859_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4019851409179594_dp
           B=0.5063389934378671_dp*0.1_dp
           V=0.5829627677107342_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4635614567449800_dp
           B=0.1032422269160612_dp
           V=0.6048693376081110_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5215860931591575_dp
           B=0.1566322094006254_dp
           V=0.6202362317732461_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5758202499099271_dp
           B=0.2098082827491099_dp
           V=0.6299005328403779_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6259893683876795_dp
           B=0.2618824114553391_dp
           V=0.6347722390609353_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5313795124811891_dp
           B=0.5263245019338556_dp*0.1_dp
           V=0.6203778981238834_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5893317955931995_dp
           B=0.1061059730982005_dp
           V=0.6308414671239979_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6426246321215801_dp
           B=0.1594171564034221_dp
           V=0.6362706466959498_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6511904367376113_dp
           B=0.5354789536565540_dp*0.1_dp
           V=0.6375414170333233_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
         end subroutine LD1730

         subroutine LD2030(X,Y,Z,W,N)
           real(dp), intent(out) :: X(2030)
           real(dp), intent(out) :: Y(2030)
           real(dp), intent(out) :: Z(2030)
           real(dp), intent(out) :: W(2030)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.4656031899197431_dp*0.0001_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.5421549195295507_dp*0.001_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2540835336814348_dp*0.1_dp
           V=0.1778522133346553_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6399322800504915_dp*0.1_dp
           V=0.2811325405682796_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1088269469804125_dp
           V=0.3548896312631459_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1570670798818287_dp
           V=0.4090310897173364_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2071163932282514_dp
           V=0.4493286134169965_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2578914044450844_dp
           V=0.4793728447962723_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3085687558169623_dp
           V=0.5015415319164265_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3584719706267024_dp
           V=0.5175127372677937_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4070135594428709_dp
           V=0.5285522262081019_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4536618626222638_dp
           V=0.5356832703713962_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4979195686463577_dp
           V=0.5397914736175170_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5393075111126999_dp
           V=0.5416899441599930_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6115617676843916_dp
           V=0.5419308476889938_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6414308435160159_dp
           V=0.5416936902030596_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6664099412721607_dp
           V=0.5419544338703164_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6859161771214913_dp
           V=0.5428983656630975_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6993625593503890_dp
           V=0.5442286500098193_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7062393387719380_dp
           V=0.5452250345057301_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7479028168349763_dp*0.1_dp
           V=0.2568002497728530_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1848951153969366_dp
           V=0.3827211700292145_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3059529066581305_dp
           V=0.4579491561917824_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4285556101021362_dp
           V=0.5042003969083574_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5468758653496526_dp
           V=0.5312708889976025_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6565821978343439_dp
           V=0.5438401790747117_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1253901572367117_dp
           B=0.3681917226439641_dp*0.1_dp
           V=0.3316041873197344_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1775721510383941_dp
           B=0.7982487607213301_dp*0.1_dp
           V=0.3899113567153771_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2305693358216114_dp
           B=0.1264640966592335_dp
           V=0.4343343327201309_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2836502845992063_dp
           B=0.1751585683418957_dp
           V=0.4679415262318919_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3361794746232590_dp
           B=0.2247995907632670_dp
           V=0.4930847981631031_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3875979172264824_dp
           B=0.2745299257422246_dp
           V=0.5115031867540091_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4374019316999074_dp
           B=0.3236373482441118_dp
           V=0.5245217148457367_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4851275843340022_dp
           B=0.3714967859436741_dp
           V=0.5332041499895321_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5303391803806868_dp
           B=0.4175353646321745_dp
           V=0.5384583126021542_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5726197380596287_dp
           B=0.4612084406355461_dp
           V=0.5411067210798852_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2431520732564863_dp
           B=0.4258040133043952_dp*0.1_dp
           V=0.4259797391468714_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3002096800895869_dp
           B=0.8869424306722721_dp*0.1_dp
           V=0.4604931368460021_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3558554457457432_dp
           B=0.1368811706510655_dp
           V=0.4871814878255202_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4097782537048887_dp
           B=0.1860739985015033_dp
           V=0.5072242910074885_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4616337666067458_dp
           B=0.2354235077395853_dp
           V=0.5217069845235350_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5110707008417874_dp
           B=0.2842074921347011_dp
           V=0.5315785966280310_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5577415286163795_dp
           B=0.3317784414984102_dp
           V=0.5376833708758905_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6013060431366950_dp
           B=0.3775299002040700_dp
           V=0.5408032092069521_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3661596767261781_dp
           B=0.4599367887164592_dp*0.1_dp
           V=0.4842744917904866_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4237633153506581_dp
           B=0.9404893773654421_dp*0.1_dp
           V=0.5048926076188130_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4786328454658452_dp
           B=0.1431377109091971_dp
           V=0.5202607980478373_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5305702076789774_dp
           B=0.1924186388843570_dp
           V=0.5309932388325743_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5793436224231788_dp
           B=0.2411590944775190_dp
           V=0.5377419770895208_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6247069017094747_dp
           B=0.2886871491583605_dp
           V=0.5411696331677717_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4874315552535204_dp
           B=0.4804978774953206_dp*0.1_dp
           V=0.5197996293282420_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5427337322059053_dp
           B=0.9716857199366665_dp*0.1_dp
           V=0.5311120836622945_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5943493747246700_dp
           B=0.1465205839795055_dp
           V=0.5384309319956951_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6421314033564943_dp
           B=0.1953579449803574_dp
           V=0.5421859504051886_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6020628374713980_dp
           B=0.4916375015738108_dp*0.1_dp
           V=0.5390948355046314_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6529222529856881_dp
           B=0.9861621540127005_dp*0.1_dp
           V=0.5433312705027845_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
         end subroutine LD2030

         subroutine LD2354(X,Y,Z,W,N)
           real(dp), intent(out) :: X(2354)
           real(dp), intent(out) :: Y(2354)
           real(dp), intent(out) :: Z(2354)
           real(dp), intent(out) :: W(2354)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.3922616270665292_dp*0.0001_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.4703831750854424_dp*0.001_dp
           Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.4678202801282136_dp*0.001_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2290024646530589_dp*0.1_dp
           V=0.1437832228979900_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5779086652271284_dp*0.1_dp
           V=0.2303572493577644_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.9863103576375984_dp*0.1_dp
           V=0.2933110752447454_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1428155792982185_dp
           V=0.3402905998359838_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1888978116601463_dp
           V=0.3759138466870372_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2359091682970210_dp
           V=0.4030638447899798_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2831228833706171_dp
           V=0.4236591432242211_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3299495857966693_dp
           V=0.4390522656946746_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3758840802660796_dp
           V=0.4502523466626247_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4204751831009480_dp
           V=0.4580577727783541_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4633068518751051_dp
           V=0.4631391616615899_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5039849474507313_dp
           V=0.4660928953698676_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5421265793440747_dp
           V=0.4674751807936953_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6092660230557310_dp
           V=0.4676414903932920_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6374654204984869_dp
           V=0.4674086492347870_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6615136472609892_dp
           V=0.4674928539483207_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6809487285958127_dp
           V=0.4680748979686447_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6952980021665196_dp
           V=0.4690449806389040_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7041245497695400_dp
           V=0.4699877075860818_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6744033088306065_dp*0.1_dp
           V=0.2099942281069176_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1678684485334166_dp
           V=0.3172269150712804_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2793559049539613_dp
           V=0.3832051358546523_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3935264218057639_dp
           V=0.4252193818146985_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5052629268232558_dp
           V=0.4513807963755000_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6107905315437531_dp
           V=0.4657797469114178_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1135081039843524_dp
           B=0.3331954884662588_dp*0.1_dp
           V=0.2733362800522836_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1612866626099378_dp
           B=0.7247167465436538_dp*0.1_dp
           V=0.3235485368463559_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2100786550168205_dp
           B=0.1151539110849745_dp
           V=0.3624908726013453_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2592282009459942_dp
           B=0.1599491097143677_dp
           V=0.3925540070712828_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3081740561320203_dp
           B=0.2058699956028027_dp
           V=0.4156129781116235_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3564289781578164_dp
           B=0.2521624953502911_dp
           V=0.4330644984623263_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4035587288240703_dp
           B=0.2982090785797674_dp
           V=0.4459677725921312_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4491671196373903_dp
           B=0.3434762087235733_dp
           V=0.4551593004456795_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4928854782917489_dp
           B=0.3874831357203437_dp
           V=0.4613341462749918_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5343646791958988_dp
           B=0.4297814821746926_dp
           V=0.4651019618269806_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5732683216530990_dp
           B=0.4699402260943537_dp
           V=0.4670249536100625_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2214131583218986_dp
           B=0.3873602040643895_dp*0.1_dp
           V=0.3549555576441708_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2741796504750071_dp
           B=0.8089496256902013_dp*0.1_dp
           V=0.3856108245249010_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3259797439149485_dp
           B=0.1251732177620872_dp
           V=0.4098622845756882_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3765441148826891_dp
           B=0.1706260286403185_dp
           V=0.4286328604268950_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4255773574530558_dp
           B=0.2165115147300408_dp
           V=0.4427802198993945_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4727795117058430_dp
           B=0.2622089812225259_dp
           V=0.4530473511488561_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5178546895819012_dp
           B=0.3071721431296201_dp
           V=0.4600805475703138_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5605141192097460_dp
           B=0.3508998998801138_dp
           V=0.4644599059958017_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6004763319352512_dp
           B=0.3929160876166931_dp
           V=0.4667274455712508_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3352842634946949_dp
           B=0.4202563457288019_dp*0.1_dp
           V=0.4069360518020356_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3891971629814670_dp
           B=0.8614309758870850_dp*0.1_dp
           V=0.4260442819919195_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4409875565542281_dp
           B=0.1314500879380001_dp
           V=0.4408678508029063_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4904893058592484_dp
           B=0.1772189657383859_dp
           V=0.4518748115548597_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5375056138769549_dp
           B=0.2228277110050294_dp
           V=0.4595564875375116_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5818255708669969_dp
           B=0.2677179935014386_dp
           V=0.4643988774315846_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6232334858144959_dp
           B=0.3113675035544165_dp
           V=0.4668827491646946_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4489485354492058_dp
           B=0.4409162378368174_dp*0.1_dp
           V=0.4400541823741973_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5015136875933150_dp
           B=0.8939009917748489_dp*0.1_dp
           V=0.4514512890193797_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5511300550512623_dp
           B=0.1351806029383365_dp
           V=0.4596198627347549_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5976720409858000_dp
           B=0.1808370355053196_dp
           V=0.4648659016801781_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6409956378989354_dp
           B=0.2257852192301602_dp
           V=0.4675502017157673_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5581222330827514_dp
           B=0.4532173421637160_dp*0.1_dp
           V=0.4598494476455523_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6074705984161695_dp
           B=0.9117488031840314_dp*0.1_dp
           V=0.4654916955152048_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6532272537379033_dp
           B=0.1369294213140155_dp
           V=0.4684709779505137_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6594761494500487_dp
           B=0.4589901487275583_dp*0.1_dp
           V=0.4691445539106986_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
         end subroutine LD2354

         subroutine LD2702(X,Y,Z,W,N)
           real(dp), intent(out) :: X(2702)
           real(dp), intent(out) :: Y(2702)
           real(dp), intent(out) :: Z(2702)
           real(dp), intent(out) :: W(2702)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.2998675149888161_dp*0.0001_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.4077860529495355_dp*0.001_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2065562538818703_dp*0.1_dp
           V=0.1185349192520667_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5250918173022379_dp*0.1_dp
           V=0.1913408643425751_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.8993480082038376_dp*0.1_dp
           V=0.2452886577209897_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1306023924436019_dp
           V=0.2862408183288702_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1732060388531418_dp
           V=0.3178032258257357_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2168727084820249_dp
           V=0.3422945667633690_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2609528309173586_dp
           V=0.3612790520235922_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3049252927938952_dp
           V=0.3758638229818521_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3483484138084404_dp
           V=0.3868711798859953_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3908321549106406_dp
           V=0.3949429933189938_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4320210071894814_dp
           V=0.4006068107541156_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4715824795890053_dp
           V=0.4043192149672723_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5091984794078453_dp
           V=0.4064947495808078_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5445580145650803_dp
           V=0.4075245619813152_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6072575796841768_dp
           V=0.4076423540893566_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6339484505755803_dp
           V=0.4074280862251555_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6570718257486958_dp
           V=0.4074163756012244_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6762557330090709_dp
           V=0.4077647795071246_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6911161696923790_dp
           V=0.4084517552782530_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7012841911659961_dp
           V=0.4092468459224052_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7064559272410020_dp
           V=0.4097872687240906_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6123554989894765_dp*0.1_dp
           V=0.1738986811745028_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1533070348312393_dp
           V=0.2659616045280191_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2563902605244206_dp
           V=0.3240596008171533_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3629346991663361_dp
           V=0.3621195964432943_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4683949968987538_dp
           V=0.3868838330760539_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5694479240657952_dp
           V=0.4018911532693111_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6634465430993955_dp
           V=0.4089929432983252_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1033958573552305_dp
           B=0.3034544009063584_dp*0.1_dp
           V=0.2279907527706409_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1473521412414395_dp
           B=0.6618803044247135_dp*0.1_dp
           V=0.2715205490578897_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1924552158705967_dp
           B=0.1054431128987715_dp
           V=0.3057917896703976_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2381094362890328_dp
           B=0.1468263551238858_dp
           V=0.3326913052452555_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2838121707936760_dp
           B=0.1894486108187886_dp
           V=0.3537334711890037_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3291323133373415_dp
           B=0.2326374238761579_dp
           V=0.3700567500783129_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3736896978741460_dp
           B=0.2758485808485768_dp
           V=0.3825245372589122_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4171406040760013_dp
           B=0.3186179331996921_dp
           V=0.3918125171518296_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4591677985256915_dp
           B=0.3605329796303794_dp
           V=0.3984720419937579_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4994733831718418_dp
           B=0.4012147253586509_dp
           V=0.4029746003338211_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5377731830445096_dp
           B=0.4403050025570692_dp
           V=0.4057428632156627_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5737917830001331_dp
           B=0.4774565904277483_dp
           V=0.4071719274114857_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2027323586271389_dp
           B=0.3544122504976147_dp*0.1_dp
           V=0.2990236950664119_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2516942375187273_dp
           B=0.7418304388646328_dp*0.1_dp
           V=0.3262951734212878_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3000227995257181_dp
           B=0.1150502745727186_dp
           V=0.3482634608242413_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3474806691046342_dp
           B=0.1571963371209364_dp
           V=0.3656596681700892_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3938103180359209_dp
           B=0.1999631877247100_dp
           V=0.3791740467794218_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4387519590455703_dp
           B=0.2428073457846535_dp
           V=0.3894034450156905_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4820503960077787_dp
           B=0.2852575132906155_dp
           V=0.3968600245508371_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5234573778475101_dp
           B=0.3268884208674639_dp
           V=0.4019931351420050_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5627318647235282_dp
           B=0.3673033321675939_dp
           V=0.4052108801278599_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5996390607156954_dp
           B=0.4061211551830290_dp
           V=0.4068978613940934_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3084780753791947_dp
           B=0.3860125523100059_dp*0.1_dp
           V=0.3454275351319704_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3589988275920223_dp
           B=0.7928938987104867_dp*0.1_dp
           V=0.3629963537007920_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4078628415881973_dp
           B=0.1212614643030087_dp
           V=0.3770187233889873_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4549287258889735_dp
           B=0.1638770827382693_dp
           V=0.3878608613694378_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5000278512957279_dp
           B=0.2065965798260176_dp
           V=0.3959065270221274_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5429785044928199_dp
           B=0.2489436378852235_dp
           V=0.4015286975463570_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5835939850491711_dp
           B=0.2904811368946891_dp
           V=0.4050866785614717_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6216870353444856_dp
           B=0.3307941957666609_dp
           V=0.4069320185051913_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4151104662709091_dp
           B=0.4064829146052554_dp*0.1_dp
           V=0.3760120964062763_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4649804275009218_dp
           B=0.8258424547294755_dp*0.1_dp
           V=0.3870969564418064_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5124695757009662_dp
           B=0.1251841962027289_dp
           V=0.3955287790534055_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5574711100606224_dp
           B=0.1679107505976331_dp
           V=0.4015361911302668_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5998597333287227_dp
           B=0.2102805057358715_dp
           V=0.4053836986719548_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6395007148516600_dp
           B=0.2518418087774107_dp
           V=0.4073578673299117_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5188456224746252_dp
           B=0.4194321676077518_dp*0.1_dp
           V=0.3954628379231406_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5664190707942778_dp
           B=0.8457661551921499_dp*0.1_dp
           V=0.4017645508847530_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6110464353283153_dp
           B=0.1273652932519396_dp
           V=0.4059030348651293_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6526430302051563_dp
           B=0.1698173239076354_dp
           V=0.4080565809484880_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6167551880377548_dp
           B=0.4266398851548864_dp*0.1_dp
           V=0.4063018753664651_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6607195418355383_dp
           B=0.8551925814238349_dp*0.1_dp
           V=0.4087191292799671_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
         end subroutine LD2702

         subroutine LD3074(X,Y,Z,W,N)
           real(dp), intent(out) :: X(3074)
           real(dp), intent(out) :: Y(3074)
           real(dp), intent(out) :: Z(3074)
           real(dp), intent(out) :: W(3074)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.2599095953754734_dp*0.0001_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.3603134089687541_dp*0.001_dp
           Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.3586067974412447_dp*0.001_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1886108518723392_dp*0.1_dp
           V=0.9831528474385880_dp*0.0001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4800217244625303_dp*0.1_dp
           V=0.1605023107954450_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.8244922058397242_dp*0.1_dp
           V=0.2072200131464099_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1200408362484023_dp
           V=0.2431297618814187_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1595773530809965_dp
           V=0.2711819064496707_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2002635973434064_dp
           V=0.2932762038321116_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2415127590139982_dp
           V=0.3107032514197368_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2828584158458477_dp
           V=0.3243808058921213_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3239091015338138_dp
           V=0.3349899091374030_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3643225097962194_dp
           V=0.3430580688505218_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4037897083691802_dp
           V=0.3490124109290343_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4420247515194127_dp
           V=0.3532148948561955_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4787572538464938_dp
           V=0.3559862669062833_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5137265251275234_dp
           V=0.3576224317551411_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5466764056654611_dp
           V=0.3584050533086076_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6054859420813535_dp
           V=0.3584903581373224_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6308106701764562_dp
           V=0.3582991879040586_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6530369230179584_dp
           V=0.3582371187963125_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6718609524611158_dp
           V=0.3584353631122350_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6869676499894013_dp
           V=0.3589120166517785_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6980467077240748_dp
           V=0.3595445704531601_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7048241721250522_dp
           V=0.3600943557111074_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5591105222058232_dp*0.1_dp
           V=0.1456447096742039_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1407384078513916_dp
           V=0.2252370188283782_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2364035438976309_dp
           V=0.2766135443474897_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3360602737818170_dp
           V=0.3110729491500851_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4356292630054665_dp
           V=0.3342506712303391_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5321569415256174_dp
           V=0.3491981834026860_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6232956305040554_dp
           V=0.3576003604348932_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.9469870086838469_dp*0.1_dp
           B=0.2778748387309470_dp*0.1_dp
           V=0.1921921305788564_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1353170300568141_dp
           B=0.6076569878628364_dp*0.1_dp
           V=0.2301458216495632_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1771679481726077_dp
           B=0.9703072762711040_dp*0.1_dp
           V=0.2604248549522893_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2197066664231751_dp
           B=0.1354112458524762_dp
           V=0.2845275425870697_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2624783557374927_dp
           B=0.1750996479744100_dp
           V=0.3036870897974840_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3050969521214442_dp
           B=0.2154896907449802_dp
           V=0.3188414832298066_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3472252637196021_dp
           B=0.2560954625740152_dp
           V=0.3307046414722089_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3885610219026360_dp
           B=0.2965070050624096_dp
           V=0.3398330969031360_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4288273776062765_dp
           B=0.3363641488734497_dp
           V=0.3466757899705373_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4677662471302948_dp
           B=0.3753400029836788_dp
           V=0.3516095923230054_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5051333589553359_dp
           B=0.4131297522144286_dp
           V=0.3549645184048486_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5406942145810492_dp
           B=0.4494423776081795_dp
           V=0.3570415969441392_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5742204122576457_dp
           B=0.4839938958841502_dp
           V=0.3581251798496118_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1865407027225188_dp
           B=0.3259144851070796_dp*0.1_dp
           V=0.2543491329913348_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2321186453689432_dp
           B=0.6835679505297343_dp*0.1_dp
           V=0.2786711051330776_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2773159142523882_dp
           B=0.1062284864451989_dp
           V=0.2985552361083679_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3219200192237254_dp
           B=0.1454404409323047_dp
           V=0.3145867929154039_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3657032593944029_dp
           B=0.1854018282582510_dp
           V=0.3273290662067609_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4084376778363622_dp
           B=0.2256297412014750_dp
           V=0.3372705511943501_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4499004945751427_dp
           B=0.2657104425000896_dp
           V=0.3448274437851510_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4898758141326335_dp
           B=0.3052755487631557_dp
           V=0.3503592783048583_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5281547442266309_dp
           B=0.3439863920645423_dp
           V=0.3541854792663162_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5645346989813992_dp
           B=0.3815229456121914_dp
           V=0.3565995517909428_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5988181252159848_dp
           B=0.4175752420966734_dp
           V=0.3578802078302898_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2850425424471603_dp
           B=0.3562149509862536_dp*0.1_dp
           V=0.2958644592860982_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3324619433027876_dp
           B=0.7330318886871096_dp*0.1_dp
           V=0.3119548129116835_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3785848333076282_dp
           B=0.1123226296008472_dp
           V=0.3250745225005984_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4232891028562115_dp
           B=0.1521084193337708_dp
           V=0.3355153415935208_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4664287050829722_dp
           B=0.1921844459223610_dp
           V=0.3435847568549328_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5078458493735726_dp
           B=0.2321360989678303_dp
           V=0.3495786831622488_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5473779816204180_dp
           B=0.2715886486360520_dp
           V=0.3537767805534621_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5848617133811376_dp
           B=0.3101924707571355_dp
           V=0.3564459815421428_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6201348281584888_dp
           B=0.3476121052890973_dp
           V=0.3578464061225468_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3852191185387871_dp
           B=0.3763224880035108_dp*0.1_dp
           V=0.3239748762836212_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4325025061073423_dp
           B=0.7659581935637135_dp*0.1_dp
           V=0.3345491784174287_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4778486229734490_dp
           B=0.1163381306083900_dp
           V=0.3429126177301782_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5211663693009000_dp
           B=0.1563890598752899_dp
           V=0.3492420343097421_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5623469504853703_dp
           B=0.1963320810149200_dp
           V=0.3537399050235257_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6012718188659246_dp
           B=0.2357847407258738_dp
           V=0.3566209152659172_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6378179206390117_dp
           B=0.2743846121244060_dp
           V=0.3581084321919782_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4836936460214534_dp
           B=0.3895902610739024_dp*0.1_dp
           V=0.3426522117591512_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5293792562683797_dp
           B=0.7871246819312640_dp*0.1_dp
           V=0.3491848770121379_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5726281253100033_dp
           B=0.1187963808202981_dp
           V=0.3539318235231476_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6133658776169068_dp
           B=0.1587914708061787_dp
           V=0.3570231438458694_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6515085491865307_dp
           B=0.1983058575227646_dp
           V=0.3586207335051714_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5778692716064976_dp
           B=0.3977209689791542_dp*0.1_dp
           V=0.3541196205164025_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6207904288086192_dp
           B=0.7990157592981152_dp*0.1_dp
           V=0.3574296911573953_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6608688171046802_dp
           B=0.1199671308754309_dp
           V=0.3591993279818963_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6656263089489130_dp
           B=0.4015955957805969_dp*0.1_dp
           V=0.3595855034661997_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
         end subroutine LD3074

         subroutine LD3470(X,Y,Z,W,N)
           real(dp), intent(out) :: X(3470)
           real(dp), intent(out) :: Y(3470)
           real(dp), intent(out) :: Z(3470)
           real(dp), intent(out) :: W(3470)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.2040382730826330_dp*0.0001_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.3178149703889544_dp*0.001_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1721420832906233_dp*0.1_dp
           V=0.8288115128076110_dp*0.0001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4408875374981770_dp*0.1_dp
           V=0.1360883192522954_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7594680813878681_dp*0.1_dp
           V=0.1766854454542662_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1108335359204799_dp
           V=0.2083153161230153_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1476517054388567_dp
           V=0.2333279544657158_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1856731870860615_dp
           V=0.2532809539930247_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2243634099428821_dp
           V=0.2692472184211158_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2633006881662727_dp
           V=0.2819949946811885_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3021340904916283_dp
           V=0.2920953593973030_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3405594048030089_dp
           V=0.2999889782948352_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3783044434007372_dp
           V=0.3060292120496902_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4151194767407910_dp
           V=0.3105109167522192_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4507705766443257_dp
           V=0.3136902387550312_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4850346056573187_dp
           V=0.3157984652454632_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5176950817792470_dp
           V=0.3170516518425422_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5485384240820989_dp
           V=0.3176568425633755_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6039117238943308_dp
           V=0.3177198411207062_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6279956655573113_dp
           V=0.3175519492394733_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6493636169568952_dp
           V=0.3174654952634756_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6677644117704504_dp
           V=0.3175676415467654_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6829368572115624_dp
           V=0.3178923417835410_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6946195818184121_dp
           V=0.3183788287531909_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7025711542057026_dp
           V=0.3188755151918807_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7066004767140119_dp
           V=0.3191916889313849_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5132537689946062_dp*0.1_dp
           V=0.1231779611744508_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1297994661331225_dp
           V=0.1924661373839880_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2188852049401307_dp
           V=0.2380881867403424_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3123174824903457_dp
           V=0.2693100663037885_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4064037620738195_dp
           V=0.2908673382834366_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4984958396944782_dp
           V=0.3053914619381535_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5864975046021365_dp
           V=0.3143916684147777_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6686711634580175_dp
           V=0.3187042244055363_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.8715738780835950_dp*0.1_dp
           B=0.2557175233367578_dp*0.1_dp
           V=0.1635219535869790_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1248383123134007_dp
           B=0.5604823383376681_dp*0.1_dp
           V=0.1968109917696070_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1638062693383378_dp
           B=0.8968568601900765_dp*0.1_dp
           V=0.2236754342249974_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2035586203373176_dp
           B=0.1254086651976279_dp
           V=0.2453186687017181_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2436798975293774_dp
           B=0.1624780150162012_dp
           V=0.2627551791580541_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2838207507773806_dp
           B=0.2003422342683208_dp
           V=0.2767654860152220_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3236787502217692_dp
           B=0.2385628026255263_dp
           V=0.2879467027765895_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3629849554840691_dp
           B=0.2767731148783578_dp
           V=0.2967639918918702_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4014948081992087_dp
           B=0.3146542308245309_dp
           V=0.3035900684660351_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4389818379260225_dp
           B=0.3519196415895088_dp
           V=0.3087338237298308_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4752331143674377_dp
           B=0.3883050984023654_dp
           V=0.3124608838860167_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5100457318374018_dp
           B=0.4235613423908649_dp
           V=0.3150084294226743_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5432238388954868_dp
           B=0.4574484717196220_dp
           V=0.3165958398598402_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5745758685072442_dp
           B=0.4897311639255524_dp
           V=0.3174320440957372_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1723981437592809_dp
           B=0.3010630597881105_dp*0.1_dp
           V=0.2182188909812599_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2149553257844597_dp
           B=0.6326031554204694_dp*0.1_dp
           V=0.2399727933921445_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2573256081247422_dp
           B=0.9848566980258631_dp*0.1_dp
           V=0.2579796133514652_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2993163751238106_dp
           B=0.1350835952384266_dp
           V=0.2727114052623535_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3407238005148000_dp
           B=0.1725184055442181_dp
           V=0.2846327656281355_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3813454978483264_dp
           B=0.2103559279730725_dp
           V=0.2941491102051334_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4209848104423343_dp
           B=0.2482278774554860_dp
           V=0.3016049492136107_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4594519699996300_dp
           B=0.2858099509982883_dp
           V=0.3072949726175648_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4965640166185930_dp
           B=0.3228075659915428_dp
           V=0.3114768142886460_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5321441655571562_dp
           B=0.3589459907204151_dp
           V=0.3143823673666223_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5660208438582166_dp
           B=0.3939630088864310_dp
           V=0.3162269764661535_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5980264315964364_dp
           B=0.4276029922949089_dp
           V=0.3172164663759821_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2644215852350733_dp
           B=0.3300939429072552_dp*0.1_dp
           V=0.2554575398967435_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3090113743443063_dp
           B=0.6803887650078501_dp*0.1_dp
           V=0.2701704069135677_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3525871079197808_dp
           B=0.1044326136206709_dp
           V=0.2823693413468940_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3950418005354029_dp
           B=0.1416751597517679_dp
           V=0.2922898463214289_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4362475663430163_dp
           B=0.1793408610504821_dp
           V=0.3001829062162428_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4760661812145854_dp
           B=0.2170630750175722_dp
           V=0.3062890864542953_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5143551042512103_dp
           B=0.2545145157815807_dp
           V=0.3108328279264746_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5509709026935597_dp
           B=0.2913940101706601_dp
           V=0.3140243146201245_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5857711030329428_dp
           B=0.3274169910910705_dp
           V=0.3160638030977130_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6186149917404392_dp
           B=0.3623081329317265_dp
           V=0.3171462882206275_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3586894569557064_dp
           B=0.3497354386450040_dp*0.1_dp
           V=0.2812388416031796_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4035266610019441_dp
           B=0.7129736739757095_dp*0.1_dp
           V=0.2912137500288045_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4467775312332510_dp
           B=0.1084758620193165_dp
           V=0.2993241256502206_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4883638346608543_dp
           B=0.1460915689241772_dp
           V=0.3057101738983822_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5281908348434601_dp
           B=0.1837790832369980_dp
           V=0.3105319326251432_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5661542687149311_dp
           B=0.2212075390874021_dp
           V=0.3139565514428167_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6021450102031452_dp
           B=0.2580682841160985_dp
           V=0.3161543006806366_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6360520783610050_dp
           B=0.2940656362094121_dp
           V=0.3172985960613294_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4521611065087196_dp
           B=0.3631055365867002_dp*0.1_dp
           V=0.2989400336901431_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4959365651560963_dp
           B=0.7348318468484350_dp*0.1_dp
           V=0.3054555883947677_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5376815804038283_dp
           B=0.1111087643812648_dp
           V=0.3104764960807702_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5773314480243768_dp
           B=0.1488226085145408_dp
           V=0.3141015825977616_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6148113245575056_dp
           B=0.1862892274135151_dp
           V=0.3164520621159896_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6500407462842380_dp
           B=0.2231909701714456_dp
           V=0.3176652305912204_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5425151448707213_dp
           B=0.3718201306118944_dp*0.1_dp
           V=0.3105097161023939_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5841860556907931_dp
           B=0.7483616335067346_dp*0.1_dp
           V=0.3143014117890550_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6234632186851500_dp
           B=0.1125990834266120_dp
           V=0.3168172866287200_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6602934551848843_dp
           B=0.1501303813157619_dp
           V=0.3181401865570968_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6278573968375105_dp
           B=0.3767559930245720_dp*0.1_dp
           V=0.3170663659156037_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6665611711264577_dp
           B=0.7548443301360158_dp*0.1_dp
           V=0.3185447944625510_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
         end subroutine LD3470

         subroutine LD3890(X,Y,Z,W,N)
           real(dp), intent(out) :: X(3890)
           real(dp), intent(out) :: Y(3890)
           real(dp), intent(out) :: Z(3890)
           real(dp), intent(out) :: W(3890)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.1807395252196920_dp*0.0001_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.2848008782238827_dp*0.001_dp
           Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.2836065837530581_dp*0.001_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1587876419858352_dp*0.1_dp
           V=0.7013149266673816_dp*0.0001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4069193593751206_dp*0.1_dp
           V=0.1162798021956766_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7025888115257997_dp*0.1_dp
           V=0.1518728583972105_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1027495450028704_dp
           V=0.1798796108216934_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1371457730893426_dp
           V=0.2022593385972785_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1727758532671953_dp
           V=0.2203093105575464_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2091492038929037_dp
           V=0.2349294234299855_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2458813281751915_dp
           V=0.2467682058747003_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2826545859450066_dp
           V=0.2563092683572224_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3191957291799622_dp
           V=0.2639253896763318_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3552621469299578_dp
           V=0.2699137479265108_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3906329503406230_dp
           V=0.2745196420166739_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4251028614093031_dp
           V=0.2779529197397593_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4584777520111870_dp
           V=0.2803996086684265_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4905711358710193_dp
           V=0.2820302356715842_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5212011669847385_dp
           V=0.2830056747491068_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5501878488737995_dp
           V=0.2834808950776839_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6025037877479342_dp
           V=0.2835282339078929_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6254572689549016_dp
           V=0.2833819267065800_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6460107179528248_dp
           V=0.2832858336906784_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6639541138154251_dp
           V=0.2833268235451244_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6790688515667495_dp
           V=0.2835432677029253_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6911338580371512_dp
           V=0.2839091722743049_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6999385956126490_dp
           V=0.2843308178875841_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7053037748656896_dp
           V=0.2846703550533846_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4732224387180115_dp*0.1_dp
           V=0.1051193406971900_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1202100529326803_dp
           V=0.1657871838796974_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2034304820664855_dp
           V=0.2064648113714232_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2912285643573002_dp
           V=0.2347942745819741_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3802361792726768_dp
           V=0.2547775326597726_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4680598511056146_dp
           V=0.2686876684847025_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5528151052155599_dp
           V=0.2778665755515867_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6329386307803041_dp
           V=0.2830996616782929_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.8056516651369069_dp*0.1_dp
           B=0.2363454684003124_dp*0.1_dp
           V=0.1403063340168372_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1156476077139389_dp
           B=0.5191291632545936_dp*0.1_dp
           V=0.1696504125939477_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1520473382760421_dp
           B=0.8322715736994519_dp*0.1_dp
           V=0.1935787242745390_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1892986699745931_dp
           B=0.1165855667993712_dp
           V=0.2130614510521968_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2270194446777792_dp
           B=0.1513077167409504_dp
           V=0.2289381265931048_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2648908185093273_dp
           B=0.1868882025807859_dp
           V=0.2418630292816186_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3026389259574136_dp
           B=0.2229277629776224_dp
           V=0.2523400495631193_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3400220296151384_dp
           B=0.2590951840746235_dp
           V=0.2607623973449605_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3768217953335510_dp
           B=0.2951047291750847_dp
           V=0.2674441032689209_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4128372900921884_dp
           B=0.3307019714169930_dp
           V=0.2726432360343356_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4478807131815630_dp
           B=0.3656544101087634_dp
           V=0.2765787685924545_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4817742034089257_dp
           B=0.3997448951939695_dp
           V=0.2794428690642224_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5143472814653344_dp
           B=0.4327667110812024_dp
           V=0.2814099002062895_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5454346213905650_dp
           B=0.4645196123532293_dp
           V=0.2826429531578994_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5748739313170252_dp
           B=0.4948063555703345_dp
           V=0.2832983542550884_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1599598738286342_dp
           B=0.2792357590048985_dp*0.1_dp
           V=0.1886695565284976_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1998097412500951_dp
           B=0.5877141038139065_dp*0.1_dp
           V=0.2081867882748234_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2396228952566202_dp
           B=0.9164573914691377_dp*0.1_dp
           V=0.2245148680600796_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2792228341097746_dp
           B=0.1259049641962687_dp
           V=0.2380370491511872_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3184251107546741_dp
           B=0.1610594823400863_dp
           V=0.2491398041852455_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3570481164426244_dp
           B=0.1967151653460898_dp
           V=0.2581632405881230_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3949164710492144_dp
           B=0.2325404606175168_dp
           V=0.2653965506227417_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4318617293970503_dp
           B=0.2682461141151439_dp
           V=0.2710857216747087_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4677221009931678_dp
           B=0.3035720116011973_dp
           V=0.2754434093903659_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5023417939270955_dp
           B=0.3382781859197439_dp
           V=0.2786579932519380_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5355701836636128_dp
           B=0.3721383065625942_dp
           V=0.2809011080679474_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5672608451328771_dp
           B=0.4049346360466055_dp
           V=0.2823336184560987_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5972704202540162_dp
           B=0.4364538098633802_dp
           V=0.2831101175806309_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2461687022333596_dp
           B=0.3070423166833368_dp*0.1_dp
           V=0.2221679970354546_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2881774566286831_dp
           B=0.6338034669281885_dp*0.1_dp
           V=0.2356185734270703_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3293963604116978_dp
           B=0.9742862487067941_dp*0.1_dp
           V=0.2469228344805590_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3697303822241377_dp
           B=0.1323799532282290_dp
           V=0.2562726348642046_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4090663023135127_dp
           B=0.1678497018129336_dp
           V=0.2638756726753028_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4472819355411712_dp
           B=0.2035095105326114_dp
           V=0.2699311157390862_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4842513377231437_dp
           B=0.2390692566672091_dp
           V=0.2746233268403837_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5198477629962928_dp
           B=0.2742649818076149_dp
           V=0.2781225674454771_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5539453011883145_dp
           B=0.3088503806580094_dp
           V=0.2805881254045684_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5864196762401251_dp
           B=0.3425904245906614_dp
           V=0.2821719877004913_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6171484466668390_dp
           B=0.3752562294789468_dp
           V=0.2830222502333124_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3350337830565727_dp
           B=0.3261589934634747_dp*0.1_dp
           V=0.2457995956744870_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3775773224758284_dp
           B=0.6658438928081572_dp*0.1_dp
           V=0.2551474407503706_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4188155229848973_dp
           B=0.1014565797157954_dp
           V=0.2629065335195311_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4586805892009344_dp
           B=0.1368573320843822_dp
           V=0.2691900449925075_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4970895714224235_dp
           B=0.1724614851951608_dp
           V=0.2741275485754276_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5339505133960747_dp
           B=0.2079779381416412_dp
           V=0.2778530970122595_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5691665792531440_dp
           B=0.2431385788322288_dp
           V=0.2805010567646741_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6026387682680377_dp
           B=0.2776901883049853_dp
           V=0.2822055834031040_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6342676150163307_dp
           B=0.3113881356386632_dp
           V=0.2831016901243473_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4237951119537067_dp
           B=0.3394877848664351_dp*0.1_dp
           V=0.2624474901131803_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4656918683234929_dp
           B=0.6880219556291447_dp*0.1_dp
           V=0.2688034163039377_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5058857069185980_dp
           B=0.1041946859721635_dp
           V=0.2738932751287636_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5443204666713996_dp
           B=0.1398039738736393_dp
           V=0.2777944791242523_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5809298813759742_dp
           B=0.1753373381196155_dp
           V=0.2806011661660987_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6156416039447128_dp
           B=0.2105215793514010_dp
           V=0.2824181456597460_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6483801351066604_dp
           B=0.2450953312157051_dp
           V=0.2833585216577828_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5103616577251688_dp
           B=0.3485560643800719_dp*0.1_dp
           V=0.2738165236962878_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5506738792580681_dp
           B=0.7026308631512033_dp*0.1_dp
           V=0.2778365208203180_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5889573040995292_dp
           B=0.1059035061296403_dp
           V=0.2807852940418966_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6251641589516930_dp
           B=0.1414823925236026_dp
           V=0.2827245949674705_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6592414921570178_dp
           B=0.1767207908214530_dp
           V=0.2837342344829828_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5930314017533384_dp
           B=0.3542189339561672_dp*0.1_dp
           V=0.2809233907610981_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6309812253390175_dp
           B=0.7109574040369549_dp*0.1_dp
           V=0.2829930809742694_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6666296011353230_dp
           B=0.1067259792282730_dp
           V=0.2841097874111479_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6703715271049922_dp
           B=0.3569455268820809_dp*0.1_dp
           V=0.2843455206008783_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
         end subroutine LD3890

         subroutine LD4334(X,Y,Z,W,N)
           real(dp), intent(out) :: X(4334)
           real(dp), intent(out) :: Y(4334)
           real(dp), intent(out) :: Z(4334)
           real(dp), intent(out) :: W(4334)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.1449063022537883_dp*0.0001_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.2546377329828424_dp*0.001_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1462896151831013_dp*0.1_dp
           V=0.6018432961087496_dp*0.0001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3769840812493139_dp*0.1_dp
           V=0.1002286583263673_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6524701904096891_dp*0.1_dp
           V=0.1315222931028093_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.9560543416134648_dp*0.1_dp
           V=0.1564213746876724_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1278335898929198_dp
           V=0.1765118841507736_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1613096104466031_dp
           V=0.1928737099311080_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1955806225745371_dp
           V=0.2062658534263270_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2302935218498028_dp
           V=0.2172395445953787_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2651584344113027_dp
           V=0.2262076188876047_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2999276825183209_dp
           V=0.2334885699462397_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3343828669718798_dp
           V=0.2393355273179203_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3683265013750518_dp
           V=0.2439559200468863_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4015763206518108_dp
           V=0.2475251866060002_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4339612026399770_dp
           V=0.2501965558158773_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4653180651114582_dp
           V=0.2521081407925925_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4954893331080803_dp
           V=0.2533881002388081_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5243207068924930_dp
           V=0.2541582900848261_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5516590479041704_dp
           V=0.2545365737525860_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6012371927804176_dp
           V=0.2545726993066799_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6231574466449819_dp
           V=0.2544456197465555_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6429416514181271_dp
           V=0.2543481596881064_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6604124272943595_dp
           V=0.2543506451429194_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6753851470408250_dp
           V=0.2544905675493763_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6876717970626160_dp
           V=0.2547611407344429_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6970895061319234_dp
           V=0.2551060375448869_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7034746912553310_dp
           V=0.2554291933816039_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7067017217542295_dp
           V=0.2556255710686343_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4382223501131123_dp*0.1_dp
           V=0.9041339695118195_dp*0.0001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1117474077400006_dp
           V=0.1438426330079022_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1897153252911440_dp
           V=0.1802523089820518_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2724023009910331_dp
           V=0.2060052290565496_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3567163308709902_dp
           V=0.2245002248967466_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4404784483028087_dp
           V=0.2377059847731150_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5219833154161411_dp
           V=0.2468118955882525_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5998179868977553_dp
           V=0.2525410872966528_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6727803154548222_dp
           V=0.2553101409933397_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7476563943166086_dp*0.1_dp
           B=0.2193168509461185_dp*0.1_dp
           V=0.1212879733668632_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1075341482001416_dp
           B=0.4826419281533887_dp*0.1_dp
           V=0.1472872881270931_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1416344885203259_dp
           B=0.7751191883575742_dp*0.1_dp
           V=0.1686846601010828_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1766325315388586_dp
           B=0.1087558139247680_dp
           V=0.1862698414660208_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2121744174481514_dp
           B=0.1413661374253096_dp
           V=0.2007430956991861_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2479669443408145_dp
           B=0.1748768214258880_dp
           V=0.2126568125394796_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2837600452294113_dp
           B=0.2089216406612073_dp
           V=0.2224394603372113_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3193344933193984_dp
           B=0.2431987685545972_dp
           V=0.2304264522673135_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3544935442438745_dp
           B=0.2774497054377770_dp
           V=0.2368854288424087_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3890571932288154_dp
           B=0.3114460356156915_dp
           V=0.2420352089461772_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4228581214259090_dp
           B=0.3449806851913012_dp
           V=0.2460597113081295_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4557387211304052_dp
           B=0.3778618641248256_dp
           V=0.2491181912257687_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4875487950541643_dp
           B=0.4099086391698978_dp
           V=0.2513528194205857_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5181436529962997_dp
           B=0.4409474925853973_dp
           V=0.2528943096693220_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5473824095600661_dp
           B=0.4708094517711291_dp
           V=0.2538660368488136_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5751263398976174_dp
           B=0.4993275140354637_dp
           V=0.2543868648299022_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1489515746840028_dp
           B=0.2599381993267017_dp*0.1_dp
           V=0.1642595537825183_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1863656444351767_dp
           B=0.5479286532462190_dp*0.1_dp
           V=0.1818246659849308_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2238602880356348_dp
           B=0.8556763251425254_dp*0.1_dp
           V=0.1966565649492420_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2612723375728160_dp
           B=0.1177257802267011_dp
           V=0.2090677905657991_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2984332990206190_dp
           B=0.1508168456192700_dp
           V=0.2193820409510504_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3351786584663333_dp
           B=0.1844801892177727_dp
           V=0.2278870827661928_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3713505522209120_dp
           B=0.2184145236087598_dp
           V=0.2348283192282090_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4067981098954663_dp
           B=0.2523590641486229_dp
           V=0.2404139755581477_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4413769993687534_dp
           B=0.2860812976901373_dp
           V=0.2448227407760734_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4749487182516394_dp
           B=0.3193686757808996_dp
           V=0.2482110455592573_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5073798105075426_dp
           B=0.3520226949547602_dp
           V=0.2507192397774103_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5385410448878654_dp
           B=0.3838544395667890_dp
           V=0.2524765968534880_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5683065353670530_dp
           B=0.4146810037640963_dp
           V=0.2536052388539425_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5965527620663510_dp
           B=0.4443224094681121_dp
           V=0.2542230588033068_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2299227700856157_dp
           B=0.2865757664057584_dp*0.1_dp
           V=0.1944817013047896_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2695752998553267_dp
           B=0.5923421684485993_dp*0.1_dp
           V=0.2067862362746635_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3086178716611389_dp
           B=0.9117817776057715_dp*0.1_dp
           V=0.2172440734649114_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3469649871659077_dp
           B=0.1240593814082605_dp
           V=0.2260125991723423_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3845153566319655_dp
           B=0.1575272058259175_dp
           V=0.2332655008689523_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4211600033403215_dp
           B=0.1912845163525413_dp
           V=0.2391699681532458_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4567867834329882_dp
           B=0.2250710177858171_dp
           V=0.2438801528273928_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4912829319232061_dp
           B=0.2586521303440910_dp
           V=0.2475370504260665_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5245364793303812_dp
           B=0.2918112242865407_dp
           V=0.2502707235640574_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5564369788915756_dp
           B=0.3243439239067890_dp
           V=0.2522031701054241_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5868757697775287_dp
           B=0.3560536787835351_dp
           V=0.2534511269978784_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6157458853519617_dp
           B=0.3867480821242581_dp
           V=0.2541284914955151_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3138461110672113_dp
           B=0.3051374637507278_dp*0.1_dp
           V=0.2161509250688394_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3542495872050569_dp
           B=0.6237111233730755_dp*0.1_dp
           V=0.2248778513437852_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3935751553120181_dp
           B=0.9516223952401907_dp*0.1_dp
           V=0.2322388803404617_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4317634668111147_dp
           B=0.1285467341508517_dp
           V=0.2383265471001355_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4687413842250821_dp
           B=0.1622318931656033_dp
           V=0.2432476675019525_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5044274237060283_dp
           B=0.1959581153836453_dp
           V=0.2471122223750674_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5387354077925727_dp
           B=0.2294888081183837_dp
           V=0.2500291752486870_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5715768898356105_dp
           B=0.2626031152713945_dp
           V=0.2521055942764682_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6028627200136111_dp
           B=0.2950904075286713_dp
           V=0.2534472785575503_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6325039812653463_dp
           B=0.3267458451113286_dp
           V=0.2541599713080121_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3981986708423407_dp
           B=0.3183291458749821_dp*0.1_dp
           V=0.2317380975862936_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4382791182133300_dp
           B=0.6459548193880908_dp*0.1_dp
           V=0.2378550733719775_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4769233057218166_dp
           B=0.9795757037087952_dp*0.1_dp
           V=0.2428884456739118_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5140823911194238_dp
           B=0.1316307235126655_dp
           V=0.2469002655757292_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5496977833862983_dp
           B=0.1653556486358704_dp
           V=0.2499657574265851_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5837047306512727_dp
           B=0.1988931724126510_dp
           V=0.2521676168486082_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6160349566926879_dp
           B=0.2320174581438950_dp
           V=0.2535935662645334_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6466185353209440_dp
           B=0.2645106562168662_dp
           V=0.2543356743363214_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4810835158795404_dp
           B=0.3275917807743992_dp*0.1_dp
           V=0.2427353285201535_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5199925041324341_dp
           B=0.6612546183967181_dp*0.1_dp
           V=0.2468258039744386_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5571717692207494_dp
           B=0.9981498331474143_dp*0.1_dp
           V=0.2500060956440310_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5925789250836378_dp
           B=0.1335687001410374_dp
           V=0.2523238365420979_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6261658523859670_dp
           B=0.1671444402896463_dp
           V=0.2538399260252846_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6578811126669331_dp
           B=0.2003106382156076_dp
           V=0.2546255927268069_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5609624612998100_dp
           B=0.3337500940231335_dp*0.1_dp
           V=0.2500583360048449_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5979959659984670_dp
           B=0.6708750335901803_dp*0.1_dp
           V=0.2524777638260203_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6330523711054002_dp
           B=0.1008792126424850_dp
           V=0.2540951193860656_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6660960998103972_dp
           B=0.1345050343171794_dp
           V=0.2549524085027472_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6365384364585819_dp
           B=0.3372799460737052_dp*0.1_dp
           V=0.2542569507009158_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6710994302899275_dp
           B=0.6755249309678028_dp*0.1_dp
           V=0.2552114127580376_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
         end subroutine LD4334

         subroutine LD4802(X,Y,Z,W,N)
           real(dp), intent(out) :: X(4802)
           real(dp), intent(out) :: Y(4802)
           real(dp), intent(out) :: Z(4802)
           real(dp), intent(out) :: W(4802)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.9687521879420705_dp*0.0001_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.2307897895367918_dp*0.001_dp
           Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.2297310852498558_dp*0.001_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2335728608887064_dp*0.1_dp
           V=0.7386265944001919_dp*0.0001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4352987836550653_dp*0.1_dp
           V=0.8257977698542210_dp*0.0001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6439200521088801_dp*0.1_dp
           V=0.9706044762057630_dp*0.0001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.9003943631993181_dp*0.1_dp
           V=0.1302393847117003_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1196706615548473_dp
           V=0.1541957004600968_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1511715412838134_dp
           V=0.1704459770092199_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1835982828503801_dp
           V=0.1827374890942906_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2165081259155405_dp
           V=0.1926360817436107_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2496208720417563_dp
           V=0.2008010239494833_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2827200673567900_dp
           V=0.2075635983209175_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3156190823994346_dp
           V=0.2131306638690909_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3481476793749115_dp
           V=0.2176562329937335_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3801466086947226_dp
           V=0.2212682262991018_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4114652119634011_dp
           V=0.2240799515668565_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4419598786519751_dp
           V=0.2261959816187525_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4714925949329543_dp
           V=0.2277156368808855_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4999293972879466_dp
           V=0.2287351772128336_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5271387221431248_dp
           V=0.2293490814084085_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5529896780837761_dp
           V=0.2296505312376273_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6000856099481712_dp
           V=0.2296793832318756_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6210562192785175_dp
           V=0.2295785443842974_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6401165879934240_dp
           V=0.2295017931529102_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6571144029244334_dp
           V=0.2295059638184868_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6718910821718863_dp
           V=0.2296232343237362_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6842845591099010_dp
           V=0.2298530178740771_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6941353476269816_dp
           V=0.2301579790280501_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7012965242212991_dp
           V=0.2304690404996513_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7056471428242644_dp
           V=0.2307027995907102_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4595557643585895_dp*0.1_dp
           V=0.9312274696671092_dp*0.0001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1049316742435023_dp
           V=0.1199919385876926_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1773548879549274_dp
           V=0.1598039138877690_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2559071411236127_dp
           V=0.1822253763574900_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3358156837985898_dp
           V=0.1988579593655040_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4155835743763893_dp
           V=0.2112620102533307_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4937894296167472_dp
           V=0.2201594887699007_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5691569694793316_dp
           V=0.2261622590895036_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6405840854894251_dp
           V=0.2296458453435705_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7345133894143348_dp*0.1_dp
           B=0.2177844081486067_dp*0.1_dp
           V=0.1006006990267000_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1009859834044931_dp
           B=0.4590362185775188_dp*0.1_dp
           V=0.1227676689635876_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1324289619748758_dp
           B=0.7255063095690877_dp*0.1_dp
           V=0.1467864280270117_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1654272109607127_dp
           B=0.1017825451960684_dp
           V=0.1644178912101232_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1990767186776461_dp
           B=0.1325652320980364_dp
           V=0.1777664890718961_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2330125945523278_dp
           B=0.1642765374496765_dp
           V=0.1884825664516690_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2670080611108287_dp
           B=0.1965360374337889_dp
           V=0.1973269246453848_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3008753376294316_dp
           B=0.2290726770542238_dp
           V=0.2046767775855328_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3344475596167860_dp
           B=0.2616645495370823_dp
           V=0.2107600125918040_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3675709724070786_dp
           B=0.2941150728843141_dp
           V=0.2157416362266829_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4001000887587812_dp
           B=0.3262440400919066_dp
           V=0.2197557816920721_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4318956350436028_dp
           B=0.3578835350611916_dp
           V=0.2229192611835437_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4628239056795531_dp
           B=0.3888751854043678_dp
           V=0.2253385110212775_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4927563229773636_dp
           B=0.4190678003222840_dp
           V=0.2271137107548774_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5215687136707969_dp
           B=0.4483151836883852_dp
           V=0.2283414092917525_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5491402346984905_dp
           B=0.4764740676087880_dp
           V=0.2291161673130077_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5753520160126075_dp
           B=0.5034021310998277_dp
           V=0.2295313908576598_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1388326356417754_dp
           B=0.2435436510372806_dp*0.1_dp
           V=0.1438204721359031_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1743686900537244_dp
           B=0.5118897057342652_dp*0.1_dp
           V=0.1607738025495257_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2099737037950268_dp
           B=0.8014695048539634_dp*0.1_dp
           V=0.1741483853528379_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2454492590908548_dp
           B=0.1105117874155699_dp
           V=0.1851918467519151_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2807219257864278_dp
           B=0.1417950531570966_dp
           V=0.1944628638070613_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3156842271975842_dp
           B=0.1736604945719597_dp
           V=0.2022495446275152_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3502090945177752_dp
           B=0.2058466324693981_dp
           V=0.2087462382438514_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3841684849519686_dp
           B=0.2381284261195919_dp
           V=0.2141074754818308_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4174372367906016_dp
           B=0.2703031270422569_dp
           V=0.2184640913748162_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4498926465011892_dp
           B=0.3021845683091309_dp
           V=0.2219309165220329_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4814146229807701_dp
           B=0.3335993355165720_dp
           V=0.2246123118340624_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5118863625734701_dp
           B=0.3643833735518232_dp
           V=0.2266062766915125_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5411947455119144_dp
           B=0.3943789541958179_dp
           V=0.2280072952230796_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5692301500357246_dp
           B=0.4234320144403542_dp
           V=0.2289082025202583_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5958857204139576_dp
           B=0.4513897947419260_dp
           V=0.2294012695120025_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2156270284785766_dp
           B=0.2681225755444491_dp*0.1_dp
           V=0.1722434488736947_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2532385054909710_dp
           B=0.5557495747805614_dp*0.1_dp
           V=0.1830237421455091_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2902564617771537_dp
           B=0.8569368062950249_dp*0.1_dp
           V=0.1923855349997633_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3266979823143256_dp
           B=0.1167367450324135_dp
           V=0.2004067861936271_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3625039627493614_dp
           B=0.1483861994003304_dp
           V=0.2071817297354263_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3975838937548699_dp
           B=0.1803821503011405_dp
           V=0.2128250834102103_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4318396099009774_dp
           B=0.2124962965666424_dp
           V=0.2174513719440102_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4651706555732742_dp
           B=0.2445221837805913_dp
           V=0.2211661839150214_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4974752649620969_dp
           B=0.2762701224322987_dp
           V=0.2240665257813102_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5286517579627517_dp
           B=0.3075627775211328_dp
           V=0.2262439516632620_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5586001195731895_dp
           B=0.3382311089826877_dp
           V=0.2277874557231869_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5872229902021319_dp
           B=0.3681108834741399_dp
           V=0.2287854314454994_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6144258616235123_dp
           B=0.3970397446872839_dp
           V=0.2293268499615575_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2951676508064861_dp
           B=0.2867499538750441_dp*0.1_dp
           V=0.1912628201529828_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3335085485472725_dp
           B=0.5867879341903510_dp*0.1_dp
           V=0.1992499672238701_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3709561760636381_dp
           B=0.8961099205022284_dp*0.1_dp
           V=0.2061275533454027_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4074722861667498_dp
           B=0.1211627927626297_dp
           V=0.2119318215968572_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4429923648839117_dp
           B=0.1530748903554898_dp
           V=0.2167416581882652_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4774428052721736_dp
           B=0.1851176436721877_dp
           V=0.2206430730516600_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5107446539535904_dp
           B=0.2170829107658179_dp
           V=0.2237186938699523_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5428151370542935_dp
           B=0.2487786689026271_dp
           V=0.2260480075032884_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5735699292556964_dp
           B=0.2800239952795016_dp
           V=0.2277098884558542_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6029253794562866_dp
           B=0.3106445702878119_dp
           V=0.2287845715109671_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6307998987073145_dp
           B=0.3404689500841194_dp
           V=0.2293547268236294_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3752652273692719_dp
           B=0.2997145098184479_dp*0.1_dp
           V=0.2056073839852528_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4135383879344028_dp
           B=0.6086725898678011_dp*0.1_dp
           V=0.2114235865831876_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4506113885153907_dp
           B=0.9238849548435643_dp*0.1_dp
           V=0.2163175629770551_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4864401554606072_dp
           B=0.1242786603851851_dp
           V=0.2203392158111650_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5209708076611709_dp
           B=0.1563086731483386_dp
           V=0.2235473176847839_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5541422135830122_dp
           B=0.1882696509388506_dp
           V=0.2260024141501235_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5858880915113817_dp
           B=0.2199672979126059_dp
           V=0.2277675929329182_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6161399390603444_dp
           B=0.2512165482924867_dp
           V=0.2289102112284834_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6448296482255090_dp
           B=0.2818368701871888_dp
           V=0.2295027954625118_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4544796274917948_dp
           B=0.3088970405060312_dp*0.1_dp
           V=0.2161281589879992_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4919389072146628_dp
           B=0.6240947677636835_dp*0.1_dp
           V=0.2201980477395102_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5279313026985183_dp
           B=0.9430706144280313_dp*0.1_dp
           V=0.2234952066593166_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5624169925571135_dp
           B=0.1263547818770374_dp
           V=0.2260540098520838_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5953484627093287_dp
           B=0.1583430788822594_dp
           V=0.2279157981899988_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6266730715339185_dp
           B=0.1900748462555988_dp
           V=0.2291296918565571_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6563363204278871_dp
           B=0.2213599519592567_dp
           V=0.2297533752536649_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5314574716585696_dp
           B=0.3152508811515374_dp*0.1_dp
           V=0.2234927356465995_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5674614932298185_dp
           B=0.6343865291465561_dp*0.1_dp
           V=0.2261288012985219_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6017706004970264_dp
           B=0.9551503504223951_dp*0.1_dp
           V=0.2280818160923688_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6343471270264178_dp
           B=0.1275440099801196_dp
           V=0.2293773295180159_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6651494599127802_dp
           B=0.1593252037671960_dp
           V=0.2300528767338634_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6050184986005704_dp
           B=0.3192538338496105_dp*0.1_dp
           V=0.2281893855065666_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6390163550880400_dp
           B=0.6402824353962306_dp*0.1_dp
           V=0.2295720444840727_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6711199107088448_dp
           B=0.9609805077002909_dp*0.1_dp
           V=0.2303227649026753_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6741354429572275_dp
           B=0.3211853196273233_dp*0.1_dp
           V=0.2304831913227114_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
         end subroutine LD4802

         subroutine LD5294(X,Y,Z,W,N)
           real(dp), intent(out) :: X(5294)
           real(dp), intent(out) :: Y(5294)
           real(dp), intent(out) :: Z(5294)
           real(dp), intent(out) :: W(5294)
           integer, intent(out) :: N
           real(dp) :: A,B,V

           N=1
           V=0.9080510764308163_dp*0.0001_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.2084824361987793_dp*0.001_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2303261686261450_dp*0.1_dp
           V=0.5011105657239616_dp*0.0001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3757208620162394_dp*0.1_dp
           V=0.5942520409683854_dp*0.0001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5821912033821852_dp*0.1_dp
           V=0.9564394826109721_dp*0.0001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.8403127529194872_dp*0.1_dp
           V=0.1185530657126338_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1122927798060578_dp
           V=0.1364510114230331_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1420125319192987_dp
           V=0.1505828825605415_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1726396437341978_dp
           V=0.1619298749867023_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2038170058115696_dp
           V=0.1712450504267789_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2352849892876508_dp
           V=0.1789891098164999_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2668363354312461_dp
           V=0.1854474955629795_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2982941279900452_dp
           V=0.1908148636673661_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3295002922087076_dp
           V=0.1952377405281833_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3603094918363593_dp
           V=0.1988349254282232_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3905857895173920_dp
           V=0.2017079807160050_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4202005758160837_dp
           V=0.2039473082709094_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4490310061597227_dp
           V=0.2056360279288953_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4769586160311491_dp
           V=0.2068525823066865_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5038679887049750_dp
           V=0.2076724877534488_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5296454286519961_dp
           V=0.2081694278237885_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5541776207164850_dp
           V=0.2084157631219326_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5990467321921213_dp
           V=0.2084381531128593_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6191467096294587_dp
           V=0.2083476277129307_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6375251212901849_dp
           V=0.2082686194459732_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6540514381131168_dp
           V=0.2082475686112415_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6685899064391510_dp
           V=0.2083139860289915_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6810013009681648_dp
           V=0.2084745561831237_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6911469578730340_dp
           V=0.2087091313375890_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6988956915141736_dp
           V=0.2089718413297697_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7041335794868720_dp
           V=0.2092003303479793_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7067754398018567_dp
           V=0.2093336148263241_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3840368707853623_dp*0.1_dp
           V=0.7591708117365267_dp*0.0001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.9835485954117399_dp*0.1_dp
           V=0.1083383968169186_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1665774947612998_dp
           V=0.1403019395292510_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2405702335362910_dp
           V=0.1615970179286436_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3165270770189046_dp
           V=0.1771144187504911_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3927386145645443_dp
           V=0.1887760022988168_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4678825918374656_dp
           V=0.1973474670768214_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5408022024266935_dp
           V=0.2033787661234659_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6104967445752438_dp
           V=0.2072343626517331_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6760910702685738_dp
           V=0.2091177834226918_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6655644120217392_dp*0.1_dp
           B=0.1936508874588424_dp*0.1_dp
           V=0.9316684484675566_dp*0.0001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.9446246161270182_dp*0.1_dp
           B=0.4252442002115869_dp*0.1_dp
           V=0.1116193688682976_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1242651925452509_dp
           B=0.6806529315354374_dp*0.1_dp
           V=0.1298623551559414_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1553438064846751_dp
           B=0.9560957491205369_dp*0.1_dp
           V=0.1450236832456426_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1871137110542670_dp
           B=0.1245931657452888_dp
           V=0.1572719958149914_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2192612628836257_dp
           B=0.1545385828778978_dp
           V=0.1673234785867195_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2515682807206955_dp
           B=0.1851004249723368_dp
           V=0.1756860118725188_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2838535866287290_dp
           B=0.2160182608272384_dp
           V=0.1826776290439367_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3159578817528521_dp
           B=0.2470799012277111_dp
           V=0.1885116347992865_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3477370882791392_dp
           B=0.2781014208986402_dp
           V=0.1933457860170574_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3790576960890540_dp
           B=0.3089172523515731_dp
           V=0.1973060671902064_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4097938317810200_dp
           B=0.3393750055472244_dp
           V=0.2004987099616311_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4398256572859637_dp
           B=0.3693322470987730_dp
           V=0.2030170909281499_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4690384114718480_dp
           B=0.3986541005609877_dp
           V=0.2049461460119080_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4973216048301053_dp
           B=0.4272112491408562_dp
           V=0.2063653565200186_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5245681526132446_dp
           B=0.4548781735309936_dp
           V=0.2073507927381027_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5506733911803888_dp
           B=0.4815315355023251_dp
           V=0.2079764593256122_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5755339829522475_dp
           B=0.5070486445801855_dp
           V=0.2083150534968778_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1305472386056362_dp
           B=0.2284970375722366_dp*0.1_dp
           V=0.1262715121590664_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1637327908216477_dp
           B=0.4812254338288384_dp*0.1_dp
           V=0.1414386128545972_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1972734634149637_dp
           B=0.7531734457511935_dp*0.1_dp
           V=0.1538740401313898_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2308694653110130_dp
           B=0.1039043639882017_dp
           V=0.1642434942331432_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2643899218338160_dp
           B=0.1334526587117626_dp
           V=0.1729790609237496_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2977171599622171_dp
           B=0.1636414868936382_dp
           V=0.1803505190260828_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3307293903032310_dp
           B=0.1942195406166568_dp
           V=0.1865475350079657_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3633069198219073_dp
           B=0.2249752879943753_dp
           V=0.1917182669679069_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3953346955922727_dp
           B=0.2557218821820032_dp
           V=0.1959851709034382_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4267018394184914_dp
           B=0.2862897925213193_dp
           V=0.1994529548117882_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4573009622571704_dp
           B=0.3165224536636518_dp
           V=0.2022138911146548_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4870279559856109_dp
           B=0.3462730221636496_dp
           V=0.2043518024208592_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5157819581450322_dp
           B=0.3754016870282835_dp
           V=0.2059450313018110_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5434651666465393_dp
           B=0.4037733784993613_dp
           V=0.2070685715318472_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5699823887764627_dp
           B=0.4312557784139123_dp
           V=0.2077955310694373_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5952403350947741_dp
           B=0.4577175367122110_dp
           V=0.2081980387824712_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2025152599210369_dp
           B=0.2520253617719557_dp*0.1_dp
           V=0.1521318610377956_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2381066653274425_dp
           B=0.5223254506119000_dp*0.1_dp
           V=0.1622772720185755_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2732823383651612_dp
           B=0.8060669688588620_dp*0.1_dp
           V=0.1710498139420709_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3080137692611118_dp
           B=0.1099335754081255_dp
           V=0.1785911149448736_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3422405614587601_dp
           B=0.1399120955959857_dp
           V=0.1850125313687736_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3758808773890420_dp
           B=0.1702977801651705_dp
           V=0.1904229703933298_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4088458383438932_dp
           B=0.2008799256601680_dp
           V=0.1949259956121987_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4410450550841152_dp
           B=0.2314703052180836_dp
           V=0.1986161545363960_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4723879420561312_dp
           B=0.2618972111375892_dp
           V=0.2015790585641370_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5027843561874343_dp
           B=0.2920013195600270_dp
           V=0.2038934198707418_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5321453674452458_dp
           B=0.3216322555190551_dp
           V=0.2056334060538251_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5603839113834030_dp
           B=0.3506456615934198_dp
           V=0.2068705959462289_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5874150706875146_dp
           B=0.3789007181306267_dp
           V=0.2076753906106002_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6131559381660038_dp
           B=0.4062580170572782_dp
           V=0.2081179391734803_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2778497016394506_dp
           B=0.2696271276876226_dp*0.1_dp
           V=0.1700345216228943_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3143733562261912_dp
           B=0.5523469316960465_dp*0.1_dp
           V=0.1774906779990410_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3501485810261827_dp
           B=0.8445193201626464_dp*0.1_dp
           V=0.1839659377002642_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3851430322303653_dp
           B=0.1143263119336083_dp
           V=0.1894987462975169_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4193013979470415_dp
           B=0.1446177898344475_dp
           V=0.1941548809452595_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4525585960458567_dp
           B=0.1751165438438091_dp
           V=0.1980078427252384_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4848447779622947_dp
           B=0.2056338306745660_dp
           V=0.2011296284744488_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5160871208276894_dp
           B=0.2359965487229226_dp
           V=0.2035888456966776_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5462112185696926_dp
           B=0.2660430223139146_dp
           V=0.2054516325352142_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5751425068101757_dp
           B=0.2956193664498032_dp
           V=0.2067831033092635_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6028073872853596_dp
           B=0.3245763905312779_dp
           V=0.2076485320284876_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6291338275278409_dp
           B=0.3527670026206972_dp
           V=0.2081141439525255_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3541797528439391_dp
           B=0.2823853479435550_dp*0.1_dp
           V=0.1834383015469222_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3908234972074657_dp
           B=0.5741296374713106_dp*0.1_dp
           V=0.1889540591777677_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4264408450107590_dp
           B=0.8724646633650199_dp*0.1_dp
           V=0.1936677023597375_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4609949666553286_dp
           B=0.1175034422915616_dp
           V=0.1976176495066504_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4944389496536006_dp
           B=0.1479755652628428_dp
           V=0.2008536004560983_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5267194884346086_dp
           B=0.1784740659484352_dp
           V=0.2034280351712291_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5577787810220990_dp
           B=0.2088245700431244_dp
           V=0.2053944466027758_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5875563763536670_dp
           B=0.2388628136570763_dp
           V=0.2068077642882360_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6159910016391269_dp
           B=0.2684308928769185_dp
           V=0.2077250949661599_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6430219602956268_dp
           B=0.2973740761960252_dp
           V=0.2082062440705320_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4300647036213646_dp
           B=0.2916399920493977_dp*0.1_dp
           V=0.1934374486546626_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4661486308935531_dp
           B=0.5898803024755659_dp*0.1_dp
           V=0.1974107010484300_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5009658555287261_dp
           B=0.8924162698525409_dp*0.1_dp
           V=0.2007129290388658_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5344824270447704_dp
           B=0.1197185199637321_dp
           V=0.2033736947471293_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5666575997416371_dp
           B=0.1502300756161382_dp
           V=0.2054287125902493_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5974457471404752_dp
           B=0.1806004191913564_dp
           V=0.2069184936818894_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6267984444116886_dp
           B=0.2106621764786252_dp
           V=0.2078883689808782_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6546664713575417_dp
           B=0.2402526932671914_dp
           V=0.2083886366116359_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5042711004437253_dp
           B=0.2982529203607657_dp*0.1_dp
           V=0.2006593275470817_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5392127456774380_dp
           B=0.6008728062339922_dp*0.1_dp
           V=0.2033728426135397_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5726819437668618_dp
           B=0.9058227674571398_dp*0.1_dp
           V=0.2055008781377608_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6046469254207278_dp
           B=0.1211219235803400_dp
           V=0.2070651783518502_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6350716157434952_dp
           B=0.1515286404791580_dp
           V=0.2080953335094320_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6639177679185454_dp
           B=0.1816314681255552_dp
           V=0.2086284998988521_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5757276040972253_dp
           B=0.3026991752575440_dp*0.1_dp
           V=0.2055549387644668_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6090265823139755_dp
           B=0.6078402297870770_dp*0.1_dp
           V=0.2071871850267654_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6406735344387661_dp
           B=0.9135459984176636_dp*0.1_dp
           V=0.2082856600431965_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6706397927793709_dp
           B=0.1218024155966590_dp
           V=0.2088705858819358_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6435019674426665_dp
           B=0.3052608357660639_dp*0.1_dp
           V=0.2083995867536322_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6747218676375681_dp
           B=0.6112185773983089_dp*0.1_dp
           V=0.2090509712889637_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
         end subroutine LD5294

         subroutine LD5810(X,Y,Z,W,N)
           real(dp), intent(out) :: X(5810)
           real(dp), intent(out) :: Y(5810)
           real(dp), intent(out) :: Z(5810)
           real(dp), intent(out) :: W(5810)
           integer, intent(out) :: N

           real(dp) :: A,B,V

           N=1
           V=0.9735347946175486_dp*0.00001_dp
           Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.1907581241803167_dp*0.001_dp
           Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
           V=0.1901059546737578_dp*0.001_dp
           Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1182361662400277_dp*0.1_dp
           V=0.3926424538919212_dp*0.0001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3062145009138958_dp*0.1_dp
           V=0.6667905467294382_dp*0.0001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5329794036834243_dp*0.1_dp
           V=0.8868891315019135_dp*0.0001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7848165532862220_dp*0.1_dp
           V=0.1066306000958872_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1054038157636201_dp
           V=0.1214506743336128_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1335577797766211_dp
           V=0.1338054681640871_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1625769955502252_dp
           V=0.1441677023628504_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1921787193412792_dp
           V=0.1528880200826557_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2221340534690548_dp
           V=0.1602330623773609_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2522504912791132_dp
           V=0.1664102653445244_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2823610860679697_dp
           V=0.1715845854011323_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3123173966267560_dp
           V=0.1758901000133069_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3419847036953789_dp
           V=0.1794382485256736_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3712386456999758_dp
           V=0.1823238106757407_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3999627649876828_dp
           V=0.1846293252959976_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4280466458648093_dp
           V=0.1864284079323098_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4553844360185711_dp
           V=0.1877882694626914_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4818736094437834_dp
           V=0.1887716321852025_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5074138709260629_dp
           V=0.1894381638175673_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5319061304570707_dp
           V=0.1898454899533629_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5552514978677286_dp
           V=0.1900497929577815_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5981009025246183_dp
           V=0.1900671501924092_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6173990192228116_dp
           V=0.1899837555533510_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6351365239411131_dp
           V=0.1899014113156229_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6512010228227200_dp
           V=0.1898581257705106_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6654758363948120_dp
           V=0.1898804756095753_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6778410414853370_dp
           V=0.1899793610426402_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6881760887484110_dp
           V=0.1901464554844117_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6963645267094598_dp
           V=0.1903533246259542_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7023010617153579_dp
           V=0.1905556158463228_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.7059004636628753_dp
           V=0.1907037155663528_dp*0.001_dp
           Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3552470312472575_dp*0.1_dp
           V=0.5992997844249967_dp*0.0001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.9151176620841283_dp*0.1_dp
           V=0.9749059382456978_dp*0.0001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1566197930068980_dp
           V=0.1241680804599158_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2265467599271907_dp
           V=0.1437626154299360_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2988242318581361_dp
           V=0.1584200054793902_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3717482419703886_dp
           V=0.1694436550982744_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4440094491758889_dp
           V=0.1776617014018108_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5145337096756642_dp
           V=0.1836132434440077_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5824053672860230_dp
           V=0.1876494727075983_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6468283961043370_dp
           V=0.1899906535336482_dp*0.001_dp
           Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6095964259104373_dp*0.1_dp
           B=0.1787828275342931_dp*0.1_dp
           V=0.8143252820767350_dp*0.0001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.8811962270959388_dp*0.1_dp
           B=0.3953888740792096_dp*0.1_dp
           V=0.9998859890887728_dp*0.0001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1165936722428831_dp
           B=0.6378121797722990_dp*0.1_dp
           V=0.1156199403068359_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1460232857031785_dp
           B=0.8985890813745037_dp*0.1_dp
           V=0.1287632092635513_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1761197110181755_dp
           B=0.1172606510576162_dp
           V=0.1398378643365139_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2066471190463718_dp
           B=0.1456102876970995_dp
           V=0.1491876468417391_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2374076026328152_dp
           B=0.1746153823011775_dp
           V=0.1570855679175456_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2682305474337051_dp
           B=0.2040383070295584_dp
           V=0.1637483948103775_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2989653312142369_dp
           B=0.2336788634003698_dp
           V=0.1693500566632843_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3294762752772209_dp
           B=0.2633632752654219_dp
           V=0.1740322769393633_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3596390887276086_dp
           B=0.2929369098051601_dp
           V=0.1779126637278296_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3893383046398812_dp
           B=0.3222592785275512_dp
           V=0.1810908108835412_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4184653789358347_dp
           B=0.3512004791195743_dp
           V=0.1836529132600190_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4469172319076166_dp
           B=0.3796385677684537_dp
           V=0.1856752841777379_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4745950813276976_dp
           B=0.4074575378263879_dp
           V=0.1872270566606832_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5014034601410262_dp
           B=0.4345456906027828_dp
           V=0.1883722645591307_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5272493404551239_dp
           B=0.4607942515205134_dp
           V=0.1891714324525297_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5520413051846366_dp
           B=0.4860961284181720_dp
           V=0.1896827480450146_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5756887237503077_dp
           B=0.5103447395342790_dp
           V=0.1899628417059528_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1225039430588352_dp
           B=0.2136455922655793_dp*0.1_dp
           V=0.1123301829001669_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1539113217321372_dp
           B=0.4520926166137188_dp*0.1_dp
           V=0.1253698826711277_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1856213098637712_dp
           B=0.7086468177864818_dp*0.1_dp
           V=0.1366266117678531_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2174998728035131_dp
           B=0.9785239488772918_dp*0.1_dp
           V=0.1462736856106918_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2494128336938330_dp
           B=0.1258106396267210_dp
           V=0.1545076466685412_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2812321562143480_dp
           B=0.1544529125047001_dp
           V=0.1615096280814007_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3128372276456111_dp
           B=0.1835433512202753_dp
           V=0.1674366639741759_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3441145160177973_dp
           B=0.2128813258619585_dp
           V=0.1724225002437900_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3749567714853510_dp
           B=0.2422913734880829_dp
           V=0.1765810822987288_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4052621732015610_dp
           B=0.2716163748391453_dp
           V=0.1800104126010751_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4349335453522385_dp
           B=0.3007127671240280_dp
           V=0.1827960437331284_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4638776641524965_dp
           B=0.3294470677216479_dp
           V=0.1850140300716308_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4920046410462687_dp
           B=0.3576932543699155_dp
           V=0.1867333507394938_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5192273554861704_dp
           B=0.3853307059757764_dp
           V=0.1880178688638289_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5454609081136522_dp
           B=0.4122425044452694_dp
           V=0.1889278925654758_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5706220661424140_dp
           B=0.4383139587781027_dp
           V=0.1895213832507346_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5946286755181518_dp
           B=0.4634312536300553_dp
           V=0.1898548277397420_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.1905370790924295_dp
           B=0.2371311537781979_dp*0.1_dp
           V=0.1349105935937341_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2242518717748009_dp
           B=0.4917878059254806_dp*0.1_dp
           V=0.1444060068369326_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2577190808025936_dp
           B=0.7595498960495142_dp*0.1_dp
           V=0.1526797390930008_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2908724534927187_dp
           B=0.1036991083191100_dp
           V=0.1598208771406474_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3236354020056219_dp
           B=0.1321348584450234_dp
           V=0.1659354368615331_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3559267359304543_dp
           B=0.1610316571314789_dp
           V=0.1711279910946440_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3876637123676956_dp
           B=0.1901912080395707_dp
           V=0.1754952725601440_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4187636705218842_dp
           B=0.2194384950137950_dp
           V=0.1791247850802529_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4491449019883107_dp
           B=0.2486155334763858_dp
           V=0.1820954300877716_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4787270932425445_dp
           B=0.2775768931812335_dp
           V=0.1844788524548449_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5074315153055574_dp
           B=0.3061863786591120_dp
           V=0.1863409481706220_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5351810507738336_dp
           B=0.3343144718152556_dp
           V=0.1877433008795068_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5619001025975381_dp
           B=0.3618362729028427_dp
           V=0.1887444543705232_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5875144035268046_dp
           B=0.3886297583620408_dp
           V=0.1894009829375006_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6119507308734495_dp
           B=0.4145742277792031_dp
           V=0.1897683345035198_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2619733870119463_dp
           B=0.2540047186389353_dp*0.1_dp
           V=0.1517327037467653_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.2968149743237949_dp
           B=0.5208107018543989_dp*0.1_dp
           V=0.1587740557483543_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3310451504860488_dp
           B=0.7971828470885599_dp*0.1_dp
           V=0.1649093382274097_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3646215567376676_dp
           B=0.1080465999177927_dp
           V=0.1701915216193265_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3974916785279360_dp
           B=0.1368413849366629_dp
           V=0.1746847753144065_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4295967403772029_dp
           B=0.1659073184763559_dp
           V=0.1784555512007570_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4608742854473447_dp
           B=0.1950703730454614_dp
           V=0.1815687562112174_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4912598858949903_dp
           B=0.2241721144376724_dp
           V=0.1840864370663302_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5206882758945558_dp
           B=0.2530655255406489_dp
           V=0.1860676785390006_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5490940914019819_dp
           B=0.2816118409731066_dp
           V=0.1875690583743703_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5764123302025542_dp
           B=0.3096780504593238_dp
           V=0.1886453236347225_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6025786004213506_dp
           B=0.3371348366394987_dp
           V=0.1893501123329645_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6275291964794956_dp
           B=0.3638547827694396_dp
           V=0.1897366184519868_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3348189479861771_dp
           B=0.2664841935537443_dp*0.1_dp
           V=0.1643908815152736_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.3699515545855295_dp
           B=0.5424000066843495_dp*0.1_dp
           V=0.1696300350907768_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4042003071474669_dp
           B=0.8251992715430854_dp*0.1_dp
           V=0.1741553103844483_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4375320100182624_dp
           B=0.1112695182483710_dp
           V=0.1780015282386092_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4699054490335947_dp
           B=0.1402964116467816_dp
           V=0.1812116787077125_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5012739879431952_dp
           B=0.1694275117584291_dp
           V=0.1838323158085421_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5315874883754966_dp
           B=0.1985038235312689_dp
           V=0.1859113119837737_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5607937109622117_dp
           B=0.2273765660020893_dp
           V=0.1874969220221698_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5888393223495521_dp
           B=0.2559041492849764_dp
           V=0.1886375612681076_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6156705979160163_dp
           B=0.2839497251976899_dp
           V=0.1893819575809276_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6412338809078123_dp
           B=0.3113791060500690_dp
           V=0.1897794748256767_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4076051259257167_dp
           B=0.2757792290858463_dp*0.1_dp
           V=0.1738963926584846_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4423788125791520_dp
           B=0.5584136834984293_dp*0.1_dp
           V=0.1777442359873466_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4760480917328258_dp
           B=0.8457772087727143_dp*0.1_dp
           V=0.1810010815068719_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5085838725946297_dp
           B=0.1135975846359248_dp
           V=0.1836920318248129_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5399513637391218_dp
           B=0.1427286904765053_dp
           V=0.1858489473214328_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5701118433636380_dp
           B=0.1718112740057635_dp
           V=0.1875079342496592_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5990240530606021_dp
           B=0.2006944855985351_dp
           V=0.1887080239102310_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6266452685139695_dp
           B=0.2292335090598907_dp
           V=0.1894905752176822_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6529320971415942_dp
           B=0.2572871512353714_dp
           V=0.1898991061200695_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.4791583834610126_dp
           B=0.2826094197735932_dp*0.1_dp
           V=0.1809065016458791_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5130373952796940_dp
           B=0.5699871359683649_dp*0.1_dp
           V=0.1836297121596799_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5456252429628476_dp
           B=0.8602712528554394_dp*0.1_dp
           V=0.1858426916241869_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5768956329682385_dp
           B=0.1151748137221281_dp
           V=0.1875654101134641_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6068186944699046_dp
           B=0.1442811654136362_dp
           V=0.1888240751833503_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6353622248024907_dp
           B=0.1731930321657680_dp
           V=0.1896497383866979_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6624927035731797_dp
           B=0.2017619958756061_dp
           V=0.1900775530219121_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5484933508028488_dp
           B=0.2874219755907391_dp*0.1_dp
           V=0.1858525041478814_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.5810207682142106_dp
           B=0.5778312123713695_dp*0.1_dp
           V=0.1876248690077947_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6120955197181352_dp
           B=0.8695262371439526_dp*0.1_dp
           V=0.1889404439064607_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6416944284294319_dp
           B=0.1160893767057166_dp
           V=0.1898168539265290_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6697926391731260_dp
           B=0.1450378826743251_dp
           V=0.1902779940661772_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6147594390585488_dp
           B=0.2904957622341456_dp*0.1_dp
           V=0.1890125641731815_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6455390026356783_dp
           B=0.5823809152617197_dp*0.1_dp
           V=0.1899434637795751_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6747258588365477_dp
           B=0.8740384899884715_dp*0.1_dp
           V=0.1904520856831751_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           A=0.6772135750395347_dp
           B=0.2919946135808105_dp*0.1_dp
           V=0.1905534498734563_dp*0.001_dp
           Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
           N=N-1
         end subroutine LD5810

       end module lebedev

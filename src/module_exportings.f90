!==============================================================================!
! MODULE ResultExportings                                                      !
!                                                                              !
!    This module join the subroutines to evaluate the density dependent        !
! interaction and its hamiltonian and also the diagonalizing basis for the     !
! single particle hamiltonian with the j_z observable.                         !
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
MODULE ResultExportings

use Constants
use MathMethods
use Basis
use Hamiltonian
use DensityDep
use Gradient

implicit none
PUBLIC

!! Attributes



!! Methods

CONTAINS




END MODULE ResultExportings

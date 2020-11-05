!*****************************************************************************************
! Neutron intensity as a function of distance for an extremely simplified case
! of monoenergetic neutrons where absorption is the ONLY mechanism modeled 
! Code also computes the resultant energy release. 
!
! This file contains the following program units
! MODULE     neutron0_mod
! PROGRAM    ss_neutron0
! SUBROUTINE setup_neutron0
! SUBROUTINE solve_neutron0
! SUBROUTINE finish_neutron0
!*****************************************************************************************
   

!*****************************************************************************************
MODULE neutron0_mod
!*****************************************************************************************
  IMPLICIT NONE
  PUBLIC

  REAL    :: Cnq      ! Conversion factor between absorbed Neutrons and Joules
  REAL    :: X(0:3)   ! X locations within the wall (distance from X=0)
  REAL    :: I(0:3)   ! Intensity (Neutrons/unit-area) at the four X locations
  REAL    :: Q(3)     ! Energy release within wall between point X(i) and X(i-1)
  REAL    :: T(3)     ! Temperature within wall between point X(i) and X(i-1)
  REAL    :: RHS(3)   ! Right-hand-side for assignment of Q values
  REAL    :: TCS(3)   ! Total Cross Section values
  REAL    :: R(3)     ! Residual array 
  REAL    :: R_max    ! maximum value in the Residual array
  INTEGER :: uout     ! Output unit for fortran write statements

END MODULE neutron0_mod
! ========================


!*****************************************************************************************
PROGRAM ss_neutron0
!*****************************************************************************************
  IMPLICIT NONE
!
! Problem Setup
!  CALL setup_neutron0
!
! Solve problem
!  CALL solve_neutron0
!
! Final output
!  CALL finish_neutron0
!
END PROGRAM ss_neutron0
! ========================


!*****************************************************************************************
SUBROUTINE setup_neutron0
!*****************************************************************************************
  USE neutron0_mod
  IMPLICIT NONE

! Problem Setup

  X(0)   = 0.0      ! Mesh locations chosen to align with faces in ss_1dcon code
  X(1)   = 0.05
  X(2)   = 0.15
  X(3)   = 0.20

  Cnq    = 4000.
  I(0)   = 1.0      ! Neutron Intensity at X=0
  T      = 400.0    ! Temperature array set to 400.0
  CALL solve_neutron0 ! initialize q with initial T
  uout   = 9        ! standard output goes to unit 6
  if(uout.ne.6)   OPEN(unit=uout, file='neutron0_out', status='unknown')

END SUBROUTINE setup_neutron0
! ========================



!*****************************************************************************************
SUBROUTINE solve_neutron0
!*****************************************************************************************
  USE neutron0_mod
  IMPLICIT NONE
  REAL    :: dx      ! X(j)-X(j-1)
  INTEGER :: j       ! index
  
  ! Solve problem 
  DO j=1,3   ! Loop over the three X locations in the wall

     ! Temperature dependant Total Cross Section (I made this polynomial up)
     TCS(j) = 1.042-0.00255*T(j) + 2.36e-5*T(j)**2 - 4.08e-8*T(j)**3 + 2.86e-17*T(j)**6
     IF(T(j).gt.705.) TCS(j) = 0.189

     ! Intensity
     dx = X(j)-X(j-1)
     I(j) = I(j-1) * exp(-TCS(j)*dx)

     ! Energy Release
      Q(j) = -Cnq*(I(j)-I(j-1))/dx

  END DO
    
END SUBROUTINE solve_neutron0
! ========================


!*****************************************************************************************
SUBROUTINE finish_neutron0
!*****************************************************************************************
  USE neutron0_mod
  IMPLICIT NONE
  INTEGER :: j, jmax=3
  
! Final output
  write(uout,*) "  X            I            Q            T            TCS"
  write(uout,7) X(0), I(0), T(1)
  DO j=1,jmax   ! Loop over the wall locations, dealing with staggered grid
     if(j.eq.1) write(uout,6) (X(j)+X(j-1))/2.,  Q(j), TCS(j)
     if(j.gt.1) write(uout,8)  X(j-1), I(j-1)
     if(j.gt.1 .and. j.lt.jmax) write(uout,9) (X(j)+X(j-1))/2.,  Q(j), T(j), TCS(j)
     if(j.gt.1 .and. j.eq.jmax) write(uout,6) (X(j)+X(j-1))/2.,  Q(j), TCS(j)
     if(j.eq.jmax) write(uout,7) X(j), I(j), T(j)
  END DO  

  write(uout,*) " Finished running program ss_neutron "
  if(uout.ne.6)   CLOSE(unit=uout)

6 format( es11.3, 5x,'-',9x, es11.3, 5x,"-", 9x, es11.3) 
7 format( es11.3, 2x, es11.3, 5x,"-",9x,es11.3 )
8 format( es11.3, 2x, es11.3, 5x,"-",12x,"-" )
9 format( es11.3, 5x,'-',9x, es11.3, 2x, es11.3, 2x, es11.3) 

END SUBROUTINE finish_neutron0
! ========================

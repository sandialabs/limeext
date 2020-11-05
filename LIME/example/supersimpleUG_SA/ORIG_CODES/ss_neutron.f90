PROGRAM ss_neutron
! ==============================================================================
! Neutron intensity as a function of distance for an extremely simplified case
! of monoenergitic neutrons where absorption is the ONLY mechanism modeled 
! Code also computes the resultant energy release. 
! ==============================================================================

  IMPLICIT NONE
  INTEGER :: j       ! mesh index
  INTEGER :: jmax=3  ! max index value
  REAL    :: dx      ! X(j)-X(j-1)
  REAL    :: Cnq     ! Conversion factor between absorbed Neutrons and Joules
  REAL    :: X(0:3)  ! X locations within the wall (distance from X=0)
  REAL    :: I(0:3)  ! Intensity (Neutrons/unit-area) at the four X locations
  REAL    :: Q(3)    ! Energy release within wall between point X(i) and X(i-1)
  REAL    :: T(3)    ! Temperature within wall between point X(i) and X(i-1)
  REAL    :: TCS(3)  ! Total Cross Section values
  INTEGER :: uout    ! Output unit for fortran write statements

! ==============================================================================
! Problem Setup
  Cnq    = 4000.
  I(0)   = 1.0      ! Neutron Intensity at X=0
  T      = 400.0

  X(0)   = 0.0      ! Mesh locations chosen to align with faces in ss_1dcon code
  X(1)   = 0.05
  X(2)   = 0.15
  X(3)   = 0.20

  write(6,*) "      index   X(j)           I(j)            Q(j)"
  write(6,*) 0, X(0), I(0)
! ==============================================================================
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

!9 format(i4, 2x, es11.3, 2x, es11.3, 2x, es11.3)
! ==============================================================================
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
! ==============================================================================
END PROGRAM ss_neutron


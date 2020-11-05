!*****************************************************************************************
! Transient thermal response of a conductive thin wall subject to rad. BC on
! the left side and a convective BC on the right side.(Lumped mass approximation)
! Temporal discretization is fully implicit backwords Euler
!
! This file contains the following program units
! MODULE     thinwall_mod
! PROGRAM    ss_thinwall
! SUBROUTINE setup_thinwall
! SUBROUTINE solve_thinwall
! SUBROUTINE finish_thinwall
! SUBROUTINE timestep_thinwall
! SUBROUTINE update_thinwall
! SUBROUTINE residual_thinwall
! SUBROUTINE residual2_thinwall  (thread safe version)
! ****************************************************************************************

   
!*****************************************************************************************
MODULE thinwall_mod
!*****************************************************************************************
  IMPLICIT NONE
  PUBLIC

  REAL    :: rhoCpW ! density*specific heat*width of wall
  REAL    :: Tfluid ! Fluid Temperature on the right-side boundary
  REAL    :: hconv  ! Convective heat trasnver coeficient on the right-side boundary
  REAL    :: Trad   ! Radiation Temperature seen by the left-side boundary
  REAL    :: time   ! time (sec)
  REAL    :: dt     ! time step
  REAL    :: T(1)   ! temperature (lumped mass approximation) of the wall
  REAL    :: T0     ! temperature at previous time step
  REAL    :: R(1)   ! Residual array 
  REAL    :: dRdT   ! derivative of Residual wrt Temperature 
  REAL    :: delT   ! Newton-based temperature update
  REAL    :: coef1  ! coeficient used in energy equation
  REAL    :: coef2  ! coeficient used in energy equation
  REAL    :: coef3  ! coeficient used in energy equation
  REAL    :: tol    ! convergence criteria tolerance
  REAL    :: qrad   ! radiation heat trasnfer on the left side (per unit area)
  REAL    :: qfluid ! convective heat transfer on the right side (per unit area)
  INTEGER :: n      ! time-step loop index
  INTEGER :: iter   ! iteration counter during iterative solve loop
  INTEGER :: uout   ! Output unit for fortran write statements
  INTEGER :: n_vars ! Number of unknowns and residual equations

END MODULE thinwall_mod
! ========================



!*****************************************************************************************
!PROGRAM ss_thinwall
!*****************************************************************************************
!  IMPLICIT NONE
!
! Problem Setup
!  CALL setup_thinwall
!
! Solve problem
!  CALL solve_thinwall(25)
!
! Final output
!  CALL finish_thinwall
!
!END PROGRAM ss_thinwall
! ========================



!*****************************************************************************************
SUBROUTINE setup_thinwall
!*****************************************************************************************
  USE thinwall_mod
  IMPLICIT NONE

  rhoCpW = 100.
  Tfluid = 275.0
  hconv  = 10.
  Trad   = 400.
  time   = 0.0
  dt     = 2.0
  T(1)   = 300. 
  T0     = T(1)
  tol    = 2.e-3 

  coef1 = 5.67e-8*dt/rhoCpW  ! Stefan-Boltzmann constant = 5.57e-8 W/(m^2 K^4)
  coef2 = dt*hconv/rhoCpW
  coef3 = coef1/(1.0+coef2)

  n_vars = 1        ! # unknowns (variables) and corresponding residual equations
  uout   = 11       ! standard output goes to unit 11
  if(uout.ne.6)   OPEN(unit=uout, file='thinwall_out', status='unknown')
  write(uout,'(a)') " time     qrad        T          qfluid   iter    Resid" 
  write(uout,9) time, 0.0,T(1),0.0
9 format(f8.1,2x, es9.2, es12.3, 2x, es10.2, i4, 3x, es9.2)

END SUBROUTINE setup_thinwall
! ========================


!*****************************************************************************************
SUBROUTINE solve_thinwall(nstep)
!*****************************************************************************************
  IMPLICIT NONE
  INTEGER :: n, nstep

  DO n=1,nstep   ! Begining of the time-step loop

     CALL time_step_thinwall

     CALL update_thinwall

  END DO     ! End of the time-step loop

END SUBROUTINE solve_thinwall
! ========================


!*****************************************************************************************
SUBROUTINE finish_thinwall
!*****************************************************************************************
  USE thinwall_mod, only: uout
  IMPLICIT NONE

  write(uout,*) " Finished running program ss_thinwall "

END SUBROUTINE finish_thinwall
! ========================


!*****************************************************************************************
SUBROUTINE time_step_thinwall
!*****************************************************************************************
  USE thinwall_mod
  IMPLICIT NONE

     ! ------------------------------------------------------------------------
     ! Solve using a simple Newton-method-based iteration scheme
     DO iter = 1,15

        ! Compute Residual
        ! CALL residual_thinwall
        CALL residual2_thinwall(T(1),T0,Trad,Tfluid,coef1,coef2,coef3, R)
        ! Check for convergence
        if(abs(R(1)) .lt. tol) EXIT   ! exit out of iteration loop when convergence criteria met

        ! Derivative of the Residual with respect to T ( 1x1 jacobian)
        dRdT = 1.0 + 4.0*coef3*T(1)**3
        
        ! Newton update
        delT = -R(1) / dRdT

        ! write(uout,'(a,i3,a,3es11.3)') "DEBUG: iter=", iter, "  Res, delT, T+delT= ",R(1), delT, T+delT

        ! update to temperature estimate
        T(1) = T(1) + delT
     END DO
     ! ------------------------------------------------------------------------
     IF(abs(R(1)) .gt. tol) write(uout,*) "Warning: ThinWall Iteration did not converge!"
  

END SUBROUTINE time_step_thinwall
! ========================


!*****************************************************************************************
SUBROUTINE update_thinwall
!*****************************************************************************************
  USE thinwall_mod
  IMPLICIT NONE

     T0   = T(1)          
     time = time + dt
     qrad =  5.67e-8*(Trad**4 - T(1)**4)    ! Stefan-Boltzmann constant = 5.57e-8 W/(m^2 K^4)
     qfluid = hconv*(Tfluid-T(1))
     write(uout,9) time, qrad,T(1),qfluid,iter, abs(R)

9 format(f8.1,2x, es9.2, es12.3, 2x, es10.2, i4, 3x, es9.2)

END SUBROUTINE update_thinwall
! ========================


!*****************************************************************************************
SUBROUTINE residual_thinwall
!*****************************************************************************************
  USE thinwall_mod
  IMPLICIT NONE

  R(1) = T(1) + coef3*T(1)**4 - (T0 + coef1*Trad**4 + coef2*Tfluid)/(1.0+coef2)

END SUBROUTINE residual_thinwall
! ========================



!*****************************************************************************************
SUBROUTINE residual2_thinwall(T,T0,Trad,Tfluid,coef1,coef2,coef3, R)
!*****************************************************************************************
  IMPLICIT NONE
  REAL,INTENT(in)  :: T      ! temperature (lumped mass approximation) of the wall
  REAL,INTENT(in)  :: T0     ! temperature at previous time step
  REAL,INTENT(in)  :: Trad   ! Radiation Temperature seen by the left-side boundary
  REAL,INTENT(in)  :: Tfluid ! Fluid Temperature on the right-side boundary
  REAL,INTENT(in)  :: coef1  ! coeficient used in energy equation
  REAL,INTENT(in)  :: coef2  ! coeficient used in energy equation
  REAL,INTENT(in)  :: coef3  ! coeficient used in energy equation
  REAL,INTENT(out) :: R(1)   ! Residual 
! --------------------------------------------------------------------

  R(1) = T + coef3*T**4 - (T0 + coef1*Trad**4 + coef2*Tfluid)/(1.0+coef2)

END SUBROUTINE residual2_thinwall
! ========================



PROGRAM ss_thinwall
! ==============================================================================
! Transient thermal response of a conductive thin wall subject to rad. BC on
! the left side and a convective BC on the right side.(Lumped mass approximation)
! Temporal discretization is first-order Euler implicit
! ==============================================================================

  IMPLICIT NONE

  REAL    :: rhoCpW ! density*specific heat*width of wall
  REAL    :: Tfluid ! Fluid Temperature on the right-side boundary
  REAL    :: hconv  ! Convective heat trasnver coeficient on the right-side boundary
  REAL    :: Trad   ! Radiation Temperature seen by the left-side boundary
  REAL    :: time   ! time (sec)
  REAL    :: dt     ! time step
  REAL    :: T      ! temperature (lumped mass approximation) of the wall
  REAL    :: T0     ! temperature at previous time step
  REAL    :: R      ! Residual  
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
! ==============================================================================
! Problem Setup
  rhoCpW = 100.
  Tfluid = 275.0
  hconv  = 10.
  Trad   = 400.
  time   = 0.0
  dt     = 2.0
  T      = 300. 
  T0     = T
  tol    = 2.e-3 
  coef1 = 5.67e-8*dt/rhoCpW  ! Stefan-Boltzmann constant = 5.57e-8 W/(m^2 K^4)
  coef2 = dt*hconv/rhoCpW
  coef3 = coef1/(1.0+coef2)

  write(6,'(a)') " time     qrad        T          qfluid   iter    Resid" 
  write(6,9) time, 0.0,T,0.0
! ==============================================================================
! Solve problem

  DO n=1,25   ! Top of the time-step loop

     ! ------------------------------------------------------------------------
     ! Solve using a simple Newton-method-based iteration scheme
     DO iter = 1,15

        ! Compute Residual
        R = T + coef3*T**4 - (T0 + coef1*Trad**4 + coef2*Tfluid)/(1.0+coef2)

        ! Check for convergence
        if(abs(R) .lt. tol) EXIT   ! exit out of iteration loop when convergence criteria met

        ! Derivative of the Residual with respect to T ( 1x1 jacobian)
        dRdT = 1.0 + 4.0*coef3*T**3
        
        ! Newton update
        delT = -R / dRdT

        ! write(6,'(a,i3,a,3es11.3)') "DEBUG: iter=", iter, "  Res, delT, T+delT= ",R, delT, T+delT

        ! update to temperature estimate
        T = T + delT
     END DO
     ! ------------------------------------------------------------------------
     IF(abs(R) .gt. tol) write(6,*) "Warning: ThinWall Iteration did not converge!"

     T0   = T           
     time = time + dt
     qrad =  5.67e-8*(Trad**4 - T**4)    ! Stefan-Boltzmann constant = 5.57e-8 W/(m^2 K^4)
     qfluid = hconv*(Tfluid-T)
     write(6,9) time, qrad,T,qfluid,iter, abs(R)

  END DO     ! Bottom of the time-step loop

9 format(f6.1,2x, es9.2, es12.3, 2x, es10.2, i4, 3x, es8.2)
! ==============================================================================
! Final output
  write(6,*) " Finished running program ss_thinwall "
! ==============================================================================
END PROGRAM ss_thinwall


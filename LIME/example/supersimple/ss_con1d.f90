!*****************************************************************************************
! Transient 1D conduction, with Heat Src and with specified surface BCs
! Spatial discretization fixed as a 3 node finite difference approximation
! Temporal discretization is fully implicit first-order Euler
! 
! This file contains the following program units
! MODULE     con1d_mod
! PROGRAM    ss_con1d
! SUBROUTINE setup_con1d
! SUBROUTINE solve_con1d
! SUBROUTINE finish_con1d
! SUBROUTINE take_time_step_con1d
! SUBROUTINE s_substitution_con1d
! SUBROUTINE update_con1d
! SUBROUTINE residual_con1d
! SUBROUTINE jacobian_con1d  (not currently used, has not been tested and debugged)
! ****************************************************************************************


!*****************************************************************************************
MODULE con1d_mod
!*****************************************************************************************
  IMPLICIT NONE
  PUBLIC

  ! Geometry, mesh and material property related -------------------------------
  REAL    :: W      ! total width of 1-D domain
  REAL    :: dx     ! width of finite difference discretization
  REAL    :: k      ! thermal conductivity
  REAL    :: rhoCp  ! density*specific heat
  ! Boundary Condition related -------------------------------------------------
  REAL    :: qlt    ! specified heat flux at the left-side boundary
  REAL    :: qrt    ! computed heat flux at the right-side boundary
  REAL    :: h_c    ! convective heat trasnfer coeficient at the left-side boundary
  REAL    :: Tf     ! fluid temperature "seen" by the left-side boundary
  REAL    :: Tr     ! far-field radiation temperature "seen" by the right-side boundary
  ! State variables and/or source term related ---------------------------------
  REAL    :: time   ! time (sec)
  REAL    :: tmax   ! maximum time 
  REAL    :: dt     ! time step
  REAL    :: Qs(3)  ! heat source array
  REAL    :: T(3)   ! temperature array during solution iteration
  REAL    :: T0(3)  ! temperature array (converged) at previous time step
  ! Solver and Solution --------------------------------------------------------
  REAL    :: R(3)   ! Residual array 
  REAL    :: R_max  ! maximum value in the Residual array
  REAL    :: coef1  ! coeficient used in energy equation
  REAL    :: coef2  ! coeficient used in energy equation
  REAL    :: coef3  ! coeficient used in energy equation
  REAL    :: tol    ! convergence criteria tolerance
  REAL    :: ebal   ! normalized energy balance  (=1.0 at steady state)
  REAL    :: f      ! relaxation factor for SS iteration method
  INTEGER :: iter   ! iteration counter during iterative solve loop
  INTEGER :: uout   ! Output unit for fortran write statements
  INTEGER :: n_vars ! Number of unknowns and residual equations

END MODULE con1d_mod
! ========================



!*****************************************************************************************
!PROGRAM ss_con1d
!*****************************************************************************************
!  IMPLICIT NONE
!
! Problem Setup
!  CALL setup_con1d
!
! Solve problem
!  CALL solve_con1d(25)
!
! Final output
!  CALL finish_con1d
!
!END PROGRAM ss_con1d
! ========================


!*****************************************************************************************
SUBROUTINE setup_con1d
!*****************************************************************************************

  USE con1d_mod
  IMPLICIT NONE

  W    = 0.20
  dx   = W/2.0
  k    = 0.5
  rhoCp= 100.
  h_c  = 0.1        ! setting = 0 makes left side BC adiabatic
  Tf   = 350.0      ! left side fluid temperature
  Tr   = 300.0      ! right side radiation temperature
  time = 0.0
  tmax = 50.0
  dt   = 2.0
  Qs   = 3000.      ! heat source array initialized 
  T    = 300.       ! Temperature array initialized 
  T0   = T          ! 
  tol  = 2.e-3      ! 
  ebal = 0.0
  uout   = 8        ! standard output goes to unit 6
  n_vars = 3        ! # unknowns (variables) and corresponding residual equations
  if(uout.ne.6) OPEN(unit=uout, file='con1d_out', status='unknown')

  coef1 = dt*k/(rhoCp*dx*dx)
  coef2 = dt/(rhoCp*0.5*dx)
  coef3 = dt/rhoCp

  write(uout,'(a)') " time     qlt        T(1)        T(2)        T(3)         qrt       ebal      iter    R_max" 
  write(uout,9) time, qlt,T(1),T(2),T(3),qrt,ebal,iter,R_max
9 format(f8.1,2x, es9.2, 3es12.3, 2x, 2es10.2, 2x, i4, 3x, es9.2)

END SUBROUTINE setup_con1d
! ========================


!*****************************************************************************************
SUBROUTINE solve_con1d(nstep)
!*****************************************************************************************
  IMPLICIT NONE
  INTEGER :: n, nstep

  DO n=1,nstep   ! Begining of the time-step loop

     CALL take_time_step_con1d

     CALL update_con1d

  END DO     ! End of the time-step loop

END SUBROUTINE solve_con1d
! ========================



!*****************************************************************************************
SUBROUTINE finish_con1d
!*****************************************************************************************
  USE con1d_mod
  IMPLICIT NONE

  write(uout,*) " Finished running program ss_con1d "
  if(uout.ne.6)   CLOSE(unit=uout)

END SUBROUTINE finish_con1d
! ========================



!*****************************************************************************************
SUBROUTINE take_time_step_con1d
!*****************************************************************************************
  IMPLICIT NONE

  INTEGER :: METHOD   !Solver method flag: 1 = Successive substitution, else Newtons Method
  METHOD = 1

!  IF(METHOD .eq. 1) THEN
  
     CALL s_substitution_con1d

!  ELSE

!     CALL newton_con1d
  
!  ENDIF

END SUBROUTINE take_time_step_con1d
! ========================


!*****************************************************************************************
SUBROUTINE s_substitution_con1d
!*****************************************************************************************
  USE con1d_mod
  IMPLICIT NONE
  INTEGER i

     ! ------------------------------------------------------------------------
     ! Solve using succesive substitution iteration method (requires relaxation factor f)
     DO i = 1,500

        ! Update boundary heat fluxes
        qrt = 5.67e-8*(Tr**4 - T(3)**4)    ! Stefan-Boltzmann constant = 5.57e-8 W/(m^2 K^4)
        qlt = h_c*(Tf - T(1))

        f = 0.20  ! I've hardwired this value here. Other conditions would require tuning
        T(1) = f*(T0(1) + 2.*coef1*(-T(1) +   T(2)        ) + coef2*qlt + coef3*Qs(1)) + (1.-f)*T(1)
        T(2) = f*(T0(2) +    coef1*( T(1) -2.*T(2) + T(3) ) +             coef3*Qs(2)) + (1.-f)*T(2)
        T(3) = f*(T0(3) + 2.*coef1*(          T(2) - T(3) ) + coef2*qrt + coef3*Qs(3)) + (1.-f)*T(3)

        ! Compute residual vector and check for convergence
        CALL residual_con1d
        if(R_max .lt. tol) THEN   ! exit out of iteration loop when convergence criteria met
          ! write(uout,*) "CON1D Iteration converged in ",i," iterations."
          iter = iter + i
          EXIT
        endif

     END DO
     ! ------------------------------------------------------------------------
     IF(R_max .gt. tol) write(uout,*) "Warning: CON1D Iteration did not converge (",R_max," > ",tol,") !"

END SUBROUTINE s_substitution_con1d
! ========================


!*****************************************************************************************
SUBROUTINE update_con1d
!*****************************************************************************************
  USE con1d_mod
  IMPLICIT NONE
  REAL :: Qbar  ! volume weighted average heat source

  T0   = T           ! f90 vector operation updating T0 
  time = time + dt
  Qbar =  0.25*Qs(1)+0.5*Qs(2)+0.25*Qs(3)
  ebal = (-(qlt+qrt)/W)/Qbar   ! must equal 1.0 at steady state
  write(uout,9) time, qlt,T(1),T(2),T(3),qrt,ebal,iter,R_max

  ! Reset iter counter
  iter = 0

9 format(f8.1,2x, es9.2, 3es12.3, 2x, 2es10.2, 3x, i4, 2x, es9.2)

END SUBROUTINE update_con1d
! ========================



!*****************************************************************************************
SUBROUTINE residual_con1d
!*****************************************************************************************
  USE con1d_mod
  IMPLICIT NONE

  ! Compute residual vector and check for convergence
  qrt =  5.67e-8*(Tr**4 - T(3)**4)    ! Stefan-Boltzmann constant = 5.57e-8 W/(m^2 K^4)
  qlt = h_c*(Tf - T(1))
  R(1) = T0(1) + 2.*coef1*(-T(1) +   T(2)        ) + coef2*qlt + coef3*Qs(1) - T(1)
  R(2) = T0(2) +    coef1*( T(1) -2.*T(2) + T(3) ) +             coef3*Qs(2) - T(2)
  R(3) = T0(3) + 2.*coef1*(          T(2) - T(3) ) + coef2*qrt + coef3*Qs(3) - T(3)
  R_max = MAXVAL(ABS(R))
  ! print *, time, iter, R_max,"[",Qs(1),", ",Qs(2),", ",Qs(3),"]"  ! Debug

END SUBROUTINE residual_con1d
! ========================



!*****************************************************************************************
SUBROUTINE jacobian_con1d
!*****************************************************************************************
  USE con1d_mod
  IMPLICIT NONE

  REAL :: J(3,3)  ! Jacobian matrix.  Note  J(eq_index, T_index)

  ! Compute Jacobian matrix
  J(1,1) = -(1+ 2.*coef1) - coef2*h_c
  J(1,2) =   2.*coef1
  J(1,3) =   0.0
  J(2,1) =   coef1
  J(2,2) = -(1+ 2.*coef1)
  J(2,3) =   coef1
  J(3,1) =   0.0
  J(3,2) =   2.*coef1
  J(3,3) = -(1+ 2.*coef1) - coef2*5.67e-8*4.0*T(3)**3   ! This is the only non-linear term

END SUBROUTINE jacobian_con1d
! ========================


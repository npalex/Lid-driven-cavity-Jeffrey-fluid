! =====================================================
subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =====================================================

! Roe-solver for the 2D incompressible Cauchy momentum equations
! 	and the viscoelastic model for a Maxwell Fluid

! waves: 5
! equations: 5

! Conserved quantities:
!       1 x_momentum
!       2 y_momentum
!		3 normal stress tau_xx
!		4 shear stress tau_xy
!		5 normal stress tau_yy
!		6 pressure

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell

! This data is along a slice in the x-direction if ixy=1
!                            or the y-direction if ixy=2.

! On output, wave contains the waves, s the speeds,
! and amdq, apdq the decomposition of the flux difference
!   f(qr(i-1)) - f(ql(i))
! into leftgoing and rightgoing parts respectively.
! With the Roe solver we have
!    amdq  =  A^- \Delta q    and    apdq  =  A^+ \Delta q
! where A is the Roe matrix.  An entropy fix can also be incorporated
! into the flux differences.

! Note that the i'th Riemann problem has left state qr(:,i-1)
!                                    and right state ql(:,i)
! From the basic clawpack routines, this routine is called with ql = qr


! This Riemann solver differs from the original clawpack Riemann solver
! for the interleaved indices

    implicit double precision (a-h,o-z)

    dimension   ql(meqn,           1-mbc:maxm+mbc)
    dimension   qr(meqn,           1-mbc:maxm+mbc)
    dimension    s(mwaves,         1-mbc:maxm+mbc)
    dimension wave(meqn,   mwaves, 1-mbc:maxm+mbc)
    dimension  apdq(meqn,          1-mbc:maxm+mbc)
    dimension  amdq(meqn,          1-mbc:maxm+mbc)

!   # Re set in setrun.f90 file
	double precision Re, Wi, gamma, E, Ma
    common /cparam/  Re, Wi, gamma, E, Ma

!   # Roe averages quantities of each interface
    parameter (maxm2 = 1800)
    double precision u(-6:maxm2+7),v(-6:maxm2+7)


!   local arrays
!   ------------
    dimension delta(5)
!    logical :: efix

!    data efix /.False./    !# Use entropy fix for transonic rarefactions
	
	! Compute the Roe-averaged variables needed in the Roe solver.
	do 10 i = 2-mbc, mx+mbc
		u(i) = (qr(1,i-1)+ql(1,i))*0.50d0
        v(i) = (qr(2,i-1)+ql(2,i))*0.50d0
    10 END DO 
	
!   # find a1 through a5, the coefficients of the 5 eigenvectors 
!	# and compute the waves:
	if (ixy == 1) then
		do 20 i = 2-mbc, mx+mbc
			delta(1) = ql(1,i) - qr(1,i-1)
			delta(2) = ql(2,i) - qr(2,i-1)
			delta(3) = ql(3,i) - qr(3,i-1)
			delta(4) = ql(4,i) - qr(4,i-1)
			delta(5) = ql(5,i) - qr(5,i-1)
		
			!-- eigenvalues
			s(1,i) = 0.d0 
			s(2,i) = u(i)/2.d0 - dsqrt((u(i)/2.d0)**2 + 1.d0/(Re*Wi))
			s(3,i) = u(i)/2.d0 + dsqrt((u(i)/2.d0)**2 + 1.d0/(Re*Wi))
			s(4,i) = u(i) - dsqrt(u(i)**2 + 2.d0/(Re*Wi))
			s(5,i) = u(i) + dsqrt(u(i)**2 + 2.d0/(Re*Wi))
		
			!-- coefficients (alpha's)
			a1 = delta(5) 
		
			a2 = -Re*v(i)*(s(3,i)/(s(2,i)-s(3,i)))*delta(1) &
					- delta(2)/Wi/(s(2,i)-s(3,i)) &
					- v(i)*delta(3)*(Re*Wi*s(3,i)*u(i) + 1.d0)/(s(2,i)-s(3,i)) &
					- (s(3,i)/(s(2,i)-s(3,i)))*delta(4)
				
			a3 = Re*v(i)*(s(2,i)/(s(2,i)-s(3,i)))*delta(1) &
					+ delta(2)/Wi/(s(2,i)-s(3,i)) &
					+ v(i)*delta(3)*(Re*Wi*s(2,i)*u(i) + 1.d0)/(s(2,i)-s(3,i)) &
					+ (s(2,i)/(s(2,i)-s(3,i)))*delta(4)
				
			a4 = Re*v(i)*(s(5,i)/(s(4,i)-s(5,i)))*delta(1) &
					+ v(i)*delta(3)*(Re*Wi*s(5,i)*u(i) + 1.d0)/(s(4,i)-s(5,i))
				
			a5 = -Re*v(i)*(s(4,i)/(s(4,i)-s(5,i)))*delta(1) &
					- v(i)*delta(3)*(Re*Wi*s(4,i)*u(i) + 1.d0)/(s(4,i)-s(5,i))
		
!     	 # Compute the waves.
			wave(1,1,i) = 0.d0
			wave(2,1,i) = 0.d0
			wave(3,1,i) = 0.d0
			wave(4,1,i) = 0.d0
			wave(5,1,i) = a1*1.d0
		
			wave(1,2,i) = 0.d0
			wave(2,2,i) = -s(2,i)*Wi*a2
			wave(3,2,i) = 0.d0
			wave(4,2,i) = 1.d0*a2
			wave(5,2,i) = 0.d0
		
			wave(1,3,i) = 0.d0
			wave(2,3,i) = -s(3,i)*Wi*a3
			wave(3,3,i) = 0.d0
			wave(4,3,i) = 1.d0*a3
			wave(5,3,i) = 0.d0
		
			wave(1,4,i) = (-Wi*s(4,i)*u(i) - 1.d0/Re)*(Re*(s(5,i)/(s(4,i)-s(5,i)))*delta(1) &
					+ delta(3)*(Re*Wi*s(5,i)*u(i) + 1.d0)/(s(4,i)-s(5,i)))
			wave(2,4,i) = -Wi*s(4,i)*a4
			wave(3,4,i) = s(4,i)*(Re*(s(5,i)/(s(4,i)-s(5,i)))*delta(1) &
					+ delta(3)*(Re*Wi*s(5,i)*u(i) + 1.d0)/(s(4,i)-s(5,i)))
			wave(4,4,i) = 1.d0*a4
			wave(5,4,i) = 0.d0
		
			wave(1,5,i) = -(-Wi*s(5,i)*u(i) - 1.d0/Re)*(Re*(s(4,i)/(s(4,i)-s(5,i)))*delta(1) &
					+ delta(3)*(Re*Wi*s(4,i)*u(i) + 1.d0)/(s(4,i)-s(5,i)))
			wave(2,5,i) = -Wi*s(5,i)*a5
			wave(3,5,i) = -s(5,i)*(Re*(s(4,i)/(s(4,i)-s(5,i)))*delta(1) &
					+ delta(3)*(Re*Wi*s(4,i)*u(i) + 1.d0)/(s(4,i)-s(5,i)))
			wave(4,5,i) = 1.d0*a5
			wave(5,5,i) = 0.d0

		20 END DO
		
	else
		do 30 i = 2-mbc, mx+mbc
			delta(1) = ql(1,i) - qr(1,i-1)
			delta(2) = ql(2,i) - qr(2,i-1)
			delta(3) = ql(3,i) - qr(3,i-1)
			delta(4) = ql(4,i) - qr(4,i-1)
			delta(5) = ql(5,i) - qr(5,i-1)
		
			!-- eigenvalues
			s(1,i) = 0.d0 
			s(2,i) = v(i)/2.d0 - dsqrt((v(i)/2.d0)**2 + 1.d0/(Re*Wi))
			s(3,i) = v(i)/2.d0 + dsqrt((v(i)/2.d0)**2 + 1.d0/(Re*Wi))
			s(4,i) = v(i) - dsqrt(v(i)**2 + 2.d0/(Re*Wi))
			s(5,i) = v(i) + dsqrt(v(i)**2 + 2.d0/(Re*Wi))
		
			!-- coefficients (alpha's)
			a1 = delta(3) 
		
			a2 = -Re*u(i)*(s(3,i)/(s(2,i)-s(3,i)))*delta(2) &
					- delta(1)/Wi/(s(2,i)-s(3,i)) &
					- u(i)/2.d0*delta(5)*(Re*Wi*s(3,i)*s(4,i) + Re*Wi*s(3,i)*s(5,i) + 2.d0)/(s(2,i)-s(3,i)) &
					- (s(3,i)/(s(2,i)-s(3,i)))*delta(4)
				
			a3 = Re*u(i)*(s(2,i)/(s(2,i)-s(3,i)))*delta(2) &
					+ delta(1)/Wi/(s(2,i)-s(3,i)) &
					+ u(i)/2.d0*delta(5)*(Re*Wi*s(2,i)*s(4,i) + Re*Wi*s(2,i)*s(5,i) + 2.d0)/(s(2,i)-s(3,i)) &
					+ (s(2,i)/(s(2,i)-s(3,i)))*delta(4)
				
			a4 = -2.d0*delta(2)/Wi/(s(4,i)-s(5,i)) &
					- delta(5)*s(5,i)/(s(4,i)-s(5,i))
				
			a5 = 2.d0*delta(2)/Wi/(s(4,i)-s(5,i)) &
					+ delta(5)*s(4,i)/(s(4,i)-s(5,i))
		
		
!     	 # Compute the waves.
			wave(1,1,i) = 0.d0
			wave(2,1,i) = 0.d0
			wave(3,1,i) = a1*1.d0
			wave(4,1,i) = 0.d0
			wave(5,1,i) = 0.d0
		
			wave(1,2,i) = -s(2,i)*Wi*a2
			wave(2,2,i) = 0.d0
			wave(3,2,i) = 0.d0
			wave(4,2,i) = 1.d0*a2
			wave(5,2,i) = 0.d0
		
			wave(1,3,i) = -s(3,i)*Wi*a3
			wave(2,3,i) = 0.d0
			wave(3,3,i) = 0.d0
			wave(4,3,i) = 1.d0*a3
			wave(5,3,i) = 0.d0
		
			wave(1,4,i) = -Wi*u(i)*a4
			wave(2,4,i) = -Wi*s(4,i)/2.d0*a4
			wave(3,4,i) = 0.d0
			wave(4,4,i) = (-Re*Wi*u(i)*v(i) + Re*Wi*s(4,i)*u(i)/2.d0)*a4
			wave(5,4,i) = 1.d0*a4
		
			wave(1,5,i) = -Wi*u(i)*a5
			wave(2,5,i) = -Wi*s(5,i)/2.d0*a5
			wave(3,5,i) = 0.d0
			wave(4,5,i) = (-Re*Wi*u(i)*v(i) + Re*Wi*s(5,i)*u(i)/2.d0)*a5
			wave(5,5,i) = 1.d0*a5
		
		30 END DO
	
	endif

!    # compute flux differences amdq and apdq.
!    ---------------------------------------

!     # amdq = SUM s*wave   over left-going waves
!     # apdq = SUM s*wave   over right-going waves

    do m=1,5
        do i=2-mbc, mx+mbc
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do mw=1,mwaves
                if (s(mw,i) < 0.d0) then
                    amdq(m,i) = amdq(m,i) + s(mw,i)*wave(m,mw,i)
                else
                    apdq(m,i) = apdq(m,i) + s(mw,i)*wave(m,mw,i)
                endif
            end do
        end do
    end do

    return
    end subroutine rpn2



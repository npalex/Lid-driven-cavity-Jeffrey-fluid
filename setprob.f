      subroutine setprob
      implicit double precision (a-h,o-z)
      character*25 fname
	 
	  real(kind=8) :: Re, Wi, gamma, E, Ma
      common /cparam/ Re, Wi, gamma, E, Ma
c
c     # Set the material parameters for the INSE equations
c
      iunit = 7
      fname = 'setprob.data'

      call opendatafile(iunit, fname)
                
c     # Reynolds number
      read(7,*) Re
	  
c     # Weissenburg number
	  read(7,*) Wi
	  
c	  # viscosity ratio (solvent viscosity/mixture viscosity)
	  read(7,*) gamma
	  
c     # viscoelastic mach number
	  Ma = dsqrt(Re*Wi)

c     # Elasticity number
      E = Wi/Re
	  
      return
      end

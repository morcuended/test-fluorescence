C Test fluorescence subroutines
      PROGRAM testing
      implicit none
      integer          ier, iflag, nflpht, l
      double precision zem, edep, lamb, zemision
      double precision temperature, density
C
      double precision temp, pressure, rhof
      
      external         temp, pressure, rhof
      external         BLOCK1
      
#define __AIRINC__
#define __ATMOSINC__

#include "test_fluorescence_corsika.h"
C ---------------------------------------------------

      zem = 1.0d0  ! in cm
      edep = 1.0d0  ! in MeV 
      write(6,*) " Height", zem, "cm"
      write(6,*) " Deposited energy", edep, "MeV"
      
      DO  L = 1, 5
        DATM(L) = 1.D0 / CATM(L)
      ENDDO
C 
C-First test the basic functions
!       write(6,*) " AATM ", AATM
!       write(6,*) " BATM ", BATM
!       write(6,*) " CATM ", CATM
!       write(6,*) " DATM ", DATM
C
      density = rhof(zem)
      write(6,*) " Density is ", density, "g/cm**3"   
C
      temperature = temp(zem)
      write(6,*) " Temperature is ", temperature     

      ! Add fluorescence initialization
      CALL FYINI
      ! Determine the number of fluorescence photons to be emitted in
      ! each step.
      ier = 0
      iflag = 2

      call fyield ( zemision, edep, nflpht, iflag, ier, 0, lamb )
      write(*,*)  " Number of fluorescence photons:", nflpht
      if ( nflpht .lt. 1 ) return

      if ( ier .eq. 1 ) write(*,*)
     *   " Warning!, too many fluorescence photons produced at once" 

      END PROGRAM testing


      subroutine fyield(zem,edep,nphots,iflag,ier,nbunch,lambda)

C  Determines the number of fluorescence photons and its wavelength 
c  emitted for a given height and energy deposit. It depends on the 
c  atmospheric pressure (hPa) and temperature (K).
c  It is called from fluor.
c--------------------------------------------------------------------
c     Input arguments:
c--------------------------------------------------------------------
c     iflag     =   flag indicating what the subroutine will return 
c                     if = 2 -> returns only number of f-photons
c                     if = 3 -> returns wavelength of each bunch
c     zem       =   z emission (in cm)
c     edep      =   deposited energy (in MeV)
c--------------------------------------------------------------------
c     Output parameters:
c--------------------------------------------------------------------
c     nphots    =   number of yielded photons
c     ier       =   flag indicating error (if n > nphmax)
c     nbunch    =   number of produced bunches of photons
c     lambda    =   wavelength of each produced bunch
c--------------------------------------------------------------------

      implicit none
c#define __CEREN1INC__
c#define __RUNPARINC__
c#include "corsika.h"
c
c-Arrays for fyield calculation
c       
      integer ier,iflag,nbunch,nphmax
      parameter (nphmax = 3d5) 
      double precision zem, edep, lambda
      double precision xn

      integer i, j, k, nfltab, ifwlenl, ifwlenu
      parameter (nfltab = 57)
      common /fluordat/lamb, il, pl, alp, ifwlenl, ifwlenu
      double precision lamb(nfltab),il(nfltab),pl(nfltab),alp(nfltab)

      double precision pzero, tzero
      parameter (pzero = 800.d0, Tzero = 293.d0) !hPa, K respectively
      double precision yabs
      parameter (yabs = 7.04d0)
      double precision pprimt, pprimti
      double precision temp, pressure
      external temp, pressure
      double precision phtyield, alfa
      integer nphots
      double precision xran(nphmax)

c-Big table to keep the fluor yield as a function of wavelength
c-and height in the atmosphere

      integer izem, maxiz
      parameter (maxiz = 101)
      common /fztab/ztab, yint, yield, yieldi
      double precision y, ztab(maxiz), yield(nfltab,maxiz),
     &                 yieldi(nfltab,maxiz), yint(maxiz)
      
      save
c------------------------------------------------------------     

      izem = max(int(log10(zem/1.0D4)/0.03)+1,1)
      if ( izem.lt.1 ) then
       write(*,*) "Fyield Error: izem < 0", izem
       izem = 1
      endif 
      if ( izem.gt.maxiz ) then
       write(*,*) "Fyield Error: z too high ",zem/1.0d5,"km ",izem
       izem = maxiz
      endif 
c
      if (iflag.eq.2) then
c       
c  Determine the number of photons poisson distributed
c
c   Interpolate except for the highest zem
c
        if (izem.lt.(maxiz-1)) then
          alfa = (zem-ztab(izem))/(ztab(izem+1)-ztab(izem)) 
          phtyield = (yint(izem)+alfa*(yint(izem+1)-yint(izem))) * edep
        else  
          phtyield = yint(maxiz) * edep  
        endif

        nphots = int(phtyield) ! Uncomment this for tests
C        call mpoiss(phtyield, nphots )

c  Warning in case of getting a very large number of photons 
        if (nphots .gt. nphmax) then 
            write(*,*) "zem, edep ",zem, edep
            write(*,*) "Error: ",nphmax - nphots,
     *                                 " remaining photons not stored"
            ier = 1
            nphots = nphmax
        endif
        return
      endif

c  Return only a wavelength per bunch now
      if ( iflag .eq. 3 ) then
        if ( nbunch .gt. 0 ) then
          call rmmard(xran, nbunch, 3)
          do i = 1, nbunch
            do j = ifwlenl, ifwlenu
              if (xran(i) .lt. yieldi(j,izem)) then
                lambda = lamb(j)
                exit
              endif
            enddo
          enddo
        endif
       endif

      return 
      end 


      subroutine fyini

c  Initializes the fluorescence yield data  
c  
c--------------------------------------------------------------------
c     Input arguments:
c--------------------------------------------------------------------
c
c     all input arguments enter through Corsika commons
c
c--------------------------------------------------------------------
c     Output parameters:
c--------------------------------------------------------------------
c
c     ifwlenl = lower index of the fluorescen band to use
c     ifwlenu = upper index of the fluorescen band to use
c     Are initialized from wavlgl and wavlgu
c
c     The fluorescence yield tables as a function of wavelength and 
c     height are initilized
c
c--------------------------------------------------------------------

      implicit none
c#define __CEREN1INC__
c#define __RUNPARINC__
c#include "corsika.h"

c
c-Arrays for fyield calculation
c       
      integer i, k, nfltab, ifwlenl, ifwlenu
      parameter (nfltab = 57)
      double precision zem
      common /fluordat/lamb, il, pl, alp, ifwlenl, ifwlenu
      double precision lamb(nfltab),il(nfltab),pl(nfltab),alp(nfltab)
      double precision wavlgl, wavlgu
c
      double precision pzero, tzero
      parameter (pzero = 800.d0, Tzero = 293.d0) !hPa, K respectively
      double precision yabs
      parameter (yabs = 7.04d0)
      double precision pprimt, pprimti
      double precision temp, pressure
      external temp, pressure
      double precision phtyield
c
c-Big table to keep the fluor yield as a function of wavelength
c-and height in the atosphere
c
      integer izem, maxiz
      parameter (maxiz = 101)
      common /fztab/ztab, yint, yield, yieldi
      double precision y, ztab(maxiz), yield(nfltab,maxiz),
     &                 yieldi(nfltab,maxiz), yint(maxiz)

c      
      save

      data ifwlenl/-1/, ifwlenu/-1/
      
c  Fluorescence spectrum and parameters.

      data (lamb(k), il(k), pl(k), alp(k), k=1, nfltab)          
c       WLen           Il       Pressure     Alpha !   p_w
     */ 281.20d0 ,     0.18d0,   19.00d0,   0.00d0,!  0.00d0,
     *  281.80d0 ,     0.30d0,   20.70d0,   0.00d0,!  0.00d0,
     *  296.20d0 ,     5.16d0,   18.50d0,   0.00d0,!  0.00d0,
     *  297.70d0 ,     2.77d0,   17.30d0,   0.00d0,!  0.00d0,
     *  302.00d0 ,     0.41d0,   21.00d0,   0.00d0,!  0.00d0,
     *  308.00d0 ,     1.44d0,   21.00d0,   0.00d0,!  0.00d0,
     *  311.70d0 ,     7.24d0,   18.70d0,   0.00d0,!  0.00d0,
     *  313.60d0 ,    11.05d0,   12.27d0,  -0.13d0,!  1.20d0,
     *  315.90d0 ,    39.33d0,   11.88d0,  -0.19d0,!  1.10d0,
     *  317.70d0 ,     0.46d0,   21.00d0,   0.00d0,!  0.00d0,
     *  326.80d0 ,     0.80d0,   19.00d0,   0.00d0,!  0.00d0,
     *  328.50d0 ,     3.80d0,   20.70d0,   0.00d0,!  0.00d0,
     *  330.90d0 ,     2.15d0,   16.90d0,   0.00d0,!  0.00d0,
     *  333.90d0 ,     4.02d0,   15.50d0,   0.00d0,!  0.00d0,
     *  337.10d0 ,   100.00d0,   15.89d0,  -0.35d0,!  1.28d0,
     *  346.30d0 ,     1.74d0,   21.00d0,   0.00d0,!  0.00d0,
     *  350.00d0 ,     2.79d0,   15.20d0,  -0.38d0,!  1.50d0,
     *  353.70d0 ,    21.35d0,   12.70d0,  -0.22d0,!  1.27d0,
     *  357.70d0 ,    67.41d0,   15.39d0,  -0.35d0,!  1.30d0,
     *  365.90d0 ,     1.13d0,   21.00d0,   0.00d0,!  0.00d0,
     *  367.20d0 ,     0.54d0,   19.00d0,   0.00d0,!  0.00d0,
     *  371.10d0 ,     4.97d0,   14.80d0,  -0.24d0,!  1.30d0,
     *  375.60d0 ,    17.87d0,   12.82d0,  -0.17d0,!  1.10d0,
     *  380.50d0 ,    27.20d0,   16.51d0,  -0.34d0,!  1.40d0,
     *  385.80d0 ,     0.50d0,   19.00d0,   0.00d0,!  0.00d0,
     *  387.70d0 ,     1.17d0,    7.60d0,   0.00d0,!  0.00d0,
     *  388.50d0 ,     0.83d0,    3.90d0,   0.00d0,!  0.00d0,
     *  391.40d0 ,    28.00d0,    2.94d0,  -0.79d0,!  0.33d0,
     *  394.30d0 ,     3.36d0,   13.70d0,  -0.20d0,!  1.20d0,
     *  399.80d0 ,     8.38d0,   13.60d0,  -0.20d0,!  1.10d0,
     *  405.00d0 ,     8.07d0,   17.80d0,  -0.37d0,!  1.50d0,
     *  414.10d0 ,     0.49d0,   19.00d0,   0.00d0,!  0.00d0,
     *  420.00d0 ,     1.75d0,   13.80d0,   0.00d0,!  0.00d0,
     *  423.60d0 ,     1.04d0,    3.90d0,   0.00d0,!  0.00d0,
     *  427.00d0 ,     7.08d0,    6.38d0,   0.00d0,!  0.00d0,
     *  427.80d0 ,     4.94d0,    2.89d0,  -0.54d0,!  0.60d0,
     *  434.30d0 ,     2.23d0,   15.89d0,   0.00d0,!  0.00d0,
     *  435.50d0 ,     0.17d0,   19.00d0,   0.00d0,!  0.00d0,
     *  441.50d0 ,     0.56d0,   20.70d0,   0.00d0,!  0.00d0,
     *  448.90d0 ,     0.72d0,   12.27d0,   0.00d0,!  0.00d0,
     *  457.30d0 ,     0.91d0,   11.88d0,   0.00d0,!  0.00d0,
     *  459.70d0 ,     0.02d0,    3.90d0,   0.00d0,!  0.00d0,
     *  464.90d0 ,     0.66d0,    3.90d0,   0.00d0,!  0.00d0,
     *  466.50d0 ,     0.53d0,   15.89d0,   0.00d0,!  0.00d0,
     *  470.60d0 ,     1.93d0,    2.94d0,   0.00d0,!  0.00d0,
     *  481.30d0 ,     0.24d0,   12.27d0,   0.00d0,!  0.00d0,
     *  491.70d0 ,     0.25d0,   11.88d0,   0.00d0,!  0.00d0,
     *  503.20d0 ,     0.12d0,   15.89d0,   0.00d0,!  0.00d0,
     *  514.60d0 ,     0.16d0,    3.90d0,   0.00d0,!  0.00d0,
     *  522.50d0 ,     0.33d0,    2.94d0,   0.00d0,!  0.00d0,
     *  530.90d0 ,     0.06d0,   11.88d0,   0.00d0,!  0.00d0,
     *  545.20d0 ,     0.03d0,   15.89d0,   0.00d0,!  0.00d0,
     *  570.40d0 ,     0.01d0,    3.90d0,   0.00d0,!  0.00d0,
     *  575.00d0 ,     0.03d0,    3.90d0,   0.00d0,!  0.00d0,
     *  586.10d0 ,     0.05d0,    2.94d0,   0.00d0,!  0.00d0,
     *  593.80d0 ,     0.01d0,   15.89d0,   0.00d0,!  0.00d0,
     *  665.90d0 ,     0.01d0,    2.94d0,   0.00d0/!  0.00d0/
c------------------------------------------------------------     
c
c  Select only the fluorescence bands within the range of wavelengths
c  set by the cwavlg in the input card (same as for Cherenkov).
c  The rows used are determined only once depending on the lower
c  and upper limits.
c
C-JLC: fix wavelength limits by hand
C
      wavlgl = 280.D0
      wavlgu = 700.D0  
      do i = 1, nfltab
        ifwlenl = i
        if (lamb(i) .ge. wavlgl) exit
      enddo
      ifwlenu = 1
      do i = 1, nfltab
        if (lamb(i) .gt. wavlgu) exit
        ifwlenu = i
      enddo
c
      do izem = 1, maxiz
        zem = 1.0D4*10**(0.03*(izem-1)) ! zem in cm
        ztab(izem) = zem  

c  Determine the half of the fluorescence yield as we consider 
c  only photons moving downwards.           
        yint(izem) = 0.d0
        do i = ifwlenl, ifwlenu             
           pprimt  = pl(i) * (tzero / temp(zem))**(alp(i) - 0.5d0)
           y = 0.5d0 * yabs * 1.d-2 * il(i)
           y = y * (1.d0 + pzero/pl(i))/(1.d0 + pressure(zem)/pprimt)
           yield(i,izem) = y
           yint(izem) = yint(izem) + y 
           yieldi(i,izem) = yint(izem)         
        enddo
c  Normalize the fluorescence yield in each band
        do i = ifwlenl, ifwlenu
           yieldi(i,izem) = yieldi(i,izem) / yint(izem) 
C           write(6,*) yint(izem)
        enddo
c     
      enddo
c
      return 
      end 



      function pressure( arg )

c-----------------------------------------------------------------------
c  Calculate the pressure (hPa) at a given height (cm) assuming an ideal
c  gas under a gravitational acceleration of 9.81 m/s**2
c  This function is called from temp and fyield.
c  argument:
c   arg    = height (cm)
c-----------------------------------------------------------------------

      implicit none
      double precision arg
      double precision pressure
      double precision thick
      external         thick
      save
      double precision grav
      parameter (grav = 9.81d2 )
c-----------------------------------------------------------------------
 
      pressure=grav*thick(arg)*1.d-3
      
      end function pressure
      


      double precision function temp( arg )

c-----------------------------------------------------------------------
c  temp(erature)
c
c  Calculate the local temperature (K) of atmosphere at a certain 
c  height (cm).
c  This function is called from fyield.
c  Argument:
c   arg    = height (cm)
c-----------------------------------------------------------------------

      implicit none
      double precision arg
      double precision pressure,rhof
      external         pressure,rhof
      save
      double precision Mmolar,Rconst !in cm3 kPa k−1 g mol−1
      parameter (Mmolar=28.96, Rconst=8.3144598d3 )
c----------------------------------------------------------------------- 

      temp=Mmolar*pressure(arg)*1.d-1/(rhof(arg)*Rconst)
      
      end function temp


*-- Author :    D. HECK IK FZK KARLSRUHE       15/10/1996
C=======================================================================

      SUBROUTINE MPOISS( AMEAN,NPRAN )

C-----------------------------------------------------------------------
C  M(UON COULOMB SCATTERING) POISS(ON DISTRIBUTION)
C
C  GENERATES A RANDOM NUMBER POISSON DISTRIBUTED WITH MEAN VALUE AMEAN.
C  THIS SUBROUTINE IS IN ANALOGY WITH SUBROUT. GPOISS.
C  (AUTHOR: L. URBAN) OF GEANT321
C  SEE CERN PROGRAM LIBRARY LONG WRITEUP W5013.
C  THIS SUBROUTINE IS CALLED FROM MUCOUL.
C  ARGUMENTS:
C   AMEAN  =  MEAN VALUE OF RANDOM NUMBER
C   NPRAN  =  RANDOM NUMBER POISSON DISTRIBUTED
C-----------------------------------------------------------------------

      IMPLICIT NONE
#define __CONSTAINC__
#define __RANDPAINC__
#define __RUNPARINC__
#include "test_fluorescence_corsika.h"

      DOUBLE PRECISION AMEAN,AN,HMXINT,P,PLIM,RR,S,X
      INTEGER          NPRAN
      SAVE
      DATA             PLIM / 16.D0 /, HMXINT / 2.D9 /
C-----------------------------------------------------------------------

C  PROTECTION AGAINST NEGATIVE MEAN VALUES
      AN = 0.D0
      IF ( AMEAN .GT. 0.D0 ) THEN
        IF ( AMEAN .LE. PLIM ) THEN
          CALL RMMARD( RD,1,1 )
          P  = EXP( -AMEAN )
          S  = P
          IF ( RD(1) .LE. S ) GOTO 20
 10       AN = AN + 1.D0
          P  = P * AMEAN / AN
          S  = S + P
          IF ( S .LT. RD(1)  .AND.  P .GT. 1.D-30 ) GOTO 10
        ELSE
          CALL RMMARD( RD,2,1 )
          RR = SQRT( (-2.D0)*LOG( RD(1) ) )
          X  = RR * COS( PI2 * RD(2) )
          AN = MIN( MAX( 0.D0, AMEAN+X*SQRT( AMEAN ) ), HMXINT )
        ENDIF
      ENDIF
 20   NPRAN = NINT(AN)

      RETURN
      END


*-- Author :    D. HECK IK FZK KARLSRUHE   17/03/2003
C=======================================================================

      SUBROUTINE RMMARD( RVEC,LENV,ISEQ )

C-----------------------------------------------------------------------
C  R(ANDO)M (NUMBER GENERATOR OF) MAR(SAGLIA TYPE) D(OUBLE PRECISION)
C
C  THESE ROUTINES (RMMARD,RMMAQD) ARE MODIFIED VERSIONS OF ROUTINES
C  FROM THE CERN LIBRARIES. DESCRIPTION OF ALGORITHM SEE:
C   http://wwwasdoc.web.cern.ch/wwwasdoc/cernlib.html (v113)
C  IT HAS BEEN CHECKED THAT RESULTS ARE BIT-IDENTICAL WITH CERN
C  DOUBLE PRECISION RANDOM NUMBER GENERATOR RMM48, DESCRIBED IN
C   http://wwwasdoc.web.cern.ch/wwwasdoc/cernlib.html (v116)
C  ARGUMENTS:
C   RVEC   = DOUBLE PREC. VECTOR FIELD TO BE FILLED WITH RANDOM NUMBERS
C   LENV   = LENGTH OF VECTOR (# OF RANDNUMBERS TO BE GENERATED)
C   ISEQ   = # OF RANDOM SEQUENCE
C
C  VERSION OF D. HECK FOR DOUBLE PRECISION RANDOM NUMBERS.
C-----------------------------------------------------------------------

      IMPLICIT NONE
#define __RANMA3INC__
#define __RANMA4INC__
#include "test_fluorescence_corsika.h"

      DOUBLE PRECISION RVEC(*)
      INTEGER          ISEQ,IVEC,LENV
      SAVE
C-----------------------------------------------------------------------

      IF ( ISEQ .GT. 0  .AND.  ISEQ .LE. KSEQ ) JSEQ = ISEQ

      DO  IVEC = 1, LENV
        UNI = U(I97(JSEQ),JSEQ) - U(J97(JSEQ),JSEQ)
        IF ( UNI .LT. 0.D0 ) UNI = UNI + 1.D0
        U(I97(JSEQ),JSEQ) = UNI
        I97(JSEQ)  = I97(JSEQ) - 1
        IF ( I97(JSEQ) .EQ. 0 ) I97(JSEQ) = 97
        J97(JSEQ)  = J97(JSEQ) - 1
        IF ( J97(JSEQ) .EQ. 0 ) J97(JSEQ) = 97
        C(JSEQ)    = C(JSEQ) - CD
        IF ( C(JSEQ) .LT. 0.D0 ) C(JSEQ)  = C(JSEQ) + CM
        UNI        = UNI - C(JSEQ)
        IF ( UNI .LT. 0.D0 ) UNI = UNI + 1.D0
C  AN EXACT ZERO HERE IS VERY UNLIKELY, BUT LET''S BE SAFE.
        IF ( UNI .EQ. 0.D0 ) UNI = TWOM48
        RVEC(IVEC) = UNI
      ENDDO

      NTOT(JSEQ) = NTOT(JSEQ) + LENV
      IF ( NTOT(JSEQ) .GE. MODCNS )  THEN
        NTOT2(JSEQ) = NTOT2(JSEQ) + 1
        NTOT(JSEQ)  = NTOT(JSEQ) - MODCNS
      ENDIF

      RETURN
      END


*-- Author :    The CORSIKA development group   21/04/1994
C=======================================================================

      DOUBLE PRECISION FUNCTION RHOF( ARG )

C-----------------------------------------------------------------------
C  RHO (DENSITY) F(UNCTION)
C
C  CALCULATES DENSITY (G/CM**3) OF ATMOSPHERE DEPENDING ON HEIGHT (CM)
C  THIS FUNCTION IS CALLED FROM BOX2, LPMEFFECT, ININKG, CERENK,
C  MUTRAC, AND INRTAB.
C  ARGUMENT:
C   ARG    = HEIGHT (CM)
C-----------------------------------------------------------------------

      IMPLICIT NONE
#define __ATMOSINC__
#define __ATMOS2INC__
#define __RUNPARINC__
#include "test_fluorescence_corsika.h"

      DOUBLE PRECISION ARG
      SAVE
C-----------------------------------------------------------------------

!       write(*,*) "HLAY", HLAY(2)
!       write(*,*) "BATM", BATM(1)
!       write(*,*) "DATM1", DATM(1)
!       write(*,*) "DATM2", DATM(2)
!       write(*,*) "DATM3", DATM(3)
!       write(*,*) "DATM4", DATM(4)
!       write(*,*) "DATM5", DATM(5)
      IF     ( ARG .LT. HLAY(2) ) THEN
        RHOF = BATM(1) * DATM(1) * EXP( (-ARG) * DATM(1) )
      ELSEIF ( ARG .LT. HLAY(3) ) THEN
        RHOF = BATM(2) * DATM(2) * EXP( (-ARG) * DATM(2) )
      ELSEIF ( ARG .LT. HLAY(4) ) THEN
        RHOF = BATM(3) * DATM(3) * EXP( (-ARG) * DATM(3) )
      ELSEIF ( ARG .LT. HLAY(5) ) THEN
        RHOF = BATM(4) * DATM(4) * EXP( (-ARG) * DATM(4) )
      ELSE
        RHOF = DATM(5)
      ENDIF
!       write(*,*) "RHO", rhof

      RETURN
      END


*-- Author :    The CORSIKA development group   21/04/1994
C=======================================================================

      DOUBLE PRECISION FUNCTION THICK( ARG )

C-----------------------------------------------------------------------
C  THICK(NESS OF ATMOSPHERE)
C
C  CALCULATES THICKNESS (G/CM**2) OF ATMOSPHERE DEPENDING ON HEIGHT (CM)
C  THIS FUNCTION IS CALLED FROM AAMAIN, BOX2, BOX3, EM, INPRM, MUBREM,
C  MUDECY, MUPRPR, MUTRAC, NRANGC, NUCINT, PRANGC, START, UPDATC,
C  UPDATE, EGS4, ELECTR, HOWFAR, PHOTON, ININKG, NKG, AND CERENK.
C  ARGUMENT:
C   ARG    = HEIGHT (CM)
C-----------------------------------------------------------------------

      IMPLICIT NONE
#define __ATMOSINC__
#define __ATMOS2INC__
#define __RUNPARINC__
#include "test_fluorescence_corsika.h"

      DOUBLE PRECISION ARG
      SAVE
C-----------------------------------------------------------------------

!       write(*,*) "DATM5", DATM(5)
      IF     ( ARG .LT. HLAY(2) ) THEN
        THICK = AATM(1) + BATM(1) * EXP( (-ARG) * DATM(1) )
      ELSEIF ( ARG .LT. HLAY(3) ) THEN
        THICK = AATM(2) + BATM(2) * EXP( (-ARG) * DATM(2) )
      ELSEIF ( ARG .LT. HLAY(4) ) THEN
        THICK = AATM(3) + BATM(3) * EXP( (-ARG) * DATM(3) )
      ELSEIF ( ARG .LT. HLAY(5) ) THEN
        THICK = AATM(4) + BATM(4) * EXP( (-ARG) * DATM(4) )
      ELSE
        THICK = AATM(5) - ARG * DATM(5)
      ENDIF
      
!       write(*,*) "THICK", THICK

      RETURN
      END


*-- Author :    The CORSIKA development group   21/04/1994
C=======================================================================

      BLOCK DATA BLOCK1

C-----------------------------------------------------------------------
C  INITIALIZES DATA
C  THIS ROUTINE IS CALLED FROM AAMAIN
C-----------------------------------------------------------------------

      IMPLICIT NONE
#define __AIRINC__
#define __ATMOSINC__
#define __ATMOS2INC__
#define __BUFFSINC__
#define __CONSTAINC__
#define __EDECAYINC__
#define __GNUPRINC__
#define __KAONSINC__
#define __MUPARTINC__
#define __NKGIINC__
#define __PAMINC__
#define __PARPARINC__
#define __PARPAEINC__
#define __RESTINC__
#define __RUNPARINC__
#define __STACKFINC__
#define __STRBARINC__
#define __VERSINC__


#include "test_fluorescence_corsika.h"

      INTEGER          I,J
C-----------------------------------------------------------------------

C  AIR
      DATA COMPOS / 0.78479D0, 0.21052D0, 0.00469D0 /
      DATA PROBTA / 0.78479D0, 0.99531D0, 1.D0      /
      DATA AVERAW / 14.543D0 /
C  VALUE OF AVOGADRO REVISED JAN. 2016 BY D.H.
      DATA AVOGDR / 6.02214129D-4 /

C  ATMOS (U.S.STANDARD IS DEFAULT)
      DATA AATM / -186.555306D0,-94.919D0,  0.61289D0,0.D0,.01128292D0 /
      DATA BATM / 1222.6562D0,1144.9069D0,1305.5948D0,540.1778D0,1.D0  /
      DATA CATM / 994186.38D0,878153.55D0,636143.04D0,772170.16D0,1.D9 /

      DATA ((AATM0(I,J),I=1,5),J=0,16)
     * /-186.5562D0, -94.919D0  ,.61289D0  ,        0.D0 , .01128292D0 ,
     *  -186.5562D0, -94.919D0  ,.61289D0  ,        0.D0 , .01128292D0 ,
     *  -118.1277D0,-154.258D0  ,.4191499D0, 5.4094056D-4, .01128292D0 ,
     * -195.837264D0,-50.4128778D0,.345594007D0,5.46207D-4,.01128292D0 ,
     * -253.95047D0,-128.97714D0,.353207D0 , 5.526876D-4 , .01128292D0 ,
     * -208.12899D0,-120.26179D0,.31167036D0,5.591489D-4 , .01128292D0 ,
     * -77.875723D0,-214.96818D0,.3721868D0, 5.5309816D-4, .01128292D0 ,
     * -242.56651D0,-103.21398D0,.3349752D0, 5.527485D-4 , .01128292D0 ,
     * -195.34842D0,-71.997323D0,.3378142D0, 5.48224D-4  , .01128292D0 ,
     *    0.D0     , 0.D0       ,      0.D0, 0.D-4       , .01128292D0 ,
     *    0.D0     , 0.D0       ,      0.D0, 0.D-4       , .01128292D0 ,
     *  -137.656D0,  -37.9610D0,   .222659D0, -6.16201D-4 ,.00207722D0 ,
     *  -163.331D0,  -65.3713D0,   .402903D0, -4.79198D-4 ,.00188667D0 ,
     *  -142.801D0,  -70.1538D0,   1.14855D0, -9.10269D-4 ,.00152236D0 ,
     *  -128.601D0,  -39.5548D0,   1.13088D0, -26.4960D-4 ,.00192534D0 ,
     *  -113.139D0,  -79.0635D0,  -54.3888D0,  0.D0       ,.4210330D-2 ,
     *  -59.0293D0,  -21.5794D0,  -7.14839D0,  0.D0       ,.1901750D-3 /
      DATA ((AATM0(I,J),I=1,5),J=17,25)
     */-149.801663D0,-57.932486D0,.63631894D0,4.35453690D-4,.01128292D0, !US
     * -136.72575606,-31.636643044,1.8890234035,3.92018679839D-4,
     *                                                     .01128292,  !Jan
     * -137.25655862,-31.793978896,2.0616227547,4.12430622892D-4,
     *                                                     .01128292,  !Feb
     * -132.36885162,-29.077046629,2.090501509 ,4.35343379252D-4,
     *                                                     .01128292,  !March
     * -129.9930412 ,-21.847248438,1.5211136484,3.95590551213D-4,
     *                                                     .01128292,  !Apr
     * -125.11468467,-14.591235621,.93641128677,3.24755909854D-4,
     *                                                     .01128292,  !May
     * -126.17178851,-7.7289852811,.81676828638,3.19476768915D-4,
     *                                                     .01128292,  !June
     * -126.17216789,-8.6182537514,.74177836911,2.93507020973D-4,
     *                                                     .01128292,  !July
     * -123.27936204,-10.051493041,.84187346153,3.24225467594D-4,
     *                                                     .01128292/  !Aug
      DATA ((AATM0(I,J),I=1,5),J=26,29)
     */-126.94494665,-9.5556536981,.74939405052,2.98231169610D-4,
     *                                                     .01128292,  !Sep
     * -133.13151125,-13.973209265,0.8378263431,3.11174217600D-4,
     *                                                     .01128292,  !Oct
     * -134.72208165,-18.172382908,1.1159806845,3.52170255153D-4,
     *                                                     .01128292,  !Nov
     * -135.40825209,-22.830409026,1.4223453493,3.75129217743D-4,
     *                                                     .01128292/  !Dec
      DATA ((AATM0(I,J),I=1,5),J=30,41)
     */  -91.6956D0,   7.01491D0, 0.505452D0,-1.81302D-03, 2.07722D-03,!Jan
     *   -72.1988D0,  22.7002D0 , 0.430171D0,-1.2030D-03 , 2.07722D-03,!Feb
     *   -63.7290D0,  -1.02799D0, 0.324414D0,-4.90772D-04, 2.07722D-03,!Mar
     *   -69.7259D0,  -2.79781D0, 0.262692D0,-8.41695D-05, 2.07722D-03,!Apr
     *   -78.5551D0,  -5.33239D0, 0.312889D0,-9.20472D-05, 1.52236D-03,!May
     *   -92.6125D0,  -8.56450D0, 0.363986D0, 1.65164D-05, 2.07722D-03,!Jun
     *   -89.9639D0, -13.9697D0,  0.441631D0,-1.46525D-05, 2.07722D-03,!Jul
     *   -90.4253D0, -18.7154D0,  0.513930D0,-2.15650D-04, 1.52236D-03,!Aug
     *   -91.6860D0, -23.3519D0,  0.891302D0,-7.65666D-04, 2.07722D-03,!Sep
     *   451.616D0 , -85.5456D0,  2.06082D0 ,-1.07600D-03, 2.07722D-03,!Oct
     *  -152.853D0 ,   4.22741D0, 1.38352D0 ,-1.15014D-03, 2.07722D-03,!Nov
     *  -100.386D0 ,   5.43849D0, 0.399465D0,-1.75472D-03, 2.07722D-03/!Dec
C
      DATA ((BATM0(I,J),I=1,5),J=0,16)
     * / 1222.6562D0, 1144.9069D0, 1305.5948D0, 540.1778D0, 1.D0   ,
     *   1222.6562D0, 1144.9069D0, 1305.5948D0, 540.1778D0, 1.D0   ,
     *   1173.9861D0, 1205.7625D0, 1386.7807D0, 555.8935D0, 1.D0   ,
     *   1240.48D0  , 1117.85D0  , 1210.9D0   , 608.2128D0, 1.D0   ,
     *   1285.2782D0, 1173.1616D0, 1320.4561D0, 680.6803D0, 1.D0   ,
     *   1251.474D0 , 1173.321D0 , 1307.826D0 , 763.1139D0, 1.D0   ,
     *   1103.3362D0, 1226.5761D0, 1382.6933D0, 685.6073D0, 1.D0   ,
     *   1262.7013D0, 1139.0249D0, 1270.2886D0, 681.4061D0, 1.D0   ,
     *   1210.4D0   , 1103.8629D0, 1215.3545D0, 629.7611D0, 1.D0   ,
     *      0.D0    ,    0.D0    ,    0.D0    ,   0.D0    , 1.D0   ,
     *      0.D0    ,    0.D0    ,    0.D0    ,   0.D0    , 1.D0   ,
     *     1130.74D0,   1052.05D0,   1137.21D0,   442.512D0, 1.D0  ,
     *     1183.70D0,   1108.06D0,   1424.02D0,   207.595D0, 1.D0  ,
     *     1177.19D0,   1125.11D0,   1304.77D0,   433.823D0, 1.D0  ,
     *     1139.99D0,   1073.82D0,   1052.96D0,   492.503D0, 1.D0  ,
     *     1133.10D0,   1101.20D0,   1085.00D0,   1098.00D0, 1.D0  ,
     *     1079.00D0,   1071.90D0,   1182.00D0,   1647.10D0, 1.D0  /
      DATA ((BATM0(I,J),I=1,5),J=17,29)
     * /1183.6071D0 , 1143.0425D0 ,1322.9748D0, 655.67307D0, 1.D0 , !US
     *  1174.8298334, 1204.8233453,1637.7703583,735.96095023,1.D0 , !Jan
     *  1176.0907565, 1197.8951104,1646.4616955,755.18728657,1.D0 , !Feb
     *  1172.6227784, 1215.3964677,1617.0099282,769.51991638,1.D0 , !March
     *  1172.3291878, 1250.2922774,1542.6248413,713.1008285 ,1.D0 , !Apr
     *  1169.9511302, 1277.6768488,1493.5303781,617.9660747 ,1.D0 , !May
     *  1171.0916276, 1295.3516434,1455.3009344,595.11713507,1.D0 , !June
     *  1172.7340688, 1258.9180079,1450.0537141,583.07727715,1.D0 , !July
     *  1169.763036 , 1251.0219808,1436.6499372,627.42169844,1.D0 , !Aug
     *  1174.8676453, 1251.5588529,1440.8257549,606.31473165,1.D0 , !Sep
     *  1176.9833473, 1244.234531 ,1464.0120855,622.11207419,1.D0 , !Oct
     *  1175.7737972, 1238.9538504,1505.1614366,670.64752105,1.D0 , !Nov
     *  1174.644971 , 1227.2753683,1585.7130562,691.23389637,1.D0 / !Dec
      DATA ((BATM0(I,J),I=1,5),J=30,41)
     */  1125.71D0,    1149.81D0,    1032.68D0,    490.789D0,  1.D0,   !Jan
     *   1108.19D0,    1159.77D0,    1079.25D0,    523.956D0,  1.D0,   !Feb
     *   1102.66D0,    1093.56D0,    1198.93D0,    589.827D0,  1.D0,   !Mar
     *   1111.70D0,    1128.64D0,    1413.98D0,    587.688D0,  1.D0,   !Apr
     *   1118.46D0,    1169.09D0,    1577.71D0,    452.177D0,  1.D0,   !May
     *   1129.88D0,    1191.98D0,    1619.82D0,    411.586D0,  1.D0,   !Jun
     *   1125.73D0,    1180.47D0,    1581.43D0,    373.796D0,  1.D0,   !Jul
     *   1125.01D0,    1175.60D0,    1518.03D0,    299.006D0,  1.D0,   !Aug
     *   1125.53D0,    1169.77D0,    1431.26D0,    247.030D0,  1.D0,   !Sep
     *    849.239D0,   1113.16D0,    1322.28D0,    372.242D0,  1.D0,   !Oct
     *   1174.09D0,    1272.49D0,    975.906D0,    481.615D0,  1.D0,   !Nov
     *   1128.71D0,    1198.10D0,    858.522D0,    480.142D0,  1.D0 /  !Dec
C
      DATA ((CATM0(I,J),I=1,5),J=0,16)
     * / 994186.38D0, 878153.55D0, 636143.04D0, 772170.16D0, 1.D9  ,
     *   994186.38D0, 878153.55D0, 636143.04D0, 772170.16D0, 1.D9  ,
     *   919546.D0  , 963267.92D0, 614315.D0  , 739059.6D0 , 1.D9  ,
     *   933697.D0  , 765229.D0  , 636790.D0  , 733793.8D0 , 1.D9  ,
     *  1088310.D0  , 935485.D0  , 635137.D0  , 727312.6D0 , 1.D9  ,
     *  1032310.D0  , 925528.D0  , 645330.D0  , 720851.4D0 , 1.D9  ,
     *   932077.D0  ,1109960.D0  , 630217.D0  , 726901.3D0 , 1.D9  ,
     *  1059360.D0  , 888814.D0  , 639902.D0  , 727251.8D0 , 1.D9  ,
     *   970276.D0  , 820946.D0  , 639074.D0  , 731776.5D0 , 1.D9  ,
     *        0.D0  ,      0.D0  ,      0.D0  ,       0.D0 , 1.D9  ,
     *        0.D0  ,      0.D0  ,      0.D0  ,       0.D0 , 1.D9  ,
     *   867358.D0  , 741208.D0  , 633846.D0  , 759850.D0, 5.4303203D9,
     *   875221.D0  , 753213.D0  , 545846.D0  , 793043.D0, 5.9787908D9,
     *   861745.D0  , 765925.D0  , 581351.D0  , 775155.D0, 7.4095699D9,
     *   861913.D0  , 744955.D0  , 675928.D0  , 829627.D0, 5.8587010D9,
     *   861730.D0  , 826340.D0  , 790950.D0  , 682800.D0, 2.6798156D9,
     *   764170.D0  , 699910.D0  , 635650.D0  , 551010.D0, 59.329575D9/
      DATA ((CATM0(I,J),I=1,5),J=17,29)
     * / 954248.34D0 ,800005.34D0 ,629568.93D0 ,737521.77D0 ,1.D9  , !US
     *   982815.95248,754029.87759,594416.83822,733974.36972,1.D9  , !Jan
     *   981369.6125 ,756657.65383,592969.89671,731345.88332,1.D9  , !Feb
     *   972654.0563 ,742769.2171 ,595342.19851,728921.61954,1.D9  , !March
     *   962396.5521 ,711452.06673,603480.61835,735460.83741,1.D9  , !Apr
     *   947742.88769,685089.57509,609640.01932,747555.95526,1.D9  , !May
     *   940102.98842,661697.57543,612702.0632 ,749976.26832,1.D9  , !June
     *   934649.58886,672975.82513,614888.52458,752631.28536,1.D9  , !July
     *   931569.97625,678861.75136,617363.34491,746739.16141,1.D9  , !Aug
     *   936953.91919,678906.60516,618132.60561,750154.67709,1.D9  , !Sep
     *   954151.404  ,692708.89816,615439.43936,747969.08133,1.D9  , !Oct
     *   964877.07766,706199.57502,610242.24564,741412.74548,1.D9  , !Nov
     *   973884.44361,723759.74682,600308.13983,738390.20525,1.D9  / !Dec
      DATA ((CATM0(I,J),I=1,5),J=30,41)
     * / 821621.D0,  635444.D0,  682968.D0,  807327.D0,  54303203.D02,!Jan
     *   786271.D0,  599986.D0,  667432.D0,  780919.D0,  54303203.D02,!Feb
     *   764831.D0,  660389.D0,  636118.D0,  734909.D0,  54303203.D02,!Mar
     *   766099.D0,  641716.D0,  588082.D0,  693300.D0,  54303203.D02,!Apr
     *   776648.D0,  626683.D0,  553087.D0,  696835.D0,  74095699.D02,!May
     *   791177.D0,  618840.D0,  535235.D0,  692253.D0,  54303203.D02,!Jun
     *   784553.D0,  628042.D0,  531652.D0,  703417.D0,  54303203.D02,!Jul
     *   781628.D0,  633793.D0,  533269.D0,  737794.D0,  74095699.D02,!Aug
     *   786017.D0,  645241.D0,  545022.D0,  805419.D0,  54303203.D02,!Sep
     *   225286.D0,  789340.D0,  566132.D0,  796434.D0,  54303203.D02,!Oct
     *   891602.D0,  582119.D0,  643130.D0,  783786.D0,  54303203.D02,!Nov
     *   829352.D0,  612649.D0,  706104.D0,  806875.D0,  54303203.D02/!Dec
C
      DATA (LAYNO(J), J=0,41)
     * / 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7,
     *   8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
     *  24, 25, 26, 27, 28 /
      DATA (HLAY(I),I=1,5)
     * /     -5779.5D2 ,   4.D5  ,      1.D6  ,     4.D6,   1.D7     /
      DATA ((HLAY0(I,J),I=1,5),J=0,16)
     * /     -5779.5D2 ,   4.D5  ,      1.D6  ,     4.D6 ,  1.D7     ,
     *       -5779.5D2 ,   4.D5  ,      1.D6  ,     4.D6 ,  1.D7     ,
     *       -5.D5  , 2.66667D5  , 5.33333D5  ,     8.D5 ,  1.D7     ,!SPolJan
     *       -5.D5  , 6.66667D5  ,13.33333D5  ,    20.D5 ,  1.D7     ,!SPolAug
     *       -5.D5  ,     7.0D5  ,    11.4D5  ,    37.0D5,  1.D7     ,!US
     *       -5779.5D2,   9.4D5  ,    15.3D5  ,    31.6D5,  1.D7     ,!Jan
     *       -5779.5D2,   9.2D5  ,    15.4D5  ,    31.0D5,  1.D7     ,!Feb
     *       -5779.5D2,   9.6D5  ,    15.2D5  ,    30.7D5,  1.D7     ,!March
     *       -5779.5D2,  10.0D5  ,    14.9D5  ,    32.6D5,  1.D7     ,!Apr
     *       -5779.5D2,  10.2D5  ,    15.1D5  ,    35.9D5,  1.D7     ,!May
     *       -5779.5D2,  10.1D5  ,    16.0D5  ,    36.7D5,  1.D7     ,!June
     *       -5779.5D2,   9.6D5  ,    16.5D5  ,    37.4D5,  1.D7     ,!July
     *       -5779.5D2,   9.6D5  ,    15.9D5  ,    36.3D5,  1.D7     ,!Aug
     *       -5779.5D2,   9.5D5  ,    16.2D5  ,    37.2D5,  1.D7     ,!Sep
     *       -5779.5D2,   9.5D5  ,    15.5D5  ,    36.5D5,  1.D7     ,!Oct
     *       -5779.5D2,   9.6D5  ,    15.3D5  ,    34.6D5,  1.D7     ,!Nov
     *       -5779.5D2,   9.6D5  ,    15.6D5  ,    33.3D5,  1.D7     /!Dec
      DATA ((HLAY0(I,J),I=1,5),J=17,28)
     * / -5779.5D2,  7.8D05,  16.4D05,   40.4D05,   100.D05,           !Jan
     *   -5779.5D2,  8.0D05,  10.6D05,   40.4D05,   100.D05,           !Feb
     *   -5779.5D2,  6.7D05,  22.4D05,   40.4D05,   100.D05,           !Mar
     *   -5779.5D2,  7.6D05,  22.0D05,   40.4D05,   100.D05,           !Apr
     *   -5779.5D2,  8.4D05,  20.0D05,   39.7D05,   100.D05,           !May
     *   -5779.5D2,  8.5D05,  17.9D05,   38.4D05,   100.D05,           !Jun
     *   -5779.5D2,  8.5D05,  15.9D05,   37.5D05,   100.D05,           !Jul
     *   -5779.5D2,  8.5D05,  14.4D05,   37.5D05,   100.D05,           !Aug
     *   -5779.5D2,  8.5D05,  13.0D05,   36.2D05,   100.D05,           !Sep
     *   -5779.5D2,  3.1D05,  10.1D05,   31.5D05,   100.D05,           !Oct
     *   -5779.5D2,  8.5D05,  22.4D05,   32.4D05,   100.D05,           !Nov
     *   -5779.5D2,  8.5D05,  22.0D05,   40.4D05,   100.D05  /         !Dec

C  CONSTA
      DATA PI  / 3.141592653589793D0 /
      DATA PI2 / 6.283185307179586D0 /
      DATA OB3 / 0.333333333333333D0 /
      DATA TB3 / 0.666666666666666D0 /
C  ENEPER IS CALCULATED IN START: ENEPER = EXP(1.D0)
C  DATA FOR MUPART: CUTOFF FOR BREMSSTRAHLUNG AT 3 MEV
      DATA BCUT /0.003D0/
C  DATA FOR REST: AVERAGE ATOMIC WEIGHT, NEUTRON CONTENTS OF N,O,AR
      DATA TAR / 14.6D0 /, CONTNE / 0.5D0, 0.5D0, 0.55D0 /, LIT / 1 /
C  KAON CONSTANTS
C  REVISED NOV. 2004 BY D. HECK
      DATA CKA /     0.D0 ,  0.1D0,     0.D0 ,     0.D0,     0.D0 ,
     *               0.D0 , 0.25D0,     0.5D0,   0.75D0,     1.D0 , !10
     *               0.5D0,  0.2D0,     0.D0 ,    0.D0 ,   149.6D0,
     *             149.6D0,0.236D0,   0.206D0,  0.135D0,   0.222D0, !20
     *               0.5D0,  0.D0 ,  0.6356D0,0.692762D0, 0.8742D0,
     *             0.6783D0,.4069D0,    0.D0 ,    0.D0 ,     0.D0 , !30
     *               0.D0 ,  0.D0 ,     0.D0 ,    0.D0 ,     1.D0 ,
     *               1.0D5,  0.D0 ,     0.D0 ,    0.D0 ,     0.D0 , !40
     *               0.D0 ,  0.D0 ,     0.D0 ,    0.D0 ,     0.D0 ,
     *               0.D0 ,0.8423D0, 0.8982D0, 0.9489D0,  0.9824D0, !50
     *           -0.2154D0, 0.012D0,-0.0101D0,   1.27D0,   0.638D0,
     *             0.057D0,   0.D0 ,   1.84D0,    0.D0 ,     1.D0 , !60
     *             0.678D0, 0.076D0, 0.0099D0,   2.22D0,  0.0298D0,
     *               0.D0 ,1.288D-2, 0.0296D0,  1.96D-2,  1.194D-2, !70
     *            0.0282D0,   0.D0 , 1.310D-2, 0.0282D0,  0.0138D0,
     *            1.241D-2,   0.D0 ,    0.D0 ,    0.D0 ,     0.D0  /
C  DATA FOR ETA DECAY
C  REVISED NOV. 2004 BY D. HECK
      DATA CETA /  0.3972D0,  0.7265D0,  0.9575D0,  0.09D-2,  2.07D0 /
C  DATA FOR STRANGE BARYON DECAY
C  REVISED NOV. 2004 BY D. HECK
      DATA CSTRBA / 0.D0     , 0.D0 , 0.D0 , 0.D0 , 0.6409D0,
     *              0.5163D0 , 0.D0 , 0.D0 , 0.D0 , 0.678D0 ,
     *              0.914D0 /
C  PARPAR
C  REVISED NOV. 2004 BY D. HECK
      DATA C /6371315.D2,  6.0D5,     20.0D5,    0.D0 ,         0.D0 ,
     *           0.D0 ,    0.D0 ,      0.D0 ,    2.5D0,        2.07D0, !10
     *           8.2D0,    0.1D0,      0.D0 ,    0.D0 ,         0.D0 ,
     *           0.D0 ,    0.D0 ,      0.D0 ,    0.D0 ,         0.D0 , !20
     *          37.7D0, 1.532873D-4, 9.386417D0, 2.D-3,  29.9792458D9,
     *           1.D0 ,    0.D0 ,     1.57D0,    0.D0 ,       0.021D0, !30
     *           0.D0 ,    0.D0 ,      0.D0 ,    2.0D1,         0.D0 ,
     *           0.D0 ,    0.D0 ,      0.D0 ,    0.D0 ,         0.D0 , !40
     *           0.D0 ,    0.D0 ,      0.D0 ,    0.D0 ,         0.D0 ,
     *           0.D0 ,    0.D0 ,      0.D0 ,    0.D0 , 137.0359991D0 /!50
C  RUNPAR , STACKF
      DATA MONIIN /  5 /, MONIOU /  6 /, MPATAP / 90 /, MEXST / 96 /,
     *     MDEBUG /  6 /, NUCNUC / 11 /, MATMFI / 19 /, MDBASE / 45 /,
     *     MTABOUT / 46 /,MLONGOUT / 48 /
     *    ,LSTCK / 23 /

C  UNRELEASABLE ENERGY (REST MASS) FOR ENERGY DEPOSIT
C  THE REST MASSES FOR NUCLEI ARE SET IN PAMAF
      DATA RESTMS/       0.,  -.511D-3,   .511D-3,       0.,  .105658 ,
     *             .105658 ,        0.,        0.,       0.,        0.,  ! 10
     *                   0.,        0., .9395653 , .938272 , -.938272 ,
     *                   0.,        0.,     .938 ,    .938 ,     .938 ,  ! 20
     *                .938 ,     .938 ,     .938 ,    .938 ,-.9395653 ,
     *               -.938 ,    -.938 ,    -.938 ,   -.938 ,    -.938 ,  ! 30
     *               -.938 ,    -.938 ,        0.,       0.,        0.,
     *                   0.,        0.,        0.,       0.,        0.,  ! 40
     *                  1D9,        0.,        0.,       0.,        0.,
     *                   0.,        0.,        0.,       0.,        0.,  ! 50
     *                   0.,        0.,        0.,    .938 ,     .938 ,
     *                .938 ,     .938 ,    -.938 ,   -.938 ,    -.938 ,  ! 60
     *               -.938 ,        0.,        0.,       0.,        0.,
     *                   0.,        0.,        0.,       0.,        0.,  ! 70
     *                   0.,        0.,        0.,       0.,        0.,
     *       5925*0.D0/

C  GNUPR
C  NITROGEN TARGET  14
      DATA ((SE14(I,J),I=1,3),J=1,14)
     *                       / 0.472000D+00,-0.426710D-02, 0.726439D-04,
     *                         0.230324D+00,-0.989733D-03,-0.807077D-05,
     *                         0.138623D+00, 0.609624D-03,-0.401675D-04,
     *                         0.827139D-01, 0.135103D-02,-0.360236D-04,
     *                         0.445693D-01, 0.137582D-02,-0.137674D-04,
     *                         0.206106D-01, 0.998620D-03, 0.422867D-05,
     *                         0.792756D-02, 0.559858D-03, 0.957875D-05,
     *                         0.247793D-02, 0.247480D-03, 0.701650D-05,
     *                         0.615535D-03, 0.860096D-04, 0.324410D-05,
     *                         0.118279D-03, 0.230732D-04, 0.104282D-05,
     *                         0.169210D-04, 0.461424D-05, 0.235175D-06,
     *                         0.169481D-05, 0.647634D-06, 0.358189D-07,
     *                         0.105988D-06, 0.568994D-07, 0.332920D-08,
     *                         0.311374D-08, 0.235385D-08, 0.143213D-09/
C  OXYGEN TARGET  16
      DATA ((SE16(I,J), I=1,3),J=1,16)
     *                        /0.475002D+00,-0.434401D-02, 0.734217D-04,
     *                         0.230261D+00,-0.966152D-03,-0.982228D-05,
     *                         0.137372D+00, 0.642454D-03,-0.408490D-04,
     *                         0.813380D-01, 0.135241D-02,-0.354835D-04,
     *                         0.437870D-01, 0.135776D-02,-0.134429D-04,
     *                         0.204919D-01, 0.988538D-03, 0.398723D-05,
     *                         0.812995D-02, 0.567070D-03, 0.942943D-05,
     *                         0.269031D-02, 0.263160D-03, 0.728079D-05,
     *                         0.732711D-03, 0.993722D-04, 0.366933D-05,
     *                         0.161940D-03, 0.303662D-04, 0.134776D-05,
     *                         0.285325D-04, 0.740356D-05, 0.371648D-06,
     *                         0.390910D-05, 0.140655D-05, 0.768260D-07,
     *                         0.401145D-06, 0.200620D-06, 0.116200D-07,
     *                         0.290010D-07, 0.202033D-07, 0.121929D-08,
     *                         0.131709D-08, 0.128046D-08, 0.795482D-10,
     *                         0.282645D-10, 0.384068D-10, 0.243535D-11/
C  ARGON TARGET  40
      DATA ((SE40(I,J),I=1,3),J=1,18)
     *                       / 0.318084D+00,-0.352566D-02, 0.829469D-04,
     *                         0.193581D+00,-0.238538D-02, 0.404919D-04,
     *                         0.148699D+00,-0.118791D-02,-0.130378D-04,
     *                         0.117201D+00, 0.966097D-04,-0.536044D-04,
     *                         0.876737D-01, 0.106482D-02,-0.612882D-04,
     *                         0.600279D-01, 0.150343D-02,-0.412273D-04,
     *                         0.370180D-01, 0.147347D-02,-0.130096D-04,
     *                         0.204422D-01, 0.117625D-02, 0.743960D-05,
     *                         0.101003D-01, 0.807913D-03, 0.155153D-04,
     *                         0.447163D-02, 0.489622D-03, 0.146804D-04,
     *                         0.177806D-02, 0.265260D-03, 0.102802D-04,
     *                         0.636671D-03, 0.129412D-03, 0.591434D-05,
     *                         0.205809D-03, 0.571042D-04, 0.291674D-05,
     *                         0.601981D-04, 0.228546D-04, 0.126074D-05,
     *                         0.159631D-04, 0.831226D-05, 0.484001D-06,
     *                         0.384379D-05, 0.275100D-05, 0.166440D-06,
     *                         0.841490D-06, 0.829259D-06, 0.515615D-07,
     *                         0.167633D-06, 0.227810D-06, 0.144446D-07/
      DATA((SE40(I,J),I=1,3),J=19,36)
     *                        /0.304029D-07, 0.570494D-07, 0.366843D-08,
     *                         0.502077D-08, 0.130224D-07, 0.845876D-09,
     *                         0.754786D-09, 0.270844D-08, 0.177211D-09,
     *                         0.103229D-09, 0.512862D-09, 0.337323D-10,
     *                         0.128308D-10, 0.883149D-10, 0.583066D-11,
     *                         0.144721D-11, 0.138082D-10, 0.914113D-12,
     *                         0.147837D-12, 0.195621D-11, 0.129757D-12,
     *                         0.136429D-13, 0.250465D-12, 0.166371D-13,
     *                         0.113379D-14, 0.288894D-13, 0.192092D-14,
     *                         0.845213D-16, 0.299003D-14, 0.198959D-15,
     *                         0.562496D-17, 0.276346D-15, 0.183981D-16,
     *                         0.332222D-18, 0.226723D-16, 0.151001D-17,
     *                         0.172872D-19, 0.163915D-17, 0.109200D-18,
     *                         0.785321D-21, 0.103480D-18, 0.689517D-20,
     *                         0.307886D-22, 0.563885D-20, 0.375787D-21,
     *                         0.102630D-23, 0.261299D-21, 0.174154D-22,
     *                         0.285163D-25, 0.100944D-22, 0.672832D-24,
     *                         0.642589D-27, 0.316302D-24, 0.210839D-25/
      DATA((SE40(I,J),I=1,3),J=37,40)
     *                        /0.112817D-28, 0.772286D-26, 0.514807D-27,
     *                         0.144773D-30, 0.137838D-27, 0.918858D-29,
     *                         0.120779D-32, 0.159956D-29, 0.106632D-30,
     *                         0.491605D-35, 0.905709D-32, 0.603784D-33/
      SAVE
      END

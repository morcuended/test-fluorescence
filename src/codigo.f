# 1 "test_fluorescence_corsika.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "test_fluorescence_corsika.F"
! Test fluorescence subroutines
      PROGRAM testing
      implicit none
      integer ier, iflag, nflpht
      double precision zem, edep, lamb, zemision

      zem = 5.0d5
      edep = 1.0d1


      ! Add fluorescence initialization
      CALL FYINI
      ! Determine the number of fluorescence photons to be emitted in
      ! each step.
      ier = 0
      iflag = 2

      call fyield ( zemision, edep, nflpht, iflag, ier, 0, lamb )
      write(*,*) "Number of fluorescence photons:",nflpht
      if ( nflpht .lt. 1 ) return

      if ( ier .eq. 1 ) write(*,*)
     * "warning!, too many fluorescence photons produced at once"

      END PROGRAM testing


      subroutine fyield(zem,edep,nphots,iflag,ier,nbunch,lambda)

C Determines the number of fluorescence photons and its wavelength
c emitted for a given height and energy deposit. It depends on the
c atmospheric pressure (hPa) and temperature (K).
c It is called from fluor.
c--------------------------------------------------------------------
c Input arguments:
c--------------------------------------------------------------------
c iflag = flag indicating what the subroutine will return
c if = 2 -> returns only number of f-photons
c if = 3 -> returns wavelength of each bunch
c zem = z emission (in cm)
c edep = deposited energy (in MeV)
c--------------------------------------------------------------------
c Output parameters:
c--------------------------------------------------------------------
c nphots = number of yielded photons
c ier = flag indicating error (if n > nphmax)
c nbunch = number of produced bunches of photons
c lambda = wavelength of each produced bunch
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
     & yieldi(nfltab,maxiz), yint(maxiz)

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
c Determine the number of photons poisson distributed
c
c Interpolate except for the highest zem
c
        if (izem.lt.(maxiz-1)) then
          alfa = (zem-ztab(izem))/(ztab(izem+1)-ztab(izem))
          phtyield = (yint(izem)+alfa*(yint(izem+1)-yint(izem))) * edep
        else
          phtyield = yint(maxiz) * edep
        endif

c nphots = int(phtyield) ! Uncomment this for tests
        call mpoiss(phtyield, nphots )
        nphots = phtyield

c Warning in case of getting a very large number of photons
        if (nphots .gt. nphmax) then
            write(*,*) "zem, edep ",zem, edep
            write(*,*) "Error: ",nphmax - nphots,
     * " remaining photons not stored"
            ier = 1
            nphots = nphmax
        endif
        return
      endif

c Return only a wavelength per bunch now
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

c Initializes the fluorescence yield data
c
c--------------------------------------------------------------------
c Input arguments:
c--------------------------------------------------------------------
c
c all input arguments enter through Corsika commons
c
c--------------------------------------------------------------------
c Output parameters:
c--------------------------------------------------------------------
c
c ifwlenl = lower index of the fluorescen band to use
c ifwlenu = upper index of the fluorescen band to use
c Are initialized from wavlgl and wavlgu
c
c The fluorescence yield tables as a function of wavelength and
c height are initilized
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
     & yieldi(nfltab,maxiz), yint(maxiz)

c
      save

      data ifwlenl/-1/, ifwlenu/-1/

c Fluorescence spectrum and parameters.

      data (lamb(k), il(k), pl(k), alp(k), k=1, nfltab)
c WLen Il Pressure Alpha ! p_w
     */ 281.20d0 , 0.18d0, 19.00d0, 0.00d0,! 0.00d0,
     * 281.80d0 , 0.30d0, 20.70d0, 0.00d0,! 0.00d0,
     * 296.20d0 , 5.16d0, 18.50d0, 0.00d0,! 0.00d0,
     * 297.70d0 , 2.77d0, 17.30d0, 0.00d0,! 0.00d0,
     * 302.00d0 , 0.41d0, 21.00d0, 0.00d0,! 0.00d0,
     * 308.00d0 , 1.44d0, 21.00d0, 0.00d0,! 0.00d0,
     * 311.70d0 , 7.24d0, 18.70d0, 0.00d0,! 0.00d0,
     * 313.60d0 , 11.05d0, 12.27d0, -0.13d0,! 1.20d0,
     * 315.90d0 , 39.33d0, 11.88d0, -0.19d0,! 1.10d0,
     * 317.70d0 , 0.46d0, 21.00d0, 0.00d0,! 0.00d0,
     * 326.80d0 , 0.80d0, 19.00d0, 0.00d0,! 0.00d0,
     * 328.50d0 , 3.80d0, 20.70d0, 0.00d0,! 0.00d0,
     * 330.90d0 , 2.15d0, 16.90d0, 0.00d0,! 0.00d0,
     * 333.90d0 , 4.02d0, 15.50d0, 0.00d0,! 0.00d0,
     * 337.10d0 , 100.00d0, 15.89d0, -0.35d0,! 1.28d0,
     * 346.30d0 , 1.74d0, 21.00d0, 0.00d0,! 0.00d0,
     * 350.00d0 , 2.79d0, 15.20d0, -0.38d0,! 1.50d0,
     * 353.70d0 , 21.35d0, 12.70d0, -0.22d0,! 1.27d0,
     * 357.70d0 , 67.41d0, 15.39d0, -0.35d0,! 1.30d0,
     * 365.90d0 , 1.13d0, 21.00d0, 0.00d0,! 0.00d0,
     * 367.20d0 , 0.54d0, 19.00d0, 0.00d0,! 0.00d0,
     * 371.10d0 , 4.97d0, 14.80d0, -0.24d0,! 1.30d0,
     * 375.60d0 , 17.87d0, 12.82d0, -0.17d0,! 1.10d0,
     * 380.50d0 , 27.20d0, 16.51d0, -0.34d0,! 1.40d0,
     * 385.80d0 , 0.50d0, 19.00d0, 0.00d0,! 0.00d0,
     * 387.70d0 , 1.17d0, 7.60d0, 0.00d0,! 0.00d0,
     * 388.50d0 , 0.83d0, 3.90d0, 0.00d0,! 0.00d0,
     * 391.40d0 , 28.00d0, 2.94d0, -0.79d0,! 0.33d0,
     * 394.30d0 , 3.36d0, 13.70d0, -0.20d0,! 1.20d0,
     * 399.80d0 , 8.38d0, 13.60d0, -0.20d0,! 1.10d0,
     * 405.00d0 , 8.07d0, 17.80d0, -0.37d0,! 1.50d0,
     * 414.10d0 , 0.49d0, 19.00d0, 0.00d0,! 0.00d0,
     * 420.00d0 , 1.75d0, 13.80d0, 0.00d0,! 0.00d0,
     * 423.60d0 , 1.04d0, 3.90d0, 0.00d0,! 0.00d0,
     * 427.00d0 , 7.08d0, 6.38d0, 0.00d0,! 0.00d0,
     * 427.80d0 , 4.94d0, 2.89d0, -0.54d0,! 0.60d0,
     * 434.30d0 , 2.23d0, 15.89d0, 0.00d0,! 0.00d0,
     * 435.50d0 , 0.17d0, 19.00d0, 0.00d0,! 0.00d0,
     * 441.50d0 , 0.56d0, 20.70d0, 0.00d0,! 0.00d0,
     * 448.90d0 , 0.72d0, 12.27d0, 0.00d0,! 0.00d0,
     * 457.30d0 , 0.91d0, 11.88d0, 0.00d0,! 0.00d0,
     * 459.70d0 , 0.02d0, 3.90d0, 0.00d0,! 0.00d0,
     * 464.90d0 , 0.66d0, 3.90d0, 0.00d0,! 0.00d0,
     * 466.50d0 , 0.53d0, 15.89d0, 0.00d0,! 0.00d0,
     * 470.60d0 , 1.93d0, 2.94d0, 0.00d0,! 0.00d0,
     * 481.30d0 , 0.24d0, 12.27d0, 0.00d0,! 0.00d0,
     * 491.70d0 , 0.25d0, 11.88d0, 0.00d0,! 0.00d0,
     * 503.20d0 , 0.12d0, 15.89d0, 0.00d0,! 0.00d0,
     * 514.60d0 , 0.16d0, 3.90d0, 0.00d0,! 0.00d0,
     * 522.50d0 , 0.33d0, 2.94d0, 0.00d0,! 0.00d0,
     * 530.90d0 , 0.06d0, 11.88d0, 0.00d0,! 0.00d0,
     * 545.20d0 , 0.03d0, 15.89d0, 0.00d0,! 0.00d0,
     * 570.40d0 , 0.01d0, 3.90d0, 0.00d0,! 0.00d0,
     * 575.00d0 , 0.03d0, 3.90d0, 0.00d0,! 0.00d0,
     * 586.10d0 , 0.05d0, 2.94d0, 0.00d0,! 0.00d0,
     * 593.80d0 , 0.01d0, 15.89d0, 0.00d0,! 0.00d0,
     * 665.90d0 , 0.01d0, 2.94d0, 0.00d0/! 0.00d0/
c------------------------------------------------------------
c
c Select only the fluorescence bands within the range of wavelengths
c set by the cwavlg in the input card (same as for Cherenkov).
c The rows used are determined only once depending on the lower
c and upper limits.
c
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

c Determine the half of the fluorescence yield as we consider
c only photons moving downwards.
        yint(izem) = 0.d0
        do i = ifwlenl, ifwlenu
           pprimt = pl(i) * (tzero / temp(zem))**(alp(i) - 0.5d0)
           y = 0.5d0 * yabs * 1.d-2 * il(i)
           y = y * (1.d0 + pzero/pl(i))/(1.d0 + pressure(zem)/pprimt)
           yield(i,izem) = y
           yint(izem) = yint(izem) + y
           yieldi(i,izem) = yint(izem)
        enddo
c Normalize the fluorescence yield in each band
        do i = ifwlenl, ifwlenu
           yieldi(i,izem) = yieldi(i,izem) / yint(izem)
        enddo
c
      enddo
c
      return
      end



      function pressure( arg )

c-----------------------------------------------------------------------
c Calculate the pressure (hPa) at a given height (cm) assuming an ideal
c gas under a gravitational acceleration of 9.81 m/s**2
c This function is called from temp and fyield.
c argument:
c arg = height (cm)
c-----------------------------------------------------------------------

      implicit none
      double precision arg
      double precision pressure
      double precision thick
      external thick
      save
      double precision grav
      parameter (grav = 9.81d2 )
c-----------------------------------------------------------------------

      pressure=grav*thick(arg)*1.d-3

      end function pressure



      double precision function temp( arg )

c-----------------------------------------------------------------------
c temp(erature)
c
c Calculate the local temperature (K) of atmosphere at a certain
c height (cm).
c This function is called from fyield.
c Argument:
c arg = height (cm)
c-----------------------------------------------------------------------

      implicit none
      double precision arg
      double precision pressure,rhof
      external pressure,rhof
      save
      double precision Mmolar,Rconst !in cm3 kPa k−1 g mol−1
      parameter (Mmolar=28.96, Rconst=8.3144598d3 )
c-----------------------------------------------------------------------

      temp=Mmolar*pressure(arg)*1.d-1/(rhof(arg)*Rconst)

      end function temp


*-- Author : D. HECK IK FZK KARLSRUHE 15/10/1996
C=======================================================================

      SUBROUTINE MPOISS( AMEAN,NPRAN )

C-----------------------------------------------------------------------
C M(UON COULOMB SCATTERING) POISS(ON DISTRIBUTION)
C
C GENERATES A RANDOM NUMBER POISSON DISTRIBUTED WITH MEAN VALUE AMEAN.
C THIS SUBROUTINE IS IN ANALOGY WITH SUBROUT. GPOISS.
C (AUTHOR: L. URBAN) OF GEANT321
C SEE CERN PROGRAM LIBRARY LONG WRITEUP W5013.
C THIS SUBROUTINE IS CALLED FROM MUCOUL.
C ARGUMENTS:
C AMEAN = MEAN VALUE OF RANDOM NUMBER
C NPRAN = RANDOM NUMBER POISSON DISTRIBUTED
C-----------------------------------------------------------------------

      IMPLICIT NONE



# 1 "test_fluorescence_corsika.h" 1
# 3712 "test_fluorescence_corsika.h"
      COMMON /CRCONSTA/PI,PI2,OB3,TB3,ENEPER,SQRT3
      DOUBLE PRECISION PI,PI2,OB3,TB3,ENEPER,SQRT3
# 4403 "test_fluorescence_corsika.h"
      COMMON /CRRANDPA/RD,FAC,U1,U2,NSEQ,ISEED,KNOR
      DOUBLE PRECISION RD(3000),FAC,U1,U2



      INTEGER ISEED(3,10),NSEQ

      LOGICAL KNOR
# 4489 "test_fluorescence_corsika.h"
      COMMON /CRRUNPAR/FIXHEI,THICK0,HILOECM,HILOELB,SIG1I,TARG1I,
     * STEPFC,RCUT,RCUT2,
# 4501 "test_fluorescence_corsika.h"
     * NRRUN,NSHOW,MPATAP,MONIIN,
     * MONIOU,MDEBUG,NUCNUC,MTABOUT,MLONGOUT,
     * ISEED1I,
# 4519 "test_fluorescence_corsika.h"
     * LSTCK,

     * LSTCK1,LSTCK2,

c#if __ANAHIST__||__AUGERHIST__||__MUONHIST__
c * LUNHST,
c#endif
     * ISHOWNO,ISHW,NOPART,NRECS,NBLKS,MAXPRT,NDEBDL,
     * N1STTR,MDBASE,



     * DEBDEL,DEBUG,FDECAY,FEGS,FIRSTI,FIXINC,FIXTAR,
     * FIX1I,FMUADD,FNKG,FPRINT,FDBASE,FPAROUT,FTABOUT,
     * FLONGOUT,GHEISH,GHESIG,GHEISDB,USELOW,TMARGIN




     * ,FOUTFILE,IFINAM
# 4563 "test_fluorescence_corsika.h"
      COMMON /CRRUNPAC/DATDIR,DSN,DSNTAB,DSNLONG,HOST,USER
# 4580 "test_fluorescence_corsika.h"
     * ,FILOUT
# 4591 "test_fluorescence_corsika.h"
      DOUBLE PRECISION FIXHEI,THICK0,HILOECM,HILOELB,SIG1I,TARG1I,
     * STEPFC,RCUT,RCUT2






      INTEGER NRRUN,NSHOW,MPATAP,MONIIN,MONIOU,MDEBUG,NUCNUC,
     * ISHOWNO,ISHW,NOPART,NRECS,NBLKS,MAXPRT,NDEBDL,
     * N1STTR,MDBASE,MTABOUT,MLONGOUT,ISEED1I(3)
# 4623 "test_fluorescence_corsika.h"
      INTEGER LSTCK

     * ,LSTCK1,LSTCK2
# 4634 "test_fluorescence_corsika.h"
      CHARACTER*132 FILOUT

      CHARACTER*255 DSN,DSNTAB,DSNLONG
      CHARACTER*132 DATDIR
      CHARACTER*60 HOST,USER
# 4654 "test_fluorescence_corsika.h"
      LOGICAL DEBDEL,DEBUG,FDECAY,FEGS,FIRSTI,FIXINC,FIXTAR,
     * FIX1I,FMUADD,FNKG,FPRINT,FDBASE,FPAROUT,FTABOUT,
     * FLONGOUT,GHEISH,GHESIG,GHEISDB,USELOW,TMARGIN
# 4667 "test_fluorescence_corsika.h"
      LOGICAL FOUTFILE
      INTEGER IFINAM
# 388 "test_fluorescence_corsika.F" 2

      DOUBLE PRECISION AMEAN,AN,HMXINT,P,PLIM,RR,S,X
      INTEGER NPRAN
      SAVE
      DATA PLIM / 16.D0 /, HMXINT / 2.D9 /
C-----------------------------------------------------------------------

C PROTECTION AGAINST NEGATIVE MEAN VALUES
      AN = 0.D0
      IF ( AMEAN .GT. 0.D0 ) THEN
        IF ( AMEAN .LE. PLIM ) THEN
          CALL RMMARD( RD,1,1 )
          P = EXP( -AMEAN )
          S = P
          IF ( RD(1) .LE. S ) GOTO 20
 10 AN = AN + 1.D0
          P = P * AMEAN / AN
          S = S + P
          IF ( S .LT. RD(1) .AND. P .GT. 1.D-30 ) GOTO 10
        ELSE
          CALL RMMARD( RD,2,1 )
          RR = SQRT( (-2.D0)*LOG( RD(1) ) )
          X = RR * COS( PI2 * RD(2) )
          AN = MIN( MAX( 0.D0, AMEAN+X*SQRT( AMEAN ) ), HMXINT )
        ENDIF
      ENDIF
 20 NPRAN = NINT(AN)

      RETURN
      END


*-- Author : D. HECK IK FZK KARLSRUHE 17/03/2003
C=======================================================================

      SUBROUTINE RMMARD( RVEC,LENV,ISEQ )

C-----------------------------------------------------------------------
C R(ANDO)M (NUMBER GENERATOR OF) MAR(SAGLIA TYPE) D(OUBLE PRECISION)
C
C THESE ROUTINES (RMMARD,RMMAQD) ARE MODIFIED VERSIONS OF ROUTINES
C FROM THE CERN LIBRARIES. DESCRIPTION OF ALGORITHM SEE:
C http:
C IT HAS BEEN CHECKED THAT RESULTS ARE BIT-IDENTICAL WITH CERN
C DOUBLE PRECISION RANDOM NUMBER GENERATOR RMM48, DESCRIBED IN
C http:
C ARGUMENTS:
C RVEC = DOUBLE PREC. VECTOR FIELD TO BE FILLED WITH RANDOM NUMBERS
C LENV = LENGTH OF VECTOR (# OF RANDNUMBERS TO BE GENERATED)
C ISEQ = # OF RANDOM SEQUENCE
C
C VERSION OF D. HECK FOR DOUBLE PRECISION RANDOM NUMBERS.
C-----------------------------------------------------------------------

      IMPLICIT NONE


# 1 "test_fluorescence_corsika.h" 1
# 4416 "test_fluorescence_corsika.h"
      INTEGER KSEQ



      PARAMETER (KSEQ = 9)

      COMMON /CRRANMA3/CD,CINT,CM,TWOM24,TWOM48,MODCNS
      DOUBLE PRECISION CD,CINT,CM,TWOM24,TWOM48
      INTEGER MODCNS





      COMMON /CRRANMA4/C,U,IJKL,I97,J97,NTOT,NTOT2,JSEQ
      DOUBLE PRECISION C(KSEQ),U(97,KSEQ),UNI
      INTEGER IJKL(KSEQ),I97(KSEQ),J97(KSEQ),
     * NTOT(KSEQ),NTOT2(KSEQ),JSEQ
# 446 "test_fluorescence_corsika.F" 2

      DOUBLE PRECISION RVEC(*)
      INTEGER ISEQ,IVEC,LENV
      SAVE
C-----------------------------------------------------------------------

      IF ( ISEQ .GT. 0 .AND. ISEQ .LE. KSEQ ) JSEQ = ISEQ

      DO IVEC = 1, LENV
        UNI = U(I97(JSEQ),JSEQ) - U(J97(JSEQ),JSEQ)
        IF ( UNI .LT. 0.D0 ) UNI = UNI + 1.D0
        U(I97(JSEQ),JSEQ) = UNI
        I97(JSEQ) = I97(JSEQ) - 1
        IF ( I97(JSEQ) .EQ. 0 ) I97(JSEQ) = 97
        J97(JSEQ) = J97(JSEQ) - 1
        IF ( J97(JSEQ) .EQ. 0 ) J97(JSEQ) = 97
        C(JSEQ) = C(JSEQ) - CD
        IF ( C(JSEQ) .LT. 0.D0 ) C(JSEQ) = C(JSEQ) + CM
        UNI = UNI - C(JSEQ)
        IF ( UNI .LT. 0.D0 ) UNI = UNI + 1.D0
C AN EXACT ZERO HERE IS VERY UNLIKELY, BUT LET''S BE SAFE.
        IF ( UNI .EQ. 0.D0 ) UNI = TWOM48
        RVEC(IVEC) = UNI
      ENDDO

      NTOT(JSEQ) = NTOT(JSEQ) + LENV
      IF ( NTOT(JSEQ) .GE. MODCNS ) THEN
        NTOT2(JSEQ) = NTOT2(JSEQ) + 1
        NTOT(JSEQ) = NTOT(JSEQ) - MODCNS
      ENDIF

      RETURN
      END


*-- Author : The CORSIKA development group 21/04/1994
C=======================================================================

      DOUBLE PRECISION FUNCTION RHOF( ARG )

C-----------------------------------------------------------------------
C RHO (DENSITY) F(UNCTION)
C
C CALCULATES DENSITY (G/CM**3) OF ATMOSPHERE DEPENDING ON HEIGHT (CM)
C THIS FUNCTION IS CALLED FROM BOX2, LPMEFFECT, ININKG, CERENK,
C MUTRAC, AND INRTAB.
C ARGUMENT:
C ARG = HEIGHT (CM)
C-----------------------------------------------------------------------

      IMPLICIT NONE



# 1 "test_fluorescence_corsika.h" 1
# 3568 "test_fluorescence_corsika.h"
      COMMON /CRATMOS/ AATM,AATM0,BATM,BATM0,CATM,CATM0,DATM,MODATM
     * ,MATMFI,LATMNEW
      DOUBLE PRECISION AATM(5),AATM0(5,0:42),BATM(5),BATM0(5,0:42),
     * CATM(5),CATM0(5,0:42),DATM(5)
      INTEGER MODATM,MATMFI
      LOGICAL LATMNEW





      COMMON /CRATMOS2/HLAY,HLAY0,THICKL,LAYNO,LAYNEW
      DOUBLE PRECISION HLAY(6),HLAY0(5,0:28),THICKL(5)
      INTEGER LAYNO(0:41)
      LOGICAL LAYNEW
# 4489 "test_fluorescence_corsika.h"
      COMMON /CRRUNPAR/FIXHEI,THICK0,HILOECM,HILOELB,SIG1I,TARG1I,
     * STEPFC,RCUT,RCUT2,
# 4501 "test_fluorescence_corsika.h"
     * NRRUN,NSHOW,MPATAP,MONIIN,
     * MONIOU,MDEBUG,NUCNUC,MTABOUT,MLONGOUT,
     * ISEED1I,
# 4519 "test_fluorescence_corsika.h"
     * LSTCK,

     * LSTCK1,LSTCK2,

c#if __ANAHIST__||__AUGERHIST__||__MUONHIST__
c * LUNHST,
c#endif
     * ISHOWNO,ISHW,NOPART,NRECS,NBLKS,MAXPRT,NDEBDL,
     * N1STTR,MDBASE,



     * DEBDEL,DEBUG,FDECAY,FEGS,FIRSTI,FIXINC,FIXTAR,
     * FIX1I,FMUADD,FNKG,FPRINT,FDBASE,FPAROUT,FTABOUT,
     * FLONGOUT,GHEISH,GHESIG,GHEISDB,USELOW,TMARGIN




     * ,FOUTFILE,IFINAM
# 4563 "test_fluorescence_corsika.h"
      COMMON /CRRUNPAC/DATDIR,DSN,DSNTAB,DSNLONG,HOST,USER
# 4580 "test_fluorescence_corsika.h"
     * ,FILOUT
# 4591 "test_fluorescence_corsika.h"
      DOUBLE PRECISION FIXHEI,THICK0,HILOECM,HILOELB,SIG1I,TARG1I,
     * STEPFC,RCUT,RCUT2






      INTEGER NRRUN,NSHOW,MPATAP,MONIIN,MONIOU,MDEBUG,NUCNUC,
     * ISHOWNO,ISHW,NOPART,NRECS,NBLKS,MAXPRT,NDEBDL,
     * N1STTR,MDBASE,MTABOUT,MLONGOUT,ISEED1I(3)
# 4623 "test_fluorescence_corsika.h"
      INTEGER LSTCK

     * ,LSTCK1,LSTCK2
# 4634 "test_fluorescence_corsika.h"
      CHARACTER*132 FILOUT

      CHARACTER*255 DSN,DSNTAB,DSNLONG
      CHARACTER*132 DATDIR
      CHARACTER*60 HOST,USER
# 4654 "test_fluorescence_corsika.h"
      LOGICAL DEBDEL,DEBUG,FDECAY,FEGS,FIRSTI,FIXINC,FIXTAR,
     * FIX1I,FMUADD,FNKG,FPRINT,FDBASE,FPAROUT,FTABOUT,
     * FLONGOUT,GHEISH,GHESIG,GHEISDB,USELOW,TMARGIN
# 4667 "test_fluorescence_corsika.h"
      LOGICAL FOUTFILE
      INTEGER IFINAM
# 501 "test_fluorescence_corsika.F" 2

      DOUBLE PRECISION ARG
      SAVE
C-----------------------------------------------------------------------


      IF ( ARG .LT. HLAY(2) ) THEN
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

      RETURN
      END


*-- Author : The CORSIKA development group 21/04/1994
C=======================================================================

      DOUBLE PRECISION FUNCTION THICK( ARG )

C-----------------------------------------------------------------------
C THICK(NESS OF ATMOSPHERE)
C
C CALCULATES THICKNESS (G/CM**2) OF ATMOSPHERE DEPENDING ON HEIGHT (CM)
C THIS FUNCTION IS CALLED FROM AAMAIN, BOX2, BOX3, EM, INPRM, MUBREM,
C MUDECY, MUPRPR, MUTRAC, NRANGC, NUCINT, PRANGC, START, UPDATC,
C UPDATE, EGS4, ELECTR, HOWFAR, PHOTON, ININKG, NKG, AND CERENK.
C ARGUMENT:
C ARG = HEIGHT (CM)
C-----------------------------------------------------------------------

      IMPLICIT NONE



# 1 "test_fluorescence_corsika.h" 1
# 3568 "test_fluorescence_corsika.h"
      COMMON /CRATMOS/ AATM,AATM0,BATM,BATM0,CATM,CATM0,DATM,MODATM
     * ,MATMFI,LATMNEW
      DOUBLE PRECISION AATM(5),AATM0(5,0:42),BATM(5),BATM0(5,0:42),
     * CATM(5),CATM0(5,0:42),DATM(5)
      INTEGER MODATM,MATMFI
      LOGICAL LATMNEW





      COMMON /CRATMOS2/HLAY,HLAY0,THICKL,LAYNO,LAYNEW
      DOUBLE PRECISION HLAY(6),HLAY0(5,0:28),THICKL(5)
      INTEGER LAYNO(0:41)
      LOGICAL LAYNEW
# 4489 "test_fluorescence_corsika.h"
      COMMON /CRRUNPAR/FIXHEI,THICK0,HILOECM,HILOELB,SIG1I,TARG1I,
     * STEPFC,RCUT,RCUT2,
# 4501 "test_fluorescence_corsika.h"
     * NRRUN,NSHOW,MPATAP,MONIIN,
     * MONIOU,MDEBUG,NUCNUC,MTABOUT,MLONGOUT,
     * ISEED1I,
# 4519 "test_fluorescence_corsika.h"
     * LSTCK,

     * LSTCK1,LSTCK2,

c#if __ANAHIST__||__AUGERHIST__||__MUONHIST__
c * LUNHST,
c#endif
     * ISHOWNO,ISHW,NOPART,NRECS,NBLKS,MAXPRT,NDEBDL,
     * N1STTR,MDBASE,



     * DEBDEL,DEBUG,FDECAY,FEGS,FIRSTI,FIXINC,FIXTAR,
     * FIX1I,FMUADD,FNKG,FPRINT,FDBASE,FPAROUT,FTABOUT,
     * FLONGOUT,GHEISH,GHESIG,GHEISDB,USELOW,TMARGIN




     * ,FOUTFILE,IFINAM
# 4563 "test_fluorescence_corsika.h"
      COMMON /CRRUNPAC/DATDIR,DSN,DSNTAB,DSNLONG,HOST,USER
# 4580 "test_fluorescence_corsika.h"
     * ,FILOUT
# 4591 "test_fluorescence_corsika.h"
      DOUBLE PRECISION FIXHEI,THICK0,HILOECM,HILOELB,SIG1I,TARG1I,
     * STEPFC,RCUT,RCUT2






      INTEGER NRRUN,NSHOW,MPATAP,MONIIN,MONIOU,MDEBUG,NUCNUC,
     * ISHOWNO,ISHW,NOPART,NRECS,NBLKS,MAXPRT,NDEBDL,
     * N1STTR,MDBASE,MTABOUT,MLONGOUT,ISEED1I(3)
# 4623 "test_fluorescence_corsika.h"
      INTEGER LSTCK

     * ,LSTCK1,LSTCK2
# 4634 "test_fluorescence_corsika.h"
      CHARACTER*132 FILOUT

      CHARACTER*255 DSN,DSNTAB,DSNLONG
      CHARACTER*132 DATDIR
      CHARACTER*60 HOST,USER
# 4654 "test_fluorescence_corsika.h"
      LOGICAL DEBDEL,DEBUG,FDECAY,FEGS,FIRSTI,FIXINC,FIXTAR,
     * FIX1I,FMUADD,FNKG,FPRINT,FDBASE,FPAROUT,FTABOUT,
     * FLONGOUT,GHEISH,GHESIG,GHEISDB,USELOW,TMARGIN
# 4667 "test_fluorescence_corsika.h"
      LOGICAL FOUTFILE
      INTEGER IFINAM
# 544 "test_fluorescence_corsika.F" 2

      DOUBLE PRECISION ARG
      SAVE
C-----------------------------------------------------------------------


      IF ( ARG .LT. HLAY(2) ) THEN
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

      RETURN
      END

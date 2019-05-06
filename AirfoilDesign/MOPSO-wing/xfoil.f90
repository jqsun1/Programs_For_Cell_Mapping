!-Keeping just the necessary stuff to do the airfoil calculations
SUBROUTINE XFOILfun(CONVERGEDflg, CLexport, CDexport, CMexport, Xpoints,Ypoints,Reynolds,Alphas,npoints,samples)
	INCLUDE 'XFOIL.INC'
	integer npoints, samples, nonConvCounter
	real*8 Xpoints(npoints), Ypoints(npoints), Reynolds, Alphas(samples), CLexport(samples), CDexport(samples), CMexport(samples), Alpha
	logical*1 CONVERGEDflg(samples);
	real ratio
	nonConvCounter = 0
	ratio = 0
	DBUGMOD  = .FALSE.
	CALL INIT	
	REINF1 = Reynolds
	
	
	!NB = 2*npoints-1
	NB = npoints
	
!	print*,'Ready to start the Fortran part?','npoints=',npoints,', Npoints=', NB
	do 123 icounter = 1, NB
		XB(icounter) = Real(Xpoints(icounter))
		YB(icounter) = Real(Ypoints(icounter))
		if (DBUGMOD) print*,icounter,XB(icounter),YB(icounter) 
	123 continue
	
	!print*,'What we got so far:'
	!print*,'  Reynolds= ', Reynolds
	!print*,'Alfa= ', ADEG
	!print*,'npoints= ',npoints
		
		


	
!---- max panel angle threshold for warning
    DATA ANGTOL / 40.0 /
    NSIDE = npoints

CONVERGED = .TRUE.
	CALL LOAD    
    CALL ABCOPY(.TRUE.)
if (.NOT.CONVERGED) RETURN ! return if there is a problem in LOAD and ABCOPY
	LU = 8
	!trying to input an airfoil for instance
!	CALL NACA(2104)
      DATA EPS1 / 1.0E-4 /

NITER = NITER1
if (DBUGMOD) WRITE(*,*) 'Alfa	    CL 		     CD		     CL/CD		 CM'
9911 FORMAT(F5.2,7X, F8.5,7X, F8.5,7X,F8.5,7X,F8.5,4X,L2)


!Alpha = Alphas(1)
!ADEG = Real(Alpha)
!ALFA = DTOR*ADEG
!
!CALL OPER
!if (DBUGMOD) WRITE(*,9911) ADEG, CL, CD, CL/CD, CM,CONVERGED
!CLexport(1) = CL
!CDexport(1) = CD
!CMexport(1) = CM
!CONVERGEDflg(1) = CONVERGED
!ratio = CL/CD


do 12345 icounter = 1, samples
	Alpha = Alphas(icounter)
	ADEG = Real(Alpha)
	!print *,'*********Now Dealing with Alpha= ',ADEG
	ALFA = DTOR*ADEG
	CONVERGED = .TRUE.
	CALL OPER
	if (DBUGMOD) WRITE(*,9911) ADEG, CL, CD, CL/CD, CM ,CONVERGED
	CLexport(icounter) = CL
	CDexport(icounter) = CD
	CMexport(icounter) = CM
	CONVERGEDflg(icounter) = CONVERGED
	if (CONVERGED) THEN
		IF (CL/CD<ratio) THEN
			RETURN
		ELSE
			ratio = CL/CD
		ENDIF
	ELSE
		nonConvCounter = nonConvCounter+1
		IF (nonConvCounter>3) RETURN
	ENDIF
12345 continue	
	
	
!##	USE params
	
!---- default paneling parameters
!	N = NPAN




	
!!!	print*,'For this simulation, the results are: Cl = ', CL, &
!!!		',  Cd = ',CD, &
!!!		',  Alfa = ',ADEG
	
RETURN
END






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE INIT
!---------------------------------------------------
!     Variable initialization/default routine.
!     See file XFOIL.INC for variable description.
!---------------------------------------------------
      INCLUDE 'XFOIL.INC'
!
      PI = 4.0*ATAN(1.0)
      HOPI = 0.50/PI
      QOPI = 0.25/PI
      DTOR = PI/180.0
!
!---- default Cp/Cv (air)
      GAMMA = 1.4
      GAMM1 = GAMMA - 1.0
!
!---- set unity freestream speed
      QINF = 1.0
!
!---- initialize freestream Mach number to zero
      MATYP = 1
      MINF1 = 0.
!
      ALFA = 0.0
      COSA = 1.0
      SINA = 0.0
!
      DO 10 I=1, IQX
        GAMU(I,1) = 0.
        GAMU(I,2) = 0.
        GAM(I) = 0.
        GAM_A(I) = 0.
   10 CONTINUE
      PSIO = 0.
!
      CL = 0.
      CM = 0.
      CD = 0.
!
      SIGTE = 0.0
      GAMTE = 0.0
      SIGTE_A = 0.
      GAMTE_A = 0.
!
      DO 20 I=1, IZX
        SIG(I) = 0.
   20 CONTINUE
!
      NQSP = 0
      DO 30 K=1, IPX
        ALQSP(K) = 0.
        CLQSP(K) = 0.
        CMQSP(K) = 0.
        DO 302 I=1, IBX
          QSPEC(I,K) = 0.
 302    CONTINUE
 30   CONTINUE
!
      AWAKE = 0.0
      AVISC = 0.0
!
      KIMAGE = 1
      YIMAGE = -10.0
      LIMAGE = .FALSE.
!
!---- output coordinate file delimiters
!-    KDELIM =  0  blanks
!-           =  1  commas
!-           =  2  tabs
      KDELIM = 0
!
      LGAMU  = .FALSE.
      LQINU  = .FALSE.
      LVISC  = .FALSE.
      LWAKE  = .FALSE.
      LPACC  = .FALSE.
      LBLINI = .FALSE.
      LIPAN  = .FALSE.
      LQAIJ  = .FALSE.
      LADIJ  = .FALSE.
      LWDIJ  = .FALSE.
      LCPXX  = .FALSE.
      LQVDES = .FALSE.
      LQSPEC = .FALSE.
      LQREFL = .FALSE.
      LVCONV = .FALSE.
      LCPREF = .FALSE.
      LFOREF = .FALSE.
      LPFILE = .FALSE.
      LPFILX = .FALSE.
      LPPSHO = .FALSE.
      LBFLAP = .FALSE.
      LFLAP  = .FALSE.
      LEIW   = .FALSE.
      LSCINI = .FALSE.
      LPLOT  = .FALSE.
      LCLIP  = .FALSE.
      LVLAB  = .TRUE.
      LCMINP = .FALSE.
      LHMOMP = .FALSE.
      LFREQP = .TRUE.
!
      LCURS  = .TRUE.
      LLAND  = .TRUE.
      LGSAME = .FALSE.
!
      LGPARM = .TRUE.
      LPLCAM = .FALSE.
!
!---- input airfoil will not be normalized
      LNORM = .FALSE.
!
!---- airfoil will not be forced symmetric
      LQSYM = .FALSE.
      LGSYM = .FALSE.
!
!---- endpoint slopes will be matched
      LQSLOP = .TRUE.
      LGSLOP = .TRUE.
      LCSLOP = .TRUE.
!
!---- grids on Qspec(s) and buffer airfoil geometry plots will be plotted
      LQGRID = .TRUE.
      LGGRID = .TRUE.
      LGTICK = .TRUE.
!
!---- no grid on Cp plots
      LCPGRD = .FALSE.
!
!---- overlay inviscid Cp on viscous Cp(x) plot
      LCPINV = .TRUE.

!---- grid and no symbols are to be used on BL variable plots
      LBLGRD = .TRUE.
      LBLSYM = .FALSE.
!
!---- buffer and current airfoil flap hinge coordinates
      XBF = 0.0
      YBF = 0.0
      XOF = 0.0
      YOF = 0.0
!
      NCPREF = 0
!                                               n
!---- circle plane array size (257, or largest 2  + 1 that will fit array size)
      ANN = LOG(FLOAT((2*IQX)-1))/LOG(2.0)
      NN = INT( ANN + 0.00001 )
      NC1 = 2**NN + 1
      NC1 = MIN( NC1 , 257 )
!
!---- default paneling parameters
      NPAN = 160
      CVPAR = 1.0
      CTERAT = 0.15
      CTRRAT = 0.2
!
!---- default paneling refinement zone x/c endpoints
      XSREF1 = 1.0
      XSREF2 = 1.0
      XPREF1 = 1.0
      XPREF2 = 1.0
!
!---- no polars present to begin with
      NPOL = 0
      IPACT = 0
      DO IP = 1, NPX
        PFNAME(IP) = ' '
        PFNAMX(IP) = ' '
      ENDDO
!
!---- no reference polars
      NPOLREF = 0
!
!---- plot aspect ratio, character size
      PLOTAR = 0.55
      CH     = 0.015
!
!---- airfoil node tick-mark size (as fraction of arc length)
      GTICK = 0.0005
!
!---- Cp limits in  Cp vs x  plot
      CPMAX =  1.0
      CPMIN = -2.0
      CPDEL = -0.5
      PFAC = PLOTAR/(CPMAX-CPMIN)
!
!---- Ue limits in  Ue vs x  plot
      UEMAX =  1.8
      UEMIN = -1.0
      UEDEL =  0.2
      UFAC = PLOTAR/(UEMAX-UEMIN)
!
!---- DCp limits in CAMB loading plot
      YPMIN = -0.4
      YPMAX =  0.4
!
!---- scaling factor for Cp vector plot
      VFAC = 0.25
!
!---- offsets and scale factor for airfoil in  Cp vs x  plot
      XOFAIR = 0.09
      YOFAIR = -.01
      FACAIR = 0.70
!
!---- u/Qinf scale factor for profile plotting
      UPRWT = 0.02
!
!---- polar plot options, grid, list, legend, no CDW
      LPGRID = .TRUE.
      LPCDW  = .FALSE.
      LPLIST = .TRUE.
      LPLEGN = .TRUE.
      LAECEN = .FALSE.
      LPCDH   = .FALSE.
      LPCMDOT = .FALSE.
!
!---- axis limits and annotation deltas for polar plot
      CPOLPLF(1,ICD) = 0.0
      CPOLPLF(2,ICD) = 0.04
      CPOLPLF(3,ICD) = 0.01
!        
      CPOLPLF(1,ICL) = 0.
      CPOLPLF(2,ICL) = 1.5
      CPOLPLF(3,ICL) = 0.5
!        
      CPOLPLF(1,ICM) = -0.25
      CPOLPLF(2,ICM) =  0.0
      CPOLPLF(3,ICM) =  0.05
!        
      CPOLPLF(1,IAL) = -4.0
      CPOLPLF(2,IAL) = 10.0
      CPOLPLF(3,IAL) =  2.0
!
!---- widths of plot boxes in polar plot page
      XCDWID = 0.45
      XALWID = 0.25
      XOCWID = 0.20
!
!---- line style and color index for each polar
!
!     1  *****************************  SOLID
!     2  **** **** **** **** **** ****  LONG DASHED
!     3  ** ** ** ** ** ** ** ** ** **  SHORT DASHED
!     4  * * * * * * * * * * * * * * *  DOTTED
!     5  ***** * ***** * ***** * *****  DASH-DOT
!     6  ***** * * ***** * * ***** * *  DASH-DOT-DOT
!     7  ***** * * * ***** * * * *****  DASH-DOT-DOT-DOT
!     8  **** **** * * **** **** * *    DASH-DASH-DOT-DOT
!
!     3  red
!     4  orange
!     5  yellow
!     6  green
!     7  cyan
!     8  blue
!     9  violet
!    10  magenta
!
      DO IP=1, NPX
!c        ILINP(IP) = 1 + MOD(IP-1,8)
!c        ICOLP(IP) = 3 + MOD(IP-1,8)
!
!------ normally solid, going to dashed after IP=7
        ILINP(IP) = 1 + (IP-1)/7
!
!------ skip yellow (hard to see on white background)
        ICOLP(IP) = 3 + MOD(IP-1,7)
        IF(ICOLP(IP) .GE. 5) ICOLP(IP) = ICOLP(IP) + 1
      ENDDO
!
!---- polar variables to be written to polar save file
      IPOL(1) = IAL
      IPOL(2) = ICL
      IPOL(3) = ICD
      IPOL(4) = ICP
      IPOL(5) = ICM
      NIPOL = 5
      NIPOL0 = 5
!
      JPOL(1) = JTN
      JPOL(2) = JTI
      NJPOL = 2
!
!---- default Cm reference location
      XCMREF = 0.25
      YCMREF = 0.
!
!---- default viscous parameters
      RETYP = 1
      REINF1 = 0.
      ACRIT(1) = 9.0
      ACRIT(2) = 9.0
      XSTRIP(1) = 1.0
      XSTRIP(2) = 1.0
      XOCTR(1) = 1.0
      XOCTR(2) = 1.0
      YOCTR(1) = 0.
      YOCTR(2) = 0.
      WAKLEN = 1.0
!
      IDAMP = 0
!
!---- select default BL plotting coordinate (can be changed in VPLO)
      IXBLP = 1   !  x
!cc   IXBLP = 2   !  s
!
!---- set BL calibration parameters
      CALL BLPINI
!
!---- Newton iteration limit
      ITMAX = 20
!
!---- max number of unconverged sequence points for early exit
      NSEQEX = 4
!
!---- drop tolerance for BL system solver
      VACCEL = 0.01
!
!---- inverse-mapping auto-filter level
      FFILT = 0.0
!
!---- default overlay airfoil filename
      ONAME = ' '
!
!---- default filename prefix
      PREFIX = ' '
!
!---- Plotting flag
      IDEV = 1   ! X11 window only
!     IDEV = 2                  ! B&W PostScript output file only (no color)
!     IDEV = 3   ! both X11 and B&W PostScript file
!     IDEV = 4   ! Color PostScript output file only 
!     IDEV = 5   ! both X11 and Color PostScript file 
!
!---- Re-plotting flag (for hardcopy)
!     IDEVRP = 2   ! B&W PostScript
      IDEVRP = 4   ! Color PostScript
!
!---- PostScript output logical unit and file specification
      IPSLU = 0  ! output to file  plot.ps   on LU 4    (default case)
!     IPSLU = ?  ! output to file  plot?.ps  on LU 10+?
!
!---- screen fraction taken up by plot window upon opening
      SCRNFR = 0.80
!
!---- Default plot size in inches
!-    (Default plot window is 11.0 x 8.5)
!-   (Must be smaller than XPAGE if objects are to fit on paper page)
      SIZE = 10.0

!---- plot-window dimensions in inches for plot blowup calculations
!-    currently,  11.0 x 8.5  default window is hard-wired in libPlt
      XPAGE = 11.0
      YPAGE = 8.5
!
!---- page margins in inches
      XMARG = 0.0
      YMARG = 0.0
!
!---- set top and bottom-side colors
!c      ICOLS(1) = 5
!c      ICOLS(2) = 7
      ICOLS(1) = 8
      ICOLS(2) = 3
!
!   3  red
!   4  orange
!   5  yellow
!   6  green
!   7  cyan
!   8  blue
!   9  violet
!  10  magenta
!
!
!##      CALL PLINITIALIZE
!
!---- set up color spectrum
      NCOLOR = 64
!##      CALL COLORSPECTRUMHUES(NCOLOR,'RYGCBM')
!
!
      NNAME  = 32
      NAME   = '                                '
!CC             12345678901234567890123456789012
!
!---- MSES domain parameters (not used in XFOIL)
      ISPARS = ' -2.0  3.0  -2.5  3.5'
!
!---- set MINF, REINF, based on current CL-dependence
      CALL MRCL(1.0,MINF_CL,REINF_CL)
!
!---- set various compressibility parameters from MINF
      CALL COMSET
!
      RETURN
      END ! INIT





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE CPCALC(N,Q,QINF,MINF,CP)
!---------------------------------------------
!     Sets compressible Cp from speed.
!---------------------------------------------
      DIMENSION Q(N),CP(N)
      REAL MINF
!
      LOGICAL DENNEG
!
      BETA = SQRT(1.0 - MINF**2)
      BFAC = 0.5*MINF**2 / (1.0 + BETA)
!
      DENNEG = .FALSE.
!
      DO 20 I=1, N
        CPINC = 1.0 - (Q(I)/QINF)**2
        DEN = BETA + BFAC*CPINC
        CP(I) = CPINC / DEN
        IF(DEN .LE. 0.0) DENNEG = .TRUE.
  20  CONTINUE
!
      IF(DENNEG) THEN
       WRITE(*,*)
       WRITE(*,*) 'CPCALC: Local speed too large. ',&
                 'Compressibility corrections invalid.'
      ENDIF
!
      RETURN
      END ! CPCALC

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, &
                       XREF,YREF, &
                       CL,CM,CDP, CL_ALF,CL_MSQ)
!-----------------------------------------------------------
!     Integrates surface pressures to get CL and CM.
!     Integrates skin friction to get CDF.
!     Calculates dCL/dAlpha for prescribed-CL routines.
!-----------------------------------------------------------
      DIMENSION X(N),Y(N), GAM(N), GAM_A(N)
      REAL MINF
!
!!!---- moment-reference coordinates
!!      XREF = 0.25
!!      YREF = 0.
!
      SA = SIN(ALFA)
      CA = COS(ALFA)
!
      BETA = SQRT(1.0 - MINF**2)
      BETA_MSQ = -0.5/BETA
!
      BFAC     = 0.5*MINF**2 / (1.0 + BETA)
      BFAC_MSQ = 0.5         / (1.0 + BETA) &
               - BFAC        / (1.0 + BETA) * BETA_MSQ
!
      CL = 0.0
      CM = 0.0

      CDP = 0.0
!
      CL_ALF = 0.
      CL_MSQ = 0.
!
      I = 1
      CGINC = 1.0 - (GAM(I)/QINF)**2
      CPG1     = CGINC/(BETA + BFAC*CGINC)
      CPG1_MSQ = -CPG1/(BETA + BFAC*CGINC)*(BETA_MSQ + BFAC_MSQ*CGINC)
!
      CPI_GAM = -2.0*GAM(I)/QINF**2
      CPC_CPI = (1.0 - BFAC*CPG1)/ (BETA + BFAC*CGINC)
      CPG1_ALF = CPC_CPI*CPI_GAM*GAM_A(I)
!
      DO 10 I=1, N
        IP = I+1
        IF(I.EQ.N) IP = 1
!
        CGINC = 1.0 - (GAM(IP)/QINF)**2
        CPG2     = CGINC/(BETA + BFAC*CGINC)
        CPG2_MSQ = -CPG2/(BETA + BFAC*CGINC)*(BETA_MSQ + BFAC_MSQ*CGINC)
!
        CPI_GAM = -2.0*GAM(IP)/QINF**2
        CPC_CPI = (1.0 - BFAC*CPG2)/ (BETA + BFAC*CGINC)
        CPG2_ALF = CPC_CPI*CPI_GAM*GAM_A(IP)
!
        DX = (X(IP) - X(I))*CA + (Y(IP) - Y(I))*SA
        DY = (Y(IP) - Y(I))*CA - (X(IP) - X(I))*SA
        DG = CPG2 - CPG1
!
        AX = (0.5*(X(IP)+X(I))-XREF)*CA + (0.5*(Y(IP)+Y(I))-YREF)*SA
        AY = (0.5*(Y(IP)+Y(I))-YREF)*CA - (0.5*(X(IP)+X(I))-XREF)*SA
        AG = 0.5*(CPG2 + CPG1)
!
        DX_ALF = -(X(IP) - X(I))*SA + (Y(IP) - Y(I))*CA
        AG_ALF = 0.5*(CPG2_ALF + CPG1_ALF)
        AG_MSQ = 0.5*(CPG2_MSQ + CPG1_MSQ)
!
        CL     = CL     + DX* AG
        CDP    = CDP    - DY* AG
        CM     = CM     - DX*(AG*AX + DG*DX/12.0) &
                        - DY*(AG*AY + DG*DY/12.0)
!
        CL_ALF = CL_ALF + DX*AG_ALF + AG*DX_ALF
        CL_MSQ = CL_MSQ + DX*AG_MSQ
!
        CPG1 = CPG2
        CPG1_ALF = CPG2_ALF
        CPG1_MSQ = CPG2_MSQ
   10 CONTINUE
!
      RETURN
      END ! CLCALC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE CDCALC
      INCLUDE 'XFOIL.INC'
!
      SA = SIN(ALFA)
      CA = COS(ALFA)
!
      IF(LVISC .AND. LBLINI) THEN
!
!----- set variables at the end of the wake
       THWAKE = THET(NBL(2),2)
       URAT   = UEDG(NBL(2),2)/QINF
       UEWAKE = UEDG(NBL(2),2) * (1.0-TKLAM) / (1.0 - TKLAM*URAT**2)
       SHWAKE = DSTR(NBL(2),2)/THET(NBL(2),2)
!
!----- extrapolate wake to downstream infinity using Squire-Young relation
!      (reduces errors of the wake not being long enough)
       CD = 2.0*THWAKE * (UEWAKE/QINF)**(0.5*(5.0+SHWAKE))
!
      ELSE
!
       CD = 0.0
!
      ENDIF
!
!---- calculate friction drag coefficient
      CDF = 0.0
      DO 20 IS=1, 2
        DO 205 IBL=3, IBLTE(IS)
          I  = IPAN(IBL  ,IS)
          IM = IPAN(IBL-1,IS)
          DX = (X(I) - X(IM))*CA + (Y(I) - Y(IM))*SA
          CDF = CDF + 0.5*(TAU(IBL,IS)+TAU(IBL-1,IS))*DX * 2.0/QINF**2
 205    CONTINUE
 20   CONTINUE
!
      RETURN
      END ! CDCALC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE LOAD
! ------------------------------------------------------
!      Reads airfoil file into buffer airfoil
!      and does various initial processesing on it.
! ------------------------------------------------------
      INCLUDE 'XFOIL.INC'
!      CHARACTER*(*) FILNAM
! 
!      FNAME = FILNAM
      IF(FNAME(1:1) .EQ. ' ') CALL ASKS('Enter filename^',FNAME)
! 
      LU = 9
      ITYPE = 1
      !CALL AREAD(LU,FNAME,IBX,XB,YB,NB,NAME,ISPARS,ITYPE,1)
      IF(ITYPE.EQ.0) RETURN
! 
!      IF(ITYPE.EQ.1) CALL ASKS('Enter airfoil name^',NAME)
      NAME = 'ThisIsATestName'
      CALL STRIP(NAME,NNAME)

! 
! ---- set default prefix for other filenames
      KDOT = INDEX(FNAME,'.')
      IF(KDOT.EQ.0) THEN
       PREFIX = FNAME
      ELSE
       PREFIX = FNAME(1:KDOT-1)
      ENDIF
      CALL STRIP(PREFIX,NPREFIX)
! 
! ---- calculate airfoil area assuming counterclockwise ordering
      AREA = 0.0
      DO 50 I=1, NB
        IP = I+1
        IF(I.EQ.NB) IP = 1
        AREA = AREA + 0.5*(YB(I)+YB(IP))*(XB(I)-XB(IP))
   50 CONTINUE
! 
      IF(AREA.GE.0.0) THEN
       LCLOCK = .FALSE.
       IF (DBUGMOD) WRITE(*,1010) NB
      ELSE
! ----- if area is negative (clockwise order), reverse coordinate order
       LCLOCK = .TRUE.
       IF (DBUGMOD) WRITE(*,1011) NB
       DO 55 I=1, NB/2
         XTMP = XB(NB-I+1)
         YTMP = YB(NB-I+1)
         XB(NB-I+1) = XB(I)
         YB(NB-I+1) = YB(I)
         XB(I) = XTMP
         YB(I) = YTMP
   55  CONTINUE
      ENDIF
! 
      IF(LNORM) THEN
       CALL NORM(XB,XBP,YB,YBP,SB,NB)
       IF (DBUGMOD) WRITE(*,1020)
      ENDIF
! 
      CALL SCALC(XB,YB,SB,NB)
      CALL SEGSPL(XB,XBP,SB,NB)
      CALL SEGSPL(YB,YBP,SB,NB)
! 
      CALL GEOPAR(XB,XBP,YB,YBP,SB,NB, W1,    &
                  SBLE,CHORDB,AREAB,RADBLE,ANGBTE,    &
                  EI11BA,EI22BA,APX1BA,APX2BA,    &
                  EI11BT,EI22BT,APX1BT,APX2BT,    &
                  THICKB,CAMBRB )
! 
      XBLE = SEVAL(SBLE,XB,XBP,SB,NB)
      YBLE = SEVAL(SBLE,YB,YBP,SB,NB)
      XBTE = 0.5*(XB(1) + XB(NB))
      YBTE = 0.5*(YB(1) + YB(NB))
! 
IF (DBUGMOD) THEN
      WRITE(*,1050) XBLE,YBLE, CHORDB,    &
                    XBTE,YBTE
ENDIF
! 
! ---- set reasonable MSES domain parameters for non-MSES coordinate file
      IF(ITYPE.LE.2 .AND. ISPARS.EQ.' ') THEN
        XBLE = SEVAL(SBLE,XB,XBP,SB,NB)
        YBLE = SEVAL(SBLE,YB,YBP,SB,NB)
        XINL = XBLE - 2.0*CHORDB
        XOUT = XBLE + 3.0*CHORDB
        YBOT = YBLE - 2.5*CHORDB
        YTOP = YBLE + 3.5*CHORDB
        XINL = AINT(20.0*ABS(XINL/CHORDB)+0.5)/20.0 * SIGN(CHORDB,XINL)
        XOUT = AINT(20.0*ABS(XOUT/CHORDB)+0.5)/20.0 * SIGN(CHORDB,XOUT)
        YBOT = AINT(20.0*ABS(YBOT/CHORDB)+0.5)/20.0 * SIGN(CHORDB,YBOT)
        YTOP = AINT(20.0*ABS(YTOP/CHORDB)+0.5)/20.0 * SIGN(CHORDB,YTOP)
        WRITE(ISPARS,1005) XINL, XOUT, YBOT, YTOP
 1005   FORMAT(1X, 4F8.2 )
      ENDIF
! 
! ---- wipe out old flap hinge location
      XBF = 0.0
      YBF = 0.0
      LBFLAP = .FALSE.
! 
! ---- wipe out off-design alphas, CLs
! c      NALOFF = 0
! c      NCLOFF = 0
! 
      RETURN
! ...............................................................
 1010 FORMAT(/' Number of input coordinate points:', I4    &
             /' Counterclockwise ordering')
 1011 FORMAT(/' Number of input coordinate points:', I4    &
             /' Clockwise ordering')
 1020 FORMAT(/' Airfoil has been normalized')
 1050 FORMAT(/'  LE  x,y  =', 2F10.5,'  |   Chord =',F10.5    &
             /'  TE  x,y  =', 2F10.5,'  |'                 )
             
      END ! LOAD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE MRCL(CLS,M_CLS,R_CLS)
! -------------------------------------------
!      Sets actual Mach, Reynolds numbers
!      from unit-CL values and specified CLS
!      depending on MATYP,RETYP flags.
! -------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL M_CLS
! 
      CLA = MAX( CLS , 0.000001 )
! 
      IF(RETYP.LT.1 .OR. RETYP.GT.3) THEN
        WRITE(*,*) 'MRCL:  Illegal Re(CL) dependence trigger.'
        WRITE(*,*) '       Setting fixed Re.'
        RETYP = 1
      ENDIF
      IF(MATYP.LT.1 .OR. MATYP.GT.3) THEN
        WRITE(*,*) 'MRCL:  Illegal Mach(CL) dependence trigger.'
        WRITE(*,*) '       Setting fixed Mach.'
        MATYP = 1
      ENDIF
! 
! 
      IF(MATYP.EQ.1) THEN
! 
        MINF  = MINF1
        M_CLS = 0.
! 
      ELSE IF(MATYP.EQ.2) THEN
! 
        MINF  =  MINF1/SQRT(CLA)
        M_CLS = -0.5*MINF/CLA
! 
      ELSE IF(MATYP.EQ.3) THEN
! 
        MINF  = MINF1
        M_CLS = 0.
! 
      ENDIF
! 
! 
      IF(RETYP.EQ.1) THEN
! 
        REINF = REINF1
        R_CLS = 0.
! 
      ELSE IF(RETYP.EQ.2) THEN
! 
        REINF =  REINF1/SQRT(CLA)
        R_CLS = -0.5*REINF/CLA
! 
      ELSE IF(RETYP.EQ.3) THEN
! 
        REINF =  REINF1/CLA
        R_CLS = -REINF /CLA
! 
      ENDIF
! 
! 
      IF(MINF .GE. 0.99) THEN
        WRITE(*,*)
        WRITE(*,*) 'MRCL: CL too low for chosen Mach(CL) dependence'
        WRITE(*,*) '      Aritificially limiting Mach to  0.99'
        MINF = 0.99
        M_CLS = 0.
      ENDIF
! 
      RRAT = 1.0
      IF(REINF1 .GT. 0.0) RRAT = REINF/REINF1
! 
      IF(RRAT .GT. 100.0) THEN
        WRITE(*,*)
        WRITE(*,*) 'MRCL: CL too low for chosen Re(CL) dependence'
        WRITE(*,*) '      Aritificially limiting Re to ',REINF1*100.0
        REINF = REINF1*100.0
        R_CLS = 0.
      ENDIF
! 
      RETURN
      END ! MRCL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE COMSET
      INCLUDE 'XFOIL.INC'
! 
! ---- set Karman-Tsien parameter TKLAM
      BETA = SQRT(1.0 - MINF**2)
      BETA_MSQ = -0.5/BETA
! 
      TKLAM   = MINF**2 / (1.0 + BETA)**2
      TKL_MSQ =     1.0 / (1.0 + BETA)**2    &
          - 2.0*TKLAM/ (1.0 + BETA) * BETA_MSQ
! 
! ---- set sonic Pressure coefficient and speed
      IF(MINF.EQ.0.0) THEN
       CPSTAR = -999.0
       QSTAR = 999.0
      ELSE
       CPSTAR = 2.0 / (GAMMA*MINF**2)    &
              * (( (1.0 + 0.5*GAMM1*MINF**2)    &
                  /(1.0 + 0.5*GAMM1        ))**(GAMMA/GAMM1) - 1.0)
       QSTAR = QINF/MINF    &
             * SQRT( (1.0 + 0.5*GAMM1*MINF**2)    &
                    /(1.0 + 0.5*GAMM1        ) )
      ENDIF
! 
      RETURN
      END ! COMSET
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		SUBROUTINE NACA(IDES1)
      INCLUDE 'XFOIL.INC'
!
!---- number of points per side
      NSIDE = IQX/3
!
      IF(IDES1 .LE. 0) THEN
       CALL ASKI('Enter NACA 4 or 5-digit airfoil designation^',IDES)
      ELSE
       IDES = IDES1
      ENDIF
!
      ITYPE = 0
      IF(IDES.LE.25099) ITYPE = 5
      IF(IDES.LE.9999 ) ITYPE = 4
!
      IF(ITYPE.EQ.0) THEN
       WRITE(*,*) 'This designation not implemented.'
       RETURN
      ENDIF
!
      IF(ITYPE.EQ.4) CALL NACA4(IDES,W1,W2,W3,NSIDE,XB,YB,NB,NAME)
      IF(ITYPE.EQ.5) CALL NACA5(IDES,W1,W2,W3,NSIDE,XB,YB,NB,NAME)
      CALL STRIP(NAME,NNAME)
!
!---- see if routines didn't recognize designator
      IF(IDES.EQ.0) RETURN
!
      LCLOCK = .FALSE.
!
      XBF = 0.0
      YBF = 0.0
      LBFLAP = .FALSE.
!
      CALL SCALC(XB,YB,SB,NB)
      CALL SEGSPL(XB,XBP,SB,NB)
      CALL SEGSPL(YB,YBP,SB,NB)
!
      CALL GEOPAR(XB,XBP,YB,YBP,SB,NB, W1,    &
                  SBLE,CHORDB,AREAB,RADBLE,ANGBTE,    &
                  EI11BA,EI22BA,APX1BA,APX2BA,    &
                  EI11BT,EI22BT,APX1BT,APX2BT,    &
                  THICKB,CAMBRB )
!
      WRITE(*,1200) NB
 1200 FORMAT(/' Buffer airfoil set using', I4,' points')
!
!---- set paneling
      CALL PANGEN(.TRUE.)
!cc      CALL PANPLT
!
      RETURN
      END ! NACA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE PANGEN(SHOPAR)
!---------------------------------------------------
!     Set paneling distribution from buffer airfoil
!     geometry, thus creating current airfoil.
! 
!     If REFINE=True, bunch points at x=XSREF on
!     top side and at x=XPREF on bottom side
!     by setting a fictitious local curvature of
!     CTRRAT*(LE curvature) there.
!---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      LOGICAL SHOPAR
!
      IF(NB.LT.2) THEN
       WRITE(*,*) 'PANGEN: Buffer airfoil not available.'
       N = 0
       RETURN
      ENDIF
!
!---- Number of temporary nodes for panel distribution calculation
!       exceeds the specified panel number by factor of IPFAC.
      IPFAC = 3
      IPFAC = 5
!
!---- number of airfoil panel points
      N = NPAN
!
!C---- number of wake points
!      NW = NPAN/8 + 2
!      IF(NW.GT.IWX) THEN
!       WRITE(*,*)
!     &  'Array size (IWX) too small.  Last wake point index reduced.'
!       NW = IWX
!      ENDIF
!
!---- set arc length spline parameter
      CALL SCALC(XB,YB,SB,NB)
!
!---- spline raw airfoil coordinates
      CALL SEGSPL(XB,XBP,SB,NB)
      CALL SEGSPL(YB,YBP,SB,NB)
!
!---- normalizing length (~ chord)
      SBREF = 0.5*(SB(NB)-SB(1))
!
!---- set up curvature array
      DO I = 1, NB
        W5(I) = ABS( CURV(SB(I),XB,XBP,YB,YBP,SB,NB) ) * SBREF
      ENDDO
!
!---- locate LE point arc length value and the normalized curvature there
      CALL LEFIND(SBLE,XB,XBP,YB,YBP,SB,NB)
      CVLE = ABS( CURV(SBLE,XB,XBP,YB,YBP,SB,NB) ) * SBREF
!
!---- check for doubled point (sharp corner) at LE
      IBLE = 0
      DO I = 1, NB-1
        IF(SBLE.EQ.SB(I) .AND. SBLE.EQ.SB(I+1)) THEN
         IBLE = I
         WRITE(*,*)
         WRITE(*,*) 'Sharp leading edge'
         GO TO 21
        ENDIF
      ENDDO
 21   CONTINUE
!
!---- set LE, TE points
      XBLE = SEVAL(SBLE,XB,XBP,SB,NB)
      YBLE = SEVAL(SBLE,YB,YBP,SB,NB)
      XBTE = 0.5*(XB(1)+XB(NB))
      YBTE = 0.5*(YB(1)+YB(NB))
      CHBSQ = (XBTE-XBLE)**2 + (YBTE-YBLE)**2
!
!---- set average curvature over 2*NK+1 points within Rcurv of LE point
      NK = 3
      CVSUM = 0.
      DO K = -NK, NK
        FRAC = FLOAT(K)/FLOAT(NK)
        SBK = SBLE + FRAC*SBREF/MAX(CVLE,20.0)
        CVK = ABS( CURV(SBK,XB,XBP,YB,YBP,SB,NB) ) * SBREF
        CVSUM = CVSUM + CVK
      ENDDO
      CVAVG = CVSUM/FLOAT(2*NK+1)
!
!---- dummy curvature for sharp LE
      IF(IBLE.NE.0) CVAVG = 10.0
!
!---- set curvature attraction coefficient actually used
      CC = 6.0 * CVPAR
!
!---- set artificial curvature at TE to bunch panels there
      CVTE = CVAVG * CTERAT
      W5(1)  = CVTE
      W5(NB) = CVTE
!
!
!**** smooth curvature array for smoother panel size distribution  ****
!
!CC      CALL ASKR('Enter curvature smoothing length/c^',SMOOL)
!CC      SMOOL = 0.010
!
!---- set smoothing length = 1 / averaged LE curvature, but 
!-    no more than 5% of chord and no less than 1/4 average panel spacing
      SMOOL = MAX( 1.0/MAX(CVAVG,20.0) , 0.25 /FLOAT(NPAN/2) )
!
      SMOOSQ = (SMOOL*SBREF) ** 2
!
!---- set up tri-diagonal system for smoothed curvatures
      W2(1) = 1.0
      W3(1) = 0.0
      DO I=2, NB-1
        DSM = SB(I) - SB(I-1)
        DSP = SB(I+1) - SB(I)
        DSO = 0.5*(SB(I+1) - SB(I-1))
!
        IF(DSM.EQ.0.0 .OR. DSP.EQ.0.0) THEN
!------- leave curvature at corner point unchanged
         W1(I) = 0.0
         W2(I) = 1.0
         W3(I) = 0.0
        ELSE
         W1(I) =  SMOOSQ * (         - 1.0/DSM) / DSO
         W2(I) =  SMOOSQ * ( 1.0/DSP + 1.0/DSM) / DSO  +  1.0
         W3(I) =  SMOOSQ * (-1.0/DSP          ) / DSO
        ENDIF
      ENDDO
!
      W1(NB) = 0.0
      W2(NB) = 1.0
!
!---- fix curvature at LE point by modifying equations adjacent to LE
      DO I=2, NB-1
        IF(SB(I).EQ.SBLE .OR. I.EQ.IBLE .OR. I.EQ.IBLE+1) THEN
!------- if node falls right on LE point, fix curvature there
         W1(I) = 0.
         W2(I) = 1.0
         W3(I) = 0.
         W5(I) = CVLE
        ELSE IF(SB(I-1).LT.SBLE .AND. SB(I).GT.SBLE) THEN
!------- modify equation at node just before LE point
         DSM = SB(I-1) - SB(I-2)
         DSP = SBLE    - SB(I-1)
         DSO = 0.5*(SBLE - SB(I-2))
!
         W1(I-1) =  SMOOSQ * (         - 1.0/DSM) / DSO
         W2(I-1) =  SMOOSQ * ( 1.0/DSP + 1.0/DSM) / DSO  +  1.0
         W3(I-1) =  0.
         W5(I-1) = W5(I-1) + SMOOSQ*CVLE/(DSP*DSO)
!
!------- modify equation at node just after LE point
         DSM = SB(I) - SBLE
         DSP = SB(I+1) - SB(I)
         DSO = 0.5*(SB(I+1) - SBLE)
         W1(I) =  0.
         W2(I) =  SMOOSQ * ( 1.0/DSP + 1.0/DSM) / DSO  +  1.0
         W3(I) =  SMOOSQ * (-1.0/DSP          ) / DSO
         W5(I) = W5(I) + SMOOSQ*CVLE/(DSM*DSO)
!
         GO TO 51
        ENDIF
      ENDDO
   51 CONTINUE
!
!---- set artificial curvature at bunching points and fix it there
      DO I=2, NB-1
!------ chord-based x/c coordinate
        XOC = (  (XB(I)-XBLE)*(XBTE-XBLE)    &
               + (YB(I)-YBLE)*(YBTE-YBLE) ) / CHBSQ
!
        IF(SB(I).LT.SBLE) THEN
!------- check if top side point is in refinement area
         IF(XOC.GT.XSREF1 .AND. XOC.LT.XSREF2) THEN
          W1(I) = 0.
          W2(I) = 1.0
          W3(I) = 0.
          W5(I) = CVLE*CTRRAT
         ENDIF
        ELSE
!------- check if bottom side point is in refinement area
         IF(XOC.GT.XPREF1 .AND. XOC.LT.XPREF2) THEN
          W1(I) = 0.
          W2(I) = 1.0
          W3(I) = 0.
          W5(I) = CVLE*CTRRAT
         ENDIF
        ENDIF
      ENDDO
!
!---- solve for smoothed curvature array W5
      IF(IBLE.EQ.0) THEN
       CALL TRISOL(W2,W1,W3,W5,NB)
      ELSE
       I = 1
       CALL TRISOL(W2(I),W1(I),W3(I),W5(I),IBLE)
       I = IBLE+1
       CALL TRISOL(W2(I),W1(I),W3(I),W5(I),NB-IBLE)
      ENDIF
!
!---- find max curvature
      CVMAX = 0.
      DO I=1, NB
        CVMAX = MAX( CVMAX , ABS(W5(I)) )
      ENDDO
!
!---- normalize curvature array
      DO I=1, NB
        W5(I) = W5(I) / CVMAX
      ENDDO
!
!---- spline curvature array
      CALL SEGSPL(W5,W6,SB,NB)
!
!---- Set initial guess for node positions uniform in s.
!     More nodes than specified (by factor of IPFAC) are 
!     temporarily used  for more reliable convergence.
      NN = IPFAC*(N-1)+1
!
!---- ratio of lengths of panel at TE to one away from the TE
      RDSTE = 0.667
      RTF = (RDSTE-1.0)*2.0 + 1.0
!
      IF(IBLE.EQ.0) THEN
!
       DSAVG = (SB(NB)-SB(1))/(FLOAT(NN-3) + 2.0*RTF)
       SNEW(1) = SB(1)
       DO I=2, NN-1
         SNEW(I) = SB(1) + DSAVG * (FLOAT(I-2) + RTF)
       ENDDO
       SNEW(NN) = SB(NB)
!
      ELSE
!
       NFRAC1 = (N * IBLE) / NB
!
       NN1 = IPFAC*(NFRAC1-1)+1
       DSAVG1 = (SBLE-SB(1))/(FLOAT(NN1-2) + RTF)
       SNEW(1) = SB(1)
       DO I=2, NN1
         SNEW(I) = SB(1) + DSAVG1 * (FLOAT(I-2) + RTF)
       ENDDO
!
       NN2 = NN - NN1 + 1
       DSAVG2 = (SB(NB)-SBLE)/(FLOAT(NN2-2) + RTF)
       DO I=2, NN2-1
         SNEW(I-1+NN1) = SBLE + DSAVG2 * (FLOAT(I-2) + RTF)
       ENDDO
       SNEW(NN) = SB(NB)
!
      ENDIF
!
!---- Newton iteration loop for new node positions
      DO 10 ITER=1, 20
!
!------ set up tri-diagonal system for node position deltas
        CV1  = SEVAL(SNEW(1),W5,W6,SB,NB)
        CV2  = SEVAL(SNEW(2),W5,W6,SB,NB)
        CVS1 = DEVAL(SNEW(1),W5,W6,SB,NB)
        CVS2 = DEVAL(SNEW(2),W5,W6,SB,NB)
!
        CAVM = SQRT(CV1**2 + CV2**2)
        IF(CAVM .EQ. 0.0) THEN
          CAVM_S1 = 0.
          CAVM_S2 = 0.
        ELSE
          CAVM_S1 = CVS1 * CV1/CAVM
          CAVM_S2 = CVS2 * CV2/CAVM
        ENDIF
!
        DO 110 I=2, NN-1
          DSM = SNEW(I) - SNEW(I-1)
          DSP = SNEW(I) - SNEW(I+1)
          CV3  = SEVAL(SNEW(I+1),W5,W6,SB,NB)
          CVS3 = DEVAL(SNEW(I+1),W5,W6,SB,NB)
!
          CAVP = SQRT(CV3**2 + CV2**2)
          IF(CAVP .EQ. 0.0) THEN
            CAVP_S2 = 0.
            CAVP_S3 = 0.
          ELSE
            CAVP_S2 = CVS2 * CV2/CAVP
            CAVP_S3 = CVS3 * CV3/CAVP
          ENDIF
!
          FM = CC*CAVM + 1.0
          FP = CC*CAVP + 1.0
!
          REZ = DSP*FP + DSM*FM
!
!-------- lower, main, and upper diagonals
          W1(I) =      -FM  +  CC*               DSM*CAVM_S1
          W2(I) =  FP + FM  +  CC*(DSP*CAVP_S2 + DSM*CAVM_S2)
          W3(I) = -FP       +  CC* DSP*CAVP_S3
!
!-------- residual, requiring that
!         (1 + C*curv)*deltaS is equal on both sides of node i
          W4(I) = -REZ
!
          CV1 = CV2
          CV2 = CV3
          CVS1 = CVS2
          CVS2 = CVS3
          CAVM    = CAVP
          CAVM_S1 = CAVP_S2
          CAVM_S2 = CAVP_S3
  110   CONTINUE
!
!------ fix endpoints (at TE)
        W2(1) = 1.0
        W3(1) = 0.0
        W4(1) = 0.0
        W1(NN) = 0.0
        W2(NN) = 1.0
        W4(NN) = 0.0
!
        IF(RTF .NE. 1.0) THEN
!------- fudge equations adjacent to TE to get TE panel length ratio RTF
!
         I = 2
         W4(I) = -((SNEW(I) - SNEW(I-1)) + RTF*(SNEW(I) - SNEW(I+1)))
         W1(I) = -1.0
         W2(I) =  1.0 + RTF
         W3(I) =      - RTF
!
         I = NN-1
         W4(I) = -((SNEW(I) - SNEW(I+1)) + RTF*(SNEW(I) - SNEW(I-1)))
         W3(I) = -1.0
         W2(I) =  1.0 + RTF
         W1(I) =      - RTF
        ENDIF
!
!
!------ fix sharp LE point
        IF(IBLE.NE.0) THEN
         I = NN1
         W1(I) = 0.0
         W2(I) = 1.0
         W3(I) = 0.0
         W4(I) = SBLE - SNEW(I)
        ENDIF
!
!------ solve for changes W4 in node position arc length values
        CALL TRISOL(W2,W1,W3,W4,NN)
!
!------ find under-relaxation factor to keep nodes from changing order
        RLX = 1.0
        DMAX = 0.0
        DO I=1, NN-1
          DS  = SNEW(I+1) - SNEW(I)
          DDS = W4(I+1) - W4(I)
          DSRAT = 1.0 + RLX*DDS/DS
          IF(DSRAT.GT.4.0) RLX = (4.0-1.0)*DS/DDS
          IF(DSRAT.LT.0.2) RLX = (0.2-1.0)*DS/DDS
          DMAX = MAX(ABS(W4(I)),DMAX)
        ENDDO
!
!------ update node position
        DO I=2, NN-1
          SNEW(I) = SNEW(I) + RLX*W4(I)
        ENDDO
!
!CC        IF(RLX.EQ.1.0) WRITE(*,*) DMAX
!CC        IF(RLX.NE.1.0) WRITE(*,*) DMAX,'    RLX =',RLX
        IF(ABS(DMAX).LT.1.E-3) GO TO 11
   10 CONTINUE
      WRITE(*,*) 'Paneling convergence failed.  Continuing anyway...'
!
   11 CONTINUE
!
!---- set new panel node coordinates
      DO I=1, N
        IND = IPFAC*(I-1) + 1
        S(I) = SNEW(IND)
        X(I) = SEVAL(SNEW(IND),XB,XBP,SB,NB)
        Y(I) = SEVAL(SNEW(IND),YB,YBP,SB,NB)
      ENDDO
!
!
!---- go over buffer airfoil again, checking for corners (double points)
      NCORN = 0
      DO 25 IB=1, NB-1
        IF(SB(IB) .EQ. SB(IB+1)) THEN
!------- found one !
!
         NCORN = NCORN+1
         XBCORN = XB(IB)
         YBCORN = YB(IB)
         SBCORN = SB(IB)
!
!------- find current-airfoil panel which contains corner
         DO 252 I=1, N
!
!--------- keep stepping until first node past corner
           IF(S(I) .LE. SBCORN) GO TO 252
!
!---------- move remainder of panel nodes to make room for additional node
            DO 2522 J=N, I, -1
              X(J+1) = X(J)
              Y(J+1) = Y(J)
              S(J+1) = S(J)
 2522       CONTINUE
            N = N+1
!
            IF(N .GT. IQX-1)    &
             STOP 'PANEL: Too many panels. Increase IQX in XFOIL.INC'
!
            X(I) = XBCORN
            Y(I) = YBCORN
            S(I) = SBCORN
!
!---------- shift nodes adjacent to corner to keep panel sizes comparable
            IF(I-2 .GE. 1) THEN
             S(I-1) = 0.5*(S(I) + S(I-2))
             X(I-1) = SEVAL(S(I-1),XB,XBP,SB,NB)
             Y(I-1) = SEVAL(S(I-1),YB,YBP,SB,NB)
            ENDIF
!
            IF(I+2 .LE. N) THEN
             S(I+1) = 0.5*(S(I) + S(I+2))
             X(I+1) = SEVAL(S(I+1),XB,XBP,SB,NB)
             Y(I+1) = SEVAL(S(I+1),YB,YBP,SB,NB)
            ENDIF
!
!---------- go on to next input geometry point to check for corner
            GO TO 25
!
  252    CONTINUE
        ENDIF
   25 CONTINUE
!
      CALL SCALC(X,Y,S,N)
      CALL SEGSPL(X,XP,S,N)
      CALL SEGSPL(Y,YP,S,N)
      CALL LEFIND(SLE,X,XP,Y,YP,S,N)
!
      XLE = SEVAL(SLE,X,XP,S,N)
      YLE = SEVAL(SLE,Y,YP,S,N)
      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))
      CHORD  = SQRT( (XTE-XLE)**2 + (YTE-YLE)**2 )
!
!---- calculate panel size ratios (user info)
      DSMIN =  1000.0
      DSMAX = -1000.0
      DO 40 I=1, N-1
        DS = S(I+1)-S(I)
        IF(DS .EQ. 0.0) GO TO 40
          DSMIN = MIN(DSMIN,DS)
          DSMAX = MAX(DSMAX,DS)
   40 CONTINUE
!
      DSMIN = DSMIN*FLOAT(N-1)/S(N)
      DSMAX = DSMAX*FLOAT(N-1)/S(N)
!cc      WRITE(*,*) 'DSmin/DSavg = ',DSMIN,'     DSmax/DSavg = ',DSMAX
!
!---- set various flags for new airfoil
      LGAMU = .FALSE.
      LQINU = .FALSE.
      LWAKE = .FALSE.
      LQAIJ = .FALSE.
      LADIJ = .FALSE.
      LWDIJ = .FALSE.
      LIPAN = .FALSE.
      LBLINI = .FALSE.
      LVCONV = .FALSE.
      LSCINI = .FALSE.
      LQSPEC = .FALSE.
      LGSAME = .FALSE.
!
      IF(LBFLAP) THEN
       XOF = XBF
       YOF = YBF
       LFLAP = .TRUE.
      ENDIF
!
!---- determine if TE is blunt or sharp, calculate TE geometry parameters
      CALL TECALC
!
!---- calculate normal vectors
      CALL NCALC(X,Y,S,N,NX,NY)
!
!---- calculate panel angles for panel routines
      CALL APCALC
!
      IF(SHARP) THEN
       WRITE(*,1090) 'Sharp trailing edge'
      ELSE
       GAP = SQRT((X(1)-X(N))**2 + (Y(1)-Y(N))**2)
       WRITE(*,1090) 'Blunt trailing edge.  Gap =', GAP
      ENDIF
 1090 FORMAT(/1X,A,F9.5)
!
      IF(SHOPAR) WRITE(*,1100) NPAN, CVPAR, CTERAT, CTRRAT,    &
                               XSREF1, XSREF2, XPREF1, XPREF2
 1100 FORMAT(/' Paneling parameters used...'    &
             /'   Number of panel nodes      ' , I4    &
             /'   Panel bunching parameter   ' , F6.3    &
             /'   TE/LE panel density ratio  ' , F6.3    &
             /'   Refined-area/LE panel density ratio   ' , F6.3    &
             /'   Top    side refined area x/c limits ' , 2F6.3    &
             /'   Bottom side refined area x/c limits ' , 2F6.3)
!
      RETURN
      END ! PANGEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE TECALC
!-------------------------------------------
!     Calculates total and projected TE gap 
!     areas and TE panel strengths.
!-------------------------------------------
      INCLUDE 'XFOIL.INC'
!
!---- set TE base vector and TE bisector components
      DXTE = X(1) - X(N)
      DYTE = Y(1) - Y(N)
      DXS = 0.5*(-XP(1) + XP(N))
      DYS = 0.5*(-YP(1) + YP(N))
!
!---- normal and streamwise projected TE gap areas
      ANTE = DXS*DYTE - DYS*DXTE
      ASTE = DXS*DXTE + DYS*DYTE
!
!---- total TE gap area
      DSTE = SQRT(DXTE**2 + DYTE**2)
!
      SHARP = DSTE .LT. 0.0001*CHORD
!
      IF(SHARP) THEN
       SCS = 1.0
       SDS = 0.0
      ELSE
       SCS = ANTE/DSTE
       SDS = ASTE/DSTE
      ENDIF
!
!---- TE panel source and vorticity strengths
      SIGTE = 0.5*(GAM(1) - GAM(N))*SCS
      GAMTE = -.5*(GAM(1) - GAM(N))*SDS
!
      SIGTE_A = 0.5*(GAM_A(1) - GAM_A(N))*SCS
      GAMTE_A = -.5*(GAM_A(1) - GAM_A(N))*SDS
!
      RETURN
      END ! TECALC
 

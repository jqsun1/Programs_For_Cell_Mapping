	SUBROUTINE OPER
      INCLUDE 'XFOIL.INC'
      PARAMETER (NPRX = 1001)
      DIMENSION XPR(NPRX), YPR(NPRX), FPR(NPRX)

      DIMENSION NBLP(NPX)
      DIMENSION IPPAI(NPX), NAPOLT(NPX)
!---- retain last-command info if OPER is exited and then re-entered
      SAVE COMOLD, ARGOLD
!
!---- logical units for  polar save file,  polar dump file
      LUPLR = 9
      LUPLX = 11
!

!
      IF(N.EQ.0) THEN
       WRITE(*,*)
       WRITE(*,*) '***  No airfoil available  ***'
       RETURN
      ENDIF
!C
      IF(IPACT.NE.0) THEN
       WRITE(*,5000) IPACT
 5000  FORMAT(/'  Polar', I3,'  is active')
      ENDIF
!
!cc 500  CONTINUE
      COMOLD = COMAND
      ARGOLD = COMARG
!
	LVISC = .TRUE.
! Set the Raynolds # here
	!REINF1 = 500000
!   LVCONV      .TRUE. if converged BL solution exists	
	LVCONV = .FALSE.
! Angle of attack in degrees
	LALFA = .TRUE.
    ALFA = DTOR*ADEG
    QINF = 1.0
    CALL SPECAL
    IF(ABS(ALFA-AWAKE) .GT. 1.0E-5) LWAKE  = .FALSE.
    IF(ABS(ALFA-AVISC) .GT. 1.0E-5) LVCONV = .FALSE.
    IF(ABS(MINF-MVISC) .GT. 1.0E-5) LVCONV = .FALSE. 
    IF(LVISC) CALL VISCAL(ITMAX)
!       CALL CPX used for plotting
       CALL FCPMIN
!
!cc    IF( LVISC .AND. LPACC .AND. LVCONV ) THEN
       IF( LPACC .AND. (LVCONV .OR. .NOT.LVISC)) THEN
        CALL PLRADD(LUPLR,IPACT)
        CALL PLXADD(LUPLX,IPACT)
       ENDIF
        COMOLD = COMAND
        ARGOLD = COMARG

    
	RETURN
	END !OPER
	
	
	
	
	
	
	
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE FCPMIN
!------------------------------------------------
!     Finds minimum Cp on dist for cavitation work
!------------------------------------------------
      INCLUDE 'XFOIL.INC'
!
      XCPMNI = X(1)
      XCPMNV = X(1)
      CPMNI = CPI(1)
      CPMNV = CPV(1)
!
      DO I = 2, N + NW
        IF(CPI(I) .LT. CPMNI) THEN
         XCPMNI = X(I)
         CPMNI = CPI(I)
        ENDIF
        IF(CPV(I) .LT. CPMNV) THEN
         XCPMNV = X(I)
         CPMNV = CPV(I)
        ENDIF
      ENDDO
!

      IF(LVISC)THEN
        CPMN = CPMNV
      ELSE
        CPMN = CPMNI
!
        CPMNV = CPMNI
        XCPMNV = XCPMNI
      ENDIF
!
      RETURN
      END ! FCPMIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SPECAL
!-----------------------------------
!     Converges to specified alpha.
!-----------------------------------
      INCLUDE 'XFOIL.INC'
      REAL MINF_CLM, MSQ_CLM
!
!---- calculate surface vorticity distributions for alpha = 0, 90 degrees
      IF(.NOT.LGAMU .OR. .NOT.LQAIJ) CALL GGCALC
!
      COSA = COS(ALFA)
      SINA = SIN(ALFA)
!
!---- superimpose suitably weighted  alpha = 0, 90  distributions
      DO 50 I=1, N
        GAM(I)   =  COSA*GAMU(I,1) + SINA*GAMU(I,2)
        GAM_A(I) = -SINA*GAMU(I,1) + COSA*GAMU(I,2)
   50 CONTINUE
      PSIO = COSA*GAMU(N+1,1) + SINA*GAMU(N+1,2)
!
      CALL TECALC
      CALL QISET
!
!---- set initial guess for the Newton variable CLM
      CLM = 1.0
!
!---- set corresponding  M(CLM), Re(CLM)
      CALL MRCL(CLM,MINF_CLM,REINF_CLM)
      CALL COMSET
!
!---- set corresponding CL(M)
      CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,    &
                  CL,CM,CDP, CL_ALF,CL_MSQ)
!
!---- iterate on CLM
      DO 100 ITCL=1, 20
!
        MSQ_CLM = 2.0*MINF*MINF_CLM
        DCLM = (CL - CLM)/(1.0 - CL_MSQ*MSQ_CLM)
!
        CLM1 = CLM
        RLX = 1.0
!
!------ under-relaxation loop to avoid driving M(CL) above 1
        DO 90 IRLX=1, 12
!
          CLM = CLM1 + RLX*DCLM
!
!-------- set new freestream Mach M(CLM)
          CALL MRCL(CLM,MINF_CLM,REINF_CLM)
!
!-------- if Mach is OK, go do next Newton iteration
          IF(MATYP.EQ.1 .OR. MINF.EQ.0.0 .OR. MINF_CLM.NE.0.0) GO TO 91
!
          RLX = 0.5*RLX
   90   CONTINUE
   91   CONTINUE
!
!------ set new CL(M)
        CALL COMSET
        CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,    &
                    CL,CM,CDP,CL_ALF,CL_MSQ)
!
        IF(ABS(DCLM).LE.1.0E-6) GO TO 110
!
  100 CONTINUE
      IF (DBUGMOD) WRITE(*,*) 'SPECAL:  Minf convergence failed'
      CONVERGED = .FALSE.
  110 CONTINUE
!
!---- set final Mach, CL, Cp distributions, and hinge moment
      CALL MRCL(CL,MINF_CL,REINF_CL)
      CALL COMSET
      CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,    &
                  CL,CM,CDP, CL_ALF,CL_MSQ)
      CALL CPCALC(N,QINV,QINF,MINF,CPI)
      IF(LVISC) THEN
       CALL CPCALC(N+NW,QVIS,QINF,MINF,CPV)
       CALL CPCALC(N+NW,QINV,QINF,MINF,CPI)
      ELSE
       CALL CPCALC(N,QINV,QINF,MINF,CPI)
      ENDIF
      IF(LFLAP) CALL MHINGE
!
      RETURN
      END ! SPECAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE VISCAL(NITER1)
!----------------------------------------
!     Converges viscous operating point
!----------------------------------------
      INCLUDE 'XFOIL.INC'
!
!---- convergence tolerance
      DATA EPS1 / 1.0E-4 /
!
      NITER = NITER1
!
!---- calculate wake trajectory from current inviscid solution if necessary
      IF(.NOT.LWAKE) THEN
       CALL XYWAKE
      ENDIF
!
!---- set velocities on wake from airfoil vorticity for alpha=0, 90
      CALL QWCALC
!
!---- set velocities on airfoil and wake for initial alpha
      CALL QISET
!
      IF(.NOT.LIPAN) THEN
!
       IF(LBLINI) CALL GAMQV
!
!----- locate stagnation point arc length position and panel index
       CALL STFIND
!
!----- set  BL position -> panel position  pointers
       CALL IBLPAN
!
!----- calculate surface arc length array for current stagnation point location
       CALL XICALC
!
!----- set  BL position -> system line  pointers
       CALL IBLSYS
!
      ENDIF
!
!---- set inviscid BL edge velocity UINV from QINV
      CALL UICALC
!
      IF(.NOT.LBLINI) THEN
!
!----- set initial Ue from inviscid Ue
       DO IBL=1, NBL(1)
         UEDG(IBL,1) = UINV(IBL,1)
       ENDDO
!
       DO IBL=1, NBL(2)
         UEDG(IBL,2) = UINV(IBL,2)
       ENDDO
!
      ENDIF
!
      IF(LVCONV) THEN
!----- set correct CL if converged point exists
       CALL QVFUE
       IF(LVISC) THEN
        CALL CPCALC(N+NW,QVIS,QINF,MINF,CPV)
        CALL CPCALC(N+NW,QINV,QINF,MINF,CPI)
       ELSE
        CALL CPCALC(N,QINV,QINF,MINF,CPI)
       ENDIF
       CALL GAMQV
       CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,    &
                   CL,CM,CDP, CL_ALF,CL_MSQ)
       CALL CDCALC
      ENDIF
!
!---- set up source influence matrix if it doesn't exist
      IF(.NOT.LWDIJ .OR. .NOT.LADIJ) CALL QDCALC
!
!---- Newton iteration for entire BL solution
      IF(NITER.EQ.0) CALL ASKI('Enter number of iterations^',NITER)
      IF (DBUGMOD) THEN
      WRITE(*,*)
      WRITE(*,*) 'Solving BL system ...'
      ENDIF
      DO 1000 ITER=1, NITER
!
!------ fill Newton system for BL variables
        CALL SETBL
!
!------ solve Newton system with custom solver
        CALL BLSOLV
!
!------ update BL variables
        CALL UPDATE
!
        IF(LALFA) THEN
!------- set new freestream Mach, Re from new CL
         CALL MRCL(CL,MINF_CL,REINF_CL)
         CALL COMSET
        ELSE
!------- set new inviscid speeds QINV and UINV for new alpha
         CALL QISET
         CALL UICALC
        ENDIF
!
!------ calculate edge velocities QVIS(.) from UEDG(..)
        CALL QVFUE
!
!------ set GAM distribution from QVIS
        CALL GAMQV
!
!------ relocate stagnation point
        CALL STMOVE
!
!------ set updated CL,CD
        CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,    &
                    CL,CM,CDP,CL_ALF,CL_MSQ)
        CALL CDCALC
!
!------ display changes and test for convergence
!        IF(RLX.LT.1.0)     &
!         WRITE(*,2000) ITER, RMSBL, RMXBL, VMXBL,IMXBL,ISMXBL,RLX
!        IF(RLX.EQ.1.0)     &
!         WRITE(*,2010) ITER, RMSBL, RMXBL, VMXBL,IMXBL,ISMXBL
!         CDPDIF = CD - CDF
!         WRITE(*,2020) ALFA/DTOR, CL, CM, CD, CDF, CDPDIF
!         CDSURF = CDP + CDF
!         WRITE(*,2025) CDSURF, CDF, CDP

        IF(RMSBL .LT. EPS1) THEN
         LVCONV = .TRUE.
         AVISC = ALFA
         MVISC = MINF
         GO TO 90
        ENDIF
!
 1000 CONTINUE
      if (DBUGMOD) WRITE(*,*) 'VISCAL:  Convergence failed'
      CONVERGED = .FALSE.

!
   90 CONTINUE
      CALL CPCALC(N+NW,QINV,QINF,MINF,CPI)
      CALL CPCALC(N+NW,QVIS,QINF,MINF,CPV)
      IF(LFLAP) CALL MHINGE


        is = 1
        hkmax = 0.
        hkm = 0.0
        psep = 0.
        patt = 0.
        do ibl = 2, iblte(is)
          hki = dstr(ibl,is) / thet(ibl,is)
          hkmax = max(hki,hkmax)
          if(hkm .lt. 4.0 .and.     &
             hki .ge. 4.0      ) then
           hfrac = (4.0 - hkm) / (hki - hkm )
           pdefm = uedg(ibl-1,is)**2 * thet(ibl-1,is)
           pdefi = uedg(ibl  ,is)**2 * thet(ibl  ,is)
           psep = pdefm*(1.0-hfrac) + pdefi*hfrac
          endif
          if(hkm .gt. 4.0 .and.     &
             hki .lt. 4.0      ) then
           hfrac = (4.0 - hkm) / (hki - hkm )
           pdefm = uedg(ibl-1,is)**2 * thet(ibl-1,is)
           pdefi = uedg(ibl  ,is)**2 * thet(ibl  ,is)
           patt = pdefm*(1.0-hfrac) + pdefi*hfrac
          endif
          hkm = hki
        enddo
        delp = patt - psep

IF (DBUGMOD) THEN
        write(*,9922)     &
          acrit(is), hkmax, cd, 2.0*psep, 2.0*patt, 2.0*delp,    &
          xoctr(is)
 9922   format(1x, f10.3, f10.4, f11.6, 3f11.6, f10.4, '     #')
ENDIF


        izero = ichar('0')

!c      fnum = acrit(is)
        fnum = xstrip(is)*100.0

        iten = int(  fnum                 / 9.99999 )
        ione = int( (fnum-float(10*iten)) / 0.99999 )
        idec = int( (fnum-float(10*iten)-float(ione)) / 0.09999 )

        fname = char(iten+izero)     &
             // char(ione+izero)     &
             // char(idec+izero) // '.bl'
        lu = 44
        open(lu,file=fname,status='unknown')
        rewind(lu)
        write(lu,'(a,a)')     &
      '#       s         ue          H          P         K ',    &
      '        x    -m du/dx'
!       1234567890 1234567890 1234567890 1234567890 1234567890 1234567890
        do ibl = 2, iblte(is)
          iblm = max( ibl-1 , 2 )
          iblp = min( ibl+1 , iblte(is) )
          i  = ipan(ibl ,is)
          hk = dstr(ibl,is) / thet(ibl,is)
          ddef = dstr(ibl,is)*uedg(ibl,is)
          pdef = thet(ibl,is)*uedg(ibl,is)**2
          edef = tstr(ibl,is)*uedg(ibl,is)**3 * 0.5
          duds = (uedg(iblp,is)-uedg(iblm,is))    &
               / (xssi(iblp,is)-xssi(iblm,is))
          dpds = -ddef*duds
          write(lu,9977)     &
             xssi(ibl,is), uedg(ibl,is), hk, pdef, edef, x(i), dpds
 9977     format(1x, 3f11.4, 2f11.6, f11.3, e14.6 )
        enddo
        close(lu)



      RETURN
!....................................................................
 2000   FORMAT    &
         (/1X,I3,'   rms: ',E10.4,'   max: ',E10.4,3X,A1,' at ',I4,I3,    &
           '   RLX:',F6.3)
 2010   FORMAT    &
         (/1X,I3,'   rms: ',E10.4,'   max: ',E10.4,3X,A1,' at ',I4,I3)
 2020   FORMAT    &
         ( 1X,3X,'   a =', F7.3,'      CL =',F8.4  /    &
           1X,3X,'  Cm =', F8.4, '     CD =',F9.5,    &
                 '   =>   CDf =',F9.5,'    CDp =',F9.5)
 2025   FORMAT    &
         ( 1X,3X, 6X     ,  8X , ' Int CD =',F9.5,    &
                 '   =>   CDf =',F9.5,'    CDp =',F9.5)
      END ! VISCAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE MHINGE
!  ----------------------------------------------------
!       Calculates the hinge moment of the flap about
!       (XOF,YOF) by integrating surface pressures.
!  ----------------------------------------------------
      INCLUDE 'XFOIL.INC'
!  
      IF(.NOT.LFLAP) THEN
!  
        CALL GETXYF(X,XP,Y,YP,S,N, TOPS,BOTS,XOF,YOF)
        LFLAP = .TRUE.
!  
      ELSE
!  
!  ------ find top and bottom y at hinge x location
        TOPS = XOF
        BOTS = S(N) - XOF
        CALL SINVRT(TOPS,XOF,X,XP,S,N)      
        CALL SINVRT(BOTS,XOF,X,XP,S,N)      
!  
      ENDIF
!  
      TOPX = SEVAL(TOPS,X,XP,S,N)
      TOPY = SEVAL(TOPS,Y,YP,S,N)
      BOTX = SEVAL(BOTS,X,XP,S,N)
      BOTY = SEVAL(BOTS,Y,YP,S,N)
!  
!  
      HMOM = 0.
      HFX  = 0.
      HFY  = 0.
!  
!  ---- integrate pressures on top and bottom sides of flap
      DO 20 I=2, N
        IF(S(I-1).GE.TOPS .AND. S(I).LE.BOTS) GO TO 20
!  
         DX = X(I) - X(I-1)
         DY = Y(I) - Y(I-1)
         XMID = 0.5*(X(I)+X(I-1)) - XOF
         YMID = 0.5*(Y(I)+Y(I-1)) - YOF
         IF(LVISC) THEN
          PMID = 0.5*(CPV(I) + CPV(I-1))
         ELSE
          PMID = 0.5*(CPI(I) + CPI(I-1))
         ENDIF
         HMOM = HMOM + PMID*(XMID*DX + YMID*DY)
         HFX  = HFX  - PMID* DY
         HFY  = HFY  + PMID* DX
   20 CONTINUE
!  
!  ---- find S(I)..S(I-1) interval containing s=TOPS
      DO I=2, N
        IF(S(I).GT.TOPS) GO TO 31
      ENDDO
!  
   31 CONTINUE
!  ---- add on top surface chunk TOPS..S(I-1),  missed in the DO 20 loop.
      DX = TOPX - X(I-1)
      DY = TOPY - Y(I-1)
      XMID = 0.5*(TOPX+X(I-1)) - XOF
      YMID = 0.5*(TOPY+Y(I-1)) - YOF
      IF(S(I) .NE. S(I-1)) THEN
       FRAC = (TOPS-S(I-1))/(S(I)-S(I-1))
      ELSE
       FRAC = 0.
      ENDIF
      IF(LVISC) THEN
       TOPP = CPV(I)*FRAC + CPV(I-1)*(1.0-FRAC)
       PMID = 0.5*(TOPP+CPV(I-1))
      ELSE
       TOPP = CPI(I)*FRAC + CPI(I-1)*(1.0-FRAC)
       PMID = 0.5*(TOPP+CPI(I-1))
      ENDIF
      HMOM = HMOM + PMID*(XMID*DX + YMID*DY)
      HFX  = HFX  - PMID* DY
      HFY  = HFY  + PMID* DX
!  
!  ---- add on inside flap surface contribution from hinge to top surface
      DX = XOF - TOPX
      DY = YOF - TOPY
      XMID = 0.5*(TOPX+XOF) - XOF
      YMID = 0.5*(TOPY+YOF) - YOF
      HMOM = HMOM + PMID*(XMID*DX + YMID*DY)
      HFX  = HFX  - PMID* DY
      HFY  = HFY  + PMID* DX
!  
!  ---- find S(I)..S(I-1) interval containing s=BOTS
      DO I=N, 2, -1
        IF(S(I-1).LT.BOTS) GO TO 41
      ENDDO
!  
   41 CONTINUE
!  ---- add on bottom surface chunk BOTS..S(I),  missed in the DO 20 loop.
      DX = X(I) - BOTX
      DY = Y(I) - BOTY
      XMID = 0.5*(BOTX+X(I)) - XOF
      YMID = 0.5*(BOTY+Y(I)) - YOF
      IF(S(I) .NE. S(I-1)) THEN
       FRAC = (BOTS-S(I-1))/(S(I)-S(I-1))
      ELSE
       FRAC = 0.
      ENDIF
      IF(LVISC) THEN
       BOTP = CPV(I)*FRAC + CPV(I-1)*(1.0-FRAC)
       PMID = 0.5*(BOTP+CPV(I))
      ELSE
       BOTP = CPI(I)*FRAC + CPI(I-1)*(1.0-FRAC)
       PMID = 0.5*(BOTP+CPI(I))
      ENDIF
      HMOM = HMOM + PMID*(XMID*DX + YMID*DY)
      HFX  = HFX  - PMID* DY
      HFY  = HFY  + PMID* DX
!  
!  ---- add on inside flap surface contribution from hinge to bottom surface
      DX = BOTX - XOF
      DY = BOTY - YOF
      XMID = 0.5*(BOTX+XOF) - XOF
      YMID = 0.5*(BOTY+YOF) - YOF
      HMOM = HMOM + PMID*(XMID*DX + YMID*DY)
      HFX  = HFX  - PMID* DY
      HFY  = HFY  + PMID* DX
!  
!  ---- add on TE base thickness contribution
      DX = X(1) - X(N)
      DY = Y(1) - Y(N)
      XMID = 0.5*(X(1)+X(N)) - XOF
      YMID = 0.5*(Y(1)+Y(N)) - YOF
      IF(LVISC) THEN
       PMID = 0.5*(CPV(1)+CPV(N))
      ELSE
       PMID = 0.5*(CPI(1)+CPI(N))
      ENDIF
      HMOM = HMOM + PMID*(XMID*DX + YMID*DY)
      HFX  = HFX  - PMID* DY
      HFY  = HFY  + PMID* DX
!  
      RETURN
      END ! MHINGE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

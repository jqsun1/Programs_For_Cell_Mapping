!***********************************************************************
!!    Module:  xgdes.f
!! 
!!    Copyright (C) 2000 Mark Drela 
!! 
!!    This program is free software; you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation; either version 2 of the License, or
!!    (at your option) any later version.
!!
!!    This program is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with this program; if not, write to the Free Software
!!    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!!***********************************************************************
!!
!      SUBROUTINE GDES
!      INCLUDE 'XFOIL.INC'
!      CHARACTER*4 COMAND, COMOLD
!      LOGICAL LRECALC, LMODPL, LPLNEW
!      DIMENSION XBOX(2), YBOX(2), XRF(2)
!      DIMENSION XBT(IBX), YBT(IBX)
!!
!      CHARACTER*128 COMARG, ARGOLD
!      CHARACTER*1 CHKEY
!!
!      DIMENSION IINPUT(20)
!      DIMENSION RINPUT(20)
!      LOGICAL ERROR
!!
!      EXTERNAL NEWPLOTG
!!
!      SAVE COMOLD, ARGOLD
!!
!      COMAND = '****'
!      COMARG = ' '
!      LRECALC = .FALSE.
!!
!      IF(NB.EQ.0) THEN
!       WRITE(*,*)
!       WRITE(*,*) '***  No airfoil available  ***'
!       RETURN
!      ENDIF
!!
!      LPLCAM = .FALSE.
!      LSYM = .TRUE.
!!
!      WRITE(*,*)
!      WRITE(*,*) 'You are working with the buffer airfoil'
!!
!      CALL PLTINI
!      CALL GOFINI
!      CALL PLOTG
!!
!!====================================================
!!---- start of menu loop
! 500  CONTINUE
!      COMOLD = COMAND
!      ARGOLD = COMARG
!!
! 501  IF(LGSYM) THEN
!       CALL ASKC('.GDESs^',COMAND,COMARG)
!      ELSE
!       CALL ASKC('.GDES^',COMAND,COMARG)
!      ENDIF
!!
!!--------------------------------------------------------
!!---- process previous command ?
!      IF(COMAND(1:1).EQ.'!') THEN
!        IF(COMOLD.EQ.'****') THEN
!          WRITE(*,*) 'Previous .GDES command not valid'
!          GO TO 501
!        ELSE
!          COMAND = COMOLD
!          COMARG  = ARGOLD
!          LRECALC = .TRUE.
!        ENDIF
!      ELSE
!        LRECALC = .FALSE.
!      ENDIF
!!
!      IF(COMAND.EQ.'    ') THEN
!!----- just <return> was typed... clean up plotting and exit OPER
!       IF(LPLOT) CALL PLEND
!       LPLOT = .FALSE.
!       LGSYM = .FALSE.
!       LGEOPL = .FALSE.
!       IF(.NOT.LGSAME) THEN
!        WRITE(*,*)
!        WRITE(*,*) 'Buffer airfoil is not identical to current airfoil'
!       ENDIF
!       CALL CLRZOOM
!       RETURN
!      ENDIF
!!
!!---- extract command line numeric arguments
!      DO I=1, 20
!        IINPUT(I) = 0
!        RINPUT(I) = 0.0
!      ENDDO
!      NINPUT = 20
!      CALL GETINT(COMARG,IINPUT,NINPUT,ERROR)
!      NINPUT = 20
!      CALL GETFLT(COMARG,RINPUT,NINPUT,ERROR)
!!
!!--------------------------------------------------------
!      IF(COMAND.EQ.'?   ') THEN
!       WRITE(*,1050)
! 1050  FORMAT(
!     & /'   <cr>     Return to Top Level'
!     & /'   !        Redo previous command'
!     &//'   GSET     Set buffer  airfoil <== current airfoil'
!     & /'   eXec     Set current airfoil <== buffer  airfoil'
!     & /'   SYMM     Toggle y-symmetry flag'
!     &//'   ADEG r   Rotate about origin (degrees)'
!     & /'   ARAD r   Rotate about origin (radians)'
!     & /'   Tran rr  Translate'
!     & /'   Scal r   Scale about origin'
!     & /'   LINS rr. Linearly-varying y scale'
!     & /'   DERO     Derotate (set chord line level)'
!     &//'   TGAP rr  Change trailing edge gap'
!     & /'   LERA rr  Change leading edge radius'
!     &//'   TCPL     Toggle thickness and camber plotting'
!     & /'   TFAC rr  Scale existing thickness and camber'
!     & /'   TSET rr  Set new thickness and camber'
!     & /'   HIGH rr  Move camber and thickness highpoints'
!     & /'  .CAMB     Modify camber shape directly or via loading'
!     &//'   BEND     Display structural properties of buffer airfoil'
!     &//'   Flap rrr Deflect trailing edge flap'
!     &//'   Modi     Modify contour via cursor'
!     & /'   SLOP     Toggle modified-contour slope matching flag'
!     &//'   CORN     Double point with cursor (set sharp corner)'
!     & /'   ADDP     Add    point with cursor or keyboard x,y'
!     & /'   MOVP     Move   point with cursor or keyboard x,y'
!     & /'   DELP     Delete point with cursor'
!     & /'   NMOV r   Move all points in surface-normal direction'
!     &//'   UNIT     Normalize buffer airfoil to unit chord'
!     & /'   Dist     Determine distance between 2 cursor points'
!     & /'   CLIS     List curvatures'
!     & /'   CPLO     Plot curvatures'
!     & /'   CANG     List panel corner angles'
!     & /'   CADD ri. Add points at corners exceeding angle threshold'
!     &//'   Plot     Replot buffer airfoil'
!     & /'   INPL     Replot buffer airfoil without scaling (in inches)'
!     & /'   Blow     Blowup plot region'
!     & /'   Rese     Reset plot scale and origin'
!     & /'   Wind     Plot window adjust via cursor and keys'
!     &//'   TSIZ r   Change tick-mark size'
!     & /'   TICK     Toggle node tick-mark plotting'
!     & /'   GRID     Toggle grid plotting'
!     & /'   GPAR     Toggle geometric parameter plotting'
!     & /'   Over f   Overlay disk file airfoil'
!     &//'   SIZE r   Change absolute plot-object size'
!     & /'  .ANNO     Annotate plot'
!     & /'   HARD     Hardcopy current plot'
!     &//'   NAME s   Specify new airfoil name'
!     & /'   NINC     Increment name version number')
!
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'Z   ') THEN
!       CALL USETZOOM(.TRUE.,.TRUE.)
!       CALL REPLOT(IDEV)
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'U   ') THEN
!       CALL CLRZOOM
!       CALL REPLOT(IDEV)
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'GSET') THEN
!       NB = N
!       DO I=1, NB
!         XB(I) = X(I)
!         YB(I) = Y(I)
!       ENDDO
!       LGSAME = .TRUE.
!       CALL SCALC(XB,YB,SB,NB)
!       CALL SEGSPL(XB,XBP,SB,NB)
!       CALL SEGSPL(YB,YBP,SB,NB)
!!
!       CALL GEOPAR(XB,XBP,YB,YBP,SB,NB,W1,
!     &             SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
!     &             EI11BA,EI22BA,APX1BA,APX2BA,
!     &             EI11BT,EI22BT,APX1BT,APX2BT,
!     &             THICKB,CAMBRB )
!!
!       CALL PLTINI
!       CALL PLOTG
!       IF(LGSYM) CALL ZERCAM
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'EXEC' .OR.
!     &       COMAND.EQ.'X   '      ) THEN
!       CALL ABCOPY(.TRUE.)
!!c       CALL NAMMOD(NAME,1,1)
!!c       CALL STRIP(NAME,NNAME)
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'SYMM') THEN
!       LGSYM = .NOT.LGSYM
!       IF(LGSYM) THEN
!        WRITE(*,*) 'y-symmetry forcing enabled.'
!        CALL ZERCAM
!       ELSE
!        WRITE(*,*) 'y-symmetry forcing disabled.'
!       ENDIF
!!
!!=================================================
!!---- rotate airfoil by degrees
!      ELSEIF(COMAND.EQ.'ADEG' .OR.
!     &       COMAND.EQ.'ARAD'     ) THEN
!       IF(COMAND.EQ.'ADEG') THEN
!         IF(NINPUT.GE.1) THEN
!          ADEG = RINPUT(1)
!         ELSE
!          ADEG = 0.0
!          CALL ASKR('Enter angle change (deg)^',ADEG)
!         ENDIF
!         ARAD = ADEG*PI/180.0
!       ELSE
!         IF(NINPUT.GE.1) THEN
!          ARAD = RINPUT(1)
!         ELSE
!          ARAD = 0.0
!          CALL ASKR('Enter angle change (rad)^',ARAD)
!         ENDIF
!       ENDIF
!!
!       CALL ROTATE(XB,YB,NB,ARAD)
!!CC      CALL SCALC(XB,YB,SB,NB)
!       CALL SEGSPL(XB,XBP,SB,NB)
!       CALL SEGSPL(YB,YBP,SB,NB)
!!
!       APX1BA = APX1BA - ARAD
!       APX2BA = APX2BA - ARAD
!       APX1BT = APX1BT - ARAD
!       APX2BT = APX2BT - ARAD
!!
!       CALL NEWPEN(2)
!       CALL PLTAIR(XB,XBP,YB,YBP,SB,NB, XOFF,XSF,YOFF,YSF,'magenta')
!       CALL PLNEWP('magenta')
!       LGEOPL = .FALSE.
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'TRAN' .OR.
!     &       COMAND.EQ.'T   '      ) THEN
!       IF    (NINPUT.GE.2) THEN
!        DELX = RINPUT(1)
!        DELY = RINPUT(2)
!       ELSEIF(NINPUT.GE.1) THEN
!        DELX = RINPUT(1)
!        DELY = 0.0
!        CALL ASKR('Enter delta(y)^',DELY)
!       ELSE
!        DELX = 0.0
!        CALL ASKR('Enter delta(x)^',DELX)
!        DELY = 0.0
!        CALL ASKR('Enter delta(y)^',DELY)
!       ENDIF
!       DO I=1, NB
!         XB(I) = XB(I) + DELX
!         YB(I) = YB(I) + DELY
!       ENDDO
!!
!       CALL NEWPEN(2)
!       CALL PLTAIR(XB,XBP,YB,YBP,SB,NB, XOFF,XSF,YOFF,YSF,'magenta')
!       CALL PLNEWP('magenta')
!       LGEOPL = .FALSE.
!       LGSAME = .FALSE.
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'SCAL' .OR.
!     &       COMAND.EQ.'S   '      ) THEN
!       IF(NINPUT.GE.1) THEN
!        FAC = RINPUT(1)
!        XXFAC = FAC
!        YYFAC = FAC
!       ELSE
!        FAC = 1.0
!        CALL ASKR('Enter scale factor (0 for separate x,y scales)^',FAC)
!        XXFAC = FAC
!        YYFAC = FAC
!       ENDIF
!!
!       IF(FAC .EQ. 0.0) THEN
!        IF(NINPUT.GE.3) THEN
!         XXFAC = RINPUT(2)
!         YYFAC = RINPUT(3)
!        ELSE
!         XXFAC = 1.0
!         CALL ASKR('Enter x scale factor^',XXFAC)
!         YYFAC = 1.0
!         CALL ASKR('Enter y scale factor^',YYFAC)
!        ENDIF
!       ENDIF
!!
!       DO I=1, NB
!         XB(I) = XB(I)*XXFAC
!         YB(I) = YB(I)*YYFAC
!       ENDDO
!!
!!----- re-order if necessary to maintain counterclockwise ordering 
!       IF(XXFAC*YYFAC .LT. 0.0) THEN
!         DO I=1, NB/2
!           XTMP = XB(I)
!           YTMP = YB(I)
!           XB(I) = XB(NB-I+1)
!           YB(I) = YB(NB-I+1)
!           XB(NB-I+1) = XTMP
!           YB(NB-I+1) = YTMP
!         ENDDO
!       ENDIF
!!
!!----- re-spline new geometry
!       CALL SCALC(XB,YB,SB,NB)
!       CALL SEGSPL(XB,XBP,SB,NB)
!       CALL SEGSPL(YB,YBP,SB,NB)
!!
!       CALL GEOPAR(XB,XBP,YB,YBP,SB,NB,W1,
!     &             SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
!     &             EI11BA,EI22BA,APX1BA,APX2BA,
!     &             EI11BT,EI22BT,APX1BT,APX2BT,
!     &             THICKB,CAMBRB )
!!
!       CALL NEWPEN(2)
!       CALL PLTAIR(XB,XBP,YB,YBP,SB,NB, XOFF,XSF,YOFF,YSF,'magenta')
!       CALL PLNEWP('magenta')
!       LGEOPL = .FALSE.
!       LGSAME = .FALSE.
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'LINS') THEN
! 40    CONTINUE
!       IF(NINPUT.GE.4) THEN
!         XOC1  = RINPUT(1)
!         YFAC1 = RINPUT(2)
!         XOC2  = RINPUT(3)
!         YFAC2 = RINPUT(4)
!       ELSE
! 1001    FORMAT(/1X,A,$)
! 41      WRITE(*,1001) 'Location 1...  enter  x/c, y-scale :  '
!         READ(*,*,ERR=41) XOC1, YFAC1
! 42      WRITE(*,1001) 'Location 2...  enter  x/c, y-scale :  '
!         READ(*,*,ERR=42) XOC2, YFAC2
!       ENDIF
!!
!       IF(ABS(XOC1-XOC2) .LT. 1.0E-5) THEN
!        WRITE(*,*) 'x/c locations 1 and 2 must be different'
!        NINPUT = 0
!        GO TO 40
!       ENDIF
!!
!       CALL LEFIND(SBLE,XB,XBP,YB,YBP,SB,NB)
!       XLE = SEVAL(SBLE,XB,XBP,SB,NB)
!       YLE = SEVAL(SBLE,YB,YBP,SB,NB)
!       XTE = 0.5*(XB(1) + XB(NB))
!       YTE = 0.5*(YB(1) + YB(NB))
!       DO I=1, NB
!         XOC = (XB(I)-XLE) / (XTE-XLE)
!         FR1 = (XOC2-XOC )/(XOC2-XOC1)
!         FR2 = (XOC -XOC1)/(XOC2-XOC1)
!         YYFAC = FR1*YFAC1 + FR2*YFAC2
!         YB(I) = YB(I)*YYFAC
!       ENDDO
!!
!!----- re-spline new geometry
!       CALL SCALC(XB,YB,SB,NB)
!       CALL SEGSPL(XB,XBP,SB,NB)
!       CALL SEGSPL(YB,YBP,SB,NB)
!!
!       CALL GEOPAR(XB,XBP,YB,YBP,SB,NB,W1,
!     &             SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
!     &             EI11BA,EI22BA,APX1BA,APX2BA,
!     &             EI11BT,EI22BT,APX1BT,APX2BT,
!     &             THICKB,CAMBRB )
!!
!       CALL NEWPEN(2)
!       CALL PLTAIR(XB,XBP,YB,YBP,SB,NB, XOFF,XSF,YOFF,YSF,'magenta')
!       CALL PLNEWP('magenta')
!       LGEOPL = .FALSE.
!       LGSAME = .FALSE.
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'DERO') THEN
!       CALL LEFIND(SBLE,XB,XBP,YB,YBP,SB,NB)
!       XLE = SEVAL(SBLE,XB,XBP,SB,NB)
!       YLE = SEVAL(SBLE,YB,YBP,SB,NB)
!       XTE = 0.5*(XB(1) + XB(NB))
!       YTE = 0.5*(YB(1) + YB(NB))
!!
!       ARAD = ATAN2(YTE-YLE,XTE-XLE)
!       CALL ROTATE(XB,YB,NB,ARAD)
!       WRITE(*,1080) ARAD / DTOR
! 1080  FORMAT(/'Rotating buffer airfoil by ',F8.3,' deg.')
!!
!       CALL SCALC(XB,YB,SB,NB)
!       CALL SEGSPL(XB,XBP,SB,NB)
!       CALL SEGSPL(YB,YBP,SB,NB)
!       CALL GEOPAR(XB,XBP,YB,YBP,SB,NB,W1,
!     &             SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
!     &             EI11BA,EI22BA,APX1BA,APX2BA,
!     &             EI11BT,EI22BT,APX1BT,APX2BT,
!     &             THICKB,CAMBRB )
!!
!       CALL NEWPEN(2)
!       CALL PLTAIR(XB,XBP,YB,YBP,SB,NB, XOFF,XSF,YOFF,YSF,'magenta')
!       CALL PLNEWP('magenta')
!       LGEOPL = .FALSE.
!       LGSAME = .FALSE.
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'TGAP') THEN
!       CALL TGAP(RINPUT,NINPUT)
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'LERA') THEN
!       CALL LERAD(RINPUT,NINPUT)
!!
!!--------------------------------------------------------
!!c      ELSEIF(COMAND.EQ.'TC  ') THEN
!!c       CALL TCBUF
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'TCPL') THEN
!       LPLCAM = .NOT.LPLCAM
!       CALL PLTINI
!       CALL GOFINI
!       CALL PLOTG
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'TFAC') THEN
!       IF(.NOT.LPLCAM) THEN
!        WRITE(*,*) 'Enabling camber,thickness plotting'
!        LPLCAM = .TRUE.
!        CALL PLTINI
!        CALL GOFINI
!        CALL PLOTG
!       ENDIF
!       CALL TCSCAL(RINPUT,NINPUT)
!       CALL NEWPEN(2)
!       CALL PLTAIR(XB,XBP,YB,YBP,SB,NB, XOFF,XSF,YOFF,YSF,'magenta')
!       CALL PLNEWP('magenta')
!       CALL PLTCAM('magenta')
!       LGEOPL = .FALSE.
!       LGSAME = .FALSE.
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'TSET') THEN
!       IF(.NOT.LPLCAM) THEN
!        WRITE(*,*) 'Enabling camber,thickness plotting'
!        LPLCAM = .TRUE.
!        CALL PLTINI
!        CALL GOFINI
!        CALL PLOTG
!       ENDIF
!       CALL TCSET(RINPUT,NINPUT)
!       CALL NEWPEN(2)
!       CALL PLTAIR(XB,XBP,YB,YBP,SB,NB, XOFF,XSF,YOFF,YSF,'magenta')
!       CALL PLNEWP('magenta')
!       CALL PLTCAM('magenta')
!       LGEOPL = .FALSE.
!       LGSAME = .FALSE.
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'HIGH') THEN
!       IF(.NOT.LPLCAM) THEN
!        WRITE(*,*) 'Enabling camber,thickness plotting'
!        LPLCAM = .TRUE.
!        CALL PLTINI
!        CALL GOFINI
!        CALL PLOTG
!       ENDIF
!       CALL HIPNT(RINPUT,NINPUT)
!       CALL NEWPEN(2)
!       CALL PLTAIR(XB,XBP,YB,YBP,SB,NB, XOFF,XSF,YOFF,YSF,'magenta')
!       CALL PLNEWP('magenta')
!       CALL PLTCAM('magenta')
!       LGEOPL = .FALSE.
!       LGSAME = .FALSE.
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'CAMB') THEN
!       IF(LGSYM) THEN
!        WRITE(*,*) 'Disabling symmetry enforcement.'
!        LGSYM = .FALSE.
!       ENDIF
!       CALL CAMB
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'BEND') THEN
!       CALL BENDUMP(NB,XB,YB)
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'CANG') THEN
!       CALL CANG(XB,YB,NB,2, IMAX,AMAX)
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'CADD') THEN
!       CALL CANG(XB,YB,NB,2, IMAX,AMAX)
!       WRITE(*,*)
!!
!       XBMIN = XB(1)
!       XBMAX = XB(1)
!       DO I=1, NB
!         XBMIN = MIN(XBMIN,XB(I))
!         XBMAX = MAX(XBMAX,XB(I))
!       ENDDO
!!
!!----- default inputs
!       ATOL = 0.5*AMAX
!       ISPL = 1
!       XRF(1) = XBMIN - 0.1*(XBMAX-XBMIN)
!       XRF(2) = XBMAX + 0.1*(XBMAX-XBMIN)
!!
!       IF    (NINPUT.LE.0) THEN
!        GO TO 70
!       ELSEIF(NINPUT.LE.1) THEN
!        ATOL = RINPUT(1)
!        GO TO 71
!       ELSEIF(NINPUT.LE.2) THEN
!        ATOL = RINPUT(1)
!        ISPL = IINPUT(2)
!        GO TO 72
!       ELSEIF(NINPUT.LE.4) THEN
!        ATOL = RINPUT(1)
!        ISPL = IINPUT(2)
!        XRF(1) = RINPUT(3)
!        XRF(2) = RINPUT(4)
!        GO TO 74
!       ENDIF
!!
! 70    WRITE(*,1090) ATOL
! 1090  FORMAT(1X,
!     &   'Enter corner angle criterion for refinement (deg):', F8.3)
!       CALL READR(1,ATOL,ERROR)
!       IF(ERROR) GO TO 70
!!
! 71    WRITE(*,1091) ISPL
! 1091  FORMAT(1X,
!     &   'Enter type of spline parameter (1=uniform, 2=arclength):', I4)
!       CALL READI(1,ISPL,ERROR)
!       IF(ERROR) GO TO 71
!       IF(ISPL.LE.0) GO TO 500
!       IF(ISPL.GT.2) GO TO 71
!!
! 72    WRITE(*,1092) XRF(1), XRF(2)
! 1092  FORMAT(1X,
!     &   'Enter refinement x limits:', 2F10.5)
!       CALL READR(2,XRF,ERROR)
!       IF(ERROR) GO TO 72
!!
! 74    CONTINUE
!       IF(ISPL.EQ.1) THEN
!        SB(1) = 0.0
!        DO I = 2, NB
!          IF(XB(I).EQ.XB(I-1) .AND. YB(I).EQ.YB(I-1)) THEN
!           SB(I) = SB(I-1)
!          ELSE
!           SB(I) = SB(I-1) + 1.0
!          ENDIF
!        ENDDO
!        CALL SEGSPL(XB,XBP,SB,NB)
!        CALL SEGSPL(YB,YBP,SB,NB)
!       ENDIF
!!
!       CALL AREFINE(XB,YB,SB,XBP,YBP,NB, ATOL, 
!     &             IBX,NNEW,W1,W2,XRF(1),XRF(2))
!!
!       NBADD = NNEW - NB
!       WRITE(*,*) 'Number of points added: ', NBADD
!!
!       NB = NNEW
!       DO I = 1, NB
!         XB(I) = W1(I)
!         YB(I) = W2(I)
!       ENDDO
!       LGSAME = .FALSE.
!!
!       CALL SCALC(XB,YB,SB,NB)
!       CALL SEGSPL(XB,XBP,SB,NB)
!       CALL SEGSPL(YB,YBP,SB,NB)
!!
!       CALL GEOPAR(XB,XBP,YB,YBP,SB,NB,W1,
!     &             SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
!     &             EI11BA,EI22BA,APX1BA,APX2BA,
!     &             EI11BT,EI22BT,APX1BT,APX2BT,
!     &             THICKB,CAMBRB )
!!
!       CALL NEWPEN(2)
!       CALL PLTAIR(XB,XBP,YB,YBP,SB,NB, XOFF,XSF,YOFF,YSF,'magenta')
!       CALL PLNEWP('magenta')
!!
!       CALL CANG(XB,YB,NB,1, IMAX,AMAX)
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'CLIS') THEN
!       CALL CLIS(XB,XBP,YB,YBP,SB,NB)
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'CPLO') THEN
!       CALL PLTCRV(SBLE,XB,XBP,YB,YBP,SB,NB,W1)
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'FLAP' .OR.
!     &       COMAND.EQ.'F   '      ) THEN
!       CALL FLAP(RINPUT,NINPUT)
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'MODI' .OR.
!     &       COMAND.EQ.'M   '      ) THEN
!!----- plot current geometry if it's not on the screen
!       IF(.NOT.LGEOPL) THEN
!        CALL PLTINI
!        CALL PLOTG
!       ENDIF
!!
!       IF(LGSYM) THEN
!         DO I = 1, NB
!           XBT(I) = XB(I)
!           YBT(I) = YB(I)
!         ENDDO
!       ENDIF
!!
!       IBFRST = 1
!       IBLAST = NB
!       NSIDE = 1
!       XBOX(1) = XMARG
!       XBOX(2) = XPAGE-XMARG
!       YBOX(1) = YMARG
!       YBOX(2) = YPAGE-YMARG
!       LMODPL = .FALSE.
!       CALL MODIXY(IBX,IBFRST,IBLAST,NSIDE,
!     &             XB,YB,XBP,YBP,SB, LGSLOP,
!     &             IGMOD1,IGMOD2,ISMOD,
!     &             XBOX,YBOX, XBOX,YBOX,SIZE,
!     &             XOFF,YOFF,XSF,YSF, LMODPL,
!     &             NEWPLOTG)
!!
!       IF(LGSYM) THEN
!         DO I = 1, NB
!           XBDEL = XB(I) - XBT(I)
!           YBDEL = YB(I) - YBT(I)
!           XB(I) = XB(I) + XBDEL
!           YB(I) = YB(I) + YBDEL
!         ENDDO
!         CALL ZERCAM
!       ENDIF
!       LGSAME = .FALSE.
!!
!       CALL SCALC(XB,YB,SB,NB)
!       CALL SEGSPL(XB,XBP,SB,NB)
!       CALL SEGSPL(YB,YBP,SB,NB)
!!
!       CALL GEOPAR(XB,XBP,YB,YBP,SB,NB,W1,
!     &             SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
!     &             EI11BA,EI22BA,APX1BA,APX2BA,
!     &             EI11BT,EI22BT,APX1BT,APX2BT,
!     &             THICKB,CAMBRB )
!!
!      CALL NEWPEN(2)
!      CALL PLTAIR(XB,XBP,YB,YBP,SB,NB, XOFF,XSF,YOFF,YSF,'magenta')
!      CALL PLNEWP('magenta')
!      LGEOPL = .FALSE.
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'SLOP') THEN
!       LGSLOP = .NOT.LGSLOP
!       IF(LGSLOP) THEN
!        WRITE(*,*) 'Modified segment will be',
!     &             ' made tangent at endpoints'
!       ELSE
!        WRITE(*,*) 'Modified segment will not be',
!     &             ' made tangent at endpoints'
!       ENDIF
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'CORN') THEN
!       IF(NB.EQ.2*IQX) THEN
!        WRITE(*,*)
!     &    'Buffer airfoil arrays will overflow.  No action taken.'
!         GO TO 500
!       ENDIF
!!
!       XWS = XWIND/SIZE
!       YWS = YWIND/SIZE
!       CALL POINTF(XB,XBP,YB,YBP,SB,NB, XWS,YWS, XOFF,YOFF,XSF,YSF,
!     &             IPNT,XC,YC)
!       IF(IPNT.EQ.0) GO TO 500
!       IF(IPNT.EQ.1 .OR. IPNT.EQ.NB) THEN
!        WRITE(*,*) 'Cannot double trailing edge point. No action taken.'
!        GO TO 500
!       ENDIF
!!
!!----- add doubled point
!       DO I=NB, IPNT, -1
!         XB(I+1) = XB(I)
!         YB(I+1) = YB(I)
!       ENDDO
!       NB = NB+1
!       LGSAME = .FALSE.
!!
!!----- spline new geometry
!       CALL SCALC(XB,YB,SB,NB)
!       CALL SEGSPL(XB,XBP,SB,NB)
!       CALL SEGSPL(YB,YBP,SB,NB)
!       CALL PLTINI
!       CALL PLOTG
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'ADDP') THEN
!       CALL ADDP
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'DELP') THEN
!       CALL DELP
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'NMOV') THEN
!       IF(NINPUT.GE.1) THEN
!        DELN = RINPUT(1)
!       ELSE
!        DELN = 0.0
!        CALL ASKR('Enter normal movement (+ outward)^',DELN)
!       ENDIF
!
!       DO IB = 1, NB
!         ENX =  YBP(IB)
!         ENY = -XBP(IB)
!         ENS = SQRT(ENX**2 + ENY**2)
!         ENX = ENX/ENS
!         ENY = ENY/ENS
!         XB(IB) = XB(IB) + ENX*DELN
!         YB(IB) = YB(IB) + ENY*DELN
!       ENDDO
!
!!----- re-spline new geometry
!       CALL SCALC(XB,YB,SB,NB)
!       CALL SEGSPL(XB,XBP,SB,NB)
!       CALL SEGSPL(YB,YBP,SB,NB)
!!
!       CALL GEOPAR(XB,XBP,YB,YBP,SB,NB,W1,
!     &             SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
!     &             EI11BA,EI22BA,APX1BA,APX2BA,
!     &             EI11BT,EI22BT,APX1BT,APX2BT,
!     &             THICKB,CAMBRB )
!!
!       CALL NEWPEN(2)
!       CALL PLTAIR(XB,XBP,YB,YBP,SB,NB, XOFF,XSF,YOFF,YSF,'magenta')
!       CALL PLNEWP('magenta')
!       LGEOPL = .FALSE.
!       LGSAME = .FALSE.
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'MOVP') THEN
!       CALL MOVP(NEWPLOTG)
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'UNIT') THEN
!       CALL NORM(XB,XBP,YB,YBP,SB,NB)
!       LGSAME = .FALSE.
!!
!!----- re-spline new geometry
!       CALL SCALC(XB,YB,SB,NB)
!       CALL SEGSPL(XB,XBP,SB,NB)
!       CALL SEGSPL(YB,YBP,SB,NB)
!!
!       CALL GEOPAR(XB,XBP,YB,YBP,SB,NB,W1,
!     &             SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
!     &             EI11BA,EI22BA,APX1BA,APX2BA,
!     &             EI11BT,EI22BT,APX1BT,APX2BT,
!     &             THICKB,CAMBRB )
!!
!       CALL NEWPEN(2)
!       CALL PLTAIR(XB,XBP,YB,YBP,SB,NB, XOFF,XSF,YOFF,YSF,'magenta')
!       CALL PLNEWP('magenta')
!       LGEOPL = .FALSE.
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'DIST' .OR.
!     &       COMAND.EQ.'D   '      ) THEN
!       CALL DIST
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'HARD') THEN
!       IF(LPLOT) CALL PLEND
!       LPLOT = .FALSE.
!       CALL REPLOT(IDEVRP)
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'PLOT' .OR.
!     &       COMAND.EQ.'P   '      ) THEN
!       CALL PLTINI
!       CALL PLOTG
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'INPL') THEN
!       CALL PLTINI
!       XOFF0 = XOFF
!       YOFF0 = YOFF
!       XSF0 = XSF
!       YSF0 = YSF
!!
!       XSF = 1.0/SIZE
!       YSF = 1.0/SIZE
!!       write(*,*) 'Enter Xoff, Yoff'
!!       read (*,*) xoff, yoff
!!       xoff = -xoff
!!       yoff = -yoff
!!
!       CALL PLOTG
!       XOFF = XOFF0
!       YOFF = YOFF0
!       XSF = XSF0
!       YSF = YSF0
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'BLOW' .OR.
!     &       COMAND.EQ.'B   '      ) THEN
!       XWS = XWIND/SIZE
!       YWS = YWIND/SIZE
!       CALL OFFGET(XOFF,YOFF,XSF,YSF,XWS,YWS, .TRUE. , .TRUE. )
!       CALL GOFSET
!       CALL PLTINI
!       CALL PLOTG
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'RESE' .OR.
!     &       COMAND.EQ.'R   '      ) THEN
!       CALL PLTINI
!       CALL GOFINI
!       CALL PLOTG
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'WIND' .OR.
!     &       COMAND.EQ.'W   '      ) THEN
!       XWS = XWIND/SIZE
!       YWS = YWIND/SIZE
!!
!       WRITE(*,*) ' '
!       WRITE(*,*) 'Type I,O,P to In,Out,Pan with cursor...'
!!
! 80    CALL PLTINI
!       CALL PLOTG
!!
!       CALL GETCURSORXY(XCRS,YCRS,CHKEY)
!!
!!----- do possible pan,zoom operations based on CHKEY
!       CALL KEYOFF(XCRS,YCRS,CHKEY, XWS,YWS, XOFF,YOFF,XSF,YSF, LPLNEW)
!!
!       IF(LPLNEW) THEN
!        CALL GOFSET
!        GO TO 80
!       ENDIF
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'TSIZ') THEN
!       IF(NINPUT.GE.1) THEN
!        GTICK = RINPUT(1)
!       ELSE
!        WRITE(*,*)
!     &    'Current tick-mark size (as fraction of perimeter) =', GTICK
!        CALL ASKR('Enter new tick-mark size^',GTICK)
!       ENDIF
!       CALL PLTINI
!       CALL PLOTG
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'TICK') THEN
!       LGTICK = .NOT.LGTICK
!       CALL PLTINI
!       CALL PLOTG
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'GRID') THEN
!       LGGRID = .NOT.LGGRID
!       CALL PLTINI
!       CALL PLOTG
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'GPAR') THEN
!       LGPARM = .NOT.LGPARM
!       CALL PLTINI
!       CALL PLOTG
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'SIZE') THEN
!       IF(NINPUT.GE.1) THEN
!        SIZE = RINPUT(1)
!       ELSE
!        WRITE(*,*) 'Current plot-object size =', SIZE
!        CALL ASKR('Enter new plot-object size^',SIZE)
!       ENDIF
!       CALL PLTINI
!       CALL PLOTG
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'OVER' .OR.
!     &       COMAND.EQ.'O   '      ) THEN
!       CALL OVER(COMARG)
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'ANNO') THEN
!       IF(LPLOT) THEN
!        CALL ANNOT(CH)
!       ELSE
!        WRITE(*,*) 'No active plot to annotate'
!       ENDIF
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'NAME') THEN
!       IF(COMARG.EQ.' ') THEN
!        CALL NAMMOD(NAME,0,-1)
!       ELSE
!        NAME = COMARG
!       ENDIF
!       CALL STRIP(NAME,NNAME)
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'NINC') THEN
!       CALL NAMMOD(NAME,1,1)
!       CALL STRIP(NAME,NNAME)
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'NDEC') THEN
!       CALL NAMMOD(NAME,-1,1)
!       CALL STRIP(NAME,NNAME)
!!
!!--------------------------------------------------------
!      ELSEIF(COMAND.EQ.'SINT') THEN
!       CALL SPLNXY(XB,XBP,YB,YBP,SB,NB)
!       CALL PLTAIR(XB,XBP,YB,YBP,SB,NB, XOFF,XSF,YOFF,YSF,'cyan')
!       CALL NEWPEN(2)
!       CALL PLTAIR(XB,XBP,YB,YBP,SB,NB, XOFF,XSF,YOFF,YSF,'magenta')
!       CALL PLNEWP('magenta')
!       LGEOPL = .FALSE.
!!
!!--------------------------------------------------------
!      ELSE
!       WRITE(*,1100) COMAND
! 1100  FORMAT(' Command ',A4,' not recognized.  Type a " ? " for list.')
!       COMAND = '****'
!      ENDIF
!!
!      GO TO 500
!      END ! GDES
!
!
!      SUBROUTINE NEWPLOTG
!      CALL GOFSET
!      CALL PLTINI
!      CALL PLOTG
!      RETURN
!      END
!
!
        SUBROUTINE ABCOPY(LCONF)
        INCLUDE 'XFOIL.INC'
        LOGICAL LCONF
  !
        IF(NB.LE.1) THEN
         WRITE(*,*) 'ABCOPY: Buffer airfoil not available.'
         RETURN
        ELSEIF(NB.GT.IQX-5) THEN
         WRITE(*,*) 'Maximum number of panel nodes  : ',IQX-5
         WRITE(*,*) 'Number of buffer airfoil points: ',NB
         WRITE(*,*) 'Current airfoil cannot be set.'
         WRITE(*,*) 'Try executing PANE at Top Level instead.'
         RETURN
        ENDIF
        IF(N.NE.NB) LBLINI = .FALSE.
  !
        N = NB
        DO 101 I=1, N
          X(I) = XB(I)
          Y(I) = YB(I)
    101 CONTINUE
        LGSAME = .TRUE.
  !
        IF(LBFLAP) THEN
         XOF = XBF
         YOF = YBF
         LFLAP = .TRUE.
        ENDIF
  !
  !---- strip out doubled points
        I = 1
   102  CONTINUE
        I = I+1
        IF(X(I-1).EQ.X(I) .AND. Y(I-1).EQ.Y(I)) THEN
          DO 104 J=I, N-1
            X(J) = X(J+1)
            Y(J) = Y(J+1)
   104    CONTINUE
          N = N-1
        ENDIF
        IF(I.LT.N) GO TO 102
  !
        CALL SCALC(X,Y,S,N)
        CALL SEGSPL(X,XP,S,N)
        CALL SEGSPL(Y,YP,S,N)
  
        CALL NCALC(X,Y,S,N,NX,NY)
  
        CALL LEFIND(SLE,X,XP,Y,YP,S,N)
        XLE = SEVAL(SLE,X,XP,S,N)
        YLE = SEVAL(SLE,Y,YP,S,N)
        XTE = 0.5*(X(1)+X(N))
        YTE = 0.5*(Y(1)+Y(N))
        CHORD  = SQRT( (XTE-XLE)**2 + (YTE-YLE)**2 )
  
        CALL TECALC
        CALL APCALC
  !
        LGAMU = .FALSE.
        LQINU = .FALSE.
        LWAKE = .FALSE.
        LQAIJ = .FALSE.
        LADIJ = .FALSE.
        LWDIJ = .FALSE.
        LIPAN = .FALSE.
        LVCONV = .FALSE.
        LSCINI = .FALSE.
  !CC      LBLINI = .FALSE.
  !
        IF (DBUGMOD) THEN
        IF(LCONF) WRITE(*,1200) N
   1200 FORMAT(/' Current airfoil nodes set from buffer airfoil nodes (',    &
                I4,' )')
        ENDIF
  !
        RETURN
        END ! ABCOPY
  
!
!      SUBROUTINE GOFINI
!!----------------------------------------------------------
!!     Sets initial airfoil scaling and offset parameters   
!!----------------------------------------------------------
!      INCLUDE 'XFOIL.INC'
!!
!!---- get airfoil bounding box
!      XBMIN = XB(1)
!      YBMIN = YB(1)
!      XBMAX = XB(1)
!      YBMAX = YB(1)
!      DO I=1, NB
!        XBMIN = MIN(XBMIN,XB(I))
!        YBMIN = MIN(YBMIN,YB(I))
!        XBMAX = MAX(XBMAX,XB(I))
!        YBMAX = MAX(YBMAX,YB(I))
!      ENDDO
!!
!!---- set camber and thickness distributions
!      CALL GETCAM(XCM,YCM,NCM,XTK,YTK,NTK,
!     &             XB,XBP,YB,YBP,SB,NB )
!!      
!!---- get camber,thickness y bounds
!      CMMIN = 0.
!      CMMAX = 0.
!      DO I=1, NCM
!        CMMIN = MIN(CMMIN,YCM(I))
!        CMMAX = MAX(CMMAX,YCM(I))
!      ENDDO
!      TKMIN = 0.
!      TKMAX = 0.
!      DO I=1, NTK
!        TKMIN = MIN(TKMIN,YTK(I))
!        TKMAX = MAX(TKMAX,YTK(I))
!      ENDDO
!!
!      XRANGE = XBMAX - XBMIN
!      YRANGE = YBMAX - YBMIN
!!
!!---- set x,y scaling factors needed for O(1) size plot with "nice" limits
!      CALL SCALIT(1,0.95*XRANGE,0.0,XSF)
!      CALL SCALIT(1,0.95*YRANGE,0.0,YSF)
!!  
!!---- grid increment as a fraction of a nice upper bound on delta x
!!c      DXYG = 0.1 / XSF
!      DXYG = 0.1 / MIN(XSF,YSF)
!!  
!!---- set "nice" grid limits as integer multiples of DXYG
!!      XGMAX = DXYG*(INT(XBMAX/DXYG+1000.05) - 999)
!!      XGMIN = DXYG*(INT(XBMIN/DXYG-1000.05) + 999)
!!      YGMAX = DXYG*(INT(YBMAX/DXYG+1000.25) - 999)
!!      YGMIN = DXYG*(INT(YBMIN/DXYG-1000.25) + 999)
!!  
!!---- set "nice" grid limits as integer multiples of DXYG
!      XGMAX = DXYG*(INT(XBMAX/DXYG+1001.01) - 1000)
!      XGMIN = DXYG*(INT(XBMIN/DXYG-1001.01) + 1000)
!      YGMAX = DXYG*(INT(YBMAX/DXYG+1001.01) - 1000)
!      YGMIN = DXYG*(INT(YBMIN/DXYG-1001.01) + 1000)
!!
!!---- set bounding box for thickness/camber plot
!      DXYC = DXYG
!      XCMIN = XGMIN
!      XCMAX = XGMAX
!      YCMIN = MIN(CMMIN,-TKMAX)
!      YCMAX = MAX(CMMAX, TKMAX)
!      YCMAX = DXYC*(INT(YCMAX/DXYC+1000.25) - 999)
!      YCMIN = DXYC*(INT(YCMIN/DXYC-1000.25) + 999)
!      YCMAX = MAX(YCMAX,YCMIN+DXYC)
!!
!!---- set minimum scaling factor to fit airfoil or grid
!      IF(LGGRID) THEN
!        XRANGE = XGMAX - XGMIN
!        YRANGE = YGMAX - YGMIN
!      ELSE
!        XRANGE = XBMAX - XBMIN
!        YRANGE = YBMAX - YBMIN
!      ENDIF
!!
!!---- include y range from thickness/camber plot if present
!      IF(LPLCAM) THEN
!        YRANGE = YRANGE + (YCMAX - YCMIN)
!      ENDIF
!!
!      RANGE = MAX(XRANGE,YRANGE)
!!
!      SF = MIN( 1.0/XRANGE , PLOTAR/YRANGE )
!      XSF = SF
!      YSF = SF
!      CHG = 0.75*CH * RANGE*SF
!!--- HHY 4/24/01 keep the character size from getting too low
!
!      CHG = MAX(CHG,0.0075)
!!
!      IF(LGGRID) THEN
!!------ set offsets to position grid, with space for numerical axis annotations
!        XOFF = XGMIN - 0.05*RANGE - 3.0*CHG/SF
!        YOFF = YGMIN - 0.05*RANGE - 2.0*CHG/SF
!      ELSE
!!------ set offsets to position airfoil
!        XOFF = XBMIN - 0.05*RANGE
!        YOFF = YBMIN - 0.05*RANGE
!      ENDIF
!!
!!---- set plot limits for DCp plot (y-axis limit defaults set in INIT)
!      XPMIN = XGMIN
!      XPMAX = XGMAX
!!cc      DXYP = DXYG
!      CALL AXISADJ(YPMIN,YPMAX,PSPAN,DXYP,NTICS)
!!
!!---- set Yoffset for camber plot in scale factor YSF for geom plots
!      DYOFFC = - YGMAX + YCMIN - 2.2*CHG/YSF
!!
!!---- set the Cp scale factor for DCp plots
!      PAR = (YPAGE-2.0*YMARG)/(XPAGE-2.0*XMARG)
!      DPRANGE = YPMAX-YPMIN
!      DYPLT = MAX(0.1,PAR-PLOTAR)
!      YSFP = 0.8*DYPLT/DPRANGE
!      YSFP = YSFP/YSF
!!
!!---- set shifts to YOFF for DCp plots in scale factor YSF for geom plots
!      DYOFFP = -YCMAX+DYOFFC + YPMIN*YSFP - 2.2*CHG/YSF
!!
!      RETURN
!      END ! GOFINI  
!
!
!
!      SUBROUTINE GOFSET
!!----------------------------------------------------------
!!     Sets grid-overlay parameters
!!----------------------------------------------------------
!      INCLUDE 'XFOIL.INC'
!!
!!---- airfoil extent
!      XBMIN = XB(1)
!      YBMIN = YB(1)
!      XBMAX = XB(1)
!      YBMAX = YB(1)
!      DO I=1, NB
!        XBMIN = MIN(XBMIN,XB(I))
!        YBMIN = MIN(YBMIN,YB(I))
!        XBMAX = MAX(XBMAX,XB(I))
!        YBMAX = MAX(YBMAX,YB(I))
!      ENDDO
!!
!      RANGE = MAX( (XWIND/SIZE)/XSF , (YWIND/SIZE)/YSF )
!!
!!---- set bounding-box corner locations in user coordinates
!      XG1 = XOFF + 0.1*RANGE + 4.0*CHG/XSF
!      YG1 = YOFF + 0.1*RANGE + 2.0*CHG/YSF
!      XG2 = XOFF - 0.1*RANGE + (XWIND/SIZE)/XSF
!      YG2 = YOFF - 0.1*RANGE + (YWIND/SIZE)/YSF
!!
!!---- crunch down onto airfoil limits
!      XG1 = MAX(XG1,XBMIN)
!      XG2 = MIN(XG2,XBMAX)
!      YG1 = MAX(YG1,YBMIN)
!      YG2 = MIN(YG2,YBMAX)
!!
!!---- set x,y scaling factors needed for O(1) size plot with "nice" limits
!      CALL SCALIT(1,0.95*(XG2-XG1),0.0,GXSF)
!      CALL SCALIT(1,0.95*(YG2-YG1),0.0,GYSF)
!!  
!      GSF = GXSF
!!cc   GSF = MIN(GXSF,GYSF)
!!
!!---- grid increment as a fraction of a nice upper bound on delta x
!      DXYG = 0.1 / GSF
!!  
!!---- set "nice" grid limits as integer multiples of DXYG
!      XGMAX = DXYG*(INT(XG2/DXYG+1001.01) - 1000)
!      XGMIN = DXYG*(INT(XG1/DXYG-1001.01) + 1000)
!      YGMAX = DXYG*(INT(YG2/DXYG+1001.01) - 1000)
!      YGMIN = DXYG*(INT(YG1/DXYG-1001.01) + 1000)
!!
!      RETURN
!      END ! GOFSET
!
!
!
!      SUBROUTINE TGAP(RINPUT,NINPUT)
!!----------------------------------
!!     Used to set buffer airfoil 
!!     trailing edge gap
!!----------------------------------
!      INCLUDE 'XFOIL.INC'
!      DIMENSION RINPUT(*)
!!
!      CALL LEFIND(SBLE,XB,XBP,YB,YBP,SB,NB)
!      XBLE = SEVAL(SBLE,XB,XBP,SB,NB)
!      YBLE = SEVAL(SBLE,YB,YBP,SB,NB)
!      XBTE = 0.5*(XB(1)+XB(NB))
!      YBTE = 0.5*(YB(1)+YB(NB))
!      CHBSQ = (XBTE-XBLE)**2 + (YBTE-YBLE)**2
!!
!      DXN = XB(1) - XB(NB)
!      DYN = YB(1) - YB(NB)
!      GAP = SQRT(DXN**2 + DYN**2)
!!
!!---- components of unit vector parallel to TE gap
!      IF(GAP.GT.0.0) THEN
!       DXU = DXN / GAP
!       DYU = DYN / GAP
!      ELSE
!       DXU = -.5*(YBP(NB) - YBP(1))
!       DYU = 0.5*(XBP(NB) - XBP(1))
!      ENDIF
!!
!      IF    (NINPUT .GE. 2) THEN
!       GAPNEW = RINPUT(1)
!       DOC    = RINPUT(2)
!      ELSEIF(NINPUT .GE. 1) THEN
!       GAPNEW = RINPUT(1)
!       DOC = 1.0
!       CALL ASKR('Enter blending distance/c (0..1)^',DOC)
!      ELSE
!       WRITE(*,1000) GAP
! 1000  FORMAT(/' Current gap =',F9.5)
!       GAPNEW = 0.0
!       CALL ASKR('Enter new gap^',GAPNEW)
!       DOC = 1.0
!       CALL ASKR('Enter blending distance/c (0..1)^',DOC)
!      ENDIF
!!
!      DOC = MIN( MAX( DOC , 0.0 ) , 1.0 )
!!
!      DGAP = GAPNEW - GAP
!!
!!---- go over each point, changing the y-thickness appropriately
!      DO 30 I=1, NB
!!
!!------ chord-based x/c
!        XOC = (  (XB(I)-XBLE)*(XBTE-XBLE)
!     &         + (YB(I)-YBLE)*(YBTE-YBLE) ) / CHBSQ
!!
!!------ thickness factor tails off exponentially away from trailing edge
!        IF(DOC .EQ. 0.0) THEN
!          TFAC = 0.0
!          IF(I.EQ.1 .OR. I.EQ.NB) TFAC = 1.0
!        ELSE
!          ARG = MIN( (1.0-XOC)*(1.0/DOC-1.0) , 15.0 )
!          TFAC = EXP(-ARG)
!        ENDIF
!!
!        IF(SB(I).LE.SBLE) THEN
!         XB(I) = XB(I) + 0.5*DGAP*XOC*TFAC*DXU
!         YB(I) = YB(I) + 0.5*DGAP*XOC*TFAC*DYU
!        ELSE
!         XB(I) = XB(I) - 0.5*DGAP*XOC*TFAC*DXU
!         YB(I) = YB(I) - 0.5*DGAP*XOC*TFAC*DYU
!        ENDIF
!   30 CONTINUE
!      LGSAME = .FALSE.
!!
!      CALL SCALC(XB,YB,SB,NB)
!      CALL SEGSPL(XB,XBP,SB,NB)
!      CALL SEGSPL(YB,YBP,SB,NB)
!!
!      CALL GEOPAR(XB,XBP,YB,YBP,SB,NB,W1,
!     &            SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
!     &            EI11BA,EI22BA,APX1BA,APX2BA,
!     &            EI11BT,EI22BT,APX1BT,APX2BT,
!     &            THICKB,CAMBRB )
!!
!      CALL PLTAIR(XB,XBP,YB,YBP,SB,NB, XOFF,XSF,YOFF,YSF,'magenta')
!      CALL PLNEWP('magenta')
!!
!      LGEOPL = .FALSE.
!!
!      RETURN
!      END ! TGAP
!
!
!
!      SUBROUTINE LERAD(RINPUT,NINPUT)
!!----------------------------
!!     Changes buffer airfoil 
!!     leading edge radius.
!!----------------------------
!      INCLUDE 'XFOIL.INC'
!      DIMENSION RINPUT(*)
!!
!      IF    (NINPUT .GE. 2) THEN
!       RFAC = RINPUT(1)
!       DOC  = RINPUT(2)
!      ELSEIF(NINPUT .GE. 1) THEN
!       RFAC = RINPUT(1)
!       DOC = 1.0
!       CALL ASKR('Enter blending distance/c from LE^',DOC)
!      ELSE
!       RFAC = 1.0
!       CALL ASKR('Enter approx. new/old LE radius scaling ratio^',RFAC)
!       DOC = 1.0
!       CALL ASKR('Enter blending distance/c from LE^',DOC)
!      ENDIF
!!
!      DOC = MAX( DOC , 0.001 )
!!
!      CALL LERSCL(XB,XBP,YB,YBP,SB,NB, DOC,RFAC, W1,W2)
!!
!      DO 40 I=1, NB
!        XB(I) = W1(I)
!        YB(I) = W2(I)
!   40 CONTINUE
!      LGSAME = .FALSE.
!!
!!---- spline new coordinates
!      CALL SCALC(XB,YB,SB,NB)
!      CALL SEGSPL(XB,XBP,SB,NB)
!      CALL SEGSPL(YB,YBP,SB,NB)
!!
!      CALL GEOPAR(XB,XBP,YB,YBP,SB,NB,W1,
!     &            SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
!     &            EI11BA,EI22BA,APX1BA,APX2BA,
!     &            EI11BT,EI22BT,APX1BT,APX2BT,
!     &            THICKB,CAMBRB )
!!
!!---- find max curvature
!      CVMAX = 0.
!      DO 6 I=NB/4, (3*NB)/4
!        CV = CURV(SB(I),XB,XBP,YB,YBP,SB,NB)
!        CVMAX = MAX( ABS(CV) , CVMAX )
!    6 CONTINUE
!!
!      RADIUS = 1.0/CVMAX
!!
!      WRITE(*,1000) RADIUS
! 1000 FORMAT(/' New LE radius = ',F7.5)
!!
!      CALL PLTAIR(XB,XBP,YB,YBP,SB,NB, XOFF,XSF,YOFF,YSF,'magenta')
!      CALL PLNEWP('magenta')
!!
!      LGEOPL = .FALSE.
!!
!      RETURN
!      END ! LERAD
!
!
!
!      SUBROUTINE SCLXY
!!---------------------------------------------------
!!     Scale airfoil about LE, TE, or selected point 
!!---------------------------------------------------
!      INCLUDE 'XFOIL.INC'
!      CHARACTER*1 VAR
!!
!      CALL LEFIND(SBLE,XB,XBP,YB,YBP,SB,NB)
!      XLE = SEVAL(SBLE,XB,XBP,SB,NB)
!      YLE = SEVAL(SBLE,YB,YBP,SB,NB)
!      XTE = 0.5*(XB(1) + XB(NB))
!      YTE = 0.5*(YB(1) + YB(NB))
!!
!      WRITE(*,*) 'Enter origin for airfoil scaling:'
!      WRITE(*,*) '  L  scales about LE'
!      WRITE(*,*) '  T  scales about TE'
!      WRITE(*,*) '  P  scales about input point'
!!      
!      CALL ASKS('Select origin for scaling^',VAR)
!      IF (VAR.EQ.'L') THEN
!        XORG = XLE
!        YORG = YLE
!       ELSE IF (VAR.EQ.'T') THEN
!        XORG = XTE
!        YORG = YTE
!       ELSE 
!        XORG = 0.25
!        YORG = 0.0
!        CALL ASKR('Enter X origin for scaling^',XORG)
!        CALL ASKR('Enter Y origin for scaling^',YORG)
!      ENDIF       
!!      
!      SCL = 1.0
!      CALL ASKR('Enter scaling factor about selected point^',SCL)
!!
!      DO 10 I=1, NB
!        XB(I) = SCL*(XB(I) - XORG) + XORG 
!        YB(I) = SCL*(YB(I) - YORG) + YORG 
!   10 CONTINUE
!      LGSAME = .FALSE.
!!
!      CALL SCALC(XB,YB,SB,NB)
!      CALL SEGSPL(XB,XBP,SB,NB)
!      CALL SEGSPL(YB,YBP,SB,NB)
!!
!      CALL GEOPAR(XB,XBP,YB,YBP,SB,NB,W1,
!     &            SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
!     &            EI11BA,EI22BA,APX1BA,APX2BA,
!     &            EI11BT,EI22BT,APX1BT,APX2BT,
!     &            THICKB,CAMBRB )
!!
!      RETURN
!      END ! SCLXY
!
!
!
!      SUBROUTINE FLAP(RINPUT,NINPUT)
!!----------------------------------------------------
!!     Modifies buffer airfoil for a deflected flap.
!!     Points may be added/subtracted in the flap
!!     break vicinity to clean things up.
!!----------------------------------------------------
!      INCLUDE 'XFOIL.INC'
!      LOGICAL LCHANGE
!      DIMENSION RINPUT(*)
!!
!      LOGICAL INSID
!      LOGICAL INSIDE
!      LOGICAL LT1NEW,LT2NEW,LB1NEW,LB2NEW
!!
!      SHT = CH * MAX(XSF,YSF)
!!
!      IF(NINPUT.GE.2) THEN
!       XBF = RINPUT(1)
!       YBF = RINPUT(2)
!      ELSE
!       XBF = -999.0
!       YBF = -999.0
!      ENDIF
!!
!      CALL GETXYF(XB,XBP,YB,YBP,SB,NB, TOPS,BOTS,XBF,YBF)
!      INSID = INSIDE(XB,YB,NB,XBF,YBF)
!!
!      WRITE(*,1050) XBF, YBF
! 1050 FORMAT(/' Flap hinge: x,y =', 2F9.5 )
!!
!      IF(NINPUT.GE.3) THEN
!       DDEF = RINPUT(3)
!      ELSE
!       DDEF = 0.0
!       CALL ASKR('Enter flap deflection in degrees (+ down)^',DDEF)
!      ENDIF
!      RDEF = DDEF*PI/180.0
!      IF(RDEF .EQ. 0.0) RETURN
!!
!!
!      IF(INSID) THEN
!        ATOP = MAX( 0.0 , -RDEF )
!        ABOT = MAX( 0.0 ,  RDEF )
!      ELSE
!        CHX = DEVAL(BOTS,XB,XBP,SB,NB) - DEVAL(TOPS,XB,XBP,SB,NB)
!        CHY = DEVAL(BOTS,YB,YBP,SB,NB) - DEVAL(TOPS,YB,YBP,SB,NB)
!        FVX = SEVAL(BOTS,XB,XBP,SB,NB) + SEVAL(TOPS,XB,XBP,SB,NB)
!        FVY = SEVAL(BOTS,YB,YBP,SB,NB) + SEVAL(TOPS,YB,YBP,SB,NB)
!        CRSP = CHX*(YBF-0.5*FVY) - CHY*(XBF-0.5*FVX)
!        IF(CRSP .GT. 0.0) THEN
!!-------- flap hinge is above airfoil
!          ATOP = MAX( 0.0 ,  RDEF )
!          ABOT = MAX( 0.0 ,  RDEF )
!        ELSE
!!-------- flap hinge is below airfoil
!          ATOP = MAX( 0.0 , -RDEF )
!          ABOT = MAX( 0.0 , -RDEF )
!        ENDIF
!      ENDIF
!!
!!---- find upper and lower surface break arc length values...
!      CALL SSS(TOPS,ST1,ST2,ATOP,XBF,YBF,XB,XBP,YB,YBP,SB,NB,1)
!      CALL SSS(BOTS,SB1,SB2,ABOT,XBF,YBF,XB,XBP,YB,YBP,SB,NB,2)
!!
!!---- ... and x,y coordinates
!      XT1 = SEVAL(ST1,XB,XBP,SB,NB)
!      YT1 = SEVAL(ST1,YB,YBP,SB,NB)
!      XT2 = SEVAL(ST2,XB,XBP,SB,NB)
!      YT2 = SEVAL(ST2,YB,YBP,SB,NB)
!      XB1 = SEVAL(SB1,XB,XBP,SB,NB)
!      YB1 = SEVAL(SB1,YB,YBP,SB,NB)
!      XB2 = SEVAL(SB2,XB,XBP,SB,NB)
!      YB2 = SEVAL(SB2,YB,YBP,SB,NB)
!!
!!
!      WRITE(*,1100) XT1, YT1, XT2, YT2,
!     &              XB1, YB1, XB2, YB2
! 1100 FORMAT(/' Top breaks: x,y =  ', 2F9.5, 4X, 2F9.5
!     &       /' Bot breaks: x,y =  ', 2F9.5, 4X, 2F9.5)
!!
!!---- find points adjacent to breaks
!      DO 5 I=1, NB-1
!        IF(SB(I).LE.ST1 .AND. SB(I+1).GT.ST1) IT1 = I+1
!        IF(SB(I).LT.ST2 .AND. SB(I+1).GE.ST2) IT2 = I
!        IF(SB(I).LE.SB1 .AND. SB(I+1).GT.SB1) IB1 = I
!        IF(SB(I).LT.SB2 .AND. SB(I+1).GE.SB2) IB2 = I+1
!    5 CONTINUE
!!
!      DSAVG = (SB(NB)-SB(1))/FLOAT(NB-1)
!!
!!---- smallest fraction of s increments i+1 and i+2 away from break point
!      SFRAC = 0.33333
!!
!      IF(ATOP .NE. 0.0) THEN
!        ST1P = ST1 + SFRAC*(SB(IT1  )-ST1)
!        ST1Q = ST1 + SFRAC*(SB(IT1+1)-ST1)
!        IF(SB(IT1) .LT. ST1Q) THEN
!!-------- simply move adjacent point to ideal SFRAC location
!          XT1NEW = SEVAL(ST1Q,XB,XBP,SB,NB)
!          YT1NEW = SEVAL(ST1Q,YB,YBP,SB,NB)
!          LT1NEW = .FALSE.
!        ELSE
!!-------- make new point at SFRAC location
!          XT1NEW = SEVAL(ST1P,XB,XBP,SB,NB)
!          YT1NEW = SEVAL(ST1P,YB,YBP,SB,NB)
!          LT1NEW = .TRUE.
!        ENDIF
!!
!        ST2P = ST2 + SFRAC*(SB(IT2 )-ST2)
!        IT2Q = MAX(IT2-1,1)
!        ST2Q = ST2 + SFRAC*(SB(IT2Q)-ST2)
!        IF(SB(IT2) .GT. ST2Q) THEN
!!-------- simply move adjacent point
!          XT2NEW = SEVAL(ST2Q,XB,XBP,SB,NB)
!          YT2NEW = SEVAL(ST2Q,YB,YBP,SB,NB)
!          LT2NEW = .FALSE.
!        ELSE
!!-------- make new point
!          XT2NEW = SEVAL(ST2P,XB,XBP,SB,NB)
!          YT2NEW = SEVAL(ST2P,YB,YBP,SB,NB)
!          LT2NEW = .TRUE.
!        ENDIF
!      ENDIF
!!
!      IF(ABOT .NE. 0.0) THEN
!        SB1P = SB1 + SFRAC*(SB(IB1  )-SB1)
!        SB1Q = SB1 + SFRAC*(SB(IB1-1)-SB1)
!        IF(SB(IB1) .GT. SB1Q) THEN
!!-------- simply move adjacent point
!          XB1NEW = SEVAL(SB1Q,XB,XBP,SB,NB)
!          YB1NEW = SEVAL(SB1Q,YB,YBP,SB,NB)
!          LB1NEW = .FALSE.
!        ELSE
!!-------- make new point
!          XB1NEW = SEVAL(SB1P,XB,XBP,SB,NB)
!          YB1NEW = SEVAL(SB1P,YB,YBP,SB,NB)
!          LB1NEW = .TRUE.
!        ENDIF
!!
!        SB2P = SB2 + SFRAC*(SB(IB2 )-SB2)
!        IB2Q = MIN(IB2+1,NB)
!        SB2Q = SB2 + SFRAC*(SB(IB2Q)-SB2)
!        IF(SB(IB2) .LT. SB2Q) THEN
!!-------- simply move adjacent point
!          XB2NEW = SEVAL(SB2Q,XB,XBP,SB,NB)
!          YB2NEW = SEVAL(SB2Q,YB,YBP,SB,NB)
!          LB2NEW = .FALSE.
!        ELSE
!!-------- make new point
!          XB2NEW = SEVAL(SB2P,XB,XBP,SB,NB)
!          YB2NEW = SEVAL(SB2P,YB,YBP,SB,NB)
!          LB2NEW = .TRUE.
!        ENDIF
!      ENDIF
!!
!!c      DSTOP = ABS(SB(IT2)-SB(IT1))
!!c      DSBOT = ABS(SB(IB2)-SB(IB1))
!!
!      SIND = SIN(RDEF)
!      COSD = COS(RDEF)
!!
!!---- rotate flap points about the hinge point (XBF,YBF)
!      DO 10 I=1, NB
!        IF(I.GE.IT1 .AND. I.LE.IB1) GO TO 10
!!
!        XBAR = XB(I) - XBF
!        YBAR = YB(I) - YBF
!!
!        XB(I) = XBF  +  XBAR*COSD  +  YBAR*SIND
!        YB(I) = YBF  -  XBAR*SIND  +  YBAR*COSD
!   10 CONTINUE
!!
!      IDIF = IT1-IT2-1
!      IF(IDIF.GT.0) THEN
!!----- delete points on upper airfoil surface which "disappeared".
!       NB  = NB -IDIF
!       IT1 = IT1-IDIF
!       IB1 = IB1-IDIF
!       IB2 = IB2-IDIF
!       DO 21 I=IT2+1, NB
!         SB(I) = SB(I+IDIF)
!         XB(I) = XB(I+IDIF)
!         YB(I) = YB(I+IDIF)
!   21  CONTINUE
!      ENDIF
!!
!      IDIF = IB2-IB1-1
!      IF(IDIF.GT.0) THEN
!!----- delete points on lower airfoil surface which "disappeared".
!       NB  = NB -IDIF
!       IB2 = IB2-IDIF
!       DO 22 I=IB1+1, NB
!         SB(I) = SB(I+IDIF)
!         XB(I) = XB(I+IDIF)
!         YB(I) = YB(I+IDIF)
!   22  CONTINUE
!      ENDIF
!!
!!
!      IF(ATOP .EQ. 0.0) THEN
!!
!!------ arc length of newly created surface on top of airfoil
!        DSNEW = ABS(RDEF)*SQRT((XT1-XBF)**2 + (YT1-YBF)**2)
!!
!!------ number of points to be added to define newly created surface
!        NPADD = INT(1.5*DSNEW/DSAVG + 1.0)
!!cc     NPADD = INT(1.5*DSNEW/DSTOP + 1.0)
!!
!!------ skip everything if no points are to be added
!        IF(NPADD.EQ.0) GO TO 35
!!
!!------ increase coordinate array length to make room for the new point(s)
!        NB  = NB +NPADD
!        IT1 = IT1+NPADD
!        IB1 = IB1+NPADD
!        IB2 = IB2+NPADD
!        DO 30 I=NB, IT1, -1
!          XB(I) = XB(I-NPADD)
!          YB(I) = YB(I-NPADD)
!   30   CONTINUE
!!
!!------ add new points along the new surface circular arc segment
!        DANG = RDEF / FLOAT(NPADD)
!        XBAR = XT1 - XBF
!        YBAR = YT1 - YBF
!        DO 31 IP=1, NPADD
!          ANG = DANG*(FLOAT(IP) - 0.5)
!          CA = COS(ANG)
!          SA = SIN(ANG)
!!
!          XB(IT1-IP) = XBF  +  XBAR*CA + YBAR*SA
!          YB(IT1-IP) = YBF  -  XBAR*SA + YBAR*CA
!   31   CONTINUE
!!
!      ELSE
!!
!!------ set point in the corner and possibly two adjacent points
!        NPADD = 1
!        IF(LT2NEW) NPADD = NPADD+1
!        IF(LT1NEW) NPADD = NPADD+1
!!
!        NB  = NB +NPADD
!        IT1 = IT1+NPADD
!        IB1 = IB1+NPADD
!        IB2 = IB2+NPADD
!        DO 33 I=NB, IT1, -1
!          XB(I) = XB(I-NPADD)
!          YB(I) = YB(I-NPADD)
!   33   CONTINUE
!!
!        IF(LT1NEW) THEN
!         XB(IT1-1) = XT1NEW
!         YB(IT1-1) = YT1NEW
!         XB(IT1-2) = XT1
!         YB(IT1-2) = YT1
!        ELSE
!         XB(IT1  ) = XT1NEW
!         YB(IT1  ) = YT1NEW
!         XB(IT1-1) = XT1
!         YB(IT1-1) = YT1
!        ENDIF
!!
!        XBAR = XT2NEW - XBF
!        YBAR = YT2NEW - YBF
!        IF(LT2NEW) THEN
!          XB(IT2+1) = XBF  +  XBAR*COSD + YBAR*SIND
!          YB(IT2+1) = YBF  -  XBAR*SIND + YBAR*COSD
!        ELSE
!          XB(IT2  ) = XBF  +  XBAR*COSD + YBAR*SIND
!          YB(IT2  ) = YBF  -  XBAR*SIND + YBAR*COSD
!        ENDIF
!!
!      ENDIF
!   35 CONTINUE
!!
!!
!      IF(ABOT .EQ. 0.0) THEN
!!
!!------ arc length of newly created surface on top of airfoil
!        DSNEW = ABS(RDEF)*SQRT((XB1-XBF)**2 + (YB1-YBF)**2)
!!
!!------ number of points to be added to define newly created surface
!        NPADD = INT(1.5*DSNEW/DSAVG + 1.0)
!!cc     NPADD = INT(1.5*DSNEW/DSBOT + 1.0)
!!
!!------ skip everything if no points are to be added
!        IF(NPADD.EQ.0) GO TO 45
!!
!!------ increase coordinate array length to make room for the new point(s)
!        NB  = NB +NPADD
!        IB2 = IB2+NPADD
!        DO 40 I=NB, IB2, -1
!          XB(I) = XB(I-NPADD)
!          YB(I) = YB(I-NPADD)
!   40   CONTINUE
!!
!!------ add new points along the new surface circular arc segment
!        DANG = RDEF / FLOAT(NPADD)
!        XBAR = XB1 - XBF
!        YBAR = YB1 - YBF
!        DO 41 IP=1, NPADD
!          ANG = DANG*(FLOAT(IP) - 0.5)
!          CA = COS(ANG)
!          SA = SIN(ANG)
!!
!          XB(IB1+IP) = XBF  +  XBAR*CA + YBAR*SA
!          YB(IB1+IP) = YBF  -  XBAR*SA + YBAR*CA
!   41   CONTINUE
!!
!      ELSE
!
!!------ set point in the corner and possibly two adjacent points
!        NPADD = 1
!        IF(LB2NEW) NPADD = NPADD+1
!        IF(LB1NEW) NPADD = NPADD+1
!!
!        NB  = NB +NPADD
!        IB2 = IB2+NPADD
!        DO 43 I=NB, IB2, -1
!          XB(I) = XB(I-NPADD)
!          YB(I) = YB(I-NPADD)
!   43   CONTINUE
!!
!        IF(LB1NEW) THEN
!         XB(IB1+1) = XB1NEW
!         YB(IB1+1) = YB1NEW
!         XB(IB1+2) = XB1
!         YB(IB1+2) = YB1
!        ELSE
!         XB(IB1  ) = XB1NEW
!         YB(IB1  ) = YB1NEW
!         XB(IB1+1) = XB1
!         YB(IB1+1) = YB1
!        ENDIF
!!
!        XBAR = XB2NEW - XBF
!        YBAR = YB2NEW - YBF
!        IF(LB2NEW) THEN
!          XB(IB2-1) = XBF  +  XBAR*COSD + YBAR*SIND
!          YB(IB2-1) = YBF  -  XBAR*SIND + YBAR*COSD
!        ELSE
!          XB(IB2  ) = XBF  +  XBAR*COSD + YBAR*SIND
!          YB(IB2  ) = YBF  -  XBAR*SIND + YBAR*COSD
!        ENDIF
!!
!      ENDIF
!   45 CONTINUE
!!
!      LGSAME = .FALSE.
!!
!!
!!---- check new geometry for splinter segments 
!      STOL = 0.2
!      CALL SCHECK(XB,YB,NB, STOL, LCHANGE)
!!
!!---- spline new geometry
!      CALL SCALC(XB,YB,SB,NB)
!      CALL SEGSPL(XB,XBP,SB,NB)
!      CALL SEGSPL(YB,YBP,SB,NB)
!!
!      CALL GEOPAR(XB,XBP,YB,YBP,SB,NB,W1,
!     &            SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
!     &            EI11BA,EI22BA,APX1BA,APX2BA,
!     &            EI11BT,EI22BT,APX1BT,APX2BT,
!     &            THICKB,CAMBRB )
!!
!      LBFLAP = .TRUE.
!!
!      IF(LGSYM) THEN
!       WRITE(*,*)
!       WRITE(*,*) 'Disabling symmetry enforcement'
!       LGSYM = .FALSE.
!      ENDIF
!!
!!
!      IF(.NOT.LPLOT) THEN
!       CALL PLTINI
!      ENDIF
!!
!!---- save current color and set new color
!      CALL GETCOLOR(ICOL0)
!!
!      CALL NEWCOLORNAME('green')
!      CALL PLOT((XBF-XOFF)*XSF,(YBF-YOFF)*YSF,3)
!      CALL PLOT((XT1-XOFF)*XSF,(YT1-YOFF)*YSF,2)
!      CALL PLOT((XBF-XOFF)*XSF,(YBF-YOFF)*YSF,3)
!      CALL PLOT((XB1-XOFF)*XSF,(YB1-YOFF)*YSF,2)
!!
!      IF(ATOP .EQ. 0.0) THEN
!        XBAR = XT1 - XBF
!        YBAR = YT1 - YBF
!        XT1C = XBF  +  XBAR*COSD + YBAR*SIND
!        YT1C = YBF  -  XBAR*SIND + YBAR*COSD
!        CALL PLOT((XBF -XOFF)*XSF,(YBF -YOFF)*YSF,3)
!        CALL PLOT((XT1C-XOFF)*XSF,(YT1C-YOFF)*YSF,2)
!      ENDIF
!!
!      IF(ABOT .EQ. 0.0) THEN
!        XBAR = XB1 - XBF
!        YBAR = YB1 - YBF
!        XB1C = XBF  +  XBAR*COSD + YBAR*SIND
!        YB1C = YBF  -  XBAR*SIND + YBAR*COSD
!        CALL PLOT((XBF -XOFF)*XSF,(YBF -YOFF)*YSF,3)
!        CALL PLOT((XB1C-XOFF)*XSF,(YB1C-YOFF)*YSF,2)
!      ENDIF
!!
!      CALL NEWCOLORNAME('red')
!      CALL PLSYMB((XBF-XOFF)*XSF,(YBF-YOFF)*YSF,0.5*SHT,1,0.0,0)
!!
!      CALL PLTAIR(XB,XBP,YB,YBP,SB,NB, XOFF,XSF,YOFF,YSF,'magenta')
!      CALL PLNEWP('magenta')
!!
!      LGEOPL = .FALSE.
!!
!      CALL NEWCOLOR(ICOL0)
!      RETURN
!      END ! FLAP
!
!
!      LOGICAL FUNCTION INSIDE(X,Y,N, XF,YF)
!      DIMENSION X(N),Y(N)
!!-------------------------------------
!!     Returns .TRUE. if point XF,YF 
!!     is inside contour X(i),Y(i).
!!-------------------------------------
!!
!!---- integrate subtended angle around airfoil perimeter
!      ANGLE = 0.0
!      DO 10 I=1, N
!        IP = I+1
!        IF(I.EQ.N) IP = 1
!        XB1 = X(I)  - XF
!        YB1 = Y(I)  - YF
!        XB2 = X(IP) - XF
!        YB2 = Y(IP) - YF
!        ANGLE = ANGLE + (XB1*YB2 - YB1*XB2)
!     &                   / SQRT((XB1**2 + YB1**2)*(XB2**2 + YB2**2))
! 10   CONTINUE
!!
!!---- angle = 0 if XF,YF is outside, angle = +/- 2 pi  if XF,YF is inside
!      INSIDE = ABS(ANGLE) .GT. 1.0
!!
!      RETURN
!      END ! INSIDE
!




      SUBROUTINE GETXYF(X,XP,Y,YP,S,N, TOPS,BOTS,XF,YF)
      DIMENSION X(N),XP(N),Y(N),YP(N),S(N)
!
      IF(XF .EQ. -999.0)    &
        CALL ASKR('Enter flap hinge x location^',XF)
!
!---- find top and bottom y at hinge x location
      TOPS = S(1) + (X(1) - XF)
      BOTS = S(N) - (X(N) - XF)
      CALL SINVRT(TOPS,XF,X,XP,S,N)      
      CALL SINVRT(BOTS,XF,X,XP,S,N)      
      TOPY = SEVAL(TOPS,Y,YP,S,N)
      BOTY = SEVAL(BOTS,Y,YP,S,N)
!
      WRITE(*,1000) TOPY, BOTY
 1000 FORMAT(/'  Top    surface:  y =', F8.4,'     y/t = 1.0'    &
             /'  Bottom surface:  y =', F8.4,'     y/t = 0.0')
!
      IF(YF .EQ. -999.0)    &
       CALL ASKR(    &
        'Enter flap hinge y location (or 999 to specify y/t)^',YF)
!
      IF(YF .EQ. 999.0) THEN
        CALL ASKR('Enter flap hinge relative y/t location^',YREL)
        YF = TOPY*YREL + BOTY*(1.0-YREL)
      ENDIF
!
      RETURN
      END ! GETXYF






!
!
!
!      SUBROUTINE PLOTG
!!--------------------------------------------------------------
!!     Plots buffer airfoil with ticked chord line or grid
!!--------------------------------------------------------------
!      INCLUDE 'XFOIL.INC'
!!
!      DATA LMASK1, LMASK2, LMASK3 / -32640, -30584, -21846 /
!      INCLUDE 'XDES.INC'
!!
!!---- node tick mark size and corner symbol size
!      DTICK = GTICK*(SB(NB)-SB(1))
!      SSH   = DTICK * 3.0
!!
!      CALL NCALC(XB,YB,SB,NB,W1,W2)
!!
!      IF(LGGRID) THEN
!        CALL GRDAIR(XGMIN,XGMAX,YGMIN,YGMAX,DXYG,DXYG,CHG,.TRUE.,.TRUE.,
!     &              XOFF,XSF,YOFF,YSF, LMASK2)
!        XL0 = XMOD(XGMIN)
!        YL0 = YMOD(YGMAX) + 2.0*CH
!      ELSE
!!------ plot chord line and tick marks every 10% chord
!        CALL NEWPEN(1)
!        CALL PLOT(XMOD(0.0),YMOD(0.0),3)
!        CALL PLOT(XMOD(1.0),YMOD(0.0),2)
!        DO ITICK=1, 10
!          XPLT = FLOAT(ITICK)/10.0
!          CALL PLOT(XMOD(XPLT),YMOD(0.003),3)
!          CALL PLOT(XMOD(XPLT),YMOD(-.003),2)
!        ENDDO
!!
!        XL0 = XMOD(XBMIN)
!        YL0 = YMOD(YBMAX) + 2.0*CH
!      ENDIF
!      IF(LPLCAM)  YL0 = YSF*(YCMAX-DYOFFC-YOFF) + 2.0*CH
!!
!      CALL PLFLUSH
!!
!      CALL NEWPEN(2)
!      CALL PLTAIR(XB,XBP,YB,YBP,SB,NB, XOFF,XSF,YOFF,YSF,'black')
!!
!      IF(LGTICK) THEN
!!----- draw tiny tick mark normal to airfoil surface at each panel node
!       DO I=2, NB-1
!         CALL PLOT(XMOD(XB(I)            ),YMOD(YB(I)            ),3)
!         CALL PLOT(XMOD(XB(I)-DTICK*W1(I)),YMOD(YB(I)-DTICK*W2(I)),2)
!       ENDDO
!      ENDIF
!!
!!C---- plot symbol at nose
!!      CALL NSFIND(STLE,XB,XBP,YB,YBP,SB,NB)
!!      XT = SEVAL(STLE,XB,XBP,SB,NB)
!!      YT = SEVAL(STLE,YB,YBP,SB,NB)
!!      CALL PLSYMB(XMOD(XT),YMOD(YT),0.005*XSF,5,0.0,0)
!!
!!---- put symbol at any doubled point
!      DO I=1, NB-1
!        IF(SB(I) .EQ. SB(I+1))
!     &     CALL PLSYMB(XMOD(XB(I)),YMOD(YB(I)),SSH,5,0.0,0)
!      ENDDO
!!
!      IF(LPLCAM) THEN
!       CALL PLTCAM(' ')
!      ENDIF
!!
!      IF(LGPARM) THEN
!       CALL NEWPEN(3)
!       CALL GPARPL(XL0,YL0,0.7*CH,.TRUE.,NAME,
!     &             CHORDB,AREAB,RADBLE,ANGBTE,
!     &             EI11BA,EI22BA,APX1BA,APX2BA,
!     &             EI11BT,EI22BT,APX1BT,APX2BT,
!     &             THICKB,CAMBRB)
!      ENDIF
!!
!      CALL PLFLUSH
!!
!      LGEOPL = .TRUE.
!      NOVER = 0
!!
!      RETURN
!      END ! PLOTG
!
!
!
!      SUBROUTINE PLTCAM(COLIN)
!!--------------------------------------------
!!     Plots camber & thickness distributions
!!--------------------------------------------
!      INCLUDE 'XFOIL.INC'
!      CHARACTER*(*) COLIN
!      CHARACTER*32 COLC, COLT
!      DATA LMASK1, LMASK2, LMASK3 / -32640, -30584, -21846 /
!!
!!---- plot camber/thickness only if camber/tickness plot is being shown
!      IF(.NOT.LPLCAM) RETURN
!!
!      CALL NEWPEN(1)
!      CALL GRDAIR(XGMIN,XGMAX,YCMIN,YCMAX,DXYG,DXYG,CHG,.FALSE.,.TRUE.,
!     &            XOFF,XSF,DYOFFC+YOFF,YSF, LMASK2)
!!     
!      CALL GETCAM(XCM,YCM,NCM,XTK,YTK,NTK,
!     &            XB,XBP,YB,YBP,SB,NB )
!      CALL SCALC(XCM,YCM,SCM,NCM)
!      CALL SEGSPL(XCM,XCMP,SCM,NCM)
!      CALL SEGSPL(YCM,YCMP,SCM,NCM)
!      CALL SCALC(XTK,YTK,STK,NTK)
!      CALL SEGSPL(XTK,XTKP,STK,NTK)
!      CALL SEGSPL(YTK,YTKP,STK,NTK)
!!
!      IF(COLIN(1:1) .EQ. ' ') THEN
!       COLC = 'green'
!       COLT = 'cyan'
!      ELSE
!       COLC = COLIN
!       COLT = COLIN
!      ENDIF
!!
!      CALL NEWPEN(2)
!      YOF = YOFF + DYOFFC
!      CALL PLTAIR(XTK,XTKP,YTK,YTKP,STK,NTK,XOFF,XSF, YOF, YSF,COLT)
!      CALL PLTAIR(XTK,XTKP,YTK,YTKP,STK,NTK,XOFF,XSF,-YOF,-YSF,COLT)
!!--- Offset for camber includes offset for LE camber point
!      YOFF1C = YOFF + DYOFFC + YCM(1)
!      CALL PLTAIR(XCM,XCMP,YCM,YCMP,SCM,NCM,XOFF,XSF, YOFF1C,YSF,COLC)
!!
!      RETURN
!      END ! PLTCAM
!
!
!      SUBROUTINE PLNEWP(COLOR)
!      INCLUDE 'XFOIL.INC'
!      CHARACTER*(*) COLOR
!!
!      LOGICAL LCOLOR
!      INCLUDE 'XDES.INC'
!!
!!---- don't plot geometric parameters if camber/tickness plot is being shown
!      IF(LPLCAM) RETURN
!!
!      LCOLOR = COLOR(1:1) .NE. ' '
!!
!      IF(LCOLOR) THEN
!        CALL GETCOLOR(ICOL0)
!        CALL NEWCOLORNAME(COLOR)
!      ENDIF
!!
!      CALL NEWPEN(3)
!!
!      NOVER = NOVER + 1
!      IF(LGGRID) THEN
!       XL0 = XMOD(XGMIN) +  2.0*CH + 9.0*CH*FLOAT(NOVER)
!       YL0 = YMOD(YGMAX) +  2.0*CH
!      ELSE
!       XL0 = XMOD(XBMIN) +  2.0*CH + 9.0*CH*FLOAT(NOVER)
!       YL0 = YMOD(YBMAX) +  2.0*CH
!      ENDIF
!      
!      IF(LPLCAM)  YL0 = YSF*(YCMAX-YOFF-DYOFFC) + 2.0*CH
!!
!      IF(LGPARM) THEN
!        CALL GPARPL(XL0,YL0,0.7*CH,.FALSE.,NAME,
!     &              CHORDB,AREAB,RADBLE,ANGBTE,
!     &              EI11BA,EI22BA,APX1BA,APX2BA,
!     &              EI11BT,EI22BT,APX1BT,APX2BT,
!     &              THICKB,CAMBRB)
!      ENDIF
!!
!      IF(LCOLOR) CALL NEWCOLOR(ICOL0)
!      CALL PLFLUSH
!!
!      RETURN
!      END ! PLNEWP
!
!
!
!      SUBROUTINE GPARPL(X0,Y0,CH, LABEL, NAME,
!     &                  CHORD,AREA,RADLE,ANGTE,
!     &                  EI11A,EI22A,APX1A,APX2A,
!     &                  EI11T,EI22T,APX1T,APX2T,
!     &                  THICK,CAMBR)
!      LOGICAL LABEL
!      EXTERNAL PLCHAR
!      CHARACTER NAME*(*)
!!
!      RTD = 45.0/ATAN(1.0)
!!
!      XSPACE = 30.0*CH
!      YSPACE =  2.0*CH
!!
!      X = X0
!      Y = Y0
!!
!      IF(LABEL) THEN
!       CALL PLCHAR(X,Y,CH,'       = ',0.0, 9)
!       CALL PLMATH(X,Y,CH,'  Oq     ',0.0, 9)
!       CALL PLSUBS(X+3.0*CH,Y,CH,'TE',0.0, 2, PLCHAR)
!      ENDIF
!      CALL PLNUMB(X+9.0*CH,Y,CH,ANGTE*RTD        ,0.0, 2)
!      CALL PLMATH(999.,Y,CH,'"'              ,0.0, 1)
!      Y = Y + YSPACE
!!
!      IF(LABEL) THEN
!       CALL PLCHAR(X,Y,CH,'   r   = ',0.0, 9)
!       CALL PLSUBS(X+3.0*CH,Y,CH,'LE',0.0, 2, PLCHAR)
!      ENDIF
!      CALL PLNUMB(X+9.0*CH,Y,CH,RADLE,0.0, 5)
!      Y = Y + YSPACE
!!
!      IF(LABEL) THEN
!       CALL PLCHAR(X,Y,CH,'camber = ',0.0, 9)
!      ENDIF
!      CALL PLNUMB(X+9.0*CH,Y,CH,CAMBR,0.0, 5)
!      Y = Y + YSPACE
!!
!      IF(LABEL) THEN
!       CALL PLCHAR(X,Y,CH,'thick. = ',0.0, 9)
!      ENDIF
!      CALL PLNUMB(X+9.0*CH,Y,CH,THICK,0.0, 5)
!      Y = Y + YSPACE
!!
!      IF(LABEL) THEN
!       CALL PLCHAR(X,Y,CH,' area  = ',0.0, 9)
!      ENDIF
!      CALL PLNUMB(X+9.0*CH,Y,CH, AREA,0.0, 5)
!      Y = Y + YSPACE
!!
!!
!!      X = X0  +  XSPACE
!!      Y = Y0
!!C
!!      Y = Y + YSPACE
!!C
!!      CALL PLMATH(X,Y,1.4*CH,'I',0.0,1)
!!      CALL PLMATH(X,Y,CH,'       2     ',0.0,-1)
!!      CALL PLCHAR(X,Y,CH,' (y-y ) ds = ',0.0,-1)
!!      CALL PLNUMB(999.,Y,CH, 1000.0*EI11T,0.0,4)
!!      CALL PLMATH(999.,Y,CH,'#'   ,0.0,1)
!!      CALL PLCHAR(999.,Y,CH, '10' ,0.0,2)
!!      CALL PLMATH(999.,Y,CH,   '3',0.0,1)
!!      CALL PLSUBS(X+4.0*CH,Y,CH,'o',0.0,1,PLCHAR)
!!      Y = Y + YSPACE
!!C
!!      CALL PLMATH(X,Y,1.4*CH,'I',0.0,1)
!!      CALL PLMATH(X,Y,CH,'       2     ',0.0,-1)
!!      CALL PLCHAR(X,Y,CH,' (y-y ) dA = ',0.0,-1)
!!      CALL PLNUMB(999.,Y,CH, 1000.0*EI11A,0.0,4)
!!      CALL PLMATH(999.,Y,CH,'#'   ,0.0,1)
!!      CALL PLCHAR(999.,Y,CH, '10' ,0.0,2)
!!      CALL PLMATH(999.,Y,CH,   '3',0.0,1)
!!      CALL PLSUBS(X+4.0*CH,Y,CH,'o',0.0,1,PLCHAR)
!!      Y = Y + YSPACE
!!C
!!      CALL PLMATH(X,Y,CH,'             ',0.0,-1)
!!      CALL PLCHAR(X,Y,CH,'      area = ',0.0,-1)
!!      CALL PLNUMB(999.,Y,CH, AREA,0.0, 5)
!!      Y = Y + YSPACE
!!
!!--- Plot airfoil name over data list
!      CALL PLCHAR(X+9.0*CH,Y,CH,NAME,0.0, 12)
!!
!      RETURN
!      END ! GPARPL
!
!
!
!      SUBROUTINE GRDAIR(XGMIN,XGMAX, YGMIN,YGMAX,DXGN,DYGN,CHG,
!     &                  LXAXIS,LYAXIS,
!     &                  XOFF,XSF,YOFF,YSF, LMASK)
!      LOGICAL LXAXIS,LYAXIS
!!----------------------------------------
!!     Plots grid with axes.
!!     Intended for airfoil plot.
!!----------------------------------------
!      INCLUDE 'XDES.INC'
!!
!      CALL NEWPEN(1)
!!
!!---- plot outline
!      CALL PLOT(XMOD(XGMIN),YMOD(YGMIN),3)
!      CALL PLOT(XMOD(XGMAX),YMOD(YGMIN),2)
!      CALL PLOT(XMOD(XGMAX),YMOD(YGMAX),2)
!      CALL PLOT(XMOD(XGMIN),YMOD(YGMAX),2)
!      CALL PLOT(XMOD(XGMIN),YMOD(YGMIN),2)
!!
!      IF(LXAXIS)
!     &  CALL XAXIS(XMOD(XGMIN),YMOD(YGMIN),(XGMAX-XGMIN)*XSF,
!     &             DXGN*XSF, XGMIN,DXGN,CHG,-2)
!      IF(LYAXIS)
!     &  CALL YAXIS(XMOD(XGMIN),YMOD(YGMIN),(YGMAX-YGMIN)*YSF,
!     &             DYGN*YSF, YGMIN,DYGN,CHG,-2)
!!
!!---- fine grid
!      NXG = INT((XGMAX-XGMIN)/DXGN + 0.1)
!      NYG = INT((YGMAX-YGMIN)/DYGN + 0.1)
!      NXG = MAX(1,NXG)
!      NYG = MAX(1,NYG)
!!
!      X0 = XMOD(XGMIN)
!      Y0 = YMOD(YGMIN)
!      DXG = (XMOD(XGMAX)-X0)/NXG
!      DYG = (YMOD(YGMAX)-Y0)/NYG
!      CALL PLGRID(X0,Y0,NXG,DXG,NYG,DYG, LMASK)
!!
!      RETURN
!      END ! GRDAIR
!
!
!
!      SUBROUTINE PLTAIR(XX,XXP,YY,YYP,SS,NN, XOFF,XSF,YOFF,YSF,COLOR)
!      DIMENSION XX(NN), XXP(NN), YY(NN), YYP(NN), SS(NN)
!      CHARACTER*(*) COLOR
!!-----------------------------
!!     Plots passed-in airfoil
!!-----------------------------
!      LOGICAL LCOLOR
!      XMOD(XTMP) = XSF * (XTMP - XOFF)
!      YMOD(YTMP) = YSF * (YTMP - YOFF)
!!
!      NT = 20
!!cc      NT = 50
!!
!      LCOLOR = COLOR(1:1) .NE. ' '
!!
!      IF(LCOLOR) THEN
!        CALL GETCOLOR(ICOL0)
!        CALL NEWCOLORNAME(COLOR)
!      ENDIF
!!
!      DO 60 I=2, NN
!        DS = SS(I) - SS(I-1)
!        CALL PLOT(XMOD(XX(I-1)),YMOD(YY(I-1)),3)
!!
!!------ subdivide current panel into NT segments for smoother airfoil plot
!        DO 610 IT=1, NT
!          ST = SS(I-1) + DS*FLOAT(IT)/FLOAT(NT)
!          XT = SEVAL(ST,XX,XXP,SS,NN)
!          YT = SEVAL(ST,YY,YYP,SS,NN)
!          CALL PLOT(XMOD(XT),YMOD(YT),2)
!  610   CONTINUE
!   60 CONTINUE
!!
!      IF(LCOLOR) CALL NEWCOLOR(ICOL0)
!!
!      CALL PLFLUSH
!!
!      RETURN
!      END ! PLTAIR
!
!
!
!      SUBROUTINE OVER(FNAME1)
!!----------------------------------------------------
!!     Overlays plot of airfoil from coordinate file.
!!----------------------------------------------------
!      INCLUDE 'XFOIL.INC'
!      CHARACTER*(*) FNAME1
!!
!      CHARACTER*32 NAME0, NAMEW
!      CHARACTER*80 ISPARS0
!!
!      IF(FNAME1(1:1).NE.' ') THEN
!       FNAME = FNAME1
!      ELSE
!!----- no argument... get it somehow
!       IF(ONAME(1:1).NE.' ') THEN
!!------ offer existing default
!        WRITE(*,1100) ONAME
! 1100   FORMAT(/' Enter filename:  ', A)
!        READ(*,1000) FNAME
! 1000   FORMAT(A)
!        CALL STRIP(FNAME,NFN)
!        IF(NFN.EQ.0) FNAME = ONAME
!       ELSE
!!------ just ask for filename
!        CALL ASKS('Enter filename^',FNAME)
!       ENDIF
!      ENDIF
!!
!      LU = 9
!      CALL AREAD(LU,FNAME,2*IQX,W1,W2,NN,NAME0,ISPARS0,ITYPE,1)
!      IF(ITYPE.EQ.0) RETURN
!!
!!---- set new default filename
!      ONAME = FNAME
!!
!      IF(LNORM) THEN
!!----- normalize to unit chord
!       CALL NORM(W1,W3,W2,W4,W5,NN)
!      ELSE
!       CALL SCALC(W1,W2,W5,NN)
!       CALL SEGSPL(W1,W3,W5,NN)
!       CALL SEGSPL(W2,W4,W5,NN)
!      ENDIF
!!
!      NAMEW  = NAME
!      SWLE   = SBLE  
!      CHORDW = CHORDB
!      AREAW  = AREAB 
!      RADWLE = RADBLE
!      ANGWTE = ANGBTE
!      EI11WA = EI11BA
!      EI22WA = EI22BA
!      APX1WA = APX1BA
!      APX2WA = APX2BA
!      EI11WT = EI11BT
!      EI22WT = EI22BT
!      APX1WT = APX1BT
!      APX2WT = APX2BT
!      THICKW = THICKB
!      CAMBRW = CAMBRB
!!
!      NAME = NAME0
!      CALL GEOPAR(W1,W3,W2,W4,W5,NN,W6,
!     &            SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
!     &            EI11BA,EI22BA,APX1BA,APX2BA,
!     &            EI11BT,EI22BT,APX1BT,APX2BT,
!     &            THICKB,CAMBRB )
!!
!      IF(.NOT.LPLOT) THEN
!       CALL PLTINI
!!cc       CALL PLOT(0.05,0.30,-3)
!      ENDIF
!!
!!
!      CALL GETCOLOR(ICOL0)
!      ICOL = 3 + MOD(NOVER,6)
!      IF(ICOL .GE. 5) ICOL = ICOL + 1
!      CALL NEWPEN(2)
!      CALL NEWCOLOR(ICOL)
!!
!      CALL PLTAIR(W1,W3,W2,W4,W5,NN, XOFF,XSF, YOFF,YSF,' ')
!      CALL PLNEWP(' ')
!!
!      CALL NEWCOLOR(ICOL0)
!!
!!
!!---- restore parameters
!      NAME   = NAMEW
!      SBLE   = SWLE  
!      CHORDB = CHORDW
!      AREAB  = AREAW 
!      RADBLE = RADWLE
!      ANGBTE = ANGWTE
!      EI11BA = EI11WA
!      EI22BA = EI22WA
!      APX1BA = APX1WA
!      APX2BA = APX2WA
!      EI11BT = EI11WT
!      EI22BT = EI22WT
!      APX1BT = APX1WT
!      APX2BT = APX2WT
!      THICKB = THICKW
!      CAMBRB = CAMBRW
!!
!      RETURN
!      END ! OVER
!
!

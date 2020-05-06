      SUBROUTINE ARCSCA(IDO,JDER,JEST,JSCA,NBJ,NEL,nl,NPL,NPNE,NVJL,
     '  DL,TOL,XP,ERROR,*)

C#### Subroutine: ARCSCA
C###  Description:
C###    ARCSCA calculates arc length DL(3) of global line segment nl,
C###    using Gaussian quadrature.

C**** Newton iteration is used to calculate arc length.
C**** TOL sets convergence tolerance.
C**** JDER=0,1 as nodal derivs are not,are updated, respec.
C**** JEST=0,1 as initial estimates of DL(3) are not,are given,respec
C**** JSCA=0 : values of global arc-length derivs DL(1) & DL(2)
C**** are not given. 
C**** Global derivs in XP are assumed to be wrt arc-length if JSCA=0.
C**** Note: SUM1 is square of deriv of arclength wrt Xi
C****       SUM2 is 0.5* second deriv of arclength wrt Xi
C****       SUM3 is arclength estimate
C****       SUM4 is used is Newton method

        SUM2=0.0d0
        DO ng=1,NGA
          XI=XIGG(IG(NGA)+ng)
          W=WG_LOCAL(IG(NGA)+ng)
          DO nj=1,3
            DO k=1,2
              XA_LOCAL(k,nj)=0.0d0
              IF(nj.LE.NJT) XA_LOCAL(k,nj)=PL1(1,k,XI)*XN_LOCAL(1,nj,1)
     '                              +PL1(2,k,XI)*XN_LOCAL(1,nj,2)
            ENDDO
          ENDDO
          IF(ITYP10(1).EQ.1) THEN
            SUM1=XA_LOCAL(2,1)**2+XA_LOCAL(2,2)**2+XA_LOCAL(2,3)**2
          ENDIF
          SUM2=SUM2+W*DSQRT(SUM1)
        ENDDO
        DL(1,nl)=SUM2
        DL(2,nl)=SUM2
        DL(3,nl)=SUM2

      DO it=1,ITMAX
        IT_count=it
        SUM3=0.0d0
        SUM4=0.0d0
        DO ng=1,NGA
          XI=XIGG(IG(NGA)+ng)
          W=WG_LOCAL(IG(NGA)+ng)
          DO nj=1,NJT
            IF(NPL(1,nj).EQ.1) THEN
            ELSE IF(NPL(1,nj).EQ.4) THEN
              DO k=1,2
                XA_LOCAL(k,nj)=0.0d0
                DO n=1,2
                  XA_LOCAL(k,nj)=XA_LOCAL(k,nj)+
     '                              PH3(n,1,k,XI)*XN_LOCAL(1,nj,n)
     '                             +PH3(n,2,k,XI)*
     '                              XN_LOCAL(2,nj,n)*DL(n,nl)
                ENDDO
              ENDDO
              XA_LOCAL(3,nj)=0.0d0
C             second deriv wrt xi and total arc length
              DO n=1,2
                XA_LOCAL(3,nj)=XA_LOCAL(3,nj)+
     '            PH3(n,2,2,XI)*XN_LOCAL(2,nj,n)
              ENDDO
C             also need deriv wrt total arc length in polar coord systems
              XA_LOCAL(4,nj)=0.0d0
              DO n=1,2
                XA_LOCAL(4,nj)=XA_LOCAL(4,nj)+
     '            PH3(n,2,1,XI)*XN_LOCAL(2,nj,n)
              ENDDO
            ENDIF
          ENDDO
          IF(ITYP10(1).EQ.1) THEN
            SUM1=XA_LOCAL(2,1)**2+XA_LOCAL(2,2)**2+XA_LOCAL(2,3)**2
            SUM2=0.0d0
            DO nj=1,NJT
              SUM2=SUM2+XA_LOCAL(2,nj)*XA_LOCAL(3,nj)
            ENDDO !nj
          ELSE IF(ITYP10(1).EQ.2) THEN
            SUM1=XA_LOCAL(2,1)**2+(XA_LOCAL(2,2)*XA_LOCAL(1,1))**2
            SUM2=XA_LOCAL(2,1)*XA_LOCAL(3,1)+XA_LOCAL(2,2)*
     '           XA_LOCAL(3,2)*XA_LOCAL(1,1)**2
     '          +XA_LOCAL(1,1)*XA_LOCAL(4,1)*XA_LOCAL(2,2)**2
            IF(NJT.EQ.3) THEN
              SUM1=SUM1+XA_LOCAL(2,3)**2
              SUM2=SUM2+XA_LOCAL(2,3)*XA_LOCAL(3,3)
            ENDIF
          ENDIF
          SUM3=SUM3+W*DSQRT(SUM1)
          IF(SUM1.GT.1.0d-6) SUM4=SUM4+W*SUM2/DSQRT(SUM1)
        ENDDO !ng

        IF(JSCA.EQ.0) THEN
          !arc-length derivs DL(1) & DL(2) not given
          !& must be calculated with Newton
          !Note: see notes in FE02 folder for maths on this
          DA=-(DL(3,nl)-SUM3)/(1.0d0-SUM4)
          IF(DABS(DA).GT.1.0d6) THEN
            WRITE(OP_STRING,'('' Length of line nl='',I4,'
     '        //''' has not converged & is set to unity'')') nl
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            DL(3,nl)=1.0d0
            GOTO 80
          ENDIF
          DL(3,nl)=DL(3,nl)+DA !is new arclength
          IF(JDER.EQ.1) THEN !update nodal derivatives
            DO n=1,2
              IF(DABS(XN_LOCAL(2,1,n)).LT.1.0d0) THEN
                XN_LOCAL(2,2,n)=DSQRT(1.0d0-XN_LOCAL(2,1,n)**2)
              ELSE
                XN_LOCAL(2,1,n)=1.0d0
                XN_LOCAL(2,2,n)=0.0d0
              ENDIF
! Note: XN_LOCAL is updated so as to keep dX_dXi fixed at its previous value
              DO nj=1,NJT
                XN_LOCAL(2,nj,n)=XN_LOCAL(2,nj,n)*DL(n,nl)/DL(3,nl) !new node deriv
              ENDDO                                     !wrt s
            ENDDO
          ENDIF
          DL(1,nl)=DL(3,nl)
          DL(2,nl)=DL(3,nl)
          IF(DABS(DA).LE.DL(3,nl)*TOL) GOTO 70
        ENDIF
      ENDDO !iteration

 80   CONTINUE

      IF(JTYP2B.EQ.1.AND.NPL(4,0).GT.0) THEN
        DL(1,NPL(4,0))=-1.0d0*DL(2,nl)
        DL(2,NPL(4,0))=-1.0d0*DL(1,nl)
        DL(3,NPL(4,0))=DL(3,nl)
      ENDIF

 9998 CALL EXITS('ARCSCA')
      RETURN
 9999 CALL ERRORS('ARCSCA',ERROR)
      CALL EXITS('ARCSCA')
      RETURN 1
      END



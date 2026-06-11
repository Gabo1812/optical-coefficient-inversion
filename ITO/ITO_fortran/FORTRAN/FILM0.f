       PROGRAM NANDK
C      ----------------------------------------------------------------
C      William E. Vargas, CICIMA AND PHYSICS DEPARTMENT, UNIVERSITY OF
C      COSTA RICA.
C      ----------------------------------------------------------------
       IMPLICIT REAL*8 (A-H,O-Z)
       PARAMETER (NOBS=601,NOBS2=2*NOBS)
       PARAMETER (NO=4)
       COMMON /OPTICS0/ HF,HS
       COMMON /OPTICS1/ WL(NOBS),NS(NOBS),KS(NOBS),T(NOBS)
       COMMON /OPTICS2/ TEXP(NOBS)
       COMMON /OPTICS3/ TF(NOBS),MNF(NOBS,NOBS),MKF(NOBS,NOBS)
       COMMON /OPTICS4/ FS(0:NO),G(0:NO),W(1:NO)
       COMMON /OPTICS5/ S(1:NO)
       REAL*8 WL,TEXP,NF(NOBS),KF(NOBS),NS,KS
       REAL*8 X(NOBS2+1),T,F,L(NOBS2+1),U(NOBS2+1),X0(NOBS2+1)
       REAL*8 MINN,MAXN,MINK,MAXK,FS,G,W
       REAL*8 TF,MNF,MKF,NEWHF,NEWNF(NOBS),NEWKF(NOBS),DHF
       COMPLEX*16 EPSB(NOBS),EPSF(NOBS)
       INTEGER FLAG
       PI=4.0D0*DATAN(1.0D0)
       HC=1239.83D0
C      --------------------------------------------------------
C      INPUT FILES: FILE 10 CONTAINS THE EXPERIMENTAL DATA (THE
C      THICKNESS OF THE FILM AND THREE COLUMNS CORRESPONDING TO
C      THE WAVELENGTH (IN THE SAME UNITS OF THE THICKNESS), AND
C      TRANSMITTANCE AND REFLECTANCE VALUES AT EACH WAVELENGTH.
C      FILE 20 HAS  THE THICKNESS AND OPTICAL  CONSTANTS OF THE
C      SUBSTRATE. FILE 30 CONTAINS THE OPTICAL CONSTANTS OF THE
C      FILM. THESE VALUES ARE USUALLY TAKEN FROM LITERATURE AND
C      THEY ARE USED AS THE FIRST APPROXIMATION IN THE SPECTRAL
C      PROJECTED GRADIENT METHOD.
C      --------------------------------------------------------
C      THICKNESS (nm) OF THE ITO-FILM
       HF=140.0D0  ! 140 nm para ITO
C      ------------------------------------------------------------
C      OUTPUT FILES:  FILE 40 [50] CONTAINS THE FIRST APPROXIMATION
C      AND THE RETRIEVED VALUES OF THE REFRACTIVE INDEX [EXTINCTION
C      COEFFICIENT] AT EACH WAVELENGTH. FILE 60 CONTAINS THE  FIRST
C      APPROXIMATION AND RETRIEVED VALUES OF THE DYNAMIC CONDUCTI-
C      VITY OF THE FILM IN  TERMS OF  ENERGY IN eV. FILE 70 HAS THE 
C      EXPERIMENTAL AND  RETRIEVED VALUES OF  TRANSMITTANCE AND RE-
C      FLECTANCE.
C      ------------------------------------------------------------
       OPEN(10,FILE='Dy_T')
       OPEN(20,FILE='DATANK0_DY')
       OPEN(40,FILE='DATAN')
       OPEN(50,FILE='DATAK') 
       OPEN(60,FILE='DATAT')
       OPEN(66,FILE='DATANK1_DY')
C      -----------------------------------------------------------------
C      THE PACKING DENSITY P IS CALCULATED IF IP=0
       PF=1.0D0
       P0=PF
       IP=1
       P=1.0D0
       IF(IP.EQ.1) THEN
          P=PF
       ELSE
          P0=P
       ENDIF
C      ------------------------------------------------------
C      FILM THICKNESS (HF nm) AND SUBSTRATE THICKNESS (HS nm)
C      ------------------------------------------------------
       HS=1.0D6
       HF=HF/P
C      -------------------------------------------------------
C      HF IS THE FILM THICKNESS, HS IS THE SUBSTRATE THICKNESS
C      -------------------------------------------------------
       WRITE(*,*) 'INPUT VALUES OF FILM AND SUBSTRATE THICKNESSES:'
       WRITE(*,*)
       WRITE(*,*) 'THICKNESS OF THE FILM:',HF
       WRITE(*,*) 'THICKNESS OF THE SUBS:',HS
       WRITE(*,*)
       WRITE(*,*) 'CALCULATING...'
       WRITE(*,*)
       DEPS=0.01D0
       SIGNO=0.0D0
       DO I=1,NOBS
          READ(10,*) WL(I),TEXP(I)
          TEXP(I)=TEXP(I)
C         TEXP(I)=TEXP(I)/100.0D0
          HW=HC/WL(I)
          WLMC=WL(I)/1000.0D0
          READ(20,*) WL(I),NF(I),KF(I) 
          CALL SODA_LIME_GLASS(WLMC,NS(I),KS(I))
          EPS1=NF(I)**2-KF(I)**2
          EPS2=2.0D0*NF(I)*KF(I)
          EPSB(I)=DCMPLX(EPS1,EPS2)
       ENDDO
       CLOSE(10)
       CLOSE(20)
       CLOSE(30)
       
       DO I=1,NOBS
          X(I)=NF(I)
          X(I+NOBS)=KF(I)
          X0(I)=X(I)
          X0(I+NOBS)=X(I+NOBS)
       ENDDO
       X(NOBS2+1)=HF
       X0(NOBS2+1)=HF
 
       MINN=+1.0D3
       MAXN=-1.0D3
       MINK=+1.0D3
       MAXK=-1.0D3 
       DO I=1,NOBS
          IF(NF(I).LT.MINN) MINN=NF(I)
          IF(NF(I).GT.MAXN) MAXN=NF(I)
          IF(KF(I).LT.MINK) MINK=KF(I)
          IF(KF(I).GT.MAXK) MAXK=KF(I)
       ENDDO
 
       MINN=0.10D0*MINN
       MAXN=10.0D0*MAXN
       MINK=0.10D0*MINK
       MAXK=10.0D0*MAXK         
       DO I=1,NOBS
          L(I)=MINN
          U(I)=MAXN
       ENDDO 
       DO I=NOBS+1,NOBS2
          L(I)=MINK
          U(I)=MAXK
       ENDDO
       L(NOBS2+1)=0.9995D0*HF
       U(NOBS2+1)=1.0005D0*HF     
C      --------------------------------------------------------        
       M=10       
       EPS=1.0D-11
       MAXIT=5000
       MAXIFCNT=20000
C      ---------------------------------------------------------
       NITER=0
       FMAX=-1.0D4 
       DO 11 IWL=1,NOBS
          DO I=1,NOBS
             X(I)=NF(I)
             X(I+NOBS)=KF(I)
          ENDDO
          X(NOBS2+1)=HF     
          CALL SPG(X,L,U,M,EPS,MAXIT,MAXIFCNT,F,CGNORM,ITER,IFCNT,
     &             IGCNT,FLAG,IWL)
C
          NITER=NITER+ITER
          TF(IWL)=X(NOBS2+1)
          IF(F.GT.FMAX) FMAX=F
          DO I=1,NOBS
             MNF(I,IWL)=X(I)
             MKF(I,IWL)=X(I+NOBS)
          ENDDO
   11  CONTINUE
       NITER=NITER/NOBS
C      ----------------------------------------------------------
       CALL AVERAGES(NEWHF,NEWNF,NEWKF,DHF,DELTANF,DELTAKF)
       WRITE(*,*) 'AVERAGE THICKNESS:',NEWHF,' +-',DHF,'nm'
       WRITE(*,*) 'LARGEST DELTA N  :',DELTANF
       WRITE(*,*) 'LARGEST DELTA K  :',DELTAKF
       WRITE(*,*) 'LARGEST FUNCTION :',FMAX
       WRITE(*,*) 'SET GRADIENT NORM:',EPS      
       WRITE(*,*) 'ITERATIONS:',NITER
       
       DO I=1,NOBS
          X(I)=NEWNF(I)
          X(I+NOBS)=NEWKF(I)
          EPS1=NEWNF(I)**2-NEWKF(I)**2
          EPS2=2.0D0*NEWNF(I)*NEWKF(I)
          EPSF(I)=DCMPLX(EPS1,EPS2)
       ENDDO
       X(NOBS2+1)=NEWHF
       
       CALL EVALF(X,F)
             
       DO I=1,NOBS
         J=I+NOBS
         ENR=HC/WL(I)
         WRITE(40,100) WL(I),NF(I),NEWNF(I)
         WRITE(50,100) WL(I),KF(I),NEWKF(I) 
         WRITE(60,100) WL(I),100.0D0*TEXP(I),100.0D0*T(I)
         WRITE(66,100) WL(I),NEWNF(I),NEWKF(I)
         EPS1F=X(I)**2-X(J)**2
         EPS2F=2.0D0*X(I)*X(J)
         EPS1F0=NF(I)**2-KF(I)**2
         EPS2F0=2.0D0*NF(I)*KF(I)
         OPEN(70,FILE='EPSILON12')
         IF(P.NE.1.0D0) THEN
            CALL MG(P,EPS1F,EPS2F,EPS1P,EPS2P)
            EELFF=EPS2F/(EPS1F**2+EPS2F**2)
            EELFP=EPS2P/(EPS1P**2+EPS2P**2)
            WRITE(70,100) ENR,EPS1F,EPS1P,EPS2F,EPS2P,EELFF,EELFP
         ELSE
            EELFF=EPS2F/(EPS1F**2+EPS2F**2)
            WRITE(70,100) ENR,EPS1F0,EPS1F,EPS2F0,EPS2F
         ENDIF
       ENDDO
C      -----------------------------------------------------------------
C      EVALUATION OF THE REFLECTION COEFFICIENT AT THE INTERFACE FILM-
C      SUBSTRATE: Dy-QUARTZ
C      -----------------------------------------------------------------
       OPEN(80,FILE='RIJ')
       OPEN(90,FILE='TIJ')
       RN1=1.0D0
       RK1=0.0D0
       RN4=RN1
       RK4=RK1
       H2=HF
       H3=HS       
       DO I=1,NOBS
          RN2=NEWNF(I)
          RK2=NEWKF(I)
          RN3=NS(I)
          RK3=KS(I)
          CALL SSRIJ(RN1,RK1,RN2,RK2,R12)
          CALL SSRIJ(RN2,RK2,RN3,RK3,R23)
          CALL SSRIJ(RN3,RK3,RN4,RK4,R34)
          CR34=R34
          TAU2=4.0D0*PI*RK2*H2/WL(I)
          TAU3=4.0D0*PI*RK3*H3/WL(I)
          CALL SCRIJ(R23,TAU3,CR34,CR23)
          CALL SCRIJ(R12,TAU2,CR23,CR12)   
          CT12=1.0D0-R12
          CALL SCTJK(CT12,R12,TAU2,CR23,CT23)
          CALL SCTJK(CT23,R23,TAU3,CR34,CT34)
 
          WRITE(80,100) WL(I),R12,R23,R34,CR23,CR12
          WRITE(90,100) WL(I),CT12,CT23,CT34
       ENDDO
C      -----------------------------------------------------------------
       CLOSE(40)
       CLOSE(50)
       CLOSE(60)
       CLOSE(66)
       CLOSE(70)
       CLOSE(80)
       CLOSE(90)
C      -----------------------------------------------------------------
C      DETERMINATION OF THE PACKING DENSITY P
       IF(IP.EQ.1) THEN
          CALL PACKING(EPSB,EPSF,P,DELTAP)
          WRITE(*,*) 'METAL PACKING DENSITY: ',P0,P,DELTAP
       ENDIF       
C      -----------------------------------------------------------------
  100  FORMAT(7(E10.4E2,2X))
       END
C      -----------------------------------------------------------------
       SUBROUTINE SSRIJ(RNI,RKI,RNJ,RKJ,RIJ)
       IMPLICIT REAL*8 (A-H,O-Z)
       DUMMY1=(RNI-RNJ)**2+(RKI-RKJ)**2
       DUMMY2=(RNI+RNJ)**2+(RKI+RKJ)**2
       RIJ=DUMMY1/DUMMY2
       RETURN
       END
C      -----------------------------------------------------------------
       SUBROUTINE SCRIJ(RIJ,TAUJ,CRJK,CRIJ)
       IMPLICIT REAL*8 (A-H,O-Z)
       RJI=RIJ
       CRIJ=(1.0D0-RIJ)*CRJK*(1.0D0-RJI)*DEXP(-2.0D0*TAUJ)
       CRIJ=CRIJ/(1.0D0-RIJ*CRJK*DEXP(-2.0D0*TAUJ))
       CRIJ=RIJ+CRIJ
       RETURN
       END
C      -----------------------------------------------------------------
       SUBROUTINE SCTJK(CTIJ,RIJ,TAUJ,CRJK,CTJK)
       IMPLICIT REAL*8 (A-H,O-Z)
       CTJK=CTIJ*(1.0D0-RIJ)*DEXP(-TAUJ)
       CTJK=CTJK/(1.0D0-RIJ*CRJK*EXP(-2.0D0*TAUJ))
       RETURN
       END
C      -----------------------------------------------------------------
       SUBROUTINE MG(P,EPS1F,EPS2F,EPS1P,EPS2P)
       IMPLICIT REAL*8 (A-H,O-Z)
       COMPLEX*16 EPSMETAL,EPSAVG,B,C,DUMMY
       EPSAVG=DCMPLX(EPS1F,EPS2F)
       EPSAIR=1.0D0
C      FRACTION OF VOIDS OR EMPTY VOLUME
       F=1.0D0-P
       A=2.0D0*(1.0D0-F)
       B=(1.0D0+2.0D0*F)*EPSAIR-(2.0D0+F)*EPSAVG
       C=-(1.0D0-F)*EPSAVG*EPSAIR
       DUMMY=CDSQRT(B*B-4.0D0*A*C)
       EPSMETAL=(-B-DUMMY)/(2.0D0*A)
       IF(DIMAG(EPSMETAL).LT.0.0D0) EPSMETAL=(-B+DUMMY)/(2.0D0*A)
       EPS1P=DREAL(EPSMETAL)
       EPS2P=DIMAG(EPSMETAL)
       RETURN
       END       
C      -----------------------------------------------------------------
       SUBROUTINE PACKING(EPSB,EPSF,P,DELTAP)
       IMPLICIT REAL*8 (A-H,O-Z)
       PARAMETER (NOBS=601)
       REAL*8 VP(NOBS)
       COMPLEX*16 EPSB(NOBS),EPSF(NOBS),Z
       DO I=1,NOBS
          Z=(EPSF(I)-1.0D0)*(EPSB(I)+2.0D0)
          Z=Z/((EPSF(I)+2.0D0)*(EPSB(I)-1.0D0))
          VP(I)=CDABS(Z)
       ENDDO
       SP=0.0D0
       DO I=1,NOBS
          SP=SP+VP(I)
       ENDDO
       P=SP/DBLE(NOBS)
       DP=0.0D0
       DO I=1,NOBS
          DP=DP+(P-VP(I))**2
       ENDDO
       DELTAP=DP/DBLE(NOBS)
       RETURN
       END
C      -----------------------------------------------------------
       SUBROUTINE AVERAGES(NEWHF,NEWNF,NEWKF,DHF,DELTANF,DELTAKF)
       IMPLICIT REAL*8 (A-H,O-Z)
       PARAMETER (NOBS=601,NOBS2=2*NOBS)	
       COMMON /OPTICS3/ TF(NOBS),MNF(NOBS,NOBS),MKF(NOBS,NOBS)
       REAL*8 TF,MNF,MKF,NEWNF(NOBS),NEWKF(NOBS),NEWHF,DHF
       REAL*8 DNF(NOBS),DKF(NOBS)
       DUMMY0=DSQRT(DBLE(NOBS*(NOBS-1)))
       NEWHF=0.0D0
       DO I=1,NOBS
          NEWHF=NEWHF+TF(I)
       ENDDO
       NEWHF=NEWHF/DBLE(NOBS)
       DO I=1,NOBS
          NEWNF(I)=0.0D0
          NEWKF(I)=0.0D0
          DO J=1,NOBS
             NEWNF(I)=NEWNF(I)+MNF(I,J)
             NEWKF(I)=NEWKF(I)+MKF(I,J)
          ENDDO
          NEWNF(I)=NEWNF(I)/DBLE(NOBS)
          NEWKF(I)=NEWKF(I)/DBLE(NOBS)
       ENDDO
       DHF=0.0D0
       DO I=1,NOBS
          DUMMY1=(TF(I)-NEWHF)**2
          DHF=DHF+DUMMY1
       ENDDO
       DHF=DHF/DUMMY0
       DELTANF=0.0D0
       DELTAKF=0.0D0
       DO I=1,NOBS
          DNF(I)=0.0D0
          DKF(I)=0.0D0
          DO J=1,NOBS
             DUMMY1=(MNF(I,J)-NEWNF(I))**2
             DUMMY2=(MKF(I,J)-NEWKF(I))**2
             DNF(I)=DNF(I)+DUMMY1
             DKF(I)=DKF(I)+DUMMY2
          ENDDO
          DNF(I)=DNF(I)/DUMMY0
          DKF(I)=DKF(I)/DUMMY0
          IF(DNF(I).GT.DELTANF) DELTANF=DNF(I)
          IF(DKF(I).GT.DELTAKF) DELTAKF=DKF(I)
       ENDDO
       RETURN
       END      
C      -----------------------------------------------------------
       SUBROUTINE SPG(X,L,U,M,EPS,MAXIT,MAXIFCNT, 
     & F,CGNORM,ITER,IFCNT,IGCNT,FLAG,IWL)
       IMPLICIT REAL*8 (A-H,O-Z)
       PARAMETER (NOBS=601,NOBS2=2*NOBS)
       COMMON /OPTICS0/ HF,HS
       COMMON /OPTICS1/ WL(NOBS),NS(NOBS),KS(NOBS),T(NOBS)
       COMMON /OPTICS2/ TEXP(NOBS)  
       INTEGER M,MAXIT,MAXIFCNT,ITER,IFCNT,IGCNT,FLAG
       REAL*8 X(NOBS2+1),L(NOBS2+1),U(NOBS2+1),EPS,F,CGNORM
       REAL*8 LMIN,LMAX,WL,TEXP,NS,KS,T,HF,HS
       PARAMETER (LMIN=1.0D-30,LMAX=1.0D+30)
       PARAMETER (MMAX=600)
       REAL*8 G(NOBS2+1),CG(NOBS2+1),XNEW(NOBS2+1),XBEST(NOBS2+1)
       REAL*8 GNEW(NOBS2+1),S(NOBS2+1)
       REAL*8 Y(NOBS2+1),D(NOBS2+1),LASTFVALUES(0:MMAX-1) 
       REAL*8 GTD,STS,STY,FBEST,FNEW,LAMBDA
       ITER=0
       IFCNT=0
       IGCNT=0
       DO I=0,M-1
         LASTFVALUES(I)=-1.0D+99
       ENDDO  
       DO I=1,NOBS2+1
          X(I)=DMIN1(U(I),DMAX1(L(I),X(I)))
          XBEST(I)= X(I)
       ENDDO
       CALL EVALFG(X,F,G,IWL)
       LASTFVALUES(0)=F
       FBEST=F
       IFCNT=IFCNT+1
       IGCNT=IGCNT+1
       CGNORM=0.0D0
       DO I=1,NOBS2+1
         CG(I)=DMIN1(U(I),DMAX1(L(I),X(I)-G(I)))-X(I)
         CGNORM=DMAX1(CGNORM,DABS(CG(I)))
       ENDDO
       IF(CGNORM.NE.0.0D0) THEN
        LAMBDA=1.0D0/CGNORM
       ENDIF
 10    IF (.NOT.CGNORM.LE.EPS) THEN  
C     & .AND..NOT.ITER.GT.MAXIT .AND. .NOT.IFCNT.GT.MAXIFCNT) THEN
          ITER=ITER+1
          GTD=0.0D0
          DO I=1,NOBS2+1
              D(I)=DMIN1(U(I),DMAX1(L(I),X(I)-LAMBDA*G(I)))-X(I)
              GTD=GTD+G(I)*D(I)
          ENDDO
          CALL LINESEARCH(X,F,D,GTD,FNEW,XNEW, 
     &    M,LASTFVALUES,IFCNT,IWL)
          F=FNEW
          LASTFVALUES(MOD(ITER,M))=F
          IF(F.LT.FBEST) THEN
              FBEST=F
              DO I=1,NOBS2+1
                  XBEST(I)= XNEW(I)
              ENDDO
          ENDIF
          CALL EVALG(XNEW,GNEW,IWL)
          IGCNT=IGCNT+1
          STS=0.0D0
          STY=0.0D0
          CGNORM=0.0D0
          DO I=1,NOBS2+1
             S(I)=XNEW(I)-X(I)
             Y(I)=GNEW(I)-G(I)
             STS=STS+S(I)*S(I)
             STY=STY+S(I)*Y(I)
             X(I)=XNEW(I)
             G(I)=GNEW(I)
             CG(I)=DMIN1(U(I),DMAX1(L(I),X(I)-G(I)))-X(I)
             CGNORM=DMAX1(CGNORM,DABS(CG(I)))
          ENDDO
          IF(STY.LE.0.0D0) THEN
             LAMBDA= LMAX
          ELSE
             LAMBDA= DMIN1(LMAX,DMAX1(LMIN,STS/STY))
          ENDIF
          GOTO 10
       ENDIF
       F=FBEST
       DO I=1,NOBS2+1
          X(I)=XBEST(I)
       ENDDO
       IF(CGNORM.LE.EPS) THEN
          FLAG=0
       ELSE IF(ITER.GT.MAXIT) THEN
          FLAG=1
       ELSE IF(IFCNT.GT.MAXIFCNT) THEN
          FLAG=2
       ELSE
          FLAG=9
       ENDIF
       RETURN
       END
C      ----------------------------------------------------------------
       SUBROUTINE LINESEARCH(X,F,D,GTD,FNEW, 
     & XNEW,M,LASTFVALUES,IFCNT,IWL)
       IMPLICIT REAL*8 (A-H,O-Z)
       PARAMETER (NOBS=601,NOBS2=2*NOBS)
       COMMON /OPTICS0/ HF,HS
       COMMON /OPTICS1/ WL(NOBS),NS(NOBS),KS(NOBS),T(NOBS)
       COMMON /OPTICS2/ TEXP(NOBS)             
       INTEGER M, IFCNT
       REAL*8 X(NOBS2+1),F,D(NOBS2+1),GTD,FNEW, XNEW(NOBS2+1),HF,HS
       REAL*8 LASTFVALUES(0:M-1),GAMMA,NS,KS,T,TEXP,WL
       PARAMETER (GAMMA=1.0D-04)
       REAL*8 ALPHA, ATEMP, FMAX
       FMAX=LASTFVALUES(0)
       DO I=1,M-1
         FMAX=DMAX1(FMAX, LASTFVALUES(I))
       ENDDO
       DO I=1,NOBS2+1
          XNEW(I)=X(I)+D(I)
       ENDDO
       CALL EVALF(XNEW,FNEW)
       IFCNT=IFCNT+1
       ALPHA=1.0D0
 100   IF (.NOT.FNEW.LE.FMAX+GAMMA*ALPHA*GTD) THEN     
         IF(ALPHA.LE.0.1D0) THEN
             ALPHA= ALPHA/2.0D0
         ELSE
             ATEMP=(-GTD*ALPHA**2)/(2.0D0*(FNEW-F-ALPHA*GTD))
             IF(ATEMP.LT.0.1D0.OR.ATEMP.GT.0.9D0*ALPHA) THEN
                 ATEMP=ALPHA/2.0D0
             ENDIF
             ALPHA=ATEMP
         ENDIF
         DO I=1,NOBS2+1
             XNEW(I)=X(I)+ALPHA*D(I)
         END DO
         CALL EVALF(XNEW,FNEW)
         IFCNT=IFCNT+1
         GOTO 100
       ENDIF
       RETURN
       END
C      -----------------------------------------------------------------
       SUBROUTINE EVALF(Y,F)
       IMPLICIT REAL*8 (A-H,O-Z)
       PARAMETER (NOBS=601,NOBS2=2*NOBS)
       COMMON /OPTICS0/ HF,HS
       COMMON /OPTICS1/ WL(NOBS),NS(NOBS),KS(NOBS),T(NOBS)
       COMMON /OPTICS2/ TEXP(NOBS)
       REAL*8 WL,TEXP,NS,KS,T,Y(NOBS2+1)
       REAL*8 NF2,KF2,NS2,KS2,N2PK2F,N2PK2S,F,HF,HS
       PI=4.0D0*DATAN(1.0D0)
       FACTOR=1.0D0/DBLE(NOBS)
       F=0.0D0
       DO I=1,NOBS
         J=I+NOBS        
         T(I)=0.0D0       
         NF2=Y(I)*Y(I) 
         KF2=Y(I+NOBS)*Y(I+NOBS) 
         NS2=NS(I)*NS(I) 
         KS2=KS(I)*KS(I)
         N2PK2F=NF2+KF2
         N2PK2S=NS2+KS2
         CONSTLF=4.0D0*PI*Y(NOBS2+1)/WL(I)
         CONSTLS=4.0D0*PI*HS/WL(I)
         PHI=CONSTLF*Y(I)        
         X=DEXP(-CONSTLF*Y(I+NOBS)) 
         Y0=DEXP(-CONSTLS*KS(I))
         Y2=Y0*Y0         
         A=64.0D0*N2PK2F*N2PK2S*Y0         
         B0=(Y(I)+1.0D0)**2+Y(I+NOBS)**2
         B1=(Y(I)+NS(I))**2+(Y(I+NOBS)+KS(I))**2
         B2=(NS(I)+1.0D0)**2+KS(I)**2
         B3=(Y(I)-NS(I))**2+(Y(I+NOBS)-KS(I))**2
         B4=(NS(I)-1.0D0)**2+KS(I)**2         
         B=B0*(B1*B2-Y2*B3*B4)         
         DUMMYC=DCOS(PHI)
         DUMMYS=DSIN(PHI)
         B24P=B2+Y2*B4
         B24M=B2-Y2*B4
         C1=Y(I)**2-1.0D0+Y(I+NOBS)**2
         C2=Y(I)**2-NS(I)**2+Y(I+NOBS)**2-KS(I)**2
         C3=Y(I+NOBS)*NS(I)-KS(I)*Y(I)         
         C=2.0D0*(C1*C2*B24M-4.0D0*Y(I+NOBS)*C3*B24P)*DUMMYC
         C=C-4.0D0*(C1*C3*B24P+Y(I+NOBS)*C2*B24M)*DUMMYS         
         D0=(Y(I)-1.0D0)**2+Y(I+NOBS)**2
         D=D0*(B2*B3-Y2*B1*B4)        
         X2=X*X 
         Q=B-C*X+D*X2 
         XDQ=X/Q    
         AXDQ=A*XDQ     
         T(I)=AXDQ
          
C         Z00=((NS(I)+1.0D0)**2+KS(I)**2)/((NS(I)-1.0D0)**2+KS(I)**2)
C         Z12=(Y(I)+1.0D0)**2+Y(J)**2
C         Z23=(Y(I)+NS(I))**2+(Y(J)+KS(I))**2
C         Z34=(NS(I)+1.0D0)**2+KS(I)**2
C         H12=(Y(I)-1.0D0)**2+Y(J)**2
C         H23=(Y(I)-NS(I))**2+(Y(J)-KS(I))**2
C         W12=1.0D0-Y(I)**2-Y(J)**2
C         W23=Y(I)**2-NS(I)**2+Y(J)**2-KS(I)**2
C         V12=-2.0D0*Y(J)
C         V23=2.0D0*(Y(J)*NS(I)-KS(I)*Y(I))
C         XP=W12*W23+V12*V23
C         XM=W12*W23-V12*V23
C         YP=W12*V23+W23*V12
C         YM=W12*V23-W23*V12
         
C         B00=Z12*Z23+H12*H23*X2+2.0D0*X*(XM*DUMMYC-YP*DUMMYS)
C                                     R13:ARTICLE:+
C         R13=H12*Z23+H23*Z12*X2+2.0D0*X*(XP*DUMMYC-YM*DUMMYS)
C                                     R31:ARTICLE:-           
C         R31=H23*Z12+H12*Z23*X2+2.0D0*X*(XP*DUMMYC+YM*DUMMYS)
C         D00=256.0D0*X2*(NS(I)**2+KS(I)**2)*(Y(I)**2+Y(J)**2)**2
C         R(I)=(R13+D00*Y2/(Z00*B00-R31*Y2))/B00         
         
         TDIFF=TEXP(I)-T(I) 
         F=F+TDIFF*TDIFF
       ENDDO
       F=FACTOR*F
       RETURN
       END
c      ------------------------------------------------------------------
       SUBROUTINE EVALG(Y,G,IWL)
       IMPLICIT REAL*8 (A-H,O-Z)
       PARAMETER (NOBS=601,NOBS2=2*NOBS)
       COMMON /OPTICS0/ HF,HS
       COMMON /OPTICS1/ WL(NOBS),NS(NOBS),KS(NOBS),T(NOBS)
       COMMON /OPTICS2/ TEXP(NOBS)
       REAL*8 WL,TEXP,NS,KS,T,Y(NOBS2+1),DFDNF(NOBS),DFDKF(NOBS)
       REAL*8 NF2,KF2,NS2,KS2,N2PK2F,N2PK2S,F,HF,HS,G(NOBS2+1)
       PI=4.0D0*DATAN(1.0D0)
       FACTOR=1.0D0/DBLE(NOBS2)
       F=0.0D0
       DO I=1,NOBS 
         J=I+NOBS       
         T(I)=0.0D0
         DFDNF(I)=0.0D0
         DFDKF(I)=0.0D0
         G(I)=0.0D0
         G(I+NOBS)=0.0D0 
             
         NF2=Y(I)*Y(I) 
         KF2=Y(J)*Y(J) 
         NS2=NS(I)*NS(I) 
         KS2=KS(I)*KS(I)
         N2PK2F=NF2+KF2
         N2PK2S=NS2+KS2
         CONSTLF=4.0D0*PI*Y(NOBS2+1)/WL(I)
         CONSTLS=4.0D0*PI*HS/WL(I)
         PHI=CONSTLF*Y(I)        
         X=DEXP(-CONSTLF*Y(J)) 
         Y0=DEXP(-CONSTLS*KS(I))
         Y2=Y0*Y0         
         A=64.0D0*N2PK2F*N2PK2S*Y0         
         B0=(Y(I)+1.0D0)**2+Y(J)**2
         B1=(Y(I)+NS(I))**2+(Y(J)+KS(I))**2
         B2=(NS(I)+1.0D0)**2+KS(I)**2
         B3=(Y(I)-NS(I))**2+(Y(J)-KS(I))**2
         B4=(NS(I)-1.0D0)**2+KS(I)**2         
         B=B0*(B1*B2-Y2*B3*B4)         
         DUMMYC=DCOS(PHI)
         DUMMYS=DSIN(PHI)
         B24P=B2+Y2*B4
         B24M=B2-Y2*B4
         C1=Y(I)**2-1.0D0+Y(J)**2
         C2=Y(I)**2-NS(I)**2+Y(J)**2-KS(I)**2
         C3=Y(J)*NS(I)-KS(I)*Y(I)         
         C=2.0D0*(C1*C2*B24M-4.0D0*Y(J)*C3*B24P)*DUMMYC
         C=C-4.0D0*(C1*C3*B24P+Y(J)*C2*B24M)*DUMMYS 
         
         IF(I.EQ.IWL) THEN
           DCDH=-2.0D0*(C1*C2*B24M-4.0D0*Y(J)*C3*B24P)*DUMMYS
           DCDH=DCDH-4.0D0*(C1*C3*B24P+Y(J)*C2*B24M)*DUMMYC
           DCDH=4.0D0*PI*Y(I)*DCDH/WL(I)
           DXDH=-4.0D0*PI*Y(I+NOBS)*X/WL(I)
         ENDIF
                
         D0=(Y(I)-1.0D0)**2+Y(J)**2
         D=D0*(B2*B3-Y2*B1*B4)        
         X2=X*X 
         Q=B-C*X+D*X2 
         XDQ=X/Q    
         AXDQ=A*XDQ     
         T(I)=AXDQ 
         TDIFF=TEXP(I)-T(I) 
         F=F+TDIFF*TDIFF
         
C        -------------
C        COMPUTE DF/DT
C        -------------
         DFDT=-2.0D0*TDIFF 
C        -----------------------------------------
C        COMPUTE DF/DA, DF/DB, DF/DC, DF/DD, DF/DX
C        -----------------------------------------
         DTDA=XDQ 
         DFDA=DFDT*DTDA 

         AXDQ2=AXDQ/Q
         AXDQ3=AXDQ2/Q 

         DTDB=-AXDQ2 
         DFDB=DFDT*DTDB 

         DTDC=X*AXDQ2 
         DFDC=DFDT*DTDC 

         DTDD=-X2*AXDQ2 
         DFDD=DFDT*DTDD 
         DXC=2.0D0*D*X-C
         DTDX=A/Q-AXDQ2*DXC 
         DFDX=DFDT*DTDX 
         
         IF(I.EQ.IWL) DFDH=DFDT*(DTDX*DXDH+DTDC*DCDH)
         
C        --------------------------------
C        COMPUTE DF/DN, DF/DK, AND DF/DH 
C        --------------------------------
         DADN=128.0D0*Y(I)*N2PK2S*Y0 
         DFDNF(I)=DFDNF(I)+DFDA*DADN 

         DADK=128.0D0*Y(J)*N2PK2S*Y0
         DFDKF(I)=DFDKF(I)+DFDA*DADK 

         DBDN=2.0D0*(Y(I)+1.0D0)*(B1*B2-Y2*B3*B4)
         DBDN=DBDN+2.0D0*B0*(B2*(Y(I)+NS(I))-Y2*B4*(Y(I)-NS(I)))
         DFDNF(I)=DFDNF(I)+DFDB*DBDN 

         DBDK=2.0D0*Y(J)*(B1*B2-Y2*B3*B4)
         DBDK=DBDK+2.0D0*B0*(B2*(Y(J)+KS(I))-Y2*B4*(Y(J)-KS(I)))
         DFDKF(I)=DFDKF(I)+DFDB*DBDK 

         DCDN1=4.0D0*(Y(I)*(C1+C2)*B24M+2.0D0*Y(J)*KS(I)*B24P)*DUMMYC
         DCDN2=-4.0D0*((2.0D0*Y(I)*C3-KS(I)*C1)*B24P+
     &        2.0D0*Y(J)*Y(I)*B24M)*DUMMYS
         DCDN3=-2.0D0*CONSTLF*(C1*C2*B24M-4.0D0*Y(J)*C3*B24P)*DUMMYS
         DCDN4=-4.0D0*CONSTLF*(C1*C3*B24P+Y(J)*C2*B24M)*DUMMYC
         DCDN=DCDN1+DCDN2+DCDN3+DCDN4
         DFDNF(I)=DFDNF(I)+DFDC*DCDN 

         DCDK1=Y(J)*(C1+C2)*B24M-2.0D0*(C3+Y(J)*NS(I))*B24P
         DCDK1=4.0D0*DCDK1*DUMMYC
         DCDK2=(2.0D0*Y(J)*C3+C1*NS(I))*B24P+(C2+2.0D0*Y(J)*Y(J))*B24M
         DCDK2=-4.0D0*DCDK2*DUMMYS
         DCDK=DCDK1+DCDK2
         DFDKF(I)=DFDKF(I)+DFDC*DCDK 

         DDDN1=2.0D0*(Y(I)-1.0D0)*(B2*B3-Y2*B1*B4)
         DDDN2=2.0D0*D0*((Y(I)-NS(I))*B2-Y2*(Y(I)+NS(I))*B4)
         DDDN=DDDN1+DDDN2
         DFDNF(I)=DFDNF(I)+DFDD*DDDN 

         DDDK1=2.0D0*Y(J)*(B2*B3-Y2*B1*B4)
         DDDK2=2.0D0*D0*((Y(J)-KS(I))*B2-Y2*(Y(J)+KS(I))*B4)
         DDDK=DDDK1+DDDK2
         DFDKF(I)=DFDKF(I)+DFDD*DDDK 

         DXDK=-CONSTLF*X 
         DFDKF(I)=DFDKF(I)+DFDX*DXDK 
           
         G(I)=DFDNF(I)
         G(I+NOBS)=DFDKF(I)         
       ENDDO
       F=FACTOR*F
       G(NOBS2+1)=DFDH
       DO I=1,NOBS2+1
          G(I)=FACTOR*G(I)
       ENDDO
       
       RETURN
       END
C      --------------------------------------------------------------
       SUBROUTINE EVALFG(Y,F,G,IWL)
       IMPLICIT REAL*8 (A-H,O-Z)
       PARAMETER (NOBS=601,NOBS2=2*NOBS)
       COMMON /OPTICS0/ HF,HS
       COMMON /OPTICS1/ WL(NOBS),NS(NOBS),KS(NOBS),T(NOBS)
       COMMON /OPTICS2/ TEXP(NOBS)
       REAL*8 WL,TEXP,NS,KS,T,Y(NOBS2+1),DFDNF(NOBS),DFDKF(NOBS)
       REAL*8 NF2,KF2,NS2,KS2,N2PK2F,N2PK2S,F,HF,HS,G(NOBS2+1)
       PI=4.0D0*DATAN(1.0D0)
       FACTOR=1.0D0/DBLE(NOBS)
       F=0.0D0
       DO I=1,NOBS 
         J=I+NOBS       
         T(I)=0.0D0
         DFDNF(I)=0.0D0
         DFDKF(I)=0.0D0
         G(I)=0.0D0
         G(I+NOBS)=0.0D0 
             
         NF2=Y(I)*Y(I) 
         KF2=Y(J)*Y(J) 
         NS2=NS(I)*NS(I) 
         KS2=KS(I)*KS(I)
         N2PK2F=NF2+KF2
         N2PK2S=NS2+KS2
         CONSTLF=4.0D0*PI*Y(NOBS2+1)/WL(I)
         CONSTLS=4.0D0*PI*HS/WL(I)
         PHI=CONSTLF*Y(I)        
         X=DEXP(-CONSTLF*Y(J)) 
         Y0=DEXP(-CONSTLS*KS(I))
         Y2=Y0*Y0         
         A=64.0D0*N2PK2F*N2PK2S*Y0         
         B0=(Y(I)+1.0D0)**2+Y(J)**2
         B1=(Y(I)+NS(I))**2+(Y(J)+KS(I))**2
         B2=(NS(I)+1.0D0)**2+KS(I)**2
         B3=(Y(I)-NS(I))**2+(Y(J)-KS(I))**2
         B4=(NS(I)-1.0D0)**2+KS(I)**2         
         B=B0*(B1*B2-Y2*B3*B4)         
         DUMMYC=DCOS(PHI)
         DUMMYS=DSIN(PHI)
         B24P=B2+Y2*B4
         B24M=B2-Y2*B4
         C1=Y(I)**2-1.0D0+Y(J)**2
         C2=Y(I)**2-NS(I)**2+Y(J)**2-KS(I)**2
         C3=Y(J)*NS(I)-KS(I)*Y(I)         
         C=2.0D0*(C1*C2*B24M-4.0D0*Y(J)*C3*B24P)*DUMMYC
         C=C-4.0D0*(C1*C3*B24P+Y(J)*C2*B24M)*DUMMYS   
         
         IF(I.EQ.IWL) THEN
           DCDH=-2.0D0*(C1*C2*B24M-4.0D0*Y(J)*C3*B24P)*DUMMYS
           DCDH=DCDH-4.0D0*(C1*C3*B24P+Y(J)*C2*B24M)*DUMMYC
           DCDH=4.0D0*PI*Y(I)*DCDH/WL(I)
           DXDH=-4.0D0*PI*Y(I+NOBS)*X/WL(I)
         ENDIF  
             
         D0=(Y(I)-1.0D0)**2+Y(J)**2
         D=D0*(B2*B3-Y2*B1*B4)        
         X2=X*X 
         Q=B-C*X+D*X2 
         XDQ=X/Q    
         AXDQ=A*XDQ     
         T(I)=AXDQ 
         TDIFF=TEXP(I)-T(I) 
         F=F+TDIFF*TDIFF
         
C        -------------
C        COMPUTE DF/DT
C        -------------
         DFDT=-2.0D0*TDIFF 
C        -----------------------------------------
C        COMPUTE DF/DA, DF/DB, DF/DC, DF/DD, DF/DX
C        -----------------------------------------
         DTDA=XDQ 
         DFDA=DFDT*DTDA 

         AXDQ2=AXDQ/Q
         AXDQ3=AXDQ2/Q 

         DTDB=-AXDQ2 
         DFDB=DFDT*DTDB 

         DTDC=X*AXDQ2 
         DFDC=DFDT*DTDC 

         DTDD=-X2*AXDQ2 
         DFDD=DFDT*DTDD 
         DXC=2.0D0*D*X-C
         DTDX=A/Q-AXDQ2*DXC 
         DFDX=DFDT*DTDX 
         
         IF(I.EQ.IWL) DFDH=DFDT*(DTDX*DXDH+DTDC*DCDH)
         
C        --------------------
C        COMPUTE DF/DN, DF/DK 
C        --------------------
         DADN=128.0D0*Y(I)*N2PK2S*Y0 
         DFDNF(I)=DFDNF(I)+DFDA*DADN 

         DADK=128.0D0*Y(J)*N2PK2S*Y0
         DFDKF(I)=DFDKF(I)+DFDA*DADK 

         DBDN=2.0D0*(Y(I)+1.0D0)*(B1*B2-Y2*B3*B4)
         DBDN=DBDN+2.0D0*B0*(B2*(Y(I)+NS(I))-Y2*B4*(Y(I)-NS(I)))
         DFDNF(I)=DFDNF(I)+DFDB*DBDN 

         DBDK=2.0D0*Y(J)*(B1*B2-Y2*B3*B4)
         DBDK=DBDK+2.0D0*B0*(B2*(Y(J)+KS(I))-Y2*B4*(Y(J)-KS(I)))
         DFDKF(I)=DFDKF(I)+DFDB*DBDK 

         DCDN1=4.0D0*(Y(I)*(C1+C2)*B24M+2.0D0*Y(J)*KS(I)*B24P)*DUMMYC
         DCDN2=-4.0D0*((2.0D0*Y(I)*C3-KS(I)*C1)*B24P+
     &        2.0D0*Y(J)*Y(I)*B24M)*DUMMYS
         DCDN3=-2.0D0*CONSTLF*(C1*C2*B24M-4.0D0*Y(J)*C3*B24P)*DUMMYS
         DCDN4=-4.0D0*CONSTLF*(C1*C3*B24P+Y(J)*C2*B24M)*DUMMYC
         DCDN=DCDN1+DCDN2+DCDN3+DCDN4
         DFDNF(I)=DFDNF(I)+DFDC*DCDN 

         DCDK1=Y(J)*(C1+C2)*B24M-2.0D0*(C3+Y(J)*NS(I))*B24P
         DCDK1=4.0D0*DCDK1*DUMMYC
         DCDK2=(2.0D0*Y(J)*C3+C1*NS(I))*B24P+(C2+2.0D0*Y(J)*Y(J))*B24M
         DCDK2=-4.0D0*DCDK2*DUMMYS
         DCDK=DCDK1+DCDK2
         DFDKF(I)=DFDKF(I)+DFDC*DCDK 

         DDDN1=2.0D0*(Y(I)-1.0D0)*(B2*B3-Y2*B1*B4)
         DDDN2=2.0D0*D0*((Y(I)-NS(I))*B2-Y2*(Y(I)+NS(I))*B4)
         DDDN=DDDN1+DDDN2
         DFDNF(I)=DFDNF(I)+DFDD*DDDN 

         DDDK1=2.0D0*Y(J)*(B2*B3-Y2*B1*B4)
         DDDK2=2.0D0*D0*((Y(J)-KS(I))*B2-Y2*(Y(J)+KS(I))*B4)
         DDDK=DDDK1+DDDK2
         DFDKF(I)=DFDKF(I)+DFDD*DDDK 

         DXDK=-CONSTLF*X 
         DFDKF(I)=DFDKF(I)+DFDX*DXDK 
           
         G(I)=DFDNF(I)
         G(I+NOBS)=DFDKF(I)
         
       ENDDO
       F=FACTOR*F
       G(NOBS2+1)=DFDH
       DO I=1,NOBS2+1
          G(I)=FACTOR*G(I)
       ENDDO
       
       RETURN
       END
C      -------------------------------------------------------------------
C      ESTA SUBRUTINA PERMITE CALCULAR EL INDICE DE REFRACCION DEL CUARZO. 
C      LA REFERENCIA BASICA ES:
C    
C      I. H. Malitson "Interspecimen comparison of the refractive index of
C      fused silica" J. Opt. Soc. Am. 55, 1205-1209 (1965). 
C
C      La subrutina es válida para longitudes de onda entre 210 y 6700 nm
C      -------------------------------------------------------------------
       SUBROUTINE FUSED_QUARTZ(WL,RN,RK)
C      WL (WAVELENGTH) MUST BE IN MICROMETERS
       IMPLICIT REAL*8 (A-H,O-Z)
       RK=0.0D0
       A1=0.6961663D0
       B1=0.0684043D0
       A2=0.4079426D0
       B2=0.1162414D0
       A3=0.8974794D0
       B3=9.8961610D0
       WL2=WL*WL
       RN=A1*WL2/(WL2-B1*B1)+A2*WL2/(WL2-B2*B2)
       RN=RN+A3*WL2/(WL2-B3*B3)
       RN=DSQRT(1.0D0+RN)
       RETURN
       END
C      -----------------------------------------------------------------
       SUBROUTINE PALADIUM(HW,RN,RK)
       IMPLICIT REAL*8 (A-H,O-Z)
       PARAMETER (NO=4)
       COMMON /OPTICS4/ FS(0:NO),G(0:NO),W(1:NO)
       REAL*8 HW,RN,RK,FS,G,W,ER,EI,SIGNO,WP
       COMPLEX*16 EPSF,II
       II=DCMPLX(0.0D0,1.0D0)
       WP=9.72D0
       EPSF=1.0D0-FS(0)*WP*WP/(HW*(HW+II*G(0)))
       DO J=1,NO
          EPSF=EPSF+FS(J)*WP*WP/(W(J)*W(J)-HW*HW-II*HW*G(J))
       ENDDO       
       ER=DREAL(EPSF)
       EI=DIMAG(EPSF)
       SIGNO=1.0D0
       IF(ER.LT.0.0D0) SIGNO=-1.0D0
       DUMMY=DSQRT(1.0D0+(EI/ER)**2)
       RN=DSQRT(0.50D0*ER*(+1.0D0+SIGNO*DUMMY))
       RK=DSQRT(0.50D0*ER*(-1.0D0+SIGNO*DUMMY))
       RETURN
       END  
C     ------------------------------------------------------------------
      SUBROUTINE PDVARGAS_DL
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NO=4)
      COMMON /OPTICS4/ FS(0:NO),G(0:NO),W(1:NO)
      REAL*8 FS,G,W
      FS(0)=0.215D+00  
      G(0)=0.112D+00
      FS(1)=0.146D+00   
      G(1)=0.366D+00   
      W(1)=0.442D+00  
      FS(2)=0.437D+00   
      G(2)=0.1606D+01   
      W(2)=0.1091D+01   
      FS(3)=0.1908D+01  
      G(3)=0.1163D+02   
      W(3)=0.353D+01
      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE PDRAKIC_DL        
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NO=4)
      COMMON /OPTICS4/ FS(0:NO),G(0:NO),W(1:NO)
      REAL*8 FS,G,W
      FS(0)=0.330D+00  
      G(0)=0.008D+00
      FS(1)=0.649D+00   
      G(1)=2.950D+00   
      W(1)=0.336D+00  
      FS(2)=0.121D+00   
      G(2)=0.555D+00   
      W(2)=0.501D+00   
      FS(3)=0.638D+00  
      G(3)=0.4621D+01   
      W(3)=1.659D+00
      FS(4)=0.453D0
      G(4)=3.236D0
      W(4)=5.715D0
      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE PDRAKIC_DLBB        
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NO=4)
      COMMON /OPTICS4/ FS(0:NO),G(0:NO),W(1:NO)
      COMMON /OPTICS5/ S(1:NO)
      REAL*8 FS,G,W
      FS(0)=0.333D0 
      G(0)=0.009D0
      FS(1)=0.769D0   
      G(1)=2.343D0   
      W(1)=0.066D0
      S(1)=0.694D0  
      FS(2)=0.093D0   
      G(2)=0.497D0   
      W(2)=0.502D0 
      S(2)=0.027D0  
      FS(3)=0.309D0  
      G(3)=2.022D0   
      W(3)=2.432D0
      S(3)=1.167D0
      FS(4)=0.409D0
      G(4)=0.119D0
      W(4)=5.987D0
      S(4)=1.331D0
      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE CHIJ(FC,HW,WP2,FJ,GJ,WJ,SJ,EPSJ)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 FJ,GJ,WJ,SJ,AJR,AJI,U,V,XI,YI,FC
      COMPLEX*16 II,AJ,ZMIN,ZMAX,WMIN,WMAX,EPSJ,CWJ
      LOGICAL FLAG         
      CWJ=DCMPLX(WJ,0.0D0)
      II=DCMPLX(0.0D0,1.0D0)
      PI=4.0D0*DATAN(1.0D0)
      SQPI=DSQRT(PI)
      SQR2=DSQRT(2.0D0)
      AJR=HW*DSQRT((DSQRT(1.0D0+(GJ/HW)**2)+1.0D0))/SQR2
      AJI=HW*DSQRT((DSQRT(1.0D0+(GJ/HW)**2)-1.0D0))/SQR2
      AJ=DCMPLX(AJR,AJI)
      ZMIN=(AJ-CWJ)/(SJ*SQR2)
      XI=DREAL(ZMIN)
      YI=AIMAG(ZMIN)
      CALL WOFZ(XI,YI,U,V,FLAG)
      WMIN=DCMPLX(U,V)
      ZMAX=(AJ+CWJ)/(SJ*SQR2)
      XI=DREAL(ZMAX)
      YI=AIMAG(ZMAX)
      CALL WOFZ(XI,YI,U,V,FLAG)
      WMAX=DCMPLX(U,V)      
      EPSJ=FC*II*SQPI*FJ*WP2*(WMIN+WMAX)/(2.0D0*SQR2*SJ*AJ)
      RETURN
      END
C     -------------------------------------------------------------------
C     SUBRUTINA PARA CALCULAR EL INDICE DE REFRACCIÓN DEL SODA LIME GLASS
C     FÓRMULA: n(?) = 1.5130 - 0.003169*(?^2) + 0.003962/(?^2)
C     DONDE ? DEBE ESTAR EN MICRÓMETROS
C     -------------------------------------------------------------------
       SUBROUTINE SODA_LIME_GLASS(WL, RN, RK)
C      WL (WAVELENGTH) MUST BE IN MICROMETERS
       IMPLICIT REAL*8 (A-H,O-Z)
       RK = 0.0D0  ! Asumimos coeficiente de extinción cero para SLG
       WL2 = WL * WL
       RN = 1.5130D0 - 0.003169D0 * WL2 + 0.003962D0 / WL2
       RETURN
       END
C     ------------------------------------------------------------------
        SUBROUTINE WOFZ(XI,YI,U,V,FLAG)
C
C  GIVEN A COMPLEX NUMBER Z = (XI,YI), THIS SUBROUTINE COMPUTES
C  THE VALUE OF THE FADDEEVA-FUNCTION W(Z) = EXP(-Z**2)*ERFC(-I*Z),
C  WHERE ERFC IS THE COMPLEX COMPLEMENTARY ERROR-FUNCTION AND I
C  MEANS SQRT(-1).
C  THE ACCURACY OF THE ALGORITHM FOR Z IN THE 1ST AND 2ND QUADRANT
C  IS 14 SIGNIFICANT DIGITS; IN THE 3RD AND 4TH IT IS 13 SIGNIFICANT
C  DIGITS OUTSIDE A CIRCULAR REGION WITH RADIUS 0.126 AROUND A ZERO
C  OF THE FUNCTION.
C  ALL REAL VARIABLES IN THE PROGRAM ARE DOUBLE PRECISION.
C
C
C  THE CODE CONTAINS A FEW COMPILER-DEPENDENT PARAMETERS :
C     RMAXREAL = THE MAXIMUM VALUE OF RMAXREAL EQUALS THE ROOT OF
C                RMAX = THE LARGEST NUMBER WHICH CAN STILL BE
C                IMPLEMENTED ON THE COMPUTER IN DOUBLE PRECISION
C                FLOATING-POINT ARITHMETIC
C     RMAXEXP  = LN(RMAX) - LN(2)
C     RMAXGONI = THE LARGEST POSSIBLE ARGUMENT OF A DOUBLE PRECISION
C                GONIOMETRIC FUNCTION (DCOS, DSIN, ...)
C  THE REASON WHY THESE PARAMETERS ARE NEEDED AS THEY ARE DEFINED WILL
C  BE EXPLAINED IN THE CODE BY MEANS OF COMMENTS
C
C
C  PARAMETER LIST
C     XI     = REAL      PART OF Z
C     YI     = IMAGINARY PART OF Z
C     U      = REAL      PART OF W(Z)
C     V      = IMAGINARY PART OF W(Z)
C     FLAG   = AN ERROR FLAG INDICATING WHETHER OVERFLOW WILL
C              OCCUR OR NOT; TYPE LOGICAL;
C              THE VALUES OF THIS VARIABLE HAVE THE FOLLOWING
C              MEANING :
C              FLAG=.FALSE. : NO ERROR CONDITION
C              FLAG=.TRUE.  : OVERFLOW WILL OCCUR, THE ROUTINE
C                             BECOMES INACTIVE
C  XI, YI      ARE THE INPUT-PARAMETERS
C  U, V, FLAG  ARE THE OUTPUT-PARAMETERS
C
C  FURTHERMORE THE PARAMETER FACTOR EQUALS 2/SQRT(PI)
C
C  THE ROUTINE IS NOT UNDERFLOW-PROTECTED BUT ANY VARIABLE CAN BE
C  PUT TO 0 UPON UNDERFLOW;
C
C  REFERENCE - GPM POPPE, CMJ WIJERS; MORE EFFICIENT COMPUTATION OF
C  THE COMPLEX ERROR-FUNCTION, ACM TRANS. MATH. SOFTWARE.
C
*
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
         LOGICAL A, B, FLAG
         PARAMETER (FACTOR   = 1.12837916709551257388D0,
     *           RMAXREAL = 0.5D+154,
     *           RMAXEXP  = 708.503061461606D0,
     *           RMAXGONI = 3.53711887601422D+15)
*
         FLAG = .FALSE.
*
         XABS = DABS(XI)
         YABS = DABS(YI)
         X    = XABS/6.3D0
         Y    = YABS/4.4D0
*
C
C     THE FOLLOWING IF-STATEMENT PROTECTS
C     QRHO = (X**2 + Y**2) AGAINST OVERFLOW
C
         IF ((XABS.GT.RMAXREAL).OR.(YABS.GT.RMAXREAL)) GOTO 100
*
         QRHO=X**2+Y**2
*
         XABSQ=XABS**2
         XQUAD=XABSQ-YABS**2
         YQUAD=2.0D0*XABS*YABS
*
         A=QRHO.LT.0.085264D0
*
         IF (A) THEN
C
C  IF (QRHO.LT.0.085264D0) THEN THE FADDEEVA-FUNCTION IS EVALUATED
C  USING A POWER-SERIES (ABRAMOWITZ/STEGUN, EQUATION (7.1.5), P.297)
C  N IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
C  ACCURACY
C
           QRHO=(1.0D0-0.85*Y)*DSQRT(QRHO)
           N=IDNINT(6.0D0+72.0D0*QRHO)
           J=2*N+1
           XSUM=1.0D0/DBLE(J)
           YSUM=0.0D0
           DO 10 I=N,1,-1
             J=J-2
             XAUX=(XSUM*XQUAD-YSUM*YQUAD)/DBLE(I)
             YSUM=(XSUM*YQUAD+YSUM*XQUAD)/DBLE(I)
             XSUM=XAUX+1.0D0/DBLE(J)
 10     CONTINUE
           U1=-FACTOR*(XSUM*YABS+YSUM*XABS)+1.0D0
           V1=FACTOR*(XSUM*XABS-YSUM*YABS)
           DAUX=DEXP(-XQUAD)
           U2=DAUX*DCOS(YQUAD)
           V2=-DAUX*DSIN(YQUAD)
*
           U=U1*U2-V1*V2
           V=U1*V2+V1*U2
*
         ELSE
C
C  IF (QRHO.GT.1.O) THEN W(Z) IS EVALUATED USING THE LAPLACE
C  CONTINUED FRACTION
C  NU IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
C  ACCURACY
C
C  IF ((QRHO.GT.0.085264D0).AND.(QRHO.LT.1.0)) THEN W(Z) IS EVALUATED
C  BY A TRUNCATED TAYLOR EXPANSION, WHERE THE LAPLACE CONTINUED FRACTION
C  IS USED TO CALCULATE THE DERIVATIVES OF W(Z)
C  KAPN IS THE MINIMUM NUMBER OF TERMS IN THE TAYLOR EXPANSION NEEDED
C  TO OBTAIN THE REQUIRED ACCURACY
C  NU IS THE MINIMUM NUMBER OF TERMS OF THE CONTINUED FRACTION NEEDED
C  TO CALCULATE THE DERIVATIVES WITH THE REQUIRED ACCURACY
C
*
          IF (QRHO.GT.1.0D0) THEN
             H=0.0D0
             KAPN=0
             QRHO=DSQRT(QRHO)
             NU=IDINT(3.0D0+(1442.0D0/(26.0D0*QRHO+77.0D0)))
          ELSE
             QRHO=(1.0D0-Y)*DSQRT(1.0D0-QRHO)
             H=1.88D0*QRHO
             H2=2.0D0*H
             KAPN=IDNINT(7.0D0+34.0D0*QRHO)
             NU= IDNINT(16.0D0+26.0D0*QRHO)
          ENDIF
*
          B=(H.GT.0.0D0)
*
          IF (B) QLAMBDA=H2**KAPN
*
           RX=0.0D0
           RY=0.0D0
           SX=0.0D0
           SY=0.0D0
*
          DO 11 N=NU, 0, -1
             NP1=N+1
             TX=YABS+H+NP1*RX
             TY=XABS-NP1*RY
             C=0.5D0/(TX**2+TY**2)
             RX=C*TX
             RY=C*TY
             IF ((B).AND.(N.LE.KAPN)) THEN
               TX=QLAMBDA+SX
               SX=RX*TX-RY*SY
               SY=RY*TX+RX*SY
               QLAMBDA=QLAMBDA/H2
             ENDIF
 11     CONTINUE
*
          IF (H.EQ.0.0D0) THEN
             U=FACTOR*RX
             V=FACTOR*RY
         ELSE
             U=FACTOR*SX
             V=FACTOR*SY
         END IF
*
          IF (YABS.EQ.0.0D0) U=DEXP(-XABS**2)
*
        END IF
*
*
C
C  EVALUATION OF W(Z) IN THE OTHER QUADRANTS
C
*
         IF (YI.LT.0.0D0) THEN
*
           IF (A) THEN
             U2=2.0D0*U2
             V2=2.0D0*V2
          ELSE
             XQUAD=-XQUAD
*
C
C         THE FOLLOWING IF-STATEMENT PROTECTS 2*EXP(-Z**2)
C         AGAINST OVERFLOW
C
          IF ((YQUAD.GT.RMAXGONI).OR.
     &        (XQUAD.GT.RMAXEXP)) GOTO 100
*
             W1=2.0D0*DEXP(XQUAD)
             U2=W1*DCOS(YQUAD)
             V2=-W1*DSIN(YQUAD)
          END IF
*
           U=U2-U
           V=V2-V
           IF (XI.GT.0.0D0) V=-V
         ELSE
           IF (XI.LT.0.0D0) V=-V
         END IF
*
         RETURN
*
100   FLAG = .TRUE.
         WRITE(*,*) 'OVERFLOW...',XQUAD,YQUAD
         RETURN
        END 

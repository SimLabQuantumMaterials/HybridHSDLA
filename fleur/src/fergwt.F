      MODULE m_fergwt
c****************************************************************
c     determines the fermi energy and weights for the k-space
c     integration using gaussing-smearing method.
c                                               c.l.fu
c*****************************************************************
      CONTAINS
      SUBROUTINE fergwt(
     >                  neigd,nkptd,jspd,nkpt,jspins,ne,
     >                  delgau,zc,eig,wtkpt,
     X                  ef,
     <                  seigv,w)
#include"cpp_double.h"

      USE m_constants, ONLY : pimach
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: neigd,nkptd,jspd
      INTEGER, INTENT (IN) :: jspins,nkpt
      REAL, INTENT (IN)    :: delgau,zc
      REAL, INTENT (OUT)   :: seigv
      REAL, INTENT (INOUT) :: ef
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: ne(nkptd,jspd)
      REAL,    INTENT (IN) :: eig(neigd,nkptd,jspd)
      REAL,    INTENT (IN) :: wtkpt(nkptd)
      REAL,    INTENT (OUT):: w(neigd,nkptd,jspd)
C     ..
C     .. Local Scalars ..
      REAL chmom,de,ef0,ef1,elow,en,eps,eup,fac,fact1,s,s0,s1,s2,
     +     workf,wt,wtk,zcdiff,zero,seigv1
      INTEGER i,ifl,it,jspin,k,nbnd
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp,real,sqrt
C     ..
C     .. External Functions ..
#ifndef CPP_APC
      REAL  CPP_ALPHA_erf
      EXTERNAL  CPP_ALPHA_erf
#else
      REAL  erf
!      EXTERNAL  erf
#endif
C     ..
C     .. Data statements ..
      DATA zero/0.e0/,eps/1.e-5/,eup/3.0e0/,elow/-3.0e0/
C     ..
      fact1 = delgau/sqrt(pimach())
c     ---> determines ef
      ifl = 0
   10 DO 50 it = 1,50
         s = 0.
         DO 40 jspin = 1,jspins
            DO 30 k = 1,nkpt
               wtk = wtkpt(k)
               nbnd = ne(k,jspin)
               DO 20 i = 1,nbnd
                  en = eig(i,k,jspin)
                  de = (en-ef)/delgau
                  wt = 2.0
                  IF (de.GT.eup) wt = 0.0
                  IF (de.GE.elow .AND. de.LE.eup) THEN
                     IF (de.LT.zero) THEN
#ifndef CPP_APC
                        wt = 1. +  CPP_ALPHA_erf(-de)
#else
                        wt = 1. + erf(-de)
#endif
                     ELSE
#ifndef CPP_APC
                        !wt = CPP_ALPHA_erfc(de)
                        wt = 1. - CPP_ALPHA_erf(de)
#else
                        wt = 1. - erf(de)
#endif
                     END IF
                  END IF
                  s = s + wt*wtk
                  w(i,k,jspin) = wt/2.
   20          CONTINUE
   30       CONTINUE
   40    CONTINUE
         s = s/real(jspins)
         zcdiff = zc - s
         IF (abs(zcdiff).LT.eps) GO TO 60
         IF (ifl.EQ.0) THEN
            ifl = 1
            ef0 = ef
            ef = ef + 0.003
            s0 = s
         ELSE
            fac = (s0-s)/ (zc-s)
            IF (abs(fac).LT.1.0e-1) THEN
               ef0 = ef
               s0 = s
               IF (zcdiff.GE.zero) THEN
                  ef = ef + 0.003
               ELSE
                  ef = ef - 0.003
               END IF
            ELSE
               ef1 = ef
               ef = ef + (ef0-ef)/fac
               ef0 = ef1
               s0 = s
            END IF
         END IF
   50 CONTINUE
      eps = 1.25*eps
      WRITE (6,FMT=8000) eps
 8000 FORMAT (10x,'warning: eps has been increased to',e10.5)
      GO TO 10
   60 workf = -27.2116*ef
      WRITE (16,FMT=8010) ef,workf,s
      WRITE (6,FMT=8010) ef,workf,s
 8010 FORMAT (/,10x,'fermi energy=',f10.5,' har',3x,'work function=',
     +       f10.5,' ev',/,10x,'number of valence electrons=',f10.5)
      IF (abs(zcdiff).GT.5.0e-4) THEN
         WRITE (*,*) 'Fermi-level determination did not converege'
         STOP 'efermi : change temperature or set gauss=F'
      ENDIF
      DO 90 jspin = 1,jspins
         WRITE (6,FMT=8020) jspin
 8020    FORMAT (/,/,5x,'band-weighting factor for spin=',i5)
         DO 80 k = 1,nkpt
            nbnd = ne(k,jspin)
            WRITE (6,FMT=8030) k
 8030       FORMAT (/,5x,'k-point=',i5,/)
            DO 70 i = 1,nbnd
               w(i,k,jspin) = wtkpt(k)*w(i,k,jspin)
   70       CONTINUE
            WRITE (6,FMT=8040) (w(i,k,jspin),i=1,nbnd)
 8040       FORMAT (5x,16f5.3)
   80    CONTINUE
   90 CONTINUE
      s1 = 0.
      s2 = 0.
      seigv = 0.
      DO 120 jspin = 1,jspins
         s = 0.
         DO 110 k = 1,nkpt
            DO 100 i = 1,ne(k,jspin)
               s = s + w(i,k,jspin)
               seigv = seigv + w(i,k,jspin)*eig(i,k,jspin)
               en = eig(i,k,jspin)
               de = (en-ef)/delgau
c     ---> correction term
               IF (abs(de).LT.3.) THEN
                  de = de*de
                  s2 = s2 + exp(-de)*wtkpt(k)
               END IF
  100       CONTINUE
  110    CONTINUE
         s1 = s1 + s
  120 CONTINUE
      seigv = (2/jspins)*seigv
      seigv1 = (1/jspins)*fact1*s2
      chmom = s1 - jspins*s
      WRITE (6,FMT=8050) seigv - seigv1,s1,chmom
      WRITE (16,FMT=8050) seigv - seigv1,s1,chmom
 8050 FORMAT (/,10x,'sum of eigenvalues-correction=',f12.5,/,10x,
     +       'sum of weight                =',f12.5,/,10x,
     +       'total moment                 =',f12.5,/)

      END SUBROUTINE fergwt
      END MODULE m_fergwt


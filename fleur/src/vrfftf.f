      SUBROUTINE vrfftf(m,n,r,rt,mdimr,wsave)
C
C***BEGIN PROLOGUE  VRFFTF
C***DATE WRITTEN   850801   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM,
C             FOURIER ANALYSIS, FORWARD TRANSFORM, MULTIPLE SEQUENCES
C***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST)
C***PURPOSE  Forward real periodic transform, M sequences.
C***DESCRIPTION
C
C  Subroutine VRFFTF computes the Fourier coefficients (forward
C  transform) of a number of real periodic sequences.  Specifically,
C  for each sequence the subroutine claculates the independent
C  Fourier coefficients described below at output parameter R.
C
C  The array WSAVE which is used by subroutine VRFFTF must be
C  initialized by calling subroutine VRFFTI(N,WSAVE).
C
C
C  Input Parameters
C
C  M       the number of sequences to be transformed.
C
C  N       the length of the sequences to be transformed.  The method
C          is most efficient when N is a product of small primes,
C          however n may be any positive integer.
C
C  R       areal two-dimensional array of size MDIMX x N containing the
C          the sequences to be transformed.  The sequences are stored
C          in the ROWS of R.  Thus, the I-th sequence to be transformed,
C          X(I,J), J=0,1,...,N-1, is stored as
C
C               R(I,J) = X(I,J-1) , J=1, 2, . . . , N.
C
C  RT      a real two-dimensional work array of size MDIMX x N.
C
C  MDIMR   the row (or first) dimension of the arrays R and RT exactly
C          as they appear in the calling program.  This parameter is
C          used to specify the variable dimension of these arrays.
C
C  WSAVE   a real one-dimensional work array which must be dimensioned
C          at least N+15.  The WSAVE array must be initialized by
C          calling subroutine VRFFTI.  A different WSAVE array must be
C          used for each different value of N.  This initialization does
C          not have to be repeated so long as N remains unchanged.  The
C          same WSAVE array may be used by VRFFTF and VRFFTB.
C
C  Output Parameters
C
C  R       contains the Fourier coefficients F(K) for each of the M
C          input sequences.  Specifically, row I of R, R(I,J),
C          J=1,2,..,N, contains the independent Fourier coefficients
C          F(I,K), for the I-th input sequence stored as
C
C             R(I,1) = REAL( F(I,0) ),
C                    = SQRT(1/N)*SUM(J=0,N-1)[ X(I,J) ],
C
C             R(I,2*K) = REAL( F(I,K) )
C                      = SQRT(1/N)*SUM(J=0,N-1)[X(I,J)*COS(2J*K*PI/N)]
C
C             R(I,2*K+1) = IMAG( F(I,K) )
C                        =-SQRT(1/N)*SUM(J=0,N-1)[X(I,J)*SIN(2J*K*PI/N)]
C
C                   for K = 1, 2, . . . , M-1,
C
C              and, when N is even,
C
C              R(I,N) = REAL( F(I,N/2) ).
C                     = SQRT(1/N)*SUM(J=0,N-1)[ (-1)**J*X(I,J) ].
C
C  WSAVE   contains results which must not be destroyed between calls
C          to VRFFTF or VRFFTB.
C
C  -----------------------------------------------------------------
C
C  NOTE  -  A call of VRFFTF followed immediately by a call of
C           of VRFFTB will return the original sequences R.  Thus,
C           VRFFTB is the correctly normalized inverse of VRFFTF.
C
C  -----------------------------------------------------------------
C
C  VRFFTF is a straightforward extension of the subprogram RFFTF to
C  handle M simultaneous sequences.  RFFTF was originally developed
C  by P. N. Swarztrauber of NCAR.
C
C
C              * * * * * * * * * * * * * * * * * * * * *
C              *                                       *
C              *         PROGRAM SPECIFICATIONS        *
C              *                                       *
C              * * * * * * * * * * * * * * * * * * * * *
C
C
C     Dimension of    R(MDIMR,N), RT(MDIMR,N), WSAVE(N+15)
C     Arguments
C
C     Latest          AUGUST 1, 1985
C     Revision
C
C     Subprograms     VRFFTI, VRFTI1, VRFFTF, VRFTF1, VRADF2, VRADF3,
C     Required        VRADF4, VRADF5, VRADFG, VRFFTB, VRFTB1, VRADB2,
C                     VRADB3, VRADB4, VRADB5, VRADBG, PIMACH
C
C     Special         NONE
C     Conditions
C
C     Common          NONE
C     blocks
C
C     I/O             NONE
C
C     Precision       SINGLE
C
C     Specialist      ROLAND SWEET
C
C     Language        FORTRAN
C
C     History         WRITTEN BY LINDA LINDGREN AND ROLAND SWEET AT THE
C                     NATIONAL BUREAU OF STANDARDS (BOULDER).
C
C     Algorithm       A REAL VARIANT OF THE STOCKHAM AUTOSORT VERSION
C                     OF THE COOLEY-TUKEY FAST FOURIER TRANSFORM.
C
C     Portability     AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN 77.
C                     THE ONLY MACHINE DEPENDENT CONSTANT IS LOCATED IN
C                     THE FUNCTION PIMACH.
C
C     Required        COS,SIN
C     resident
C     Routines
C
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  VRFTF1
C***END PROLOGUE  VRFFTF
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION r(mdimr,n),rt(mdimr,n),wsave(n+15)
C***FIRST EXECUTABLE STATEMENT  VRFFTF
      IF (n.EQ.1) RETURN
      CALL vrftf1(m,n,r,rt,mdimr,wsave(1),wsave(n+1))
      RETURN
      END
      SUBROUTINE vradf2(mp,ido,l1,cc,ch,mdimc,wa1)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION ch(mdimc,ido,2,l1),cc(mdimc,ido,l1,2),wa1(ido)

      DO 20 k = 1,l1
         DO 10 m = 1,mp
            ch(m,1,1,k) = cc(m,1,k,1) + cc(m,1,k,2)
            ch(m,ido,2,k) = cc(m,1,k,1) - cc(m,1,k,2)
   10    CONTINUE
   20 CONTINUE
      IF ( ido ==  2 ) THEN
         GOTO 70
      ELSEIF ( ido < 2 ) THEN
         GOTO 100
      ENDIF
      idp2 = ido + 2
      DO 60 k = 1,l1
         DO 50 i = 3,ido,2
            ic = idp2 - i
            DO 40 m = 1,mp
               ch(m,i,1,k) = cc(m,i,k,1) +
     +                       (wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*
     +                       cc(m,i-1,k,2))
               ch(m,ic,2,k) = (wa1(i-2)*cc(m,i,k,2)-
     +                        wa1(i-1)*cc(m,i-1,k,2)) - cc(m,i,k,1)
               ch(m,i-1,1,k) = cc(m,i-1,k,1) +
     +                         (wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*
     +                         cc(m,i,k,2))
               ch(m,ic-1,2,k) = cc(m,i-1,k,1) -
     +                          (wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*
     +                          cc(m,i,k,2))
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
      IF (mod(ido,2).EQ.1) RETURN
   70 DO 90 k = 1,l1
         DO 80 m = 1,mp
            ch(m,1,2,k) = -cc(m,ido,k,2)
            ch(m,ido,1,k) = cc(m,ido,k,1)
   80    CONTINUE
   90 CONTINUE
  100 RETURN
      END
      SUBROUTINE vradf3(mp,ido,l1,cc,ch,mdimc,wa1,wa2)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      USE m_constants, ONLY : pimach
      DIMENSION ch(mdimc,ido,3,l1),cc(mdimc,ido,l1,3),wa1(ido),wa2(ido)

      arg = 2.*pimach()/3.
      taur = cos(arg)
      taui = sin(arg)
      DO 20 k = 1,l1
         DO 10 m = 1,mp
            ch(m,1,1,k) = cc(m,1,k,1) + (cc(m,1,k,2)+cc(m,1,k,3))
            ch(m,1,3,k) = taui* (cc(m,1,k,3)-cc(m,1,k,2))
            ch(m,ido,2,k) = cc(m,1,k,1) +
     +                      taur* (cc(m,1,k,2)+cc(m,1,k,3))
   10    CONTINUE
   20 CONTINUE
      IF (ido.EQ.1) RETURN
      idp2 = ido + 2
      DO 50 k = 1,l1
         DO 40 i = 3,ido,2
            ic = idp2 - i
            DO 30 m = 1,mp
               ch(m,i-1,1,k) = cc(m,i-1,k,1) +
     +                         ((wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i,
     +                         k,2))+ (wa2(i-2)*cc(m,i-1,k,
     +                         3)+wa2(i-1)*cc(m,i,k,3)))
               ch(m,i,1,k) = cc(m,i,k,1) +
     +                       ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,
     +                       2))+ (wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*cc(m,
     +                       i-1,k,3)))
               ch(m,i-1,3,k) = (cc(m,i-1,k,1)+
     +                         taur* ((wa1(i-2)*cc(m,i-1,k,
     +                         2)+wa1(i-1)*cc(m,i,k,2))+ (wa2(i-2)*cc(m,
     +                         i-1,k,3)+wa2(i-1)*cc(m,i,k,3)))) +
     +                         (taui* ((wa1(i-2)*cc(m,i,k,
     +                         2)-wa1(i-1)*cc(m,i-1,k,
     +                         2))- (wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*cc(m,
     +                         i-1,k,3))))
               ch(m,ic-1,2,k) = (cc(m,i-1,k,1)+
     +                          taur* ((wa1(i-2)*cc(m,i-1,k,
     +                          2)+wa1(i-1)*cc(m,i,k,
     +                          2))+ (wa2(i-2)*cc(m,i-1,k,
     +                          3)+wa2(i-1)*cc(m,i,k,3)))) -
     +                          (taui* ((wa1(i-2)*cc(m,i,k,
     +                          2)-wa1(i-1)*cc(m,i-1,k,
     +                          2))- (wa2(i-2)*cc(m,i,k,
     +                          3)-wa2(i-1)*cc(m,i-1,k,3))))
               ch(m,i,3,k) = (cc(m,i,k,1)+taur*
     +                       ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,
     +                       2))+ (wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*cc(m,
     +                       i-1,k,3)))) + (taui* ((wa2(i-2)*cc(m,i-1,k,
     +                       3)+wa2(i-1)*cc(m,i,k,3))- (wa1(i-2)*cc(m,
     +                       i-1,k,2)+wa1(i-1)*cc(m,i,k,2))))
               ch(m,ic,2,k) = (taui* ((wa2(i-2)*cc(m,i-1,k,
     +                        3)+wa2(i-1)*cc(m,i,k,3))- (wa1(i-2)*cc(m,
     +                        i-1,k,2)+wa1(i-1)*cc(m,i,k,2)))) -
     +                        (cc(m,i,k,1)+taur* ((wa1(i-2)*cc(m,i,k,
     +                        2)-wa1(i-1)*cc(m,i-1,k,
     +                        2))+ (wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*cc(m,
     +                        i-1,k,3))))
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      RETURN
      END
      SUBROUTINE vradf4(mp,ido,l1,cc,ch,mdimc,wa1,wa2,wa3)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION cc(mdimc,ido,l1,4),ch(mdimc,ido,4,l1),wa1(ido),wa2(ido),
     +          wa3(ido)

      hsqt2 = sqrt(2.)/2.
      DO 20 k = 1,l1
         DO 10 m = 1,mp
            ch(m,1,1,k) = (cc(m,1,k,2)+cc(m,1,k,4)) +
     +                    (cc(m,1,k,1)+cc(m,1,k,3))
            ch(m,ido,4,k) = (cc(m,1,k,1)+cc(m,1,k,3)) -
     +                      (cc(m,1,k,2)+cc(m,1,k,4))
            ch(m,ido,2,k) = cc(m,1,k,1) - cc(m,1,k,3)
            ch(m,1,3,k) = cc(m,1,k,4) - cc(m,1,k,2)
   10    CONTINUE
   20 CONTINUE
      IF ( ido ==  2 ) THEN
         GOTO 70
      ELSEIF ( ido < 2 ) THEN
         GOTO 100
      ENDIF
      idp2 = ido + 2
      DO 60 k = 1,l1
         DO 50 i = 3,ido,2
            ic = idp2 - i
            DO 40 m = 1,mp
               ch(m,i-1,1,k) = ((wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i,
     +                         k,2))+ (wa3(i-2)*cc(m,i-1,k,
     +                         4)+wa3(i-1)*cc(m,i,k,4))) +
     +                         (cc(m,i-1,k,1)+ (wa2(i-2)*cc(m,i-1,k,
     +                         3)+wa2(i-1)*cc(m,i,k,3)))
               ch(m,ic-1,4,k) = (cc(m,i-1,k,1)+
     +                          (wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i,
     +                          k,3))) - ((wa1(i-2)*cc(m,i-1,k,
     +                          2)+wa1(i-1)*cc(m,i,k,2))+
     +                          (wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i,
     +                          k,4)))
               ch(m,i,1,k) = ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,
     +                       2))+ (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,
     +                       i-1,k,4))) + (cc(m,i,k,1)+
     +                       (wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*cc(m,i-1,k,
     +                       3)))
               ch(m,ic,4,k) = ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,
     +                        k,2))+ (wa3(i-2)*cc(m,i,k,
     +                        4)-wa3(i-1)*cc(m,i-1,k,4))) -
     +                        (cc(m,i,k,1)+ (wa2(i-2)*cc(m,i,k,
     +                        3)-wa2(i-1)*cc(m,i-1,k,3)))
               ch(m,i-1,3,k) = ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,
     +                         k,2))- (wa3(i-2)*cc(m,i,k,
     +                         4)-wa3(i-1)*cc(m,i-1,k,4))) +
     +                         (cc(m,i-1,k,1)- (wa2(i-2)*cc(m,i-1,k,
     +                         3)+wa2(i-1)*cc(m,i,k,3)))
               ch(m,ic-1,2,k) = (cc(m,i-1,k,1)-
     +                          (wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i,
     +                          k,3))) - ((wa1(i-2)*cc(m,i,k,
     +                          2)-wa1(i-1)*cc(m,i-1,k,2))-
     +                          (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,i-1,
     +                          k,4)))
               ch(m,i,3,k) = ((wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i,k,
     +                       4))- (wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,
     +                       i,k,2))) + (cc(m,i,k,1)-
     +                       (wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*cc(m,i-1,k,
     +                       3)))
               ch(m,ic,2,k) = ((wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i,
     +                        k,4))- (wa1(i-2)*cc(m,i-1,k,
     +                        2)+wa1(i-1)*cc(m,i,k,2))) -
     +                        (cc(m,i,k,1)- (wa2(i-2)*cc(m,i,k,
     +                        3)-wa2(i-1)*cc(m,i-1,k,3)))
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
      IF (mod(ido,2).EQ.1) RETURN
   70 CONTINUE
      DO 90 k = 1,l1
         DO 80 m = 1,mp
            ch(m,ido,1,k) = (hsqt2* (cc(m,ido,k,2)-cc(m,ido,k,4))) +
     +                      cc(m,ido,k,1)
            ch(m,ido,3,k) = cc(m,ido,k,1) -
     +                      (hsqt2* (cc(m,ido,k,2)-cc(m,ido,k,4)))
            ch(m,1,2,k) = (-hsqt2* (cc(m,ido,k,2)+cc(m,ido,k,4))) -
     +                    cc(m,ido,k,3)
            ch(m,1,4,k) = (-hsqt2* (cc(m,ido,k,2)+cc(m,ido,k,4))) +
     +                    cc(m,ido,k,3)
   80    CONTINUE
   90 CONTINUE
  100 RETURN
      END
      SUBROUTINE vradf5(mp,ido,l1,cc,ch,mdimc,wa1,wa2,wa3,wa4)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      USE m_constants, ONLY : pimach
      DIMENSION cc(mdimc,ido,l1,5),ch(mdimc,ido,5,l1),wa1(ido),wa2(ido),
     +          wa3(ido),wa4(ido)

      arg = 2.*pimach()/5.
      tr11 = cos(arg)
      ti11 = sin(arg)
      tr12 = cos(2.*arg)
      ti12 = sin(2.*arg)
      DO 20 k = 1,l1
         DO 10 m = 1,mp
            ch(m,1,1,k) = cc(m,1,k,1) + (cc(m,1,k,5)+cc(m,1,k,2)) +
     +                    (cc(m,1,k,4)+cc(m,1,k,3))
            ch(m,ido,2,k) = cc(m,1,k,1) +
     +                      tr11* (cc(m,1,k,5)+cc(m,1,k,2)) +
     +                      tr12* (cc(m,1,k,4)+cc(m,1,k,3))
            ch(m,1,3,k) = ti11* (cc(m,1,k,5)-cc(m,1,k,2)) +
     +                    ti12* (cc(m,1,k,4)-cc(m,1,k,3))
            ch(m,ido,4,k) = cc(m,1,k,1) +
     +                      tr12* (cc(m,1,k,5)+cc(m,1,k,2)) +
     +                      tr11* (cc(m,1,k,4)+cc(m,1,k,3))
            ch(m,1,5,k) = ti12* (cc(m,1,k,5)-cc(m,1,k,2)) -
     +                    ti11* (cc(m,1,k,4)-cc(m,1,k,3))
   10    CONTINUE
   20 CONTINUE
      IF (ido.EQ.1) RETURN
      idp2 = ido + 2
      DO 50 k = 1,l1
         DO 40 i = 3,ido,2
            ic = idp2 - i
            DO 30 m = 1,mp
               ch(m,i-1,1,k) = cc(m,i-1,k,1) +
     +                         ((wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i,
     +                         k,2))+ (wa4(i-2)*cc(m,i-1,k,
     +                         5)+wa4(i-1)*cc(m,i,k,5))) +
     +                         ((wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i,
     +                         k,3))+ (wa3(i-2)*cc(m,i-1,k,
     +                         4)+wa3(i-1)*cc(m,i,k,4)))
               ch(m,i,1,k) = cc(m,i,k,1) +
     +                       ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,
     +                       2))+ (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m,
     +                       i-1,k,5))) + ((wa2(i-2)*cc(m,i,k,
     +                       3)-wa2(i-1)*cc(m,i-1,k,3))+
     +                       (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,i-1,k,
     +                       4)))
               ch(m,i-1,3,k) = cc(m,i-1,k,1) +
     +                         tr11* (wa1(i-2)*cc(m,i-1,k,2)+
     +                         wa1(i-1)*cc(m,i,k,2)+
     +                         wa4(i-2)*cc(m,i-1,k,5)+
     +                         wa4(i-1)*cc(m,i,k,5)) +
     +                         tr12* (wa2(i-2)*cc(m,i-1,k,3)+
     +                         wa2(i-1)*cc(m,i,k,3)+
     +                         wa3(i-2)*cc(m,i-1,k,4)+
     +                         wa3(i-1)*cc(m,i,k,4)) +
     +                         ti11* (wa1(i-2)*cc(m,i,k,2)-
     +                         wa1(i-1)*cc(m,i-1,k,2)-
     +                         (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m,i-1,
     +                         k,5))) + ti12* (wa2(i-2)*cc(m,i,k,3)-
     +                         wa2(i-1)*cc(m,i-1,k,3)-
     +                         (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,i-1,
     +                         k,4)))
               ch(m,ic-1,2,k) = cc(m,i-1,k,1) +
     +                          tr11* (wa1(i-2)*cc(m,i-1,k,2)+
     +                          wa1(i-1)*cc(m,i,k,2)+
     +                          wa4(i-2)*cc(m,i-1,k,5)+
     +                          wa4(i-1)*cc(m,i,k,5)) +
     +                          tr12* (wa2(i-2)*cc(m,i-1,k,3)+
     +                          wa2(i-1)*cc(m,i,k,3)+
     +                          wa3(i-2)*cc(m,i-1,k,4)+
     +                          wa3(i-1)*cc(m,i,k,4)) -
     +                          (ti11* (wa1(i-2)*cc(m,i,k,
     +                          2)-wa1(i-1)*cc(m,i-1,k,
     +                          2)- (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m,
     +                          i-1,k,5)))+ti12* (wa2(i-2)*cc(m,i,k,
     +                          3)-wa2(i-1)*cc(m,i-1,k,
     +                          3)- (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,
     +                          i-1,k,4))))
               ch(m,i,3,k) = (cc(m,i,k,1)+tr11*
     +                       ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,
     +                       2))+ (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m,
     +                       i-1,k,5)))+tr12* ((wa2(i-2)*cc(m,i,k,
     +                       3)-wa2(i-1)*cc(m,i-1,k,3))+ (wa3(i-2)*cc(m,
     +                       i,k,4)-wa3(i-1)*cc(m,i-1,k,4)))) +
     +                       (ti11* ((wa4(i-2)*cc(m,i-1,k,
     +                       5)+wa4(i-1)*cc(m,i,k,5))- (wa1(i-2)*cc(m,
     +                       i-1,k,2)+wa1(i-1)*cc(m,i,k,2)))+
     +                       ti12* ((wa3(i-2)*cc(m,i-1,k,
     +                       4)+wa3(i-1)*cc(m,i,k,4))- (wa2(i-2)*cc(m,
     +                       i-1,k,3)+wa2(i-1)*cc(m,i,k,3))))
               ch(m,ic,2,k) = (ti11* ((wa4(i-2)*cc(m,i-1,k,
     +                        5)+wa4(i-1)*cc(m,i,k,5))- (wa1(i-2)*cc(m,
     +                        i-1,k,2)+wa1(i-1)*cc(m,i,k,2)))+
     +                        ti12* ((wa3(i-2)*cc(m,i-1,k,
     +                        4)+wa3(i-1)*cc(m,i,k,4))- (wa2(i-2)*cc(m,
     +                        i-1,k,3)+wa2(i-1)*cc(m,i,k,3)))) -
     +                        (cc(m,i,k,1)+tr11* ((wa1(i-2)*cc(m,i,k,
     +                        2)-wa1(i-1)*cc(m,i-1,k,
     +                        2))+ (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m,
     +                        i-1,k,5)))+tr12* ((wa2(i-2)*cc(m,i,k,
     +                        3)-wa2(i-1)*cc(m,i-1,k,
     +                        3))+ (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,
     +                        i-1,k,4))))
               ch(m,i-1,5,k) = (cc(m,i-1,k,1)+
     +                         tr12* ((wa1(i-2)*cc(m,i-1,k,
     +                         2)+wa1(i-1)*cc(m,i,k,2))+ (wa4(i-2)*cc(m,
     +                         i-1,k,5)+wa4(i-1)*cc(m,i,k,5)))+
     +                         tr11* ((wa2(i-2)*cc(m,i-1,k,
     +                         3)+wa2(i-1)*cc(m,i,k,3))+ (wa3(i-2)*cc(m,
     +                         i-1,k,4)+wa3(i-1)*cc(m,i,k,4)))) +
     +                         (ti12* ((wa1(i-2)*cc(m,i,k,
     +                         2)-wa1(i-1)*cc(m,i-1,k,
     +                         2))- (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m,
     +                         i-1,k,5)))-ti11* ((wa2(i-2)*cc(m,i,k,
     +                         3)-wa2(i-1)*cc(m,i-1,k,
     +                         3))- (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,
     +                         i-1,k,4))))
               ch(m,ic-1,4,k) = (cc(m,i-1,k,1)+
     +                          tr12* ((wa1(i-2)*cc(m,i-1,k,
     +                          2)+wa1(i-1)*cc(m,i,k,
     +                          2))+ (wa4(i-2)*cc(m,i-1,k,
     +                          5)+wa4(i-1)*cc(m,i,k,5)))+
     +                          tr11* ((wa2(i-2)*cc(m,i-1,k,
     +                          3)+wa2(i-1)*cc(m,i,k,
     +                          3))+ (wa3(i-2)*cc(m,i-1,k,
     +                          4)+wa3(i-1)*cc(m,i,k,4)))) -
     +                          (ti12* ((wa1(i-2)*cc(m,i,k,
     +                          2)-wa1(i-1)*cc(m,i-1,k,
     +                          2))- (wa4(i-2)*cc(m,i,k,
     +                          5)-wa4(i-1)*cc(m,i-1,k,5)))-
     +                          ti11* ((wa2(i-2)*cc(m,i,k,
     +                          3)-wa2(i-1)*cc(m,i-1,k,
     +                          3))- (wa3(i-2)*cc(m,i,k,
     +                          4)-wa3(i-1)*cc(m,i-1,k,4))))
               ch(m,i,5,k) = (cc(m,i,k,1)+tr12*
     +                       ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,
     +                       2))+ (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m,
     +                       i-1,k,5)))+tr11* ((wa2(i-2)*cc(m,i,k,
     +                       3)-wa2(i-1)*cc(m,i-1,k,3))+ (wa3(i-2)*cc(m,
     +                       i,k,4)-wa3(i-1)*cc(m,i-1,k,4)))) +
     +                       (ti12* ((wa4(i-2)*cc(m,i-1,k,
     +                       5)+wa4(i-1)*cc(m,i,k,5))- (wa1(i-2)*cc(m,
     +                       i-1,k,2)+wa1(i-1)*cc(m,i,k,2)))-
     +                       ti11* ((wa3(i-2)*cc(m,i-1,k,
     +                       4)+wa3(i-1)*cc(m,i,k,4))- (wa2(i-2)*cc(m,
     +                       i-1,k,3)+wa2(i-1)*cc(m,i,k,3))))
               ch(m,ic,4,k) = (ti12* ((wa4(i-2)*cc(m,i-1,k,
     +                        5)+wa4(i-1)*cc(m,i,k,5))- (wa1(i-2)*cc(m,
     +                        i-1,k,2)+wa1(i-1)*cc(m,i,k,2)))-
     +                        ti11* ((wa3(i-2)*cc(m,i-1,k,
     +                        4)+wa3(i-1)*cc(m,i,k,4))- (wa2(i-2)*cc(m,
     +                        i-1,k,3)+wa2(i-1)*cc(m,i,k,3)))) -
     +                        (cc(m,i,k,1)+tr12* ((wa1(i-2)*cc(m,i,k,
     +                        2)-wa1(i-1)*cc(m,i-1,k,
     +                        2))+ (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m,
     +                        i-1,k,5)))+tr11* ((wa2(i-2)*cc(m,i,k,
     +                        3)-wa2(i-1)*cc(m,i-1,k,
     +                        3))+ (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,
     +                        i-1,k,4))))
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      RETURN
      END
      SUBROUTINE vradfg(mp,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,mdimc,wa)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      USE m_constants, ONLY : pimach
      DIMENSION ch(mdimc,ido,l1,ip),cc(mdimc,ido,ip,l1),
     +          c1(mdimc,ido,l1,ip),c2(mdimc,idl1,ip),
     +          ch2(mdimc,idl1,ip),wa(ido)

      tpi = 2.*pimach()
      arg = tpi/real(ip)
      dcp = cos(arg)
      dsp = sin(arg)
      ipph = (ip+1)/2
      ipp2 = ip + 2
      idp2 = ido + 2
      nbd = (ido-1)/2
      IF (ido.EQ.1) GO TO 250
      DO 20 ik = 1,idl1
         DO 10 m = 1,mp
            ch2(m,ik,1) = c2(m,ik,1)
   10    CONTINUE
   20 CONTINUE
      DO 50 j = 2,ip
         DO 40 k = 1,l1
            DO 30 m = 1,mp
               ch(m,1,k,j) = c1(m,1,k,j)
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      IF (nbd.GT.l1) GO TO 100
      is = -ido
      DO 90 j = 2,ip
         is = is + ido
         idij = is
         DO 80 i = 3,ido,2
            idij = idij + 2
            DO 70 k = 1,l1
               DO 60 m = 1,mp
                  ch(m,i-1,k,j) = wa(idij-1)*c1(m,i-1,k,j) +
     +                            wa(idij)*c1(m,i,k,j)
                  ch(m,i,k,j) = wa(idij-1)*c1(m,i,k,j) -
     +                          wa(idij)*c1(m,i-1,k,j)
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
      GO TO 150
  100 is = -ido
      DO 140 j = 2,ip
         is = is + ido
         DO 130 k = 1,l1
            idij = is
            DO 120 i = 3,ido,2
               idij = idij + 2
               DO 110 m = 1,mp
                  ch(m,i-1,k,j) = wa(idij-1)*c1(m,i-1,k,j) +
     +                            wa(idij)*c1(m,i,k,j)
                  ch(m,i,k,j) = wa(idij-1)*c1(m,i,k,j) -
     +                          wa(idij)*c1(m,i-1,k,j)
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
  150 IF (nbd.LT.l1) GO TO 200
      DO 190 j = 2,ipph
         jc = ipp2 - j
         DO 180 k = 1,l1
            DO 170 i = 3,ido,2
               DO 160 m = 1,mp
                  c1(m,i-1,k,j) = ch(m,i-1,k,j) + ch(m,i-1,k,jc)
                  c1(m,i-1,k,jc) = ch(m,i,k,j) - ch(m,i,k,jc)
                  c1(m,i,k,j) = ch(m,i,k,j) + ch(m,i,k,jc)
                  c1(m,i,k,jc) = ch(m,i-1,k,jc) - ch(m,i-1,k,j)
  160          CONTINUE
  170       CONTINUE
  180    CONTINUE
  190 CONTINUE
      GO TO 280
  200 DO 240 j = 2,ipph
         jc = ipp2 - j
         DO 230 i = 3,ido,2
            DO 220 k = 1,l1
               DO 210 m = 1,mp
                  c1(m,i-1,k,j) = ch(m,i-1,k,j) + ch(m,i-1,k,jc)
                  c1(m,i-1,k,jc) = ch(m,i,k,j) - ch(m,i,k,jc)
                  c1(m,i,k,j) = ch(m,i,k,j) + ch(m,i,k,jc)
                  c1(m,i,k,jc) = ch(m,i-1,k,jc) - ch(m,i-1,k,j)
  210          CONTINUE
  220       CONTINUE
  230    CONTINUE
  240 CONTINUE
      GO TO 280
  250 DO 270 ik = 1,idl1
         DO 260 m = 1,mp
            c2(m,ik,1) = ch2(m,ik,1)
  260    CONTINUE
  270 CONTINUE
  280 DO 310 j = 2,ipph
         jc = ipp2 - j
         DO 300 k = 1,l1
            DO 290 m = 1,mp
               c1(m,1,k,j) = ch(m,1,k,j) + ch(m,1,k,jc)
               c1(m,1,k,jc) = ch(m,1,k,jc) - ch(m,1,k,j)
  290       CONTINUE
  300    CONTINUE
  310 CONTINUE
C
      ar1 = 1.
      ai1 = 0.
      DO 370 l = 2,ipph
         lc = ipp2 - l
         ar1h = dcp*ar1 - dsp*ai1
         ai1 = dcp*ai1 + dsp*ar1
         ar1 = ar1h
         DO 330 ik = 1,idl1
            DO 320 m = 1,mp
               ch2(m,ik,l) = c2(m,ik,1) + ar1*c2(m,ik,2)
               ch2(m,ik,lc) = ai1*c2(m,ik,ip)
  320       CONTINUE
  330    CONTINUE
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         DO 360 j = 3,ipph
            jc = ipp2 - j
            ar2h = dc2*ar2 - ds2*ai2
            ai2 = dc2*ai2 + ds2*ar2
            ar2 = ar2h
            DO 350 ik = 1,idl1
               DO 340 m = 1,mp
                  ch2(m,ik,l) = ch2(m,ik,l) + ar2*c2(m,ik,j)
                  ch2(m,ik,lc) = ch2(m,ik,lc) + ai2*c2(m,ik,jc)
  340          CONTINUE
  350       CONTINUE
  360    CONTINUE
  370 CONTINUE
      DO 400 j = 2,ipph
         DO 390 ik = 1,idl1
            DO 380 m = 1,mp
               ch2(m,ik,1) = ch2(m,ik,1) + c2(m,ik,j)
  380       CONTINUE
  390    CONTINUE
  400 CONTINUE
C
      IF (ido.LT.l1) GO TO 440
      DO 430 k = 1,l1
         DO 420 i = 1,ido
            DO 410 m = 1,mp
               cc(m,i,1,k) = ch(m,i,k,1)
  410       CONTINUE
  420    CONTINUE
  430 CONTINUE
      GO TO 480
  440 DO 470 i = 1,ido
         DO 460 k = 1,l1
            DO 450 m = 1,mp
               cc(m,i,1,k) = ch(m,i,k,1)
  450       CONTINUE
  460    CONTINUE
  470 CONTINUE
  480 DO 510 j = 2,ipph
         jc = ipp2 - j
         j2 = j + j
         DO 500 k = 1,l1
            DO 490 m = 1,mp
               cc(m,ido,j2-2,k) = ch(m,1,k,j)
               cc(m,1,j2-1,k) = ch(m,1,k,jc)
  490       CONTINUE
  500    CONTINUE
  510 CONTINUE
      IF (ido.EQ.1) RETURN
      IF (nbd.LT.l1) GO TO 560
      DO 550 j = 2,ipph
         jc = ipp2 - j
         j2 = j + j
         DO 540 k = 1,l1
            DO 530 i = 3,ido,2
               ic = idp2 - i
               DO 520 m = 1,mp
                  cc(m,i-1,j2-1,k) = ch(m,i-1,k,j) + ch(m,i-1,k,jc)
                  cc(m,ic-1,j2-2,k) = ch(m,i-1,k,j) - ch(m,i-1,k,jc)
                  cc(m,i,j2-1,k) = ch(m,i,k,j) + ch(m,i,k,jc)
                  cc(m,ic,j2-2,k) = ch(m,i,k,jc) - ch(m,i,k,j)
  520          CONTINUE
  530       CONTINUE
  540    CONTINUE
  550 CONTINUE
      RETURN
  560 DO 600 j = 2,ipph
         jc = ipp2 - j
         j2 = j + j
         DO 590 i = 3,ido,2
            ic = idp2 - i
            DO 580 k = 1,l1
               DO 570 m = 1,mp
                  cc(m,i-1,j2-1,k) = ch(m,i-1,k,j) + ch(m,i-1,k,jc)
                  cc(m,ic-1,j2-2,k) = ch(m,i-1,k,j) - ch(m,i-1,k,jc)
                  cc(m,i,j2-1,k) = ch(m,i,k,j) + ch(m,i,k,jc)
                  cc(m,ic,j2-2,k) = ch(m,i,k,jc) - ch(m,i,k,j)
  570          CONTINUE
  580       CONTINUE
  590    CONTINUE
  600 CONTINUE
      RETURN
      END
      SUBROUTINE vrftf1(m,n,c,ch,mdimc,wa,fac)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION ch(mdimc,n),c(mdimc,n),wa(n),fac(15)

      nf = fac(2)
      na = 1
      l2 = n
      iw = n
      DO 110 k1 = 1,nf
         kh = nf - k1
         ip = fac(kh+3)
         l1 = l2/ip
         ido = n/l2
         idl1 = ido*l1
         iw = iw - (ip-1)*ido
         na = 1 - na
         IF (ip.NE.4) GO TO 20
         ix2 = iw + ido
         ix3 = ix2 + ido
         IF (na.NE.0) GO TO 10
         CALL vradf4(m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2),wa(ix3))
         GO TO 100
   10    CALL vradf4(m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2),wa(ix3))
         GO TO 100
   20    IF (ip.NE.2) GO TO 40
         IF (na.NE.0) GO TO 30
         CALL vradf2(m,ido,l1,c,ch,mdimc,wa(iw))
         GO TO 100
   30    CALL vradf2(m,ido,l1,ch,c,mdimc,wa(iw))
         GO TO 100
   40    IF (ip.NE.3) GO TO 60
         ix2 = iw + ido
         IF (na.NE.0) GO TO 50
         CALL vradf3(m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2))
         GO TO 100
   50    CALL vradf3(m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2))
         GO TO 100
   60    IF (ip.NE.5) GO TO 80
         ix2 = iw + ido
         ix3 = ix2 + ido
         ix4 = ix3 + ido
         IF (na.NE.0) GO TO 70
         CALL vradf5(m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         GO TO 100
   70    CALL vradf5(m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         GO TO 100
   80    IF (ido.EQ.1) na = 1 - na
         IF (na.NE.0) GO TO 90
         CALL vradfg(m,ido,ip,l1,idl1,c,c,c,ch,ch,mdimc,wa(iw))
         na = 1
         GO TO 100
   90    CALL vradfg(m,ido,ip,l1,idl1,ch,ch,ch,c,c,mdimc,wa(iw))
         na = 0
  100    l2 = l1
  110 CONTINUE
      scale = sqrt(1./n)
      IF (na.EQ.1) GO TO 140
      DO 130 j = 1,n
         DO 120 i = 1,m
            c(i,j) = scale*ch(i,j)
  120    CONTINUE
  130 CONTINUE
      RETURN
  140 DO 160 j = 1,n
         DO 150 i = 1,m
            c(i,j) = scale*c(i,j)
  150    CONTINUE
  160 CONTINUE
      RETURN
      END

      SUBROUTINE vrfftb(m,n,r,rt,mdimr,wsave)
C***BEGIN PROLOGUE  VRFFTB
C***DATE WRITTEN   850801   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM,
C             FOURIER SYNTHESIS, BACKWARD TRANSFORM, MULTIPLE SEQUENCES
C***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST)
C***PURPOSE  Backward real periodic transform, M sequences.
C***DESCRIPTION
C
C  Subroutine VRFFTB computes the synthesis (backward transform) of a
C  number of real periodic sequences from their Fourier coefficients.
C  Specifically, for each set of independent Fourier coefficients
C  F(K), the corresponding real periodic sequence is computed.
C
C  The array WSAVE which is used by subroutine VRFFTB must be
C  initialized by calling subroutine VRFFTI(N,WSAVE).
C
C
C  Input Parameters
C
C  M       the number of sets of coefficients.
C
C  N       the length of the sequences of coefficients to be
C          transformed.  The method is most efficient when N is a
C          product of small primes, however n may be any positive
C          integer.
C
C  R       areal two-dimensional array of size MDIMX x N containing the
C          coefficients to be transformed.  Each set of coefficients
C          F(K), K\0,1,..,N-1, is stored as a ROW of R.  Specifically,
C          the I-th set of independent Fourier coefficients is stored
C
C                R(I,1) = REAL( F(I,0) ),
C
C                R(I,2*K) = REAL( F(I,K) )
C
C                R(I,2*K+1) = IMAG( F(I,K) )
C
C                   for K = 1, 2, . . . , M-1,
C
C                and, when N is even,
C
C                R(I,N) = REAL( F(I,N/2) ).
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
C          same WSAVE array may be used by VRFFTB and VRFFTB.
C
C  Output Parameters
C
C  R       contains M real periodic sequences corresponding to the given
C          coefficients.  Specifically, the I-th row of R contains the
C          real periodic sequence corresponding to the I-th set of
C          independent Fourier coefficients F(I,K) stored as
C
C               R(I,J) = X(I,J-1) ,   J = 1, 2, . . . , N, where
C
C               X(I,J) = SQRT(1/N)* F(I,0) + (-1)**J*F(I,N/2)
C                        + 2*SUM(K=1,M)[ REAL(F(I,2K))*COS(2K*J*PI/N)
C                        - IMAG(F(I,2K+1))*SIN(2K*J*PI/N) ]  ,
C
C                 when N is even, and
C
C               X(I,J) = SQRT(1/N)* F(I,0) +
C                        2*SUM(K=1,M)[ REAL(F(I,2K))*COS(2K*J*PI/N)
C                        - IMAG(F(I,2K+1))*SIN(2K*J*PI/N) ]  ,
C
C                 when N is odd.
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
C  VRFFTB is a straightforward extension of the subprogram RFFTB to
C  handle M simultaneous sequences.  RFFTB was originally developed
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
C     required        VRADF4, VRADF5, VRADFG, VRFFTB, VRFTB1, VRADB2,
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
C***ROUTINES CALLED  VRFTB1
C***END PROLOGUE  VRFFTB
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION r(mdimr,n),rt(mdimr,n),wsave(n+15)

      IF (n.EQ.1) RETURN
      CALL vrftb1(m,n,r,rt,mdimr,wsave(1),wsave(n+1))
      RETURN
      END
      SUBROUTINE vradb2(mp,ido,l1,cc,ch,mdimc,wa1)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION cc(mdimc,ido,2,l1),ch(mdimc,ido,l1,2),wa1(ido)

      DO 20 k = 1,l1
         DO 10 m = 1,mp
            ch(m,1,k,1) = cc(m,1,1,k) + cc(m,ido,2,k)
            ch(m,1,k,2) = cc(m,1,1,k) - cc(m,ido,2,k)
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
               ch(m,i-1,k,1) = cc(m,i-1,1,k) + cc(m,ic-1,2,k)
               ch(m,i,k,1) = cc(m,i,1,k) - cc(m,ic,2,k)
               ch(m,i-1,k,2) = wa1(i-2)* (cc(m,i-1,1,k)-
     +                         cc(m,ic-1,2,k)) - wa1(i-1)*
     +                         (cc(m,i,1,k)+cc(m,ic,2,k))
               ch(m,i,k,2) = wa1(i-2)* (cc(m,i,1,k)+cc(m,ic,2,k)) +
     +                       wa1(i-1)* (cc(m,i-1,1,k)-cc(m,ic-1,2,k))
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
      IF (mod(ido,2).EQ.1) RETURN
   70 DO 90 k = 1,l1
         DO 80 m = 1,mp
            ch(m,ido,k,1) = cc(m,ido,1,k) + cc(m,ido,1,k)
            ch(m,ido,k,2) = - (cc(m,1,2,k)+cc(m,1,2,k))
   80    CONTINUE
   90 CONTINUE
  100 RETURN
      END
      SUBROUTINE vradb3(mp,ido,l1,cc,ch,mdimc,wa1,wa2)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      USE m_constants, ONLY : pimach
      DIMENSION cc(mdimc,ido,3,l1),ch(mdimc,ido,l1,3),wa1(ido),wa2(ido)

      arg = 2.*pimach()/3.
      taur = cos(arg)
      taui = sin(arg)
      DO 20 k = 1,l1
         DO 10 m = 1,mp
            ch(m,1,k,1) = cc(m,1,1,k) + 2.*cc(m,ido,2,k)
            ch(m,1,k,2) = cc(m,1,1,k) + (2.*taur)*cc(m,ido,2,k) -
     +                    (2.*taui)*cc(m,1,3,k)
            ch(m,1,k,3) = cc(m,1,1,k) + (2.*taur)*cc(m,ido,2,k) +
     +                    2.*taui*cc(m,1,3,k)
   10    CONTINUE
   20 CONTINUE
      IF (ido.EQ.1) RETURN
      idp2 = ido + 2
      DO 50 k = 1,l1
         DO 40 i = 3,ido,2
            ic = idp2 - i
            DO 30 m = 1,mp
               ch(m,i-1,k,1) = cc(m,i-1,1,k) +
     +                         (cc(m,i-1,3,k)+cc(m,ic-1,2,k))
               ch(m,i,k,1) = cc(m,i,1,k) + (cc(m,i,3,k)-cc(m,ic,2,k))
               ch(m,i-1,k,2) = wa1(i-2)* ((cc(m,i-1,1,k)+taur* (cc(m,
     +                         i-1,3,k)+cc(m,ic-1,2,k)))-
     +                         (taui* (cc(m,i,3,k)+cc(m,ic,2,k)))) -
     +                         wa1(i-1)* ((cc(m,i,1,k)+taur* (cc(m,i,3,
     +                         k)-cc(m,ic,2,k)))+ (taui* (cc(m,i-1,3,
     +                         k)-cc(m,ic-1,2,k))))
               ch(m,i,k,2) = wa1(i-2)* ((cc(m,i,1,k)+taur* (cc(m,i,3,
     +                       k)-cc(m,ic,2,k)))+ (taui* (cc(m,i-1,3,
     +                       k)-cc(m,ic-1,2,k)))) +
     +                       wa1(i-1)* ((cc(m,i-1,1,k)+taur* (cc(m,i-1,
     +                       3,k)+cc(m,ic-1,2,k)))-
     +                       (taui* (cc(m,i,3,k)+cc(m,ic,2,k))))
               ch(m,i-1,k,3) = wa2(i-2)* ((cc(m,i-1,1,k)+taur* (cc(m,
     +                         i-1,3,k)+cc(m,ic-1,2,k)))+
     +                         (taui* (cc(m,i,3,k)+cc(m,ic,2,k)))) -
     +                         wa2(i-1)* ((cc(m,i,1,k)+taur* (cc(m,i,3,
     +                         k)-cc(m,ic,2,k)))- (taui* (cc(m,i-1,3,
     +                         k)-cc(m,ic-1,2,k))))
               ch(m,i,k,3) = wa2(i-2)* ((cc(m,i,1,k)+taur* (cc(m,i,3,
     +                       k)-cc(m,ic,2,k)))- (taui* (cc(m,i-1,3,
     +                       k)-cc(m,ic-1,2,k)))) +
     +                       wa2(i-1)* ((cc(m,i-1,1,k)+taur* (cc(m,i-1,
     +                       3,k)+cc(m,ic-1,2,k)))+
     +                       (taui* (cc(m,i,3,k)+cc(m,ic,2,k))))
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      RETURN
      END
      SUBROUTINE vradb4(mp,ido,l1,cc,ch,mdimc,wa1,wa2,wa3)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION cc(mdimc,ido,4,l1),ch(mdimc,ido,l1,4),wa1(ido),wa2(ido),
     +          wa3(ido)

      sqrt2 = sqrt(2.)
      DO 20 k = 1,l1
         DO 10 m = 1,mp
            ch(m,1,k,3) = (cc(m,1,1,k)+cc(m,ido,4,k)) -
     +                    (cc(m,ido,2,k)+cc(m,ido,2,k))
            ch(m,1,k,1) = (cc(m,1,1,k)+cc(m,ido,4,k)) +
     +                    (cc(m,ido,2,k)+cc(m,ido,2,k))
            ch(m,1,k,4) = (cc(m,1,1,k)-cc(m,ido,4,k)) +
     +                    (cc(m,1,3,k)+cc(m,1,3,k))
            ch(m,1,k,2) = (cc(m,1,1,k)-cc(m,ido,4,k)) -
     +                    (cc(m,1,3,k)+cc(m,1,3,k))
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
               ch(m,i-1,k,1) = (cc(m,i-1,1,k)+cc(m,ic-1,4,k)) +
     +                         (cc(m,i-1,3,k)+cc(m,ic-1,2,k))
               ch(m,i,k,1) = (cc(m,i,1,k)-cc(m,ic,4,k)) +
     +                       (cc(m,i,3,k)-cc(m,ic,2,k))
               ch(m,i-1,k,2) = wa1(i-2)* ((cc(m,i-1,1,k)-cc(m,ic-1,4,
     +                         k))- (cc(m,i,3,k)+cc(m,ic,2,k))) -
     +                         wa1(i-1)* ((cc(m,i,1,k)+cc(m,ic,4,k))+
     +                         (cc(m,i-1,3,k)-cc(m,ic-1,2,k)))
               ch(m,i,k,2) = wa1(i-2)* ((cc(m,i,1,k)+cc(m,ic,4,k))+
     +                       (cc(m,i-1,3,k)-cc(m,ic-1,2,k))) +
     +                       wa1(i-1)* ((cc(m,i-1,1,k)-cc(m,ic-1,4,k))-
     +                       (cc(m,i,3,k)+cc(m,ic,2,k)))
               ch(m,i-1,k,3) = wa2(i-2)* ((cc(m,i-1,1,k)+cc(m,ic-1,4,
     +                         k))- (cc(m,i-1,3,k)+cc(m,ic-1,2,k))) -
     +                         wa2(i-1)* ((cc(m,i,1,k)-cc(m,ic,4,k))-
     +                         (cc(m,i,3,k)-cc(m,ic,2,k)))
               ch(m,i,k,3) = wa2(i-2)* ((cc(m,i,1,k)-cc(m,ic,4,k))-
     +                       (cc(m,i,3,k)-cc(m,ic,2,k))) +
     +                       wa2(i-1)* ((cc(m,i-1,1,k)+cc(m,ic-1,4,k))-
     +                       (cc(m,i-1,3,k)+cc(m,ic-1,2,k)))
               ch(m,i-1,k,4) = wa3(i-2)* ((cc(m,i-1,1,k)-cc(m,ic-1,4,
     +                         k))+ (cc(m,i,3,k)+cc(m,ic,2,k))) -
     +                         wa3(i-1)* ((cc(m,i,1,k)+cc(m,ic,4,k))-
     +                         (cc(m,i-1,3,k)-cc(m,ic-1,2,k)))
               ch(m,i,k,4) = wa3(i-2)* ((cc(m,i,1,k)+cc(m,ic,4,k))-
     +                       (cc(m,i-1,3,k)-cc(m,ic-1,2,k))) +
     +                       wa3(i-1)* ((cc(m,i-1,1,k)-cc(m,ic-1,4,k))+
     +                       (cc(m,i,3,k)+cc(m,ic,2,k)))
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
      IF (mod(ido,2).EQ.1) RETURN
   70 CONTINUE
      DO 90 k = 1,l1
         DO 80 m = 1,mp
            ch(m,ido,k,1) = (cc(m,ido,1,k)+cc(m,ido,3,k)) +
     +                      (cc(m,ido,1,k)+cc(m,ido,3,k))
            ch(m,ido,k,2) = sqrt2* ((cc(m,ido,1,k)-cc(m,ido,3,k))-
     +                      (cc(m,1,2,k)+cc(m,1,4,k)))
            ch(m,ido,k,3) = (cc(m,1,4,k)-cc(m,1,2,k)) +
     +                      (cc(m,1,4,k)-cc(m,1,2,k))
            ch(m,ido,k,4) = -sqrt2* ((cc(m,ido,1,k)-cc(m,ido,3,k))+
     +                      (cc(m,1,2,k)+cc(m,1,4,k)))
   80    CONTINUE
   90 CONTINUE
  100 RETURN
      END
      SUBROUTINE vradb5(mp,ido,l1,cc,ch,mdimc,wa1,wa2,wa3,wa4)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      USE m_constants, ONLY : pimach
      DIMENSION cc(mdimc,ido,5,l1),ch(mdimc,ido,l1,5),wa1(ido),wa2(ido),
     +          wa3(ido),wa4(ido)

      arg = 2.*pimach()/5.
      tr11 = cos(arg)
      ti11 = sin(arg)
      tr12 = cos(2.*arg)
      ti12 = sin(2.*arg)
      DO 20 k = 1,l1
         DO 10 m = 1,mp
            ch(m,1,k,1) = cc(m,1,1,k) + 2.*cc(m,ido,2,k) +
     +                    2.*cc(m,ido,4,k)
            ch(m,1,k,2) = (cc(m,1,1,k)+tr11*2.*cc(m,ido,2,k)+
     +                    tr12*2.*cc(m,ido,4,k)) -
     +                    (ti11*2.*cc(m,1,3,k)+ti12*2.*cc(m,1,5,k))
            ch(m,1,k,3) = (cc(m,1,1,k)+tr12*2.*cc(m,ido,2,k)+
     +                    tr11*2.*cc(m,ido,4,k)) -
     +                    (ti12*2.*cc(m,1,3,k)-ti11*2.*cc(m,1,5,k))
            ch(m,1,k,4) = (cc(m,1,1,k)+tr12*2.*cc(m,ido,2,k)+
     +                    tr11*2.*cc(m,ido,4,k)) +
     +                    (ti12*2.*cc(m,1,3,k)-ti11*2.*cc(m,1,5,k))
            ch(m,1,k,5) = (cc(m,1,1,k)+tr11*2.*cc(m,ido,2,k)+
     +                    tr12*2.*cc(m,ido,4,k)) +
     +                    (ti11*2.*cc(m,1,3,k)+ti12*2.*cc(m,1,5,k))
   10    CONTINUE
   20 CONTINUE
      IF (ido.EQ.1) RETURN
      idp2 = ido + 2
      DO 50 k = 1,l1
         DO 40 i = 3,ido,2
            ic = idp2 - i
            DO 30 m = 1,mp
               ch(m,i-1,k,1) = cc(m,i-1,1,k) +
     +                         (cc(m,i-1,3,k)+cc(m,ic-1,2,k)) +
     +                         (cc(m,i-1,5,k)+cc(m,ic-1,4,k))
               ch(m,i,k,1) = cc(m,i,1,k) + (cc(m,i,3,k)-cc(m,ic,2,k)) +
     +                       (cc(m,i,5,k)-cc(m,ic,4,k))
               ch(m,i-1,k,2) = wa1(i-2)* ((cc(m,i-1,1,k)+tr11* (cc(m,
     +                         i-1,3,k)+cc(m,ic-1,2,k))+tr12* (cc(m,i-1,
     +                         5,k)+cc(m,ic-1,4,k)))-
     +                         (ti11* (cc(m,i,3,k)+cc(m,ic,2,
     +                         k))+ti12* (cc(m,i,5,k)+cc(m,ic,4,k)))) -
     +                         wa1(i-1)* ((cc(m,i,1,k)+tr11* (cc(m,i,3,
     +                         k)-cc(m,ic,2,k))+tr12* (cc(m,i,5,k)-cc(m,
     +                         ic,4,k)))+ (ti11* (cc(m,i-1,3,k)-cc(m,
     +                         ic-1,2,k))+ti12* (cc(m,i-1,5,k)-cc(m,
     +                         ic-1,4,k))))
               ch(m,i,k,2) = wa1(i-2)* ((cc(m,i,1,k)+tr11* (cc(m,i,3,
     +                       k)-cc(m,ic,2,k))+tr12* (cc(m,i,5,k)-cc(m,
     +                       ic,4,k)))+ (ti11* (cc(m,i-1,3,k)-cc(m,ic-1,
     +                       2,k))+ti12* (cc(m,i-1,5,k)-cc(m,ic-1,4,
     +                       k)))) + wa1(i-1)* ((cc(m,i-1,1,
     +                       k)+tr11* (cc(m,i-1,3,k)+cc(m,ic-1,2,
     +                       k))+tr12* (cc(m,i-1,5,k)+cc(m,ic-1,4,k)))-
     +                       (ti11* (cc(m,i,3,k)+cc(m,ic,2,
     +                       k))+ti12* (cc(m,i,5,k)+cc(m,ic,4,k))))
               ch(m,i-1,k,3) = wa2(i-2)* ((cc(m,i-1,1,k)+tr12* (cc(m,
     +                         i-1,3,k)+cc(m,ic-1,2,k))+tr11* (cc(m,i-1,
     +                         5,k)+cc(m,ic-1,4,k)))-
     +                         (ti12* (cc(m,i,3,k)+cc(m,ic,2,
     +                         k))-ti11* (cc(m,i,5,k)+cc(m,ic,4,k)))) -
     +                         wa2(i-1)* ((cc(m,i,1,k)+tr12* (cc(m,i,3,
     +                         k)-cc(m,ic,2,k))+tr11* (cc(m,i,5,k)-cc(m,
     +                         ic,4,k)))+ (ti12* (cc(m,i-1,3,k)-cc(m,
     +                         ic-1,2,k))-ti11* (cc(m,i-1,5,k)-cc(m,
     +                         ic-1,4,k))))
               ch(m,i,k,3) = wa2(i-2)* ((cc(m,i,1,k)+tr12* (cc(m,i,3,
     +                       k)-cc(m,ic,2,k))+tr11* (cc(m,i,5,k)-cc(m,
     +                       ic,4,k)))+ (ti12* (cc(m,i-1,3,k)-cc(m,ic-1,
     +                       2,k))-ti11* (cc(m,i-1,5,k)-cc(m,ic-1,4,
     +                       k)))) + wa2(i-1)* ((cc(m,i-1,1,
     +                       k)+tr12* (cc(m,i-1,3,k)+cc(m,ic-1,2,
     +                       k))+tr11* (cc(m,i-1,5,k)+cc(m,ic-1,4,k)))-
     +                       (ti12* (cc(m,i,3,k)+cc(m,ic,2,
     +                       k))-ti11* (cc(m,i,5,k)+cc(m,ic,4,k))))
               ch(m,i-1,k,4) = wa3(i-2)* ((cc(m,i-1,1,k)+tr12* (cc(m,
     +                         i-1,3,k)+cc(m,ic-1,2,k))+tr11* (cc(m,i-1,
     +                         5,k)+cc(m,ic-1,4,k)))+
     +                         (ti12* (cc(m,i,3,k)+cc(m,ic,2,
     +                         k))-ti11* (cc(m,i,5,k)+cc(m,ic,4,k)))) -
     +                         wa3(i-1)* ((cc(m,i,1,k)+tr12* (cc(m,i,3,
     +                         k)-cc(m,ic,2,k))+tr11* (cc(m,i,5,k)-cc(m,
     +                         ic,4,k)))- (ti12* (cc(m,i-1,3,k)-cc(m,
     +                         ic-1,2,k))-ti11* (cc(m,i-1,5,k)-cc(m,
     +                         ic-1,4,k))))
               ch(m,i,k,4) = wa3(i-2)* ((cc(m,i,1,k)+tr12* (cc(m,i,3,
     +                       k)-cc(m,ic,2,k))+tr11* (cc(m,i,5,k)-cc(m,
     +                       ic,4,k)))- (ti12* (cc(m,i-1,3,k)-cc(m,ic-1,
     +                       2,k))-ti11* (cc(m,i-1,5,k)-cc(m,ic-1,4,
     +                       k)))) + wa3(i-1)* ((cc(m,i-1,1,
     +                       k)+tr12* (cc(m,i-1,3,k)+cc(m,ic-1,2,
     +                       k))+tr11* (cc(m,i-1,5,k)+cc(m,ic-1,4,k)))+
     +                       (ti12* (cc(m,i,3,k)+cc(m,ic,2,
     +                       k))-ti11* (cc(m,i,5,k)+cc(m,ic,4,k))))
               ch(m,i-1,k,5) = wa4(i-2)* ((cc(m,i-1,1,k)+tr11* (cc(m,
     +                         i-1,3,k)+cc(m,ic-1,2,k))+tr12* (cc(m,i-1,
     +                         5,k)+cc(m,ic-1,4,k)))+
     +                         (ti11* (cc(m,i,3,k)+cc(m,ic,2,
     +                         k))+ti12* (cc(m,i,5,k)+cc(m,ic,4,k)))) -
     +                         wa4(i-1)* ((cc(m,i,1,k)+tr11* (cc(m,i,3,
     +                         k)-cc(m,ic,2,k))+tr12* (cc(m,i,5,k)-cc(m,
     +                         ic,4,k)))- (ti11* (cc(m,i-1,3,k)-cc(m,
     +                         ic-1,2,k))+ti12* (cc(m,i-1,5,k)-cc(m,
     +                         ic-1,4,k))))
               ch(m,i,k,5) = wa4(i-2)* ((cc(m,i,1,k)+tr11* (cc(m,i,3,
     +                       k)-cc(m,ic,2,k))+tr12* (cc(m,i,5,k)-cc(m,
     +                       ic,4,k)))- (ti11* (cc(m,i-1,3,k)-cc(m,ic-1,
     +                       2,k))+ti12* (cc(m,i-1,5,k)-cc(m,ic-1,4,
     +                       k)))) + wa4(i-1)* ((cc(m,i-1,1,
     +                       k)+tr11* (cc(m,i-1,3,k)+cc(m,ic-1,2,
     +                       k))+tr12* (cc(m,i-1,5,k)+cc(m,ic-1,4,k)))+
     +                       (ti11* (cc(m,i,3,k)+cc(m,ic,2,
     +                       k))+ti12* (cc(m,i,5,k)+cc(m,ic,4,k))))
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      RETURN
      END
      SUBROUTINE vradbg(mp,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
     +                  mdimc,wa)
      USE m_constants, ONLY : pimach
      DIMENSION ch(mdimc,ido,l1,ip),cc(mdimc,ido,ip,l1),
     +          c1(mdimc,ido,l1,ip),c2(mdimc,idl1,ip),
     +          ch2(mdimc,idl1,ip),wa(ido)

      tpi = 2.*pimach()
      arg = tpi/real(ip)
      dcp = cos(arg)
      dsp = sin(arg)
      idp2 = ido + 2
      nbd = (ido-1)/2
      ipp2 = ip + 2
      ipph = (ip+1)/2
      IF (ido.LT.l1) GO TO 40
      DO 30 k = 1,l1
         DO 20 i = 1,ido
            DO 10 m = 1,mp
               ch(m,i,k,1) = cc(m,i,1,k)
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
      GO TO 80
   40 DO 70 i = 1,ido
         DO 60 k = 1,l1
            DO 50 m = 1,mp
               ch(m,i,k,1) = cc(m,i,1,k)
   50       CONTINUE
   60    CONTINUE
   70 CONTINUE
   80 DO 110 j = 2,ipph
         jc = ipp2 - j
         j2 = j + j
         DO 100 k = 1,l1
            DO 90 m = 1,mp
               ch(m,1,k,j) = cc(m,ido,j2-2,k) + cc(m,ido,j2-2,k)
               ch(m,1,k,jc) = cc(m,1,j2-1,k) + cc(m,1,j2-1,k)
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE
      IF (ido.EQ.1) GO TO 210
      IF (nbd.LT.l1) GO TO 160
      DO 150 j = 2,ipph
         jc = ipp2 - j
         DO 140 k = 1,l1
            DO 130 i = 3,ido,2
               ic = idp2 - i
               DO 120 m = 1,mp
                  ch(m,i-1,k,j) = cc(m,i-1,2*j-1,k) + cc(m,ic-1,2*j-2,k)
                  ch(m,i-1,k,jc) = cc(m,i-1,2*j-1,k) -
     +                             cc(m,ic-1,2*j-2,k)
                  ch(m,i,k,j) = cc(m,i,2*j-1,k) - cc(m,ic,2*j-2,k)
                  ch(m,i,k,jc) = cc(m,i,2*j-1,k) + cc(m,ic,2*j-2,k)
  120          CONTINUE
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE
      GO TO 210
  160 DO 200 j = 2,ipph
         jc = ipp2 - j
         DO 190 i = 3,ido,2
            ic = idp2 - i
            DO 180 k = 1,l1
               DO 170 m = 1,mp
                  ch(m,i-1,k,j) = cc(m,i-1,2*j-1,k) + cc(m,ic-1,2*j-2,k)
                  ch(m,i-1,k,jc) = cc(m,i-1,2*j-1,k) -
     +                             cc(m,ic-1,2*j-2,k)
                  ch(m,i,k,j) = cc(m,i,2*j-1,k) - cc(m,ic,2*j-2,k)
                  ch(m,i,k,jc) = cc(m,i,2*j-1,k) + cc(m,ic,2*j-2,k)
  170          CONTINUE
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE
  210 ar1 = 1.
      ai1 = 0.
      DO 270 l = 2,ipph
         lc = ipp2 - l
         ar1h = dcp*ar1 - dsp*ai1
         ai1 = dcp*ai1 + dsp*ar1
         ar1 = ar1h
         DO 230 ik = 1,idl1
            DO 220 m = 1,mp
               c2(m,ik,l) = ch2(m,ik,1) + ar1*ch2(m,ik,2)
               c2(m,ik,lc) = ai1*ch2(m,ik,ip)
  220       CONTINUE
  230    CONTINUE
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         DO 260 j = 3,ipph
            jc = ipp2 - j
            ar2h = dc2*ar2 - ds2*ai2
            ai2 = dc2*ai2 + ds2*ar2
            ar2 = ar2h
            DO 250 ik = 1,idl1
               DO 240 m = 1,mp
                  c2(m,ik,l) = c2(m,ik,l) + ar2*ch2(m,ik,j)
                  c2(m,ik,lc) = c2(m,ik,lc) + ai2*ch2(m,ik,jc)
  240          CONTINUE
  250       CONTINUE
  260    CONTINUE
  270 CONTINUE
      DO 300 j = 2,ipph
         DO 290 ik = 1,idl1
            DO 280 m = 1,mp
               ch2(m,ik,1) = ch2(m,ik,1) + ch2(m,ik,j)
  280       CONTINUE
  290    CONTINUE
  300 CONTINUE
      DO 330 j = 2,ipph
         jc = ipp2 - j
         DO 320 k = 1,l1
            DO 310 m = 1,mp
               ch(m,1,k,j) = c1(m,1,k,j) - c1(m,1,k,jc)
               ch(m,1,k,jc) = c1(m,1,k,j) + c1(m,1,k,jc)
  310       CONTINUE
  320    CONTINUE
  330 CONTINUE
      IF (ido.EQ.1) GO TO 430
      IF (nbd.LT.l1) GO TO 380
      DO 370 j = 2,ipph
         jc = ipp2 - j
         DO 360 k = 1,l1
            DO 350 i = 3,ido,2
               DO 340 m = 1,mp
                  ch(m,i-1,k,j) = c1(m,i-1,k,j) - c1(m,i,k,jc)
                  ch(m,i-1,k,jc) = c1(m,i-1,k,j) + c1(m,i,k,jc)
                  ch(m,i,k,j) = c1(m,i,k,j) + c1(m,i-1,k,jc)
                  ch(m,i,k,jc) = c1(m,i,k,j) - c1(m,i-1,k,jc)
  340          CONTINUE
  350       CONTINUE
  360    CONTINUE
  370 CONTINUE
      GO TO 430
  380 DO 420 j = 2,ipph
         jc = ipp2 - j
         DO 410 i = 3,ido,2
            DO 400 k = 1,l1
               DO 390 m = 1,mp
                  ch(m,i-1,k,j) = c1(m,i-1,k,j) - c1(m,i,k,jc)
                  ch(m,i-1,k,jc) = c1(m,i-1,k,j) + c1(m,i,k,jc)
                  ch(m,i,k,j) = c1(m,i,k,j) + c1(m,i-1,k,jc)
                  ch(m,i,k,jc) = c1(m,i,k,j) - c1(m,i-1,k,jc)
  390          CONTINUE
  400       CONTINUE
  410    CONTINUE
  420 CONTINUE
  430 CONTINUE
      IF (ido.EQ.1) RETURN
      DO 450 ik = 1,idl1
         DO 440 m = 1,mp
            c2(m,ik,1) = ch2(m,ik,1)
  440    CONTINUE
  450 CONTINUE
      DO 480 j = 2,ip
         DO 470 k = 1,l1
            DO 460 m = 1,mp
               c1(m,1,k,j) = ch(m,1,k,j)
  460       CONTINUE
  470    CONTINUE
  480 CONTINUE
      IF (nbd.GT.l1) GO TO 530
      is = -ido
      DO 520 j = 2,ip
         is = is + ido
         idij = is
         DO 510 i = 3,ido,2
            idij = idij + 2
            DO 500 k = 1,l1
               DO 490 m = 1,mp
                  c1(m,i-1,k,j) = wa(idij-1)*ch(m,i-1,k,j) -
     +                            wa(idij)*ch(m,i,k,j)
                  c1(m,i,k,j) = wa(idij-1)*ch(m,i,k,j) +
     +                          wa(idij)*ch(m,i-1,k,j)
  490          CONTINUE
  500       CONTINUE
  510    CONTINUE
  520 CONTINUE
      GO TO 580
  530 is = -ido
      DO 570 j = 2,ip
         is = is + ido
         DO 560 k = 1,l1
            idij = is
            DO 550 i = 3,ido,2
               idij = idij + 2
               DO 540 m = 1,mp
                  c1(m,i-1,k,j) = wa(idij-1)*ch(m,i-1,k,j) -
     +                            wa(idij)*ch(m,i,k,j)
                  c1(m,i,k,j) = wa(idij-1)*ch(m,i,k,j) +
     +                          wa(idij)*ch(m,i-1,k,j)
  540          CONTINUE
  550       CONTINUE
  560    CONTINUE
  570 CONTINUE
  580 RETURN
      END
      SUBROUTINE vrftb1(m,n,c,ch,mdimc,wa,fac)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION ch(mdimc,n),c(mdimc,n),wa(n),fac(15)

      nf = fac(2)
      na = 0
      l1 = 1
      iw = 1
      DO 160 k1 = 1,nf
         ip = fac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idl1 = ido*l1
         IF (ip.NE.4) GO TO 30
         ix2 = iw + ido
         ix3 = ix2 + ido
         IF (na.NE.0) GO TO 10
         CALL vradb4(m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2),wa(ix3))
         GO TO 20
   10    CALL vradb4(m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2),wa(ix3))
   20    na = 1 - na
         GO TO 150
   30    IF (ip.NE.2) GO TO 60
         IF (na.NE.0) GO TO 40
         CALL vradb2(m,ido,l1,c,ch,mdimc,wa(iw))
         GO TO 50
   40    CALL vradb2(m,ido,l1,ch,c,mdimc,wa(iw))
   50    na = 1 - na
         GO TO 150
   60    IF (ip.NE.3) GO TO 90
         ix2 = iw + ido
         IF (na.NE.0) GO TO 70
         CALL vradb3(m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2))
         GO TO 80
   70    CALL vradb3(m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2))
   80    na = 1 - na
         GO TO 150
   90    IF (ip.NE.5) GO TO 120
         ix2 = iw + ido
         ix3 = ix2 + ido
         ix4 = ix3 + ido
         IF (na.NE.0) GO TO 100
         CALL vradb5(m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         GO TO 110
  100    CALL vradb5(m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  110    na = 1 - na
         GO TO 150
  120    IF (na.NE.0) GO TO 130
         CALL vradbg(m,ido,ip,l1,idl1,c,c,c,ch,ch,mdimc,wa(iw))
         GO TO 140
  130    CALL vradbg(m,ido,ip,l1,idl1,ch,ch,ch,c,c,mdimc,wa(iw))
  140    IF (ido.EQ.1) na = 1 - na
  150    l1 = l2
         iw = iw + (ip-1)*ido
  160 CONTINUE
      scale = sqrt(1./n)
      IF (na.EQ.0) GO TO 190
      DO 180 j = 1,n
         DO 170 i = 1,m
            c(i,j) = scale*ch(i,j)
  170    CONTINUE
  180 CONTINUE
      RETURN
  190 DO 210 j = 1,n
         DO 200 i = 1,m
            c(i,j) = scale*c(i,j)
  200    CONTINUE
  210 CONTINUE
      RETURN
      END

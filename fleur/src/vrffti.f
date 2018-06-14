      SUBROUTINE vrffti(n,wsave)
C***BEGIN PROLOGUE  VRFFTI
C***DATE WRITTEN   860701   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM,
C             MULTIPLE SEQUENCES
C***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST)
C***PURPOSE  Initialization for VRFFTF and VRFFTB.
C***DESCRIPTION
C
C  Subroutine VRFFTI initializes the array WSAVE which is used in
C  both VRFFTF and VRFFTB.  The prime factorization of N together with
C  a tabulation of certain trigonometric functions are computed and
C  stored in the array WSAVE.
C
C  Input Parameter
C
C  N       the length of the sequence to be transformed.  There is no
C          restriction on N.
C
C  Output Parameter
C
C  WSAVE   a work array which must be dimensioned at least N+15.
C          The same work array can be used for both VRFFTF and VRFFTB
C          as long as N remains unchanged.  Different WSAVE arrays
C          are required for different values of N.  The contents of
C          WSAVE must not be changed between calls of VRFFTF or VRFFTB.
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
C***ROUTINES CALLED  VRFTI1
C***END PROLOGUE  VRFFTI
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION wsave(n+15)
C***FIRST EXECUTABLE STATEMENT  VRFFTI
      IF (n.EQ.1) RETURN
      CALL vrfti1(n,wsave(1),wsave(n+1))
      RETURN
      END
      SUBROUTINE vrfti1(n,wa,fac)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      USE m_constants, ONLY : pimach
      DIMENSION wa(n),fac(15),ntryh(4)
      DATA ntryh(1),ntryh(2),ntryh(3),ntryh(4)/4,2,3,5/

      nl = n
      nf = 0
      j = 0
   10 j = j + 1
      IF ( j <=4 ) THEN
        ntry = ntryh(j)
      ELSE
        ntry = ntry + 2
      ENDIF
   40 nq = nl/ntry
      nr = nl - ntry*nq
      IF ( nr /= 0 ) THEN
        GOTO 10
      ENDIF
      nf = nf + 1
      fac(nf+2) = ntry
      nl = nq
      IF (ntry.NE.2) GO TO 70
      IF (nf.EQ.1) GO TO 70
      DO 60 i = 2,nf
         ib = nf - i + 2
         fac(ib+2) = fac(ib+1)
   60 CONTINUE
      fac(3) = 2
   70 IF (nl.NE.1) GO TO 40
      fac(1) = n
      fac(2) = nf
      tpi = 2.*pimach()
      argh = tpi/real(n)
      is = 0
      nfm1 = nf - 1
      l1 = 1
      IF (nfm1.EQ.0) RETURN
      DO 100 k1 = 1,nfm1
         ip = fac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         ipm = ip - 1
         DO 90 j = 1,ipm
            ld = ld + l1
            i = is
            argld = real(ld)*argh
            fi = 0.
            DO 80 ii = 3,ido,2
               i = i + 2
               fi = fi + 1.
               arg = fi*argld
               wa(i-1) = cos(arg)
               wa(i) = sin(arg)
   80       CONTINUE
            is = is + ido
   90    CONTINUE
         l1 = l2
  100 CONTINUE
      RETURN
      END

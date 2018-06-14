      SUBROUTINE cfft(a,b,ntot,n,nspan,isn)

*-------------------------------------------------------------*
* driver routine for ccfft subroutine instead of cfft on cray *
*              and dcft, dcft2 and dcft3 essl routines on IBM *
*-------------------------------------------------------------*

      IMPLICIT NONE

C ... input variables
      INTEGER ntot,n,nspan,isn
      REAL a(ntot),b(ntot)
      
C ... local variables
      INTEGER i,ld1,ld2,n1,n2,n3,dimfft,idim,s(4)

      LOGICAL calc

      REAL, DIMENSION(:), ALLOCATABLE :: table, aux
      REAL, DIMENSION(:), ALLOCATABLE :: work, aux1, aux2
      COMPLEX, DIMENSION(:), ALLOCATABLE :: x

      INTEGER naux,naux1,naux2,lam,la1,la2
      REAL, PARAMETER :: scale = 1.0
C ... save variables
      SAVE n1,n2,n3

C ... data statements
      DATA s/1,1,1,1/

C ... now, what do we have to do ?

      IF ((ntot.eq.n).and.(n.eq.nspan)) THEN
C  ...                                          1D-FFT
          dimfft=1
          n1=n
          n2=2
          n3=2
          calc=.true. 
      ELSE 
          IF (n.eq.nspan) THEN
C  ...                                          2D or 3D first step
              n1=n
              n2=0
              calc=.false.
          ELSE 
              IF (ntot.eq.nspan) THEN
C  ...                                          2D second step or 3D third step
                  IF (n2.eq.0) THEN
                      dimfft=2
                      n2=n
                      n3=1
                      calc=.true.
                  ELSE
                      dimfft=3
                      n3=n
                      calc=.true.
                  ENDIF
              ELSE
C  ...                                          3D second step.
                  n2=n
                  calc=.false.
              ENDIF
          ENDIF
      ENDIF

      IF (calc) THEN

C ... build x from a and b
         
         ALLOCATE (x(ntot))
         x = (0.0,0.0)
         DO i=1,ntot
           x(i)=cmplx(a(i),b(i))
         ENDDO
C ... do the FFT

#ifdef CPP_AIX

         ld1 = n1
         ld2 = n1*n2
         IF (dimfft.eq.1) THEN
           naux1 = 20000
           IF (n1 > 2048) naux1 = naux1 + CEILING( 2.28*n1)
           ALLOCATE ( aux1(naux1), aux2(naux1) )
         ELSEIF (dimfft.eq.2) THEN
           naux1 = 40000 + CEILING( 2.28*(n1+n2) )
           IF ( max(n1, n2) <= 2048 ) naux1 = 40000
           naux2 = 20000 + CEILING( (2*max(n1, n2)+256) *
     +                         (2.28 + min(64, n1, n2)) )
           IF ( max(n1, n2) < 256 ) naux2 = 20000 
           ALLOCATE ( aux1(naux1), aux2(naux2) )
         ELSE IF (dimfft.eq.3) THEN
           IF  ( max(n2, n3) < 252 ) THEN
             IF (n1 <= 2048)  THEN
                 naux = 60000
             ELSE
                 naux = 60000+CEILING(4.56*n1)
             ENDIF
           ELSE
             la1 = CEILING( (2*n2+256)*(min(64,n1)   +4.56) )
             la2 = CEILING( (2*n3+256)*(min(64,n1*n2)+4.56) )
             lam = max( la2 , la1 ) 
             IF ( (n2 >= 252) .AND. (n3 <  252) ) lam = la1
             IF ( (n2 <  252) .AND. (n3 >= 252) ) lam = la2
             IF (n1 <= 2048) THEN
                 naux = 60000 + lam
             ELSE
                 naux = 60000+CEILING(4.56*n1)+lam
             ENDIF
           ENDIF
           ALLOCATE ( aux(naux) )
         ENDIF
#else

         ld1=n1
         ld2=n2
         s(1)=dimfft

#ifndef CPP_MPI
c t,j90:
         idim = 1024*n
         ALLOCATE (table(16*n+100),work(idim))
#else
c t3e:
         idim = 2*n1*max(n2,1)*max(n3,1)
         ALLOCATE (table(12*(n1+n2+n3)),work(idim))
#endif
#endif
         IF (dimfft.eq.1) THEN
#ifdef CPP_AIX
             CALL dcft(-1, x, 1, 1, x, 1, 1, n, 1,-isn, 1.0,
     +                              aux1, naux1, aux2, naux1)
             CALL dcft( 0, x, 1, 1, x, 1, 1, n, 1,-isn, 1.0,
     +                              aux1, naux1, aux2, naux1)
#else
             CALL ccfft(0,n,1.0,x,x,table,work,s)
             CALL ccfft(isn,n,1.0,x,x,table,work,s)
#endif
         ENDIF
         IF (dimfft.eq.2) THEN
#ifdef CPP_AIX
             CALL dcft2(-1, x, 1, n1, x, 1, n1, n1, n2,-isn, 1.0,
     +                                  aux1, naux1, aux2, naux2)
             CALL dcft2( 0, x, 1, n1, x, 1, n1, n1, n2,-isn, 1.0,
     +                                  aux1, naux1, aux2, naux2)
#else
             CALL ccfft2d(0,n1,n2,1.0,x,ld1,x,ld1,table,work,s)
             CALL ccfft2d(isn,n1,n2,1.0,x,ld1,x,ld1,table,work,s)
#endif
         ENDIF
         IF (dimfft.eq.3) THEN
#ifdef CPP_AIX
             CALL dcft3(x, ld1, ld2, x, ld1, ld2, n1, n2, n3,
     +                                -isn, scale, aux, naux)
#else
             CALL ccfft3d(0,n1,n2,n3,1.0,x,ld1,ld2,x,ld1,ld2,
     +                                            table,work,s)
             CALL ccfft3d(isn,n1,n2,n3,1.0,x,ld1,ld2,x,ld1,ld2,
     +                                              table,work,s)
#endif
         ENDIF

#ifdef CPP_AIX
         IF (dimfft.eq.3) THEN
           DEALLOCATE ( aux )
         ELSE
            DEALLOCATE ( aux1, aux2 )
         ENDIF
#else
         DEALLOCATE (table,work)
#endif

C ... backup a and b

         DO i=1,ntot
           a(i)=real(x(i))
           b(i)=aimag(x(i))
         ENDDO

         DEALLOCATE (x)
      ENDIF

      END

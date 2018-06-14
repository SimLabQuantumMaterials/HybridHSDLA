      COMPLEX FUNCTION anglso(theta,phi,l1,m1,is1,l2,m2,is2)
C
c calculates spin-orbit matrix for theta,phi =/= 0
c
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT(IN) :: is1,is2,l1,l2,m1,m2
      REAL,    INTENT(IN) :: theta,phi             
C     ..
C     .. Local Scalars ..
      REAL sgm1,sgm2,xlz,xlpl,xlmn,angl_r,angl_i
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,real,sqrt,isign
c
      anglso = cmplx(0.0,0.0)
      IF (l1.NE.l2) THEN
          RETURN
      ENDIF
C
      sgm1 = is1
      sgm2 = is2
      IF (l1.LT.0) THEN
          WRITE (6,FMT=*)
     +      ' PROGRAM STOPS IN ANGLSO ( L < 0 ) .'
          WRITE (6,FMT=*) ' L1 =',l1,'    L2 =',l2
          STOP 'ANGLSO'
      ELSE IF ((abs(m1).GT.l1) .OR. (abs(m2).GT.l2)) THEN
          WRITE (6,FMT=*)
     +      ' PROGRAM STOPS IN ANGLSO ( M < L OR L < M )'
          WRITE (6,FMT=*) ' L1 =',l1,'    L2 =',l2
          WRITE (6,FMT=*) ' M1 =',m1,'    M2 =',m2
          STOP 'ANGLSO'
      ELSE IF ((is1.NE.-1.AND.is1.NE.1) .OR.
     +         (is2.NE.-1.AND.is2.NE.1)) THEN
          WRITE (6,FMT=*)
     +      ' PROGRAM STOPS IN ANGLSO ( S >< +-1/2 ) .'
          WRITE (6,FMT=*) ' S1 =',0.5*sgm1,'    S2 =',0.5*sgm2
          STOP 'ANGLSO'
      END IF
C
c lz,l+,l-
C
      xlz = 0.0
      xlpl= 0.0
      xlmn= 0.0
c is1.eq.is2-2 -> <-| |+> => l+   
      IF (m1.EQ.m2+1) THEN
          xlpl = sqrt(real((l2-m2)* (l2+m2+1)))
c is1.eq.is2+2 -> <+| |-> => l-   
      ELSE IF (m1.EQ.m2-1) THEN
          xlmn = sqrt(real((l2+m2)* (l2-m2+1)))
c is1.eq.is2 -> <+| |+> => lz   
      ELSE IF (m1.EQ.m2  ) THEN
          xlz  = m2
      END IF
c
c rotated spin-orbit angular matrix
c <1| |1> or <2| |2>          
c
      IF (is1.EQ.is2) THEN
          angl_r = isign(1,is1) * ( cos(theta)*xlz +
     *                          0.5*sin(theta)*cos(phi)*(xlmn + xlpl) )
          angl_i = isign(1,is1)*0.5*sin(theta)*sin(phi)*(xlmn - xlpl)
c <1| |2>
      ELSEIF (is1.EQ.is2+2) THEN
          angl_r =  sin(theta)*xlz +  cos(phi)*(
     +            - cos(theta/2.)**2 * xlmn + sin(theta/2.)**2 * xlpl )
          angl_i = -sin(phi)*( 
     +              cos(theta/2.)**2 * xlmn + sin(theta/2.)**2 * xlpl )
c <2| |1>
      ELSEIF (is1.EQ.is2-2) THEN
          angl_r =  sin(theta)*xlz +  cos(phi)*(
     +            - cos(theta/2.)**2 * xlpl + sin(theta/2.)**2 * xlmn )
          angl_i =  sin(phi)*( 
     +              cos(theta/2.)**2 * xlpl + sin(theta/2.)**2 * xlmn )
      ELSE
          angl_r = 0.0
          angl_i = 0.0
      ENDIF
C
      anglso = cmplx(angl_r,angl_i)

      RETURN
      END

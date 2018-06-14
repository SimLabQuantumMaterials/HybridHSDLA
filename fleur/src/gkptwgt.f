      SUBROUTINE gkptwgt(
     >                   nkptd,amat,bmat,nkpt,latnam,
     X                   bk,
     <                   w,wt)
c                                                             
c     this subroutine generates the weight factor of a k-point 
c     in the irreducible wedge of a 2d-brillouin zone         
c     ver are vertex coordinates in counter clockwise order and 
c     in units of pi/a(1) and pi/a(2)                              
c     it is assumed that each k-point has the same 'area'. D.S.Wang     
c     
c     changed by                        Stefan Bl"ugel, IFF, Jan.96
c                                                           
      USE m_constants, ONLY : pimach
      IMPLICIT NONE
c                                                        
c     .. Scalar Arguments
      INTEGER,    INTENT (IN) :: nkpt,nkptd
      REAL,       INTENT (OUT):: wt
      CHARACTER*3,INTENT (IN) :: latnam
c     ..
c     .. Array Arguments ..
      REAL, INTENT (INOUT):: bk(3,nkptd)
      REAL, INTENT (IN)   :: amat(3,3),bmat(3,3)
      REAL, INTENT (OUT)  :: w(nkptd)
c     ..
c     .. Local Scalars ..
      REAL cross,dot,eps,tpii,x1,x2,y1,y2
      REAL tpi,pi,s1,s2
      INTEGER i,j,ikpt,nver
c     ..
c     .. Local Arrays ..
      REAL ver(2,4),dummy(2,nkptd)
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC abs,atan2
c     ..
c     .. Data Statements ..
      DATA ver/0.e0,0.e0,1.e0,0.e0,1.e0,1.e0,0.e0,1.e0/
      DATA eps/1.e-6/
c     ..
      pi = pimach()
      tpi = 2 * pi
      tpii =  1.0 / tpi
c                                                           
C      nver = 3
      IF ( latnam.eq.'squ' ) THEN
         nver = 3
      ELSEIF ( latnam.eq.'p-r' .or. latnam.eq.'c-r' ) THEN
         nver = 4
      ELSEIF ( latnam.eq.'hex' ) THEN
         nver = 3
         ver(2,3) = 1./3.
      ELSEIF ( latnam.eq.'hx3' .or. latnam.eq. 'obl' ) THEN
         STOP "gkptwgt:  weights for hx3 or obl not defined"
      ENDIF
c                                                          
c     transform from internal coordinates to xy-coordinates
c
c                            changed by shz Feb.96
      DO 10 ikpt = 1 , nkpt
         w(ikpt) = 0
         dummy(1,ikpt)=bk(1,ikpt)
         dummy(2,ikpt)=bk(2,ikpt)
         s1 = 0.0
         s2 = 0.0
         DO i = 1,2 
            s1 = s1+bmat(i,1)*bk(i,ikpt)
            s2 = s2+bmat(i,2)*bk(i,ikpt)
         ENDDO   
         IF (latnam.eq.'hex') THEN
           bk(1,ikpt) = s1*amat(2,2)/tpi
           bk(2,ikpt) = s2*amat(1,1)/pi
         ELSE
           bk(1,ikpt) = s1*amat(1,1)/pi
           bk(2,ikpt) = s2*amat(2,2)/pi
         ENDIF  
   10 ENDDO 

      DO 20 ikpt = 1 , nkpt
         DO 30 i = 1,nver
            x1 = ( ver(1,i)-bk(1,ikpt) ) / amat(1,1)
            y1 = ( ver(2,i)-bk(2,ikpt) ) / amat(2,2)
            j  = i + 1
            if ( j.gt.nver ) j = 1
            x2 = ( ver(1,j)-bk(1,ikpt) ) / amat(1,1)
            y2 = ( ver(2,j)-bk(2,ikpt) ) / amat(2,2)
            dot = x1*x2 + y1*y2
            cross = x1*y2 - y1*x2
            IF ( abs(cross).ge.eps ) THEN
              w(ikpt) = w(ikpt) + atan2(cross,dot)
            ENDIF 
   30    ENDDO
   20 ENDDO
c
      DO ikpt = 1 , nkpt
         w(ikpt) = w(ikpt) * tpii
      ENDDO
c   
      wt = 0.0
      DO ikpt = 1,nkpt
         wt = wt + w(ikpt)
         bk(1,ikpt)=dummy(1,ikpt)
         bk(2,ikpt)=dummy(2,ikpt)
      ENDDO

      RETURN
      END

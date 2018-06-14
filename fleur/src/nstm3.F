      MODULE m_nstm3
c***********************************************************************
c     included writing to vacwave!
c     set up mapping array to general G_parallel(j)=(gvac1(j),gvac2(j))
c             for vacuum density in order to write out information
c             on electronic structure for calculation of tunneling current    
c                            change by shz, Jan.99
c
c***********************************************************************
      CONTAINS
      SUBROUTINE nstm3(
     >                 nv2d,nop,natd,nvd,nmzd,k1d,k2d,k3d,n3d,ikpt,nv,
     >                 jspins,nop2,mx1,mx2,n2d,jspin,nmz,nvac,nkpt,nw,
     >                 tworkf,area,delz,z1,wk,k1,k2,jspd,nwdd,nat,ig2,
     >                 ig,nstr2,mrot,rgphs,amat,bmat,evac,vz,vz0,taual,
     X                 gvac1d,gvac2d)
c
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nv2d,nop,natd,nvd,nmzd,k1d,k2d,k3d
      INTEGER, INTENT (IN) :: ikpt,nv,jspins,nop2,mx1,mx2,n2d,n3d
      INTEGER, INTENT (IN) :: jspin,nmz,nvac,nkpt,nw,jspd,nwdd,nat
      REAL,    INTENT (IN) :: tworkf,area,delz,wk,z1
C     ..
C     .. Array  Arguments ..
      INTEGER, INTENT (IN) :: k1(nvd),k2(nvd)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),ig2(n3d)
      INTEGER, INTENT (IN) :: nstr2(n2d),mrot(3,3,nop)
      REAL,    INTENT (IN) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      REAL,    INTENT (IN) :: amat(3,3),bmat(3,3),evac(2)
      REAL,    INTENT (IN) :: vz(nmzd,2),vz0(2),taual(3,natd)
      INTEGER, INTENT (OUT) :: gvac1d(nv2d),gvac2d(nv2d)
C     ..
C     .. Local Scalars
      INTEGER n2,k,j,i,ivac
      REAL    dz0
C     ..
C     .. Local Arrays ..
      INTEGER gvac1(nv2d),gvac2(nv2d),gindex(nv2d)
      REAL gvacl(nv2d),gvac(2)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC sqrt
C     ..
C     .. External Subroutines ..
      EXTERNAL sort
c
      IF (ikpt.EQ.1) THEN
         n2 = 0
         DO 47 k = 1,nv
            DO 45 j = 1,n2
               IF (k1(k).EQ.gvac1(j).
     +              AND.k2(k).EQ.gvac2(j)) THEN
                  GO TO 47
               END IF
 45         CONTINUE
            n2 = n2 + 1
            IF (n2.GT.nv2d) STOP 'vacden0'
            gvac1(n2) = k1(k)
            gvac2(n2) = k2(k)
            DO i=1,2
               gvac(i)=k1(k)*bmat(1,i)+k2(k)*bmat(2,i)
            END DO
            gvacl(n2) = sqrt(real(gvac(1)**2+gvac(2)**2))
 47      CONTINUE
         CALL sort(n2,gvacl,gindex)
         DO j = 1,n2
c  gvac1d, gvac2d are now ordered by increasing length
            gvac1d(j)=gvac1(gindex(j))
            gvac2d(j)=gvac2(gindex(j))
         END DO
c 
         IF (jspin.EQ.1) THEN
            WRITE (87,'(f10.6,1x,i1,1x,f10.6)') tworkf,jspins,area
            WRITE (87,'(2(f10.6,1x))') amat(1,1), amat(2,1)
            WRITE (87,'(2(f10.6,1x))') amat(1,2), amat(2,2)
            WRITE (87,'(2(f10.6,1x))') bmat(1,1), bmat(2,1)
            WRITE (87,'(2(f10.6,1x))') bmat(1,2), bmat(2,2)
            WRITE (87,'(i2)') nop2
            DO j = 1, nop2
               WRITE (87,'(i2,1x,i2)') mrot(1,1,j), mrot(1,2,j)
               WRITE (87,'(i2,1x,i2)') mrot(2,1,j), mrot(2,2,j)
            END DO
            WRITE (87,'(i3)') n2
            DO j = 1,n2
               WRITE (87,'(3(i3,1x),f10.6)') j, gvac1(gindex(j)), 
     +              gvac2(gindex(j)),gvacl(gindex(j))
            END DO
c
c     Write info on 2D-starfunctions

            WRITE (87,'(i2,1x,i2,1x,i2)') mx1,mx2, n2d
            DO i=1, n2d
               WRITE (87,'(i2)') nstr2(i)
            END DO
            DO i=-mx1, mx1
               DO j=-mx2,mx2
                  WRITE (87,'(i2,1x,e12.4)') ig2(ig(i,j,0)),rgphs(i,j,0)
               END DO
            END DO
         END IF
         WRITE (87,'(i1,1x,i1)') jspin, nvac
         WRITE (87,'(2(e16.8,1x))') (evac(i), i=1,nvac)
         WRITE (87,'(2(e16.8,1x))') (vz0(i), i=1,nvac)
         dz0=0.0
         DO i=1, nat
            IF (abs(taual(3,i)).GT.dz0) dz0=abs(taual(3,i))
         END DO
         dz0=z1-dz0*amat(3,3)
         WRITE (87,'(i3,1x,f6.4,1x,f12.6)') nmz,delz,dz0   
         DO ivac=1,nvac
            DO i=1, nmz
               WRITE (87,'(e16.8)') vz(i,ivac)
            END DO
         END DO
         WRITE (87,'(i4)') nkpt
      END IF

c  only write here if not on T3E

#ifndef CPP_MPI

      WRITE (87,'(i3,1x,f12.6)') ikpt,wk
#endif

      END SUBROUTINE nstm3
      END MODULE m_nstm3

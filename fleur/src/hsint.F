      MODULE m_hsint
      CONTAINS
      SUBROUTINE hsint(
     >                 npotmatfile,l_noco,l_J,matsize,
     >                 jspd,nvd,n3d,k1d,k2d,k3d,
     >                 ng3,ig,rgphs,ustep,vpw,k1,k2,k3,nv,jspin,
     >                 n_size,n_rank,rk,bkpt,bbmat,nlotot,
     <                 nmat,aa,bb)
c*********************************************************************
c     initializes and sets up the hamiltonian and overlap matrices
c     for the interstitial. only the lower triangle of the hermitian
c     matrices are stored in compact real mode such that if h(i,j),
c     i.ge.j, is hermitian and a is real, then
c       a(i,j)=real( h(i,j) )  and  a(j,i)=imag( h(i,j) )
c                    m. weinert  1986
c
c For the eigenvector parallelization each pe calculates an equal share
c of columns labeled nc. Then the starting element of a columns nc is
c
c     ii = (nc-1)*( n_rank - n_size + 1 ) + n_size*(nc-1)*nc/2
c
c and, if a non-collinear matrix has to be set up, the starting column
c for the second spin-direction is
c
c     nc = int( 1. + (nv - n_rank - 1)/n_size ) + 1 .
c
c For this direction, the outer loop starts at
c
c istart = n_rank + (nc - 1)*n_size - nv .                        gb99
c
c for a lo-calculation nv has to be replaced by nv+nlotot         gb01
c
c*********************************************************************
       USE m_dotir, ONLY : dotirp
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER nmat,nlotot
      INTEGER, INTENT (IN) :: npotmatfile,ng3,n_size,n_rank,jspin
      INTEGER, INTENT (IN) :: matsize,jspd,nvd,n3d,k1d,k2d,k3d
      LOGICAL, INTENT (IN) :: l_noco,l_J
C     ..
C     .. Array Arguments ..
      COMPLEX, INTENT (IN) :: ustep(n3d)
      COMPLEX, INTENT (INOUT) :: vpw(n3d)
      INTEGER, INTENT (IN) :: nv(jspd)
      INTEGER, INTENT (IN) :: k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      REAL,    INTENT (IN) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      REAL,    INTENT (IN) :: rk(nvd,jspd),bkpt(3),bbmat(3,3)
#ifdef CPP_INVERSION
      REAL,    INTENT (OUT):: aa(matsize),bb(matsize)
#else
      COMPLEX, INTENT (OUT):: aa(matsize),bb(matsize)
#endif
C     ..
C     .. Local Scalars ..
      COMPLEX th,ts
      REAL phase,b1(3),b2(3),r2
      INTEGER i,i1,i2,i3,ii,in,j,ig3,ispin,l
      INTEGER istart,nc

      COMPLEX ust1,vp1
      COMPLEX, ALLOCATABLE :: vpw1(:)  ! for J constants
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC real
C     ..
      DO 10 i = 1,matsize
#ifdef CPP_INVERSION
         aa(i) = 0.0
         bb(i) = 0.0
#else
         aa(i) = cmplx(0.0,0.0)
         bb(i) = cmplx(0.0,0.0)
#endif
 10   CONTINUE
      ust1 = ustep(1)
      ispin = jspin
      nmat = nv(ispin)

c---> pk non-collinear
      IF (l_noco) THEN
c---> determine spin-up spin-up part of Hamiltonian- and overlapp-matrix
c---> reload V_11
         READ (npotmatfile) (vpw(ig3),ig3=1,ng3)

c--- J const
        IF( l_J) THEN
         ALLOCATE ( vpw1(n3d) )
         READ (npotmatfile) (vpw1(ig3),ig3=1,ng3)
        ENDIF
c--- J const

         nmat = nv(1) + nv(jspd)
         ispin = 1

c--- J const
        IF (l_J) THEN
         DO i = 1,ng3
            vpw(i) = (vpw(i) + vpw1(i))/2.
         END DO
        ENDIF
c--- J const

         vp1 = real(vpw(1))
      ENDIF
c---> pk non-collinear

      vp1 = vpw(1)
c---> loop over (k+g')
      ii = 0
      DO 30 i = n_rank+1, nv(ispin), n_size
c--->    loop over (k+g)
         DO 20 j = 1,i - 1
            ii = ii + 1
c-->     determine index and phase factor
            i1 = k1(i,ispin) - k1(j,ispin)
            i2 = k2(i,ispin) - k2(j,ispin)
            i3 = k3(i,ispin) - k3(j,ispin)
            in = ig(i1,i2,i3)
            IF (in.EQ.0) GO TO 20
            phase = rgphs(i1,i2,i3)
c+APW_LO
#ifdef CPP_APW
            b1(1) = bkpt(1)+k1(i,ispin) ; b2(1) = bkpt(1)+k1(j,ispin)
            b1(2) = bkpt(2)+k2(i,ispin) ; b2(2) = bkpt(2)+k2(j,ispin)
            b1(3) = bkpt(3)+k3(i,ispin) ; b2(3) = bkpt(3)+k3(j,ispin)
            r2 = dotirp(b1,b2,bbmat)
            th = phase*(0.5*r2*ustep(in)+vpw(in))
#else
            th = phase* (0.25* (rk(i,ispin)**2+rk(j,ispin)**2)*ustep(in)
     +         + vpw(in))
#endif
c-APW_LO
c--->    determine matrix element and store
            ts = phase*ustep(in)
#ifdef CPP_INVERSION
            aa(ii) = real(th)
            bb(ii) = real(ts)
#else
            aa(ii) = th
            bb(ii) = ts
#endif
   20    CONTINUE
c--->    diagonal term (g-g'=0 always first star)
         ii = ii + 1
#ifdef CPP_INVERSION
         aa(ii) = 0.5*rk(i,ispin)*rk(i,ispin)*real(ust1) + real(vp1)
         bb(ii) = real(ust1)
#else
         aa(ii) = 0.5*rk(i,ispin)*rk(i,ispin)*ust1 + vp1
         bb(ii) = ust1
#endif
   30 CONTINUE

c---> pk non-collinear
      IF (l_noco) THEN
c+gb99
      nc = int( 1. + (nv(1)+nlotot - n_rank - 1)/n_size )
      istart = n_rank + nc*n_size - (nv(1)+nlotot)
!      ii = (nv(1)+nlotot+1)*(nv(1)+nlotot+2)/2 - 1
      ii = nc*(n_rank-n_size+1) + n_size*(nc+1)*nc/2 + nv(1)+nlotot
c-gb99
      ispin = 2
c---> determine spin-down spin-down part of Hamiltonian- and ovlp-matrix
c---> reload V_22

c--- J constants 
      IF(.not.l_J) THEN
       READ (npotmatfile) (vpw(ig3),ig3=1,ng3)
       vp1 = real(vpw(1))
      ENDIF
c--- J constants

c---> loop over (k+g')
      DO i = istart+1, nv(ispin), n_size
         nc = nc + 1
c--->    loop over (k+g)
         DO j = 1,i - 1
c-gb99      ii = (nv(1)+i-1)*(nv(1)+i)/2 + nv(1) + j
            ii = (nc-1)*( n_rank - n_size + 1 ) + n_size*(nc-1)*nc/2 +
     +            nv(1)+nlotot + j
c--->       determine index and phase factor
            i1 = k1(i,ispin) - k1(j,ispin)
            i2 = k2(i,ispin) - k2(j,ispin)
            i3 = k3(i,ispin) - k3(j,ispin)
            in = ig(i1,i2,i3)
            IF (in.EQ.0) THEN
               WRITE (*,*) 'HSINT: G-G'' not in star i,j= ',i,j
            ELSE
               phase = rgphs(i1,i2,i3)
c+APW_LO
#ifdef CPP_APW
            b1(1) = bkpt(1)+k1(i,ispin) ; b2(1) = bkpt(1)+k1(j,ispin)
            b1(2) = bkpt(2)+k2(i,ispin) ; b2(2) = bkpt(2)+k2(j,ispin)
            b1(3) = bkpt(3)+k3(i,ispin) ; b2(3) = bkpt(3)+k3(j,ispin)
            r2 = dotirp(b1,b2,bbmat)
               th = phase*( 0.5*r2*ustep(in) + vpw(in) )
#else
            th = phase* (0.25* (rk(i,ispin)**2+rk(j,ispin)**2)*ustep(in)
     +         + vpw(in))
#endif
c-APW_LO
               ts = phase*ustep(in)
               aa(ii) = th
               bb(ii) = ts
            ENDIF
         ENDDO
c--->    diagonal term (g-g'=0 always first star)
c-gb99   ii = (nv(1)+i)*(nv(1)+i+1)/2
         ii = ii + 1
         aa(ii) = 0.5*rk(i,ispin)*rk(i,ispin)*ust1 + vp1
         bb(ii) = ust1
      ENDDO

c---> determine spin-down spin-up part of Hamiltonian- and ovlp-matrix
c---> reload real part of V_21
      READ (npotmatfile) (vpw(ig3),ig3=1,ng3)
      nc = int( 1. + (nv(1)+nlotot - n_rank - 1)/n_size )
c
c---> loop over (k+g')
c
      DO i = istart+1, nv(2), n_size
         nc = nc + 1
c--->    loop over (k+g)
         DO j = 1,nv(1)
c-gb99      ii = (nv(1)+i-1)*(nv(1)+i)/2 + j
            ii = (nc-1)*( n_rank - n_size + 1 ) + n_size*(nc-1)*nc/2 + j
c--->       determine index and phase factor
            i1 = k1(i,2) - k1(j,1)
            i2 = k2(i,2) - k2(j,1)
            i3 = k3(i,2) - k3(j,1)
            in = ig(i1,i2,i3)
            IF (in.EQ.0) THEN
               WRITE (*,*) 'HSINT: G-G'' not in star i,j= ',i,j
            ELSE
              aa(ii) = rgphs(i1,i2,i3)*vpw(in) 
c--- J constants 
              IF(l_J) THEN
              aa(ii) = 0
              ENDIF
c--- J constants

            ENDIF
         ENDDO
      ENDDO
c---> pk non-collinear
      ENDIF

      IF (l_J) DEALLOCATE (vpw1)

      RETURN
      END SUBROUTINE hsint
      END MODULE m_hsint

      MODULE m_apws
c*********************************************************************
c     determines the lapw list such that |k+G|<rkmax.
c     bk(i) is the nk k-point given in internal (i.e. b1,b2,b3) units.
c        m. weinert  1986
c     unit 29 removed gb 2004
c*********************************************************************
c     modified for explicit use of z-reflection symmetry in seclr4.f
c        g. bihlmayer '96
c     subroutine boxdim added to treat non-orthogonal lattice vectors
c        s.bluegel, IFF, 18.Nov.97
c*********************************************************************
      CONTAINS
      SUBROUTINE apws(
     >                jspd,nvd,nbasfcn,jspins,l_ss,qss,bk,
     >                rkmax,nk,nk0,bbmat,bmat,l_zref,n_size,jspin,
     <                bkpt,nv,rk,k1,k2,k3,matind,nred)

      USE m_dotir, ONLY : dotirp
      IMPLICIT NONE
C     .. 
C     .. Scalar Arguments ..
      INTEGER, INTENT  (IN) :: jspd,nvd,nbasfcn
      INTEGER, INTENT  (IN) :: jspins,nk,nk0,n_size,jspin
      INTEGER, INTENT (OUT) :: nred
      LOGICAL, INTENT  (IN) :: l_zref,l_ss
      REAL,    INTENT  (IN) :: rkmax
C     ..
C     .. Array Arguments ..
      REAL,    INTENT (IN)  :: qss(3),bk(3),bbmat(3,3),bmat(3,3)
      INTEGER, INTENT (OUT) :: nv(jspd)
      INTEGER, INTENT (OUT) :: matind(nbasfcn,2)
      INTEGER, INTENT (OUT) :: k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      REAL,    INTENT (OUT) :: rk(nvd,jspd),bkpt(3)
C     ..
C     .. Local Scalars ..
      REAL arltv1,arltv2,arltv3,r2,rk2,rkm,t,r2q,gla,eps
      INTEGER i,itt,j,j1,j2,j3,k,l,m,mk1,mk2,mk3,n,nk1,
     +        ispin,jsp_start,jsp_end,gmi
C     ..
C     .. Local Arrays ..
      REAL s(3),sq(3),rkq(nvd),gsk3(nvd)
      INTEGER k1rev(nvd),k2rev(nvd),k3rev(nvd),index3(nvd)
#ifdef CPP_MPI
      INTEGER              :: n_inner,n_bound
      REAL,    ALLOCATABLE :: rk_help(:)
      INTEGER, ALLOCATABLE :: k_help(:,:),pos(:)
#endif
C     ..
C     .. External Subroutines
      EXTERNAL boxdim,sort
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC int,sqrt
C     ..
c---> in a spin-spiral calculation different basis sets are used for
c---> the two spin directions, because the cutoff radius is defined
c---> by |G + k +/- qss/2| < rkmax.

      nk1 = nk0 + nk
      DO 20 i = 1,3
         bkpt(i) = bk(i)
   20 CONTINUE
#ifdef DISPLAY_INFO
      print*, 'APWS: constructing mesh for k = ', bkpt
#endif

C---> Determine rkmax box of size mk1, mk2, mk3,
c     for which |G(mk1,mk2,mk3) + (k1,k2,k3)| < rkmax
c     arltv(i) length of reciprical lattice vector along direction (i)
C
      CALL boxdim(
     >            bmat,
     <            arltv1,arltv2,arltv3)

c     (add 1+1 due to integer rounding, strange k_vector in BZ)
      mk1 = int( rkmax/arltv1 ) + 4
      mk2 = int( rkmax/arltv2 ) + 4
      mk3 = int( rkmax/arltv3 ) + 4
#ifdef DISPLAY_INFO
      print*, 'mk1,2,3 = ', mk1, mk2, mk3
#endif
      IF (l_ss) THEN
         jsp_start = 1
         jsp_end   = jspins
      ELSE
         jsp_start = jspin
         jsp_end   = jspin
      ENDIF
      rkm = rkmax
c---> if too many basis functions, reduce rkmax
      DO ispin = jsp_start,jsp_end
         GO TO 40
 30      rkm = rkm - 0.1
         WRITE (*,FMT=8000) (bkpt(i),i=1,3),rkm
 8000    FORMAT (' $$$ k=(',3f10.6,'): rkm truncated to',f12.6)
 40      rk2 = rkm*rkm
c--->    obtain vectors
         n = 0
         DO 70 j1 = -mk1,mk1
            s(1) = bkpt(1) + j1 + (2*ispin - 3)/2.0*qss(1)
            sq(1) = bkpt(1) + j1
            DO 60 j2 = -mk2,mk2
               s(2) = bkpt(2) + j2 + (2*ispin - 3)/2.0*qss(2)
               sq(2) = bkpt(2) + j2 
               DO 50 j3 = -mk3,mk3
                  s(3) = bkpt(3) + j3 + (2*ispin - 3)/2.0*qss(3)
                  sq(3) = bkpt(3) + j3
                  r2 = dotirp(s,s,bbmat)
                  r2q = dotirp(sq,sq,bbmat)
                  IF (r2.LE.rk2) THEN
                     n = n + 1
                     IF (n.GT.nvd) GO TO 30
                     k1(n,ispin) = j1
                     k2(n,ispin) = j2
                     k3(n,ispin) = j3
                     rk(n,ispin) = sqrt(r2)
                     rkq(n) = sqrt(r2q)
                  END IF
 50            CONTINUE
 60         CONTINUE
 70      CONTINUE
         nv(ispin) = n
c
c--->    sort by shell-metzner
c
c (for spin-spirals & LO's we have to sort according to the k+G's (rkq), not
c  the k+G+q's (rk). Otherwise we might couple an LO to k+G1+q and k+G2-q !)
c                                                                       gb01
         m = n
 80      m = m/2
         IF (m.LE.0) GO TO 130
         k = n - m
         j = 1
 90      i = j
 100     l = i + m
         IF (rkq(i).GT.rkq(l)) GO TO 120
 110     j = j + 1
         IF (j.GT.k) GO TO 80
         GO TO 90
 120     t = rkq(i)
         rkq(i) = rkq(l)
         rkq(l) = t
         t = rk(i,ispin)
         rk(i,ispin) = rk(l,ispin)
         rk(l,ispin) = t
         itt = k1(i,ispin)
         k1(i,ispin) = k1(l,ispin)
         k1(l,ispin) = itt
         itt = k2(i,ispin)
         k2(i,ispin) = k2(l,ispin)
         k2(l,ispin) = itt
         itt = k3(i,ispin)
         k3(i,ispin) = k3(l,ispin)
         k3(l,ispin) = itt
         i = i - m
         IF (i.LT.1) GO TO 110
         GO TO 100
 130     CONTINUE
c+gu
c--->    determine pairs of K-vectors, where K_z = K'_-z to use 
c--->    z-reflection
         IF (l_zref) THEN
            n=0
            DO i=1,nv(ispin)
               DO j=1,i
                  IF (((k1(i,ispin).EQ.k1(j,ispin)).AND.
     +                 (k2(i,ispin).EQ.k2(j,ispin))).AND.
     +                 (k3(i,ispin).EQ.-k3(j,ispin))) THEN
                     n=n+1 
                     matind(n,1)=i
                     matind(n,2)=j
                  ENDIF
               ENDDO
            ENDDO
            nred=n

#ifdef CPP_MPI
            IF (n_size.GT.1) THEN
c
c--->     order K's in sequence K_1,...K_n | K_0,... | K_-1,....K_-n
c
               ALLOCATE (pos(nv(ispin)))
               n_inner = nv(ispin) - nred
               IF (mod(nred,n_size).EQ.0) THEN
                  n_bound = nred
               ELSE
                  n_bound = (1+int( nred/n_size ))*n_size
               ENDIF
               IF (nv(ispin) - nred + n_bound.GT.nvd) THEN
                  WRITE ( 6,*) 'increase nvd by:',
     +                 nv(ispin)-nred+n_bound-nvd
                  WRITE (16,*) 'increase nvd by:',
     +                 nv(ispin)-nred+n_bound-nvd
                  STOP 'z-ref & ev || : nvd too small!'
               ENDIF

               i = 1
               j = 1
               DO n = 1, nred 
                  IF (matind(n,1).EQ.matind(n,2)) THEN
                     pos(matind(n,1)) = n_inner + i
                     i = i + 1
                  ELSE
                     pos(matind(n,1)) = j
                     pos(matind(n,2)) = j + n_bound
                     j = j + 1
                  ENDIF
               ENDDO
c--->          resort the rk,k1,k2,k3 and matind arrays:
               ALLOCATE (rk_help(nv(ispin)),k_help(3,nv(ispin)))
               DO n = 1, nv(ispin)
                  rk_help(n)  = rk(n,ispin)
                  k_help(1,n) = k1(n,ispin)
                  k_help(2,n) = k2(n,ispin)
                  k_help(3,n) = k3(n,ispin)
               ENDDO
               DO n = nv(ispin), 1, -1
                  rk(pos(n),ispin) = rk_help(n)
                  k1(pos(n),ispin) = k_help(1,n)
                  k2(pos(n),ispin) = k_help(2,n)
                  k3(pos(n),ispin) = k_help(3,n)
               ENDDO
               DO n = nred + 1, n_bound
                  rk(n,ispin) = rk(nv(ispin),ispin)
                  k1(n,ispin) = k1(nv(ispin),ispin)
                  k2(n,ispin) = k2(nv(ispin),ispin)
                  k3(n,ispin) = k3(nv(ispin),ispin)
               ENDDO
               DEALLOCATE (rk_help,k_help)
               DEALLOCATE (pos)
               nv(ispin) = nv(ispin) - nred + n_bound
            ENDIF
#endif
         ENDIF

         IF (l_ss) THEN  ! sort additionally like in strgn1... gb
           i = 1
           gla = 0.
           gsk3(1) = 0.0
           eps=1.e-10
           DO  k = 1,nv(ispin)
             IF (rkq(k)-gla.GE.eps) i=i+1
             gla = rkq(k)
             gmi = (mk1+k1(k,ispin)) +
     +             (mk2+k2(k,ispin))*(2*mk1+1) +
     +             (mk3+k3(k,ispin))*(2*mk1+1)*(2*mk2+1)
             gsk3(k) = i * (9.+(2*mk1+1)*(2*mk2+1)*(2*mk3+1)) + gmi
           ENDDO
           CALL sort(nv(ispin),gsk3,index3)
           DO  k = 1,nv(ispin)
              k1rev(k) = k1(index3(k),ispin)
              k2rev(k) = k2(index3(k),ispin)
              k3rev(k) = k3(index3(k),ispin)
              gsk3(k) =  rk(index3(k),ispin)
           ENDDO
           DO  k = 1,nv(ispin)
              k1(k,ispin) = k1rev(k)
              k2(k,ispin) = k2rev(k)
              k3(k,ispin) = k3rev(k)
              rk(k,ispin) = gsk3(k)
           ENDDO
        ENDIF
c-gu
      ENDDO

      IF ((.NOT. l_ss) .AND. (jspins.EQ.2) ) THEN
         nv(jspins-(jspin-1)) = nv(jspin)
         DO i = 1,nv(jspin)
            rk(i,jspins-(jspin-1)) = rk(i,jspin)
            k1(i,jspins-(jspin-1)) = k1(i,jspin)
            k2(i,jspins-(jspin-1)) = k2(i,jspin)
            k3(i,jspins-(jspin-1)) = k3(i,jspin)
         ENDDO
!      ELSE
!         DO i = 1, min(nv(1),nv(jspd))
!            WRITE(*,'(3i10,3x,3i10,3f12.6)') k1(i,1),k2(i,1),k3(i,1),
!     +      k1(i,jspd),k2(i,jspd),k3(i,jspd),rk(i,1),rk(i,jspd),rkq(i)
!         ENDDO
      ENDIF

      END SUBROUTINE apws
      END MODULE m_apws

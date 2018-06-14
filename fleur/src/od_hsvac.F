      MODULE m_od_hsvac
      CONTAINS
      SUBROUTINE od_hsvac(
     >     npotmatfile,nmzxyd,n2d,k1d,k2d,k3d,n3d,nmzd,jspd,
     >     nwdd,matsize,nvd,nv2d,nmzxy,nmz,ig1,nlotot,
     >     jsp,jspins,nw,vxy,vz,evac,bbmat,bmat,
     >     bkpt,k1,k2,k3,nv,ig,ig2,sk2,phi2,
     >     omtil,area,delz,z1,MM,vM,m_cyl,n2d_1,nq2,
     >     n_size,n_rank,zrfs,l_noco,l_J,qss,
     X     aa,bb,
     <     nv2)
      
c     subroutine for calculating the hamiltonian and overlap matrices in
c     the vacuum in the case of 1-dimensional calculations
c     Y. Mokrousov June 2002             
      
      USE m_cylbes
      USE m_dcylbs
      USE m_od_vacfun
      USE m_constants, only: pimach

      IMPLICIT NONE
c     ..
c     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: matsize,nmzxyd,nvd,k1d,k2d,k3d,n3d,vM
      INTEGER, INTENT (IN) :: nmzd,jspd,nwdd,nmzxy,nmz,nv2d,MM,n2d
      INTEGER, INTENT (IN) :: jsp,nw,jspins,n_size,n_rank,n2d_1,nq2
      INTEGER, INTENT (IN) :: npotmatfile,nlotot
      LOGICAL, INTENT (IN) :: zrfs,l_noco,l_J
      REAL,    INTENT (IN) :: omtil,area,delz,z1
      INTEGER, INTENT (IN) :: m_cyl
c     ..
c     .. Array Arguments ..
      COMPLEX, INTENT (INOUT) :: vxy(nmzxyd,n2d_1-1,2)
      INTEGER, INTENT (IN) :: k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      INTEGER, INTENT (IN) :: nv(jspd)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      INTEGER, INTENT (IN) :: ig2(n3d),ig1(-k3d:k3d,-MM:MM)
      INTEGER, INTENT (OUT):: nv2(jspd)
      REAL,    INTENT (INOUT) :: vz(nmzd,2,4)
      REAL,    INTENT (IN) :: evac(2,jspd,nwdd)
      REAL,    INTENT (IN) :: sk2(n2d),phi2(n2d)
      REAL,    INTENT (IN) :: bmat(3,3),bbmat(3,3),bkpt(3)
      REAL,    INTENT (IN) :: qss(3)
      
#ifdef CPP_INVERSION
      REAL,    INTENT (INOUT) :: aa(matsize),bb(matsize)
#else 
      COMPLEX, INTENT (INOUT) :: aa(matsize),bb(matsize)
#endif
c     ..
c     .. Local Scalars ..
      COMPLEX hij,sij,apw_lo,exp1,exp2,exp3,am,bm,ic
      REAL    d2,wronk,gr,gphi,qq,tpi,x,y
      INTEGER i,i2,ii,ik,j,jk,k,jspin,ipot,npot,ii0,m,l,i3,imz
      INTEGER jspin1,jspin2,jmax,irec2,irec3,ivac,ind1,gi
      INTEGER i_start,nc,nc_0,rotax,chiral,zi,m1,z,indm,indl
c     ..
c     .. Local Arrays ..

      INTEGER, ALLOCATABLE :: nvp(:,:),ind(:,:,:)
      INTEGER, ALLOCATABLE :: kvac3(:,:),map1(:,:)
      COMPLEX, ALLOCATABLE :: tddv(:,:,:,:)
      COMPLEX, ALLOCATABLE :: tduv(:,:,:,:)
      COMPLEX, ALLOCATABLE :: tudv(:,:,:,:)
      COMPLEX, ALLOCATABLE :: tuuv(:,:,:,:)
      COMPLEX, ALLOCATABLE ::  a(:,:,:),b(:,:,:)
      COMPLEX, ALLOCATABLE :: ai(:,:,:),bi(:,:,:)
      REAL, ALLOCATABLE :: bess(:),dbss(:),bess1(:)
      REAL, ALLOCATABLE :: ddnv(:,:,:),dudz(:,:,:)
      REAL, ALLOCATABLE :: duz(:,:,:)
      REAL, ALLOCATABLE :: udz(:,:,:),uz(:,:,:)
! l_J auxiliary potential array
      COMPLEX, ALLOCATABLE :: vxy1(:,:,:)
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC cmplx,conjg,exp,real,sqrt,abs,aimag
c     ..
      tpi = 2.*pimach()
      ic  = cmplx(0.,1.)
      d2 = sqrt(omtil/area)

      IF (l_J) ALLOCATE (vxy1(nmzxyd,n2d_1-1,2))

      ALLOCATE (
     &    ai(-vM:vM,nv2d,nvd),bi(-vM:vM,nv2d,nvd),
     &    nvp(nv2d,jspd),ind(n2d,nv2d,jspd),
     &    kvac3(nv2d,jspd),map1(nvd,jspd),
     &    tddv(-vM:vM,-vM:vM,nv2d,nv2d),
     &    tduv(-vM:vM,-vM:vM,nv2d,nv2d),
     &    tudv(-vM:vM,-vM:vM,nv2d,nv2d),
     &    tuuv(-vM:vM,-vM:vM,nv2d,nv2d),
     &    a(-vM:vM,nvd,jspd),b(-vM:vM,nvd,jspd),
     &    bess(-vM:vM),dbss(-vM:vM),bess1(-vM:vM),
     &    ddnv(-vM:vM,nv2d,jspd),dudz(-vM:vM,nv2d,jspd),
     &    duz(-vM:vM,nv2d,jspd),
     &    udz(-vM:vM,nv2d,jspd),uz(-vM:vM,nv2d,jspd) )

c--->     set up mapping function from 3d-->1d lapws
c--->            creating arrays ind and nvp

      DO jspin = 1,jspins

         nv2(jspin) = 0
         DO 20 k = 1,nv(jspin)
            DO 10 j = 1,nv2(jspin)
               IF (k3(k,jspin).EQ.kvac3(j,jspin)) THEN
                  map1(k,jspin) = j
                  GO TO 20
               END IF
 10         CONTINUE
            nv2(jspin) = nv2(jspin) + 1
            IF (nv2(jspin).GT.nv2d) STOP 'nv2d'
            kvac3(nv2(jspin),jspin) = k3(k,jspin)
            map1(k,jspin) = nv2(jspin)
 20      CONTINUE

         DO ik = 1,nv2d
            nvp(ik,jspin) = 0
            DO i = 1,n2d
               ind(i,ik,jspin) = 0
            END DO
         END DO

         DO k = 1,nv(jspin)
            ik = map1(k,jspin)
            nvp(ik,jspin) = nvp(ik,jspin) + 1
            ind(nvp(ik,jspin),ik,jspin) = k
         END DO
         
      ENDDO

      npot = 1      
      ivac = 1

      IF (l_noco) THEN
c--->         load the non-warping part of the potential
         READ (npotmatfile)((vz(imz,ivac,ipot),imz=1,nmzd),ipot=1,4)
         npot = 3
c--->         for J-coeff. we average the up-up and down-down parts
c--->         and off-diagonal elements of the potential matrix to zero
         IF (l_J) THEN
            vz(:,ivac,1) = (vz(:,ivac,1) + vz(:,ivac,2))/2.
            vz(:,ivac,2) =  vz(:,ivac,1)
            vz(:,ivac,3) = 0.0
            vz(:,ivac,4) = 0.0
         END IF
      ENDIF

      DO ipot = 1,npot

         IF (l_noco) THEN
            IF (.NOT.l_J) THEN
               READ (npotmatfile)((vxy(imz,k,ivac),
     +                        imz=1,nmzxy),k=1,n2d_1-1)
            END IF 
c--->  l_J we want to average the diagonal elements of the pot. matrix
            IF (l_J .AND. ipot.EQ.1) THEN
               READ (npotmatfile)((vxy(imz,k,ivac),
     +                        imz=1,nmzxy),k=1,n2d_1-1)
               READ (npotmatfile)((vxy1(imz,k,ivac),
     +                        imz=1,nmzxy),k=1,n2d_1-1)
               vxy(:,:,ivac) = (vxy(:,:,ivac)+vxy1(:,:,ivac))/2.
            END IF
                
            IF (l_J .AND. ipot.EQ.3) THEN
               READ (npotmatfile)((vxy(imz,k,ivac),
     +                        imz=1,nmzxy),k=1,n2d_1-1)
            END IF

            IF (l_J .AND. ipot.EQ.3) vxy(:,:,ivac)=cmplx(0.,0.)

         ENDIF ! loco
         
c     get the wavefunctions and set up the tuuv, etc matrices
       
         CALL od_vacfun(
     >        m_cyl,z1,nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n2d,n3d,jspd,
     >        jsp,jspins,l_noco,qss,ipot,ig,ig1,tpi,
     >        nmzxy,nmz,delz,ig2,n2d_1,
     >        bbmat,ivac,evac(1,1,nw),bkpt,MM,vM,
     >        vxy(1,1,ivac),vz,kvac3,nv2,
     <        tuuv,tddv,tudv,tduv,uz,duz,udz,dudz,ddnv)

         IF (l_noco) THEN

          DO jspin = 1,jspins

           DO k = 1,nv(jspin)
            irec3 = ig(k1(k,jspin),k2(k,jspin),k3(k,jspin))
            IF (irec3.NE.0) THEN
               irec2 = ig2(irec3)
               gr = sk2(irec2)
               gphi = phi2(irec2)
               i2 = map1(k,jspin)
               qq = gr*z1
               CALL cylbes(vM,qq,bess)
               CALL dcylbs(vM,qq,bess,dbss)
               DO m = -vM,vM
                  wronk = uz(m,i2,jspin)*dudz(m,i2,jspin) -
     -                   udz(m,i2,jspin)*duz(m,i2,jspin)
                  a(m,k,jspin)=exp(-cmplx(0.0,m*gphi))*(ic**m)*
     *                   cmplx(dudz(m,i2,jspin)*bess(m)-
     +                   udz(m,i2,jspin)*gr*dbss(m),0.0)
     *                   /(d2*wronk)
                  b(m,k,jspin)=exp(-cmplx(0.0,m*gphi))*(ic**m)*
     *                   cmplx(-duz(m,i2,jspin)*bess(m)+
     -                   uz(m,i2,jspin)*gr*dbss(m),0.0)
     *                   /(d2*wronk)
               END DO
            END IF
           ENDDO

          ENDDO  ! jspin

         ELSE 

          DO k = 1,nv(jsp)
            irec3 = ig(k1(k,jsp),k2(k,jsp),k3(k,jsp))
            IF (irec3.NE.0) THEN
               irec2 = ig2(irec3)
               gr = sk2(irec2)
               gphi = phi2(irec2)
               i2 = map1(k,jsp)
               qq = gr*z1
               CALL cylbes(vM,qq,bess) 
               CALL dcylbs(vM,qq,bess,dbss)
               DO m = -vM,vM
                  wronk = uz(m,i2,jsp)*dudz(m,i2,jsp) -
     -                   udz(m,i2,jsp)*duz(m,i2,jsp) 
                  a(m,k,1)=exp(-cmplx(0.0,m*gphi))*(ic**m)*
     *                   cmplx(dudz(m,i2,jsp)*bess(m)-
     +                   udz(m,i2,jsp)*gr*dbss(m),0.0)
     *                   /(d2*wronk)
                  
                  b(m,k,1)=exp(-cmplx(0.0,m*gphi))*(ic**m)*
     *                   cmplx(-duz(m,i2,jsp)*bess(m)+
     -                   uz(m,i2,jsp)*gr*dbss(m),0.0)
     *                   /(d2*wronk)

               END DO 
            END IF
          ENDDO

         ENDIF ! loco
c     update hamiltonian and overlap matrices

         IF (ipot.EQ.1 .OR. ipot.EQ.2) THEN
         jspin = ipot
c+gb||
         IF (ipot.EQ.1) THEN
            nc = 0
            i_start = n_rank
         ELSE
            nc = nc + nlotot
            nc_0 = nc
            i_start = mod(mod(n_rank - (nv(1)+nlotot),n_size) +
     +                                           n_size,n_size)
         ENDIF

         DO 50 i = i_start+1,nv(jspin),n_size
            ik = map1(i,jspin)
            nc = nc + 1
            IF (ipot.EQ.1) THEN
               jspin = 1
               ii0 = nc*(nc-1)/2*n_size- (nc-1)*(n_size-n_rank-1)
            ELSEIF (ipot.EQ.2) THEN
               jspin = 2
               ii0=nc*(nc-1)/2*n_size-(nc-1)*(n_size-n_rank-1)+
     +                                             nv(1)+nlotot
            ENDIF
            jspin1 = jsp
            IF (l_noco) jspin1 = jspin
            DO j = 1,i - 1
               ii = ii0 + j
c     overlap: only  (g-g') parallel=0        
               IF (map1(j,jspin).EQ.ik) THEN
                  sij = (0.0,0.0)
                  DO m = -vM,vM
                     sij = sij + conjg(a(m,i,jspin))*a(m,j,jspin) 
     +                    +conjg(b(m,i,jspin))*b(m,j,jspin)
     *                    *ddnv(m,ik,jspin1)
                  END DO

#ifdef CPP_INVERSION
                  bb(ii) = bb(ii) + real(sij)
#else 
                  bb(ii) = bb(ii) + sij
#endif
               END IF
            ENDDO
            ii = ii0 + i
            sij = (0.0,0.0)
            DO m = -vM,vM
                  sij = sij + conjg(a(m,i,jspin))*a(m,i,jspin)+ 
     +             conjg(b(m,i,jspin))*b(m,i,jspin)*ddnv(m,ik,jspin1)
            END DO 

#ifdef CPP_INVERSION
            bb(ii) = bb(ii) + real(sij)
#else 
            bb(ii) = bb(ii) + sij
#endif
 50      CONTINUE
         ENDIF ! ipot.eq.1.or.2
c   hamiltonian update 
c   for the noncylindr. contributions we use the cutoff of m_cyl        
         IF (ipot.EQ.1) THEN
            jspin1 = 1
            jspin2 = 1
            nc = 0
            i_start = n_rank
         ELSEIF (ipot.EQ.2) THEN
            jspin1 = 2
            jspin2 = 2
            nc = nc_0
            i_start = mod(mod(n_rank - (nv(1)+nlotot),n_size) +
     +                                           n_size,n_size)
         ELSEIF (ipot.EQ.3) THEN
            jspin1 = 2
            jspin2 = 1
            nc = nc_0
            i_start = mod(mod(n_rank - (nv(1)+nlotot),n_size) +
     +                                           n_size,n_size)
         ENDIF

         ai(:,:,:) = cmplx(0.,0.)
         bi(:,:,:) = cmplx(0.,0.)

         DO ik = 1,nv2(jspin1)
            DO jk = 1,nv2(jspin2)
               i3 = kvac3(ik,jspin1) - kvac3(jk,jspin2) 
               DO l = -vM,vM
                  DO m = -vM,vM
                   IF (l.EQ.m .OR. (iabs(m).LE.m_cyl 
     &                             .AND. iabs(l).LE.m_cyl)) THEN
                     ind1 = ig1(i3,m-l)
                     IF (ind1.NE.0) THEN
                        DO gi = 1,nvp(ik,jspin1)
                           i = ind(gi,ik,jspin1)
                           ai(l,jk,i) = ai(l,jk,i) +
     +                          conjg(a(m,i,jspin1))*tuuv(m,l,ik,jk)
     +                        + conjg(b(m,i,jspin1))*tudv(m,l,ik,jk)
                           bi(l,jk,i) = bi(l,jk,i) +
     +                          conjg(a(m,i,jspin1))*tduv(m,l,ik,jk)
     +                        + conjg(b(m,i,jspin1))*tddv(m,l,ik,jk)

                        END DO
                     END IF
                   END IF   ! noncyl. contributions
                  END DO
               END DO
            END DO
         END DO

         DO i = i_start+1, nv(jspin1), n_size
            ik = map1(i,jspin1)
            nc = nc + 1
            IF (ipot.EQ.1) THEN
               ii0 = nc*(nc-1)/2*n_size- (nc-1)*(n_size-n_rank-1)
               jmax = i
            ELSEIF (ipot.EQ.2) THEN
               ii0 = nc*(nc-1)/2*n_size- (nc-1)*(n_size-n_rank-1)
     +                                             + nv(1)+nlotot
               jmax = i
            ELSEIF (ipot.EQ.3) THEN
               ii0 = nc*(nc-1)/2*n_size- (nc-1)*(n_size-n_rank-1)
               jmax = nv(jspin2)
            ENDIF
            DO j = 1,jmax
               ii = ii0 + j
               jk = map1(j,jspin2)
               hij = cmplx(0.,0.)
               DO l = -vM,vM
                  hij = hij +
     +                 ai(l,jk,i)*a(l,j,jspin2) +
     +                 bi(l,jk,i)*b(l,j,jspin2)
               END DO   
#ifdef CPP_INVERSION
               aa(ii) = aa(ii) + real(hij)
#else 
               aa(ii) = aa(ii) + hij
#endif
            END DO
         END DO

      ENDDO !ipot

      IF (l_J) DEALLOCATE (vxy1)

      DEALLOCATE (ai,bi,nvp,ind,kvac3,map1,
     &    tddv,tduv,tudv,tuuv,a,b,bess,dbss,bess1,
     &    ddnv,dudz,duz,udz,uz )

      RETURN
      END SUBROUTINE od_hsvac
      END MODULE m_od_hsvac

      MODULE m_hsvac
      CONTAINS
      SUBROUTINE hsvac(
     >                 npotmatfile,nmzxyd,n2d,k1d,k2d,k3d,n3d,nmzd,jspd,
     >                 nwdd,matsize,nvd,nv2d,nvac,nmzxy,nmz,nlotot,
     >                 jsp,jspins,nw,vxy,vz,evac,bbmat,bmat,
     >                 bkpt,k1,k2,k3,kv2,nv,nq2,ig2,ig,rgphs,zrfs,
     >                 invs2,omtil,area,delz,z1,l_noco,l_J,qss,
     >                 n_size,n_rank,
     X                 aa,bb,
     <                 nv2)
c*********************************************************************
c     adds in the vacuum contributions to the the hamiltonian and
c     overlap matrices. as written, each k-point calculates the
c     vacuum functions again since there is only a single vacuum
c     parameter per vacuum.
c                m. weinert
c*********************************************************************
c     modified by R.Podloucky for speeding up and microtaskin
c*********************************************************************

      USE m_vacfun
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: matsize,nmzxyd,n2d,nvd,k1d,k2d,k3d,n3d
      INTEGER, INTENT (IN) :: nmzd,jspd,nwdd,nv2d,nmzxy,nmz,npotmatfile
      INTEGER, INTENT (IN) :: jsp,nw,jspins,nq2,nvac,n_size,n_rank
      INTEGER, INTENT (IN) :: nlotot
      LOGICAL, INTENT (IN) :: zrfs,invs2,l_noco,l_J
      REAL,    INTENT (IN) :: omtil,area,delz,z1
C     ..
C     .. Array Arguments ..
      COMPLEX, INTENT (INOUT) :: vxy(nmzxyd,n2d-1,2)
      INTEGER, INTENT (IN) :: k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      INTEGER, INTENT (IN) :: nv(jspd),kv2(2,n2d)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),ig2(n3d)
      INTEGER, INTENT (OUT):: nv2(jspd)
      REAL,    INTENT (IN) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      REAL,    INTENT (INOUT) :: vz(nmzd,2,4)
      REAL,    INTENT (IN) :: evac(2,jspd,nwdd)
      REAL,    INTENT (IN) :: bmat(3,3),bbmat(3,3),bkpt(3)
      REAL,    INTENT (IN) :: qss(3)
#ifdef CPP_INVERSION
      REAL,    INTENT (INOUT) :: aa(matsize),bb(matsize)
#else
      COMPLEX, INTENT (INOUT) :: aa(matsize),bb(matsize)
#endif
C     ..
C     .. Local Scalars ..
      COMPLEX hij,sij,apw_lo,c_1
      REAL d2,gz,sign,th,wronk,fac1
      INTEGER i,i2,ii,ik,j,jk,k,jspin,ipot,npot,ii0
      INTEGER ivac,irec,imz,igvm2,igvm2i
      INTEGER jspin1,jspin2,jmax,jsp_start,jsp_end
      INTEGER i_start,nc,nc_0
C     ..
C     .. Local Arrays ..
      INTEGER kvac1(nv2d,jspd),kvac2(nv2d,jspd),map2(nvd,jspd)
      COMPLEX tddv(nv2d,nv2d),tduv(nv2d,nv2d),
     +        tudv(nv2d,nv2d),tuuv(nv2d,nv2d),vxy_help(n2d-1)
      COMPLEX a(nvd,jspd),b(nvd,jspd)
      REAL ddnv(nv2d,jspd),dudz(nv2d,jspd),duz(nv2d,jspd),
     +     udz(nv2d,jspd),uz(nv2d,jspd)
! l_J auxiliary potential array
      COMPLEX, ALLOCATABLE :: vxy1(:,:,:)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cmplx,conjg,exp,real,sqrt
C     ..
      d2 = sqrt(omtil/area)

      IF (l_J) ALLOCATE (vxy1(nmzxyd,n2d-1,2))

c--->    set up mapping function from 3d-->2d lapws
      
      DO jspin = 1,jspins
         nv2(jspin) = 0
         DO 20 k = 1,nv(jspin)
            DO 10 j = 1,nv2(jspin)
               IF (k1(k,jspin).EQ.kvac1(j,jspin)
     +             .AND. k2(k,jspin).EQ.kvac2(j,jspin)) THEN
                  map2(k,jspin) = j
                  GO TO 20
               END IF
 10         CONTINUE
            nv2(jspin) = nv2(jspin) + 1
            IF (nv2(jspin).GT.nv2d) STOP 'hsvac: nv2d'
            kvac1(nv2(jspin),jspin) = k1(k,jspin)
            kvac2(nv2(jspin),jspin) = k2(k,jspin)
            map2(k,jspin) = nv2(jspin)
 20      CONTINUE
      ENDDO
c--->    loop over the two vacuua (1: upper; 2: lower)
      DO ivac = 1,2
         sign = 3. - 2.*ivac
         npot = 1
c---> pk non-collinear
         IF (l_noco) THEN
c--->       if the two vacuua are equivalent, the potential file has to
c--->       be backspaced, because the potential is the same at both
c--->       surfaces of the film
            IF ((ivac.EQ.2) .AND. (nvac.EQ.1)) THEN
               DO irec = 1,4
                  BACKSPACE (npotmatfile)
               ENDDO
            ENDIF
c--->       load the non-warping part of the potential
            READ (npotmatfile)((vz(imz,ivac,ipot),imz=1,nmzd),ipot=1,4)
            npot = 3
! for J-coeff. we average the up-up and down-down parts and off-diagonal elements of the
! potential matrix to zero
            IF (l_J) THEN
               vz(:,ivac,1) = (vz(:,ivac,1) + vz(:,ivac,2))/2.
               vz(:,ivac,2) =  vz(:,ivac,1)
               vz(:,ivac,3) = 0.0
               vz(:,ivac,4) = 0.0
            END IF
         ENDIF
c---> pk non-collinear

         DO ipot = 1,npot
c--->       get the wavefunctions and set up the tuuv, etc matrices
            IF (l_noco) THEN
               IF (.NOT.l_J) THEN
                  READ (npotmatfile)((vxy(imz,igvm2,ivac),
     +                           imz=1,nmzxy),igvm2=1,nq2-1)
               END IF
! l_J we want to average the diagonal elements of the potential matrix
               IF (l_J .AND. ipot.EQ.1) THEN
                  READ (npotmatfile)((vxy(imz,igvm2,ivac),
     +                           imz=1,nmzxy),igvm2=1,nq2-1)
                  READ (npotmatfile)((vxy1(imz,igvm2,ivac),
     +                           imz=1,nmzxy),igvm2=1,nq2-1)
                  vxy(:,:,ivac) = (vxy(:,:,ivac)+vxy1(:,:,ivac))/2.
               END IF
                
               IF (l_J .AND. ipot.EQ.3) THEN
                  READ (npotmatfile)((vxy(imz,igvm2,ivac),
     +                           imz=1,nmzxy),igvm2=1,nq2-1)
               END IF

               IF (nvac.EQ.1 .AND. ivac.EQ.2 .AND.(.NOT.zrfs) ) THEN
c--->          In this case (inversion, no z-reflection and thus no
c--->          2d-inversion) the coeff. of a 2d-star (containing G) of
c--->          the second vacuum is equal to the invers 2d-star
c--->          (containing -G) of the first vacuum.
! l_J no need to do this symmetrizaion twice, if the potential vxy is the same
                IF (.NOT.l_J .OR. (l_J .AND. ipot.EQ.1)) THEN
                  DO imz = 1,nmzxy
                     DO igvm2 = 2,nq2
c--->                   find the index of the invers 2d-star
                        igvm2i = ig2(ig(-kv2(1,igvm2),-kv2(2,igvm2),0))
                        vxy_help(igvm2-1) = vxy(imz,igvm2i-1,2)
                     ENDDO
                     DO igvm2 = 2,nq2
                        vxy(imz,igvm2-1,2) = vxy_help(igvm2-1)
                     ENDDO
                  ENDDO
                END IF ! l_J
               ENDIF ! ivac-nvac
! l_J we want the off-diagonal potential matrix elements to be zero
               IF (l_J .AND. ipot.EQ.3) vxy(:,:,ivac)=cmplx(0.,0.)
            ENDIF
            CALL vacfun(
     >               nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n2d,n3d,jspd,
     >               jsp,jspins,l_noco,qss,ipot,
     >               nmzxy,nmz,invs2,delz,ig2,ig,rgphs,
     >               bbmat,ivac,evac(1,1,nw),bkpt,
     >               vxy(1,1,ivac),vz,kvac1,kvac2,nv2,
     <               tuuv,tddv,tudv,tduv,uz,duz,udz,dudz,ddnv,wronk)
            fac1 = 1.0 / (d2*wronk)
c
c--->       generate a and b coeffficients
c
            IF (l_noco) THEN
               DO jspin = 1,jspins
                  DO k = 1,nv(jspin)
                     gz = sign*bmat(3,3)*k3(k,jspin)
                     i2 = map2(k,jspin)
                     th = gz*z1
                     c_1 = fac1 * cmplx( cos(th), sin(th) )
                     a(k,jspin) = - c_1 * cmplx( 
     +                              dudz(i2,jspin), gz*udz(i2,jspin) )
                     b(k,jspin) =   c_1 * cmplx( 
     +                               duz(i2,jspin), gz* uz(i2,jspin) )
                  ENDDO
               ENDDO
            ELSE
               DO k = 1,nv(jsp)
                  gz = sign*bmat(3,3)*k3(k,jsp)
                  i2 = map2(k,jsp)
                  th = gz*z1
                  c_1 = fac1 * cmplx( cos(th), sin(th) )
                  a(k,1) = - c_1 * cmplx( 
     +                           dudz(i2,jsp), gz*udz(i2,jsp) )
                  b(k,1) =   c_1 * cmplx( 
     +                            duz(i2,jsp), gz* uz(i2,jsp) )
               ENDDO
            ENDIF
c--->       update hamiltonian and overlap matrices
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
     +                                              n_size,n_size) 
            ENDIF
c-gb||
            DO 50 i = i_start+1, nv(jspin), n_size
               ik = map2(i,jspin)
               nc = nc + 1
               IF (ipot.EQ.1) THEN
                  jspin = 1
                  ii0 = nc*(nc-1)/2*n_size- (nc-1)*(n_size-n_rank-1)
               ELSEIF (ipot.EQ.2) THEN
                  jspin = 2
                  ii0=nc*(nc-1)/2*n_size-(nc-1)*(n_size-n_rank-1)+
     +                                                nv(1)+nlotot
               ENDIF
               jspin1 = jsp
               IF (l_noco) jspin1 = jspin
               DO j = 1,i - 1
                  ii = ii0 + j
c--->             overlap: only  (g-g') parallel=0       '
                  IF (map2(j,jspin).EQ.ik) THEN
                     sij = conjg(a(i,jspin))*a(j,jspin) 
     +                   + conjg(b(i,jspin))*b(j,jspin)*ddnv(ik,jspin1)
c+APW_LO
#ifdef CPP_APW
                     apw_lo = conjg(a(i,jspin)*  uz(ik,jspin1) +
     +                              b(i,jspin)* udz(ik,jspin1) ) 
     +                           * (a(j,jspin)* duz(ik,jspin1) +
     +                              b(j,jspin)*dudz(ik,jspin1) )
     +                      +      (a(j,jspin)*  uz(ik,jspin1) +
     +                              b(j,jspin)* udz(ik,jspin1) ) 
     +                      * conjg(a(i,jspin)* duz(ik,jspin1) +
     +                              b(i,jspin)*dudz(ik,jspin1) )
c            IF (i.lt.10) write (3,'(2i4,2f20.10)') i,j,apw_lo
#ifdef CPP_INVERSION
                     aa(ii) = aa(ii) + 0.25 * real(apw_lo) 
#else 
                     aa(ii) = aa(ii) + 0.25 * apw_lo
#endif
#endif
c+APW_LO
#ifdef CPP_INVERSION
                     bb(ii) = bb(ii) + real(sij)
#else 
                     bb(ii) = bb(ii) + sij
#endif
                  END IF
               ENDDO
               ii = ii0 + i
               sij = conjg(a(i,jspin))*a(i,jspin) 
     +             + conjg(b(i,jspin))*b(i,jspin)*ddnv(ik,jspin1)
#ifdef CPP_INVERSION
               bb(ii) = bb(ii) + real(sij)
#else
               bb(ii) = bb(ii) + sij
#endif
 50         CONTINUE
            ENDIF

c--->    hamiltonian update
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
     +                                              n_size,n_size) 
            ELSEIF (ipot.EQ.3) THEN
               jspin1 = 2
               jspin2 = 1
               nc = nc_0
               i_start = mod(mod(n_rank - (nv(1)+nlotot),n_size) +
     +                                              n_size,n_size) 
            ENDIF
            DO i = i_start+1, nv(jspin1), n_size
               ik = map2(i,jspin1)
               nc = nc + 1
               IF (ipot.EQ.1) THEN
                  ii0 = nc*(nc-1)/2*n_size- (nc-1)*(n_size-n_rank-1)
                  jmax = i
               ELSEIF (ipot.EQ.2) THEN
                  ii0 = nc*(nc-1)/2*n_size- (nc-1)*(n_size-n_rank-1)
     +                                                + nv(1)+nlotot
                  jmax = i
               ELSEIF (ipot.EQ.3) THEN
                  ii0 = nc*(nc-1)/2*n_size- (nc-1)*(n_size-n_rank-1)
                  jmax = nv(jspin2)
               ENDIF
               DO j = 1,jmax
                  ii = ii0 + j
                  jk = map2(j,jspin2)
                  hij = conjg(a(i,jspin1))* (tuuv(ik,jk)*a(j,jspin2)
     +                                      +tudv(ik,jk)*b(j,jspin2))
     +                + conjg(b(i,jspin1))* (tddv(ik,jk)*b(j,jspin2)
     +                                      +tduv(ik,jk)*a(j,jspin2))
#ifdef CPP_INVERSION
                  aa(ii) = aa(ii) + real(hij)
#else
                  aa(ii) = aa(ii) + hij
#endif
               ENDDO
            ENDDO

c--->    end of loop over different parts of the potential matrix
         ENDDO

c---> end of loop over vacua
      ENDDO

      IF (l_J) DEALLOCATE (vxy1)

      END SUBROUTINE hsvac
      END MODULE m_hsvac

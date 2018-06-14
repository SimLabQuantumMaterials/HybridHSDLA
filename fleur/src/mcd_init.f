      SUBROUTINE mcd_init(
     >                    lmaxd,jspins,jmtd,msh,ntypd,nstd,
     >                    jri,lmax,dx,zatom,rmsh,ncst,
     >                    vr,g,f,emcd_up,emcd_lo,itype,jspin,
     <                    ncore,e_mcd,m_mcd)

!-----------------------------------------------------------------------
!
! For a given atom-type 'itype' look, whether a core state is in the
! energy interval [emcd_lo,emcd_up] and, if found, calculate the 
! MCD-matrix elements 'm_mcd'.
!          
!-----------------------------------------------------------------------

      USE m_nabla
      USE m_dr2fdr
      USE m_constants, ONLY : c_light
      USE m_setcor
      USE m_differ
      IMPLICIT NONE

      INTEGER, PARAMETER :: l_max = 3

! Arguments ...

      INTEGER, INTENT (IN)  :: lmaxd,jspins,jmtd,msh,ntypd,nstd,itype
      INTEGER, INTENT (IN)  :: jspin
      REAL,    INTENT (IN)  :: emcd_up,emcd_lo
      INTEGER, INTENT (IN)  :: jri(ntypd),lmax(ntypd),ncst(ntypd)
      REAL,    INTENT (IN)  :: dx(ntypd),zatom(ntypd)
      REAL,    INTENT (IN)  :: rmsh(jmtd,ntypd),vr(jmtd,ntypd,jspins)
      REAL,    INTENT (IN)  :: f(jmtd,2,0:lmaxd,jspin:jspin)
      REAL,    INTENT (IN)  :: g(jmtd,2,0:lmaxd,jspin:jspin)
      INTEGER, INTENT (OUT) :: ncore(ntypd)
      REAL,    INTENT (OUT) :: e_mcd(ntypd,jspins,nstd)
      COMPLEX, INTENT (OUT) :: m_mcd(nstd,(3+1)**2,3*ntypd,2)

! Locals ...

      INTEGER kap,mue,iri,l,ispin,i,icore,korb,nst,n_core,ierr
      REAL c,t2,e,fj,fl,fn,weight,d,ms,rn,bmu
      INTEGER kappa(nstd),nprnc(nstd),l_core(nstd)
      REAL vrd(msh),occ(nstd),a(msh),b(msh),j_core(nstd),e_mcd1(nstd)
      REAL gv1(jmtd)
      REAL, ALLOCATABLE :: gc(:,:,:),fc(:,:,:)
      REAL, ALLOCATABLE :: gv(:,:,:,:),fv(:,:,:,:),dgv(:,:,:,:)

!-----------------------------------------------------------------------

      c = c_light(1.0)
      ALLOCATE ( gc(jri(itype),ncst(itype),jspins) )
      ALLOCATE ( fc(jri(itype),ncst(itype),jspins) )
      
! core setup

      ncore(itype) = 0
      bmu = 0.0
      CALL setcor(
     >            zatom(itype),nstd,1,1,bmu,
     <            nst,kappa,nprnc,occ)

      DO ispin = jspin, jspin

! extend core potential

        DO iri = 1, jri(itype)
          vrd(iri) = vr(iri,itype,ispin)
        ENDDO      
        t2 = vrd(jri(itype)) / (jri(itype) - msh)
        DO iri = jri(itype) + 1, msh
          vrd(iri) =  vrd(jri(itype))  + t2* ( iri-jri(itype) )
        ENDDO

! calculate core

        n_core = 0
        DO korb = 1, ncst(itype)
          IF (occ(korb).GT.0) THEN
            fn = nprnc(korb)
            fj = iabs(kappa(korb)) - .5e0
            weight = 2*fj + 1.e0
            fl = fj + (.5e0)*isign(1,kappa(korb))
            e = -2* (zatom(itype)/ (fn+fl))**2
            d = exp(dx(itype))
            rn = rmsh(1,itype)*( d**(msh-1) )
            CALL differ(
     >                  fn,fl,fj,c,zatom(itype),dx(itype),rmsh(1,itype),
     >                  rn,d,msh,vrd,
     X                  e,
     <                  a,b,ierr)
            IF (ierr.NE.0) STOP 'mcd_init: error in core-levels'
            IF ( (e.LE.emcd_up).AND.(e.GE.emcd_lo) ) THEN
              write(*,*) 'good    ev = ',e
              n_core = n_core + 1
              j_core(n_core) = fj
              l_core(n_core) = NINT( fl )
              e_mcd1(n_core) = e
              DO iri = 1, jri(itype)
                gc(iri,n_core,ispin) = a(iri)
                fc(iri,n_core,ispin) = b(iri)
              ENDDO
            ENDIF
          ENDIF
        ENDDO

      ENDDO

!-----------------------------------------------------------------------

      IF (n_core.GT.0) THEN
      
        ALLOCATE ( gv(jri(itype),0:l_max,jspins,2) )
        ALLOCATE (dgv(jri(itype),0:l_max,jspins,2) )
        ALLOCATE ( fv(jri(itype),0:l_max,jspins,2) )
        DO i = 1, 2
         DO iri = 3*(itype-1)+1 , 3*(itype-1)+3
           DO l = 1, (l_max+1)**2
             DO icore = 1, nstd
               m_mcd(icore,l,iri,i) = cmplx(0.0,0.0)
             ENDDO
           ENDDO
         ENDDO
        ENDDO
c
c bring LAPW wavefunctions in a proper form:
c
        DO ispin = jspin, jspin
          ms = ispin - 1.5
          DO l = 0, l_max
            DO iri = 1, jri(itype)
              gv(iri,l,ispin,1) = f(iri,1,l,ispin)   ! large component of u
              fv(iri,l,ispin,1) = f(iri,2,l,ispin)   ! small              .
              gv(iri,l,ispin,2) = g(iri,1,l,ispin)   ! large component of u
              fv(iri,l,ispin,2) = g(iri,2,l,ispin)   ! small
            ENDDO 
            gv1(:) = rmsh(:,itype) * gv(:,l,ispin,1)
            CALL dr2fdr(                                          ! deriative of u (large)
     >                  gv1,rmsh(1,itype),jri(itype),
     <                  dgv(1,l,ispin,1) )
            gv1(:) = rmsh(:,itype) * gv(:,l,ispin,2)              !              .
            CALL dr2fdr(                                          ! deriative of u (large)
     >                  gv1,rmsh(1,itype),jri(itype),
     <                  dgv(1,l,ispin,2) )
          ENDDO 
c
c
c
          DO icore = 1, n_core

            DO i = 1, 2
c              write(*,*) j_core(icore),l_core(icore),l_max,ms
              CALL nabla(
     >                 itype,icore,jri(itype),dx(itype),nstd,ntypd,
     >                 j_core(icore),l_core(icore),l_max,ms,
     >                 rmsh(1,itype),gc(1,icore,ispin),
     >                 gv(1,0,ispin,i),dgv(1,0,ispin,i),
     <                 m_mcd(1,1,1,i) )
            ENDDO

            DO i = 1, 2*icore*l_core(icore)
              ncore(itype) = ncore(itype) + 1
              IF (ncore(itype).GT.nstd) STOP 'mcd_init: nstd too small'
              e_mcd(itype,ispin,ncore(itype)) = e_mcd1(icore)
            ENDDO
          ENDDO
        ENDDO 

        DEALLOCATE (gv,fv,dgv)
      ENDIF
      DEALLOCATE (gc,fc)

       
c      DO i = 1, 2
c       DO iri = 3*(itype-1)+1 , 3*(itype-1)+3
c         write (*,*) iri
c         DO icore = 1, ncore(itype)
c           write (*,'(10f10.5)') (m_mcd(icore,l,iri,i),l=1,9)
c         ENDDO
c       ENDDO
c      ENDDO
      END SUBROUTINE mcd_init

      MODULE m_usetup
c-------------------------------------------------------------------+
c Sets up the quantities needed for the LDA+U subroutines:          |
c     radial integrals: us,dus,uds,duds                             |
c     overlap of dot u: ddn                                         |
c     potential matrix: vs_mmp                                      |
c     total energy contribution: e_ldau                             |
C                                                  G.B. Oct. 2000   |
c-------------------------------------------------------------------+
      CONTAINS
      SUBROUTINE u_setup(
     >                   jspd,ntype,lmaxd,lmaxb,jmtd,nlhd,nw,n_u,
     >                   jspins,lda_u,neq,el,vr,jri,dx,rmsh,irank,
     <                   vs_mmp,e_ldau)

      USE m_umtx
      USE m_uj2f
      USE m_nmat_rot
      USE m_vmmp
      USE m_radfun
      USE m_types, ONLY : t_utype

      IMPLICIT NONE

C ... Arguments ...
      INTEGER, INTENT (IN) :: jspd,ntype,lmaxd,lmaxb,jmtd,nlhd
      INTEGER, INTENT (IN) :: jspins,nw,irank
      INTEGER, INTENT (INOUT) :: n_u
      REAL,    INTENT (OUT):: e_ldau

      TYPE (t_utype), INTENT (IN) :: lda_u(ntype)

      INTEGER, INTENT (IN) :: jri(ntype),neq(ntype)
      REAL,    INTENT (IN) :: rmsh(jmtd,ntype),el(0:lmaxd,ntype,jspd)
      REAL,    INTENT (IN) :: dx(ntype),vr(jmtd,0:nlhd,ntype,jspd)
      COMPLEX,INTENT (OUT)::vs_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb,n_u,jspins)

C ... Local Variables ...
      INTEGER itype,ispin,j,k,l,n,jspin,urec
      INTEGER noded,nodeu
      REAL wronk,theta,phi
      LOGICAL n_mmp_exist,n_exist
      CHARACTER*8 l_type*2,l_form*9

      REAL f(jmtd,2),g(jmtd,2)
      REAL f0(n_u,jspins),f2(n_u,jspins),f4(n_u,jspins),f6(n_u,jspins)
      REAL, ALLOCATABLE :: u(:,:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: ns_mmp(:,:,:,:)
c
c look, whether density matrix exists already:
c
      INQUIRE (file='n_mmp_mat',exist=n_mmp_exist)
      IF (n_mmp_exist) THEN
c
c calculate slater integrals from u and j
c
         CALL uj2f(
     >             jspins,ntype,n_u,lda_u,
     <             f0,f2,f4,f6)
c
c set up e-e- interaction matrix
c
         ALLOCATE ( u(-lmaxb:lmaxb,-lmaxb:lmaxb,
     +                -lmaxb:lmaxb,-lmaxb:lmaxb,n_u,jspins) )
         DO ispin = 1, 1 ! jspins
         f0(:,1) = (f0(:,1) + f0(:,jspins) ) / 2
         f2(:,1) = (f2(:,1) + f2(:,jspins) ) / 2
         f4(:,1) = (f4(:,1) + f4(:,jspins) ) / 2
         f6(:,1) = (f6(:,1) + f6(:,jspins) ) / 2
           CALL umtx(
     >               lmaxb,ntype,n_u,lda_u(:)%l,f0(1,ispin),
     >               f2(1,ispin),f4(1,ispin),f6(1,ispin),
     <               u(-lmaxb,-lmaxb,-lmaxb,-lmaxb,1,ispin) )
         ENDDO
c
c read density matrix
c
         ALLOCATE (ns_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb,n_u,jspins))
         OPEN (69,file='n_mmp_mat',status='unknown',form='formatted')
         READ (69,9000) ns_mmp
 9000    FORMAT(7f20.13)
         CLOSE (69)
!
! check for possible rotation of n_mmp
!
      INQUIRE (file='n_mmp_rot',exist=n_exist)
      IF (n_exist) THEN
         OPEN (68,file='n_mmp_rot',status='old',form='formatted')
         READ(68,*) theta,phi
         CLOSE (68)
         CALL nmat_rot(phi,theta,0.0,3,n_u,jspins,
     X                 ns_mmp)
      ENDIF
c
c calculate potential matrix and total energy correction
c
         CALL v_mmp(
     >              ntype,n_u,jspins,lmaxb,lda_u,neq,ns_mmp,u,f0,f2,
     <              vs_mmp,e_ldau)
         IF (irank.eq.0) THEN
           DO jspin = 1,jspins
             n = 0
             WRITE (6,'(a7,i3)') 'spin #',jspin
             DO itype = 1,ntype
                l = lda_u(itype)%l
                IF (l.GE.0) THEN
                  n = n + 1
                  WRITE (l_type,'(i2)') 2*(2*l+1)
                  l_form = '('//l_type//'f12.7)'
                  WRITE (6,'(a20,i3)') 'n-matrix for atom # ',itype
                  WRITE (6,l_form) ((ns_mmp(k,j,n,jspin),k=-l,l),j=-l,l)
                  WRITE (6,'(a20,i3)') 'V-matrix for atom # ',itype
                  IF (lda_u(itype)%l_amf) THEN
                    WRITE (6,*) 'using the around-mean-field limit '
                  ELSE
                    WRITE (6,*) 'using the atomic limit of LDA+U '
                  ENDIF
                  WRITE (6,l_form) ((vs_mmp(k,j,n,jspin),k=-l,l),j=-l,l)
                ENDIF
             ENDDO
           ENDDO
           WRITE (6,*) e_ldau
         ENDIF
         DEALLOCATE ( u,ns_mmp )
      ELSE
         IF (irank.eq.0) 
     +      WRITE (*,*) 'no density matrix found ... skipping LDA+U'
         DO jspin = 1,jspins
           DO  n = 1, n_u
             DO j = -lmaxb,lmaxb
               DO k = -lmaxb,lmaxb
                 vs_mmp(k,j,n,jspin) = cmplx(0.0,0.0)
               ENDDO
             ENDDO
           ENDDO
         ENDDO
         e_ldau = 0.0
         n_u = 0
      ENDIF

      RETURN
      END SUBROUTINE u_setup
      END MODULE m_usetup

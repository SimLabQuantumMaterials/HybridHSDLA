      MODULE m_sympsi

c Calculates the irreducible represetantions of the wave functions.
c if k-point is in Brillouin zone boundary results are correct only for
c non-symmorphic groups (factor groups would be needed for that...). 
c jsym contains the number of irreducible rep., corresponding character
c tables are given in the file syminfo.
c
c Double groups work only with non-collinear calculations, for normal spin-orbit 
c calculations both spin up and down components would be needed...

c Jussi Enkovaara, Juelich 2004

      CONTAINS
      SUBROUTINE sympsi(
     >                  bkpt,nv,kx,ky,kz,nop,neigd,
     >                  ne,nvd,mrot,amat,nbasfcn,z,eig,l_noco,l_soc,
     <                  ksym,jsym)

      USE m_grp_k
      USE m_inv3

      IMPLICIT NONE
c
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: neigd,nvd,nbasfcn,nop
      INTEGER, INTENT (IN) :: nv,ne
      LOGICAL, INTENT(IN)  :: l_soc,l_noco
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: kx(nvd),ky(nvd),kz(nvd),mrot(3,3,nop)
      REAL,    INTENT (IN) :: bkpt(3),eig(neigd),amat(3,3)
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
      REAL,    INTENT (IN) :: z(nbasfcn,neigd)
#else
      COMPLEX, INTENT (IN) :: z(nbasfcn,neigd)
#endif
      INTEGER, INTENT (OUT):: jsym(neigd),ksym(neigd)
C     ..
C     .. Local Scalars ..
      REAL small,degthre
      INTEGER i,k,n
      INTEGER nclass,nirr,n1,n2,c,ndeg
      LOGICAL soc, char_written
C     ..
C     .. Local Arrays ..
      INTEGER mrot_k(3,3,2*nop)
      INTEGER :: mtmpinv(3,3),d
      INTEGER :: gmap(nvd,nop)
      REAL ::    kv(3),kvtest(3)
      INTEGER :: deg(ne),begdeg(0:ne),enddeg(0:ne),nndeg(ne)
      
      REAL :: norm(ne)
      LOGICAL :: symdone(ne)

#ifdef CPP_INVERSION
      REAL, ALLOCATABLE :: csum(:,:,:),overlap(:,:),chars(:,:)
      REAL, SAVE,ALLOCATABLE :: char_table(:,:)
#else
      COMPLEX, ALLOCATABLE :: csum(:,:,:),chars(:,:)
      COMPLEX, SAVE, ALLOCATABLE :: char_table(:,:)
#endif
      CHARACTER(LEN=5) :: grpname,irrname(2*nop)
      COMPLEX          :: c_table(2*nop,2*nop)
      COMPLEX, ALLOCATABLE :: su(:,:,:)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,real
C     ..
C     .. Data statements ..
c
      DATA small/1.0e-4/

      soc=l_soc.AND.l_noco
      jsym=0
      ksym=0
      IF (l_soc.AND.(.NOT.l_noco)) RETURN
      
      IF (soc) THEN
         ALLOCATE(su(2,2,2*nop))
         CALL grp_k(mrot,mrot_k,amat,bkpt,nclass,nirr,c_table,
     &     grpname,irrname,su)
      ELSE
         CALL grp_k(mrot,mrot_k,amat,bkpt,nclass,nirr,c_table,
     &        grpname,irrname)
      ENDIF
      ALLOCATE(csum(ne,ne,nclass))
      ALLOCATE(chars(ne,nclass))
      chars=0.0
      !>

      IF (ALLOCATED(char_table)) THEN
         IF (SIZE(char_table,2).NE.nclass) THEN
            DEALLOCATE(char_table)
            ALLOCATE(char_table(nirr,nclass))
            char_written=.FALSE.
         ENDIF
      ELSE
         ALLOCATE(char_table(nirr,nclass))
         char_written=.FALSE.
      ENDIF
      char_table=c_table(1:nirr,1:nclass)

      !<--map the (k+g)-vectors related by inv(rot)
      gmap=0
      DO c=1,nclass
         CALL inv3(mrot_k(:,:,c),mtmpinv,d)
         kloop: DO k=1,nv
         kv(1)=kx(k)
         kv(2)=ky(k)
         kv(3)=kz(k)
         kv=kv+bkpt
         kvtest=MATMUL(kv,mtmpinv)
c         kvtest=MATMUL(kv,mrot_k(:,:,c))
         DO i = 1,nv
            kv(1)=kx(i)
            kv(2)=ky(i)
            kv(3)=kz(i)
            kv=kv+bkpt
            IF (ABS(kvtest(1)-kv(1)).LT.small.AND.
     &           ABS(kvtest(2)-kv(2)).LT.small.AND.
     &           ABS(kvtest(3)-kv(3)).LT.small) THEN
               gmap(k,c)=i
               cycle kloop
            ENDIF
         ENDDO
         WRITE(6,*) 'Problem in symcheck, cannot find rotated kv for'
     &          , k,kx(k),ky(k),kz(k)
         RETURN
         ENDDO kloop
      ENDDO

!norms
      DO i=1,ne
         norm(i)=0.0
         IF (soc) THEN
            DO k=1,nv*2
               norm(i)=norm(i)+ABS(z(k,i))**2
            ENDDO
         ELSE
            DO k=1,nv
               norm(i)=norm(i)+ABS(z(k,i))**2
            ENDDO
         ENDIF
         norm(i)=sqrt(norm(i))
      ENDDO


      !<-- Calculate the characters
      symdone=.FALSE.
      stateloop: DO i=1,ne
            IF (symdone(i)) CYCLE stateloop
            ndeg=0
            deg=0
            degthre=0.0001
            DO n=1,ne
               IF (ABS(eig(i)-eig(n)).LT.degthre) THEN
                  ndeg=ndeg+1
                  deg(ndeg)=n
               ENDIF
            ENDDO

            csum=0.0
            DO c=1,nclass
               DO n1=1,ndeg
                  DO n2=1,ndeg
                     DO k=1,nv
#ifdef CPP_SOC
                        IF (soc) THEN
                       csum(n1,n2,c)=csum(n1,n2,c)+(conjg(z(k,deg(n1)))*
     &                   (su(1,1,c)*z(gmap(k,c),deg(n2))+
     &                    su(1,2,c)*z(gmap(k,c)+nv,deg(n2)))+
     &                    conjg(z(k+nv,deg(n1)))*
     &                   (su(2,1,c)*z(gmap(k,c),deg(n2))+
     &                    su(2,2,c)*z(gmap(k,c)+nv,deg(n2))))/
     &                      (norm(deg(n1))*norm(deg(n2)))
                        ELSE
#endif
#ifdef CPP_INVERSION
                         csum(n1,n2,c)=csum(n1,n2,c)+z(k,deg(n1))*
     &                z(gmap(k,c),deg(n2))/(norm(deg(n1))*norm(deg(n2)))
#else
                       csum(n1,n2,c)=csum(n1,n2,c)+conjg(z(k,deg(n1)))*
     &                z(gmap(k,c),deg(n2))/(norm(deg(n1))*norm(deg(n2)))
#endif
#ifdef CPP_SOC
                        ENDIF
#endif
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
               
! We might have taken degenerate states which are not degenerate due to symmetry
! so look for irreducible reps
        DO n1=1,ndeg
           chars(deg(n1),:)=0.0
           DO n2=1,ndeg
              IF (ANY(ABS(csum(n1,n2,:)).GT.0.01)) THEN
                 chars(deg(n1),:)=chars(deg(n1),:)+csum(n2,n2,:)
              ENDIF
           ENDDO
           symdone(deg(n1))=.TRUE.
        ENDDO
            

! determine the irreducible presentation
         irrloop: DO n1=1,ndeg 
!         write(*,'(i3,6(2f6.3,2x))') i,chars(deg(n1),1:nclass)
           DO c=1,nirr
              IF (ALL(ABS(chars(deg(n1),1:nclass)-
     &             char_table(c,1:nclass)).LT.0.001)) THEN
                 jsym(deg(n1))=c
                 cycle irrloop
              ELSE IF (ALL(ABS(char_table(c,1:nclass)).LT.0.001)) THEN
                 char_table(c,:)=chars(deg(n1),:)
                 jsym(deg(n1))=c
                 cycle irrloop
              ENDIF
           ENDDO
        ENDDO irrloop

      ENDDO stateloop
      !>

      IF (.NOT.char_written) THEN
         WRITE(444,124) bkpt
         WRITE(444,*) 'Group is ' ,grpname
         DO c=1,nirr
#ifdef CPP_INVERSION
            IF (ANY((char_table).GT.0.001)) THEN
#else
            IF (ANY(IMAG(char_table).GT.0.001)) THEN
#endif
          WRITE(444,123) c,irrname(c),(char_table(c,n),n=1,nclass)
          ELSE
          WRITE(444,123) c,irrname(c),(REAL(char_table(c,n)),n=1,nclass)
          ENDIF
         ENDDO
         char_written=.TRUE.
      ENDIF
 123  FORMAT(i3,1x,a5,1x,20f7.3)
 124  FORMAT('Character table for k: ',3f8.4)

      DEALLOCATE(csum)
      DEALLOCATE(chars)
      
      END SUBROUTINE

      END  MODULE m_sympsi

      MODULE m_doswrite
c
c-- now write cdninf for all kpts if on T3E
c-- now read data from tmp_dos and write to vacdos&dosinp .. dw
c
      CONTAINS
      SUBROUTINE doswrite(
     >                   neigd,nkptd,nvd,ntypd,jspd,layerd,
     >                   jspins,semic,dos,vacdos,cdinf,film,
     >                   nkpt,ntype,layers,nvac,slice,l_noco,nop2,
     >                   ndir,invs,zrfs,bmat,neq,izlay,
     >                   n2max,reclength_vw,nstm,tworkf
     >                   ,nstars,starcoeff,l_mcd,ncored,ncore,e_mcd,
     >                   emax,emin,efermi,sigma,nsld,natd,odi)

      USE m_evaldos
      USE m_od_types, ONLY : od_inp
      IMPLICIT NONE

C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: neigd,nkptd,nvd,ntypd,jspd,layerd
      INTEGER, INTENT (IN) :: jspins,layers,ndir,nvac,nkpt
      INTEGER, INTENT (IN) :: nstars,ntype,n2max,nop2,nsld,natd
      INTEGER, INTENT (IN) :: reclength_vw,nstm,ncored
      REAL,    INTENT (IN) :: tworkf,emax,emin,efermi,sigma
      LOGICAL, INTENT (IN) :: dos,film,invs,semic,l_noco,slice,zrfs
      LOGICAL, INTENT (IN) :: starcoeff,vacdos,cdinf,l_mcd
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN)  :: neq(ntypd),izlay(layerd,2),ncore(ntypd)
      REAL, INTENT(IN)      :: bmat(3,3),e_mcd(ntypd,jspins,ncored)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
c+odim

c    locals
      INTEGER :: jsym(neigd),ksym(neigd)
      REAL    :: wk,bkpt(3)
      REAL   :: eig(neigd)
      REAL   :: qal(0:3,ntypd,neigd,jspd)
      REAL   :: qis(neigd,nkptd,jspd)
      REAL   :: qvac(neigd,2,nkptd,jspd)
      REAL   :: qvlay(neigd,layerd,2)
      COMPLEX :: qstars(nstars,neigd,layerd,2)
      INTEGER :: ne,ikpt,kspin,j,i,n
      COMPLEX, ALLOCATABLE :: ac(:,:),bc(:,:)

      EXTERNAL cdninf

c     check if there is anything todo here
      IF (.NOT.(dos.OR.cdinf.OR.vacdos.OR.(nstm.EQ.3))) RETURN
c     check if settings in inp-file make any sense
      IF (vacdos.and..not.dos) THEN
         write(6,*) "STOP DOS: only set vacdos=.true. if dos=.true."
         STOP "DOS"
      ENDIF
      IF (vacdos.and.(.not.starcoeff.and.(nstars.ne.1)))THEN
          write(6,*) "STOP DOS: if stars=f set nstars=1"
         STOP "DOS"
      ENDIF

      IF (dos.AND.(ndir.GE.0)) THEN
c---  >    open files for bandstucture+ old style vacdos
         OPEN (85,file='dosinp')
         IF (vacdos) THEN
            OPEN (86,file='vacDOS')
         ENDIF
      ENDIF

      IF ((dos.AND.(ndir.GE.0)).OR.cdinf) THEN
c
c      write bandstructure or cdn-info to output-file
c         
         DO kspin = 1,jspins
            IF (dos.AND.(ndir.GE.0)) THEN
c---  >       write header information to vacdos & dosinp
c     
               IF (film) THEN
                  WRITE (85,FMT=8080) nvac,nkpt
               ELSE
                  WRITE (85,FMT=8080) jspins,nkpt
               ENDIF
 8080          FORMAT (12i6)
               WRITE (85,FMT=8080) ntype, (neq(n),n=1,ntype)
               IF (vacdos) THEN
                  WRITE (86,FMT=8080) nvac,nkpt
                  WRITE (86,FMT=8080) layers
                  WRITE (86,'(20(i3,1x))') (izlay(i,1),i=1,layers)
               ENDIF
            ENDIF

            DO ikpt=1,nkpt
               READ (84,rec=nkpt*(kspin-1)+ikpt) bkpt,wk,ne,eig,
     +              qal(:,:,:,kspin),qvac(:,:,ikpt,kspin),
     +              qis(:,ikpt,kspin),
     +              qvlay(:,:,:),qstars,ksym,jsym
               
               CALL cdninf(
     >              film,invs,zrfs,semic,l_noco,kspin,ntype,neq,
     +              nvac,slice,dos,ndir,vacdos,layers,ikpt,bkpt,
     >              wk,bmat,neigd,nkptd,nvd,ntypd,jspd,layerd,
     +              ne,eig,qal(0,1,1,kspin),qis,qvac,
     +              qvlay(:,:,:),
     >              qstars,nstars,starcoeff,ksym,jsym)
            ENDDO
            
         ENDDO                  ! end spin loop (kspin = 1,jspins)

      ENDIF

      IF (dos.AND.(ndir.GE.0)) THEN
         CLOSE(85)
         RETURN
c     ok, all done in the bandstructure/cdninf case
      ENDIF
c
c     write DOS/VACDOS
c
c     
      IF (dos.AND.(ndir.LT.0)) THEN
         CALL evaldos(
     >                film,vacdos,nstars,nvac,nkpt,nkptd,ntype,layers,
     >                jspins,neigd,ntypd,layerd,neq,emin,emax,
     >                sigma,efermi,nop2,invs,l_noco,odi,
     >                l_mcd,ncored,ncore,e_mcd,ndir,natd,nsld,bmat)
      ENDIF
c
c      Now write to vacwave if nstm=3 
c     all data
c     has been written to tmp_vacwave and must be written now
c     by PE=0 only!
c
      IF (nstm.EQ.3) THEN

         OPEN (89,file='tmp_vacwave',status='old',access='direct',
     +                                          recl=reclength_vw)
         ALLOCATE ( ac(n2max,neigd),bc(n2max,neigd) )
         DO ikpt = 1,nkpt
            WRITE(*,*) 'Read rec',ikpt,'from vacwave'
            READ(89,rec=ikpt) wk,ne,bkpt(1),bkpt(2),
     +                 eig,ac,bc
            WRITE (87,'(i3,1x,f12.6)') ikpt,wk
            i=0
            DO n = 1, ne
               IF (ABS(eig(n)-tworkf).LE.emax) i=i+1
            END DO
            WRITE (87,FMT=990) bkpt(1), bkpt(2), i, n2max
            DO n = 1, ne
               IF (ABS(eig(n)-tworkf).LE.emax) THEN
                  WRITE (87,FMT=1000) eig(n)
                  DO j=1,n2max
                     WRITE (87,FMT=1010) ac(j,n),bc(j,n)
                  END DO
               END IF
            END DO
 990        FORMAT(2(f8.4,1x),i3,1x,i3)
 1000       FORMAT(e10.4)
 1010       FORMAT(2(2e20.8,1x))
         END DO
         DEALLOCATE ( ac,bc )
c
         CLOSE(89)

      ENDIF

      RETURN
      END SUBROUTINE doswrite
      END MODULE m_doswrite

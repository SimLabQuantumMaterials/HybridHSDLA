      SUBROUTINE cdninf(
     >                  film,invs,zrfs,semic,l_noco,jspin,ntyp,neq,nvac,
     +                  slice,dos,ndir,vacdos,layers,ikpt,bkpt,wk,
     >                  bmat,neigd,nkptd,nvd,ntypd,jspd,layerd,
     +                  nbands,eig,qal,qis,qvac,qvlay,
     x                  qstars,nstars,starcoeff,jsym,ksym)
c***********************************************************************
c     this subroutine calculates the charge distribution of each state
c     and writes this information to the out file. If dos or vacdos
c     are .true. it also write the necessary information for dos or
c     bandstructure plots to the file dosinp and vacdos respectivly
c***********************************************************************
c       changed this subroutine slightly for parallisation of dosinp&
c       vacdos output (argument z replaced by ksym,jsym, removed sympsi
c       call)                                        d.wortmann 5.99
c
c******** ABBREVIATIONS ************************************************
c     qal      : l-like charge of each state
c     qvac     : vacuum charge of each state
c     qvlay    : charge in layers (z-ranges) in the vacuum of each state
c     starcoeff: T if star coefficients have been calculated
c     qstars   : star coefficients for layers (z-ranges) in vacuum
c
c***********************************************************************

      USE m_cotra, ONLY : cotra3
      IMPLICIT NONE

C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: neigd,nkptd,nvd,ntypd,jspd,layerd
      REAL wk
      INTEGER ikpt,jspin,layers,nbands,ndir,ntyp,nvac
      LOGICAL dos,film,invs,semic,l_noco,slice,vacdos,zrfs
c
c     STM Arguments
      INTEGER, INTENT (IN) ::nstars
      LOGICAL, INTENT (IN) ::starcoeff
      COMPLEX, INTENT (IN) ::qstars(nstars,neigd,layerd,2)
C     ..
C     .. Array Arguments ..
      REAL,    INTENT (IN) :: qvlay(neigd,layerd,2)
      REAL,    INTENT (IN) :: qis(neigd,nkptd,jspd),bmat(3,3)
      REAL,    INTENT (IN) :: qvac(neigd,2,nkptd,jspd)
      REAL,    INTENT (IN) :: bkpt(3),eig(neigd),qal(0:3,ntypd,neigd)
      INTEGER, INTENT (IN) :: neq(ntypd),jsym(neigd),ksym(neigd)
C     ..
C     .. Local Scalars ..
      REAL qalmax,qishlp,qvacmt,qvact
      INTEGER i,iband,ilay,iqispc,iqvacpc,ityp,itypqmax,ivac,l,lqmax
      INTEGER istar
C     ..
C     .. Local Arrays ..
      REAL cartk(3)
      INTEGER iqalpc(0:3,ntypd)
      CHARACTER chstat(0:3)
C     ..
C     .. Data statements ..
      DATA chstat/'s','p','d','f'/
C     ..


      IF (film) THEN
         WRITE (6,FMT=8000) (bkpt(i),i=1,3)
         WRITE (16,FMT=8000) (bkpt(i),i=1,3)
 8000    FORMAT (/,3x,'q(atom,l): k=',3f10.5,/,/,t8,'e',t13,'max',t18,
     +          'int',t22,'vac',t28,'spheres(s,p,d,f)')
         IF (dos .AND. (.NOT.semic)) THEN
            CALL cotra3(bkpt,cartk,bmat)
            WRITE (85,FMT=8020) cartk(1),cartk(2),cartk(3),nbands,wk
C     *************** for vacdos shz Jan.96
            IF (vacdos) THEN
               WRITE (86,FMT=8020) cartk(1),cartk(2),cartk(3),nbands,wk
            END IF
         END IF
      ELSE
         WRITE (6,FMT=8010) (bkpt(i),i=1,3)
         WRITE (16,FMT=8010) (bkpt(i),i=1,3)
 8010    FORMAT (/,3x,'q(atom,l): k=',3f10.5,/,/,t8,'e',t13,'max',t18,
     +          'int',t24,'spheres(s,p,d,f)')
         IF (dos .AND. (.NOT.semic)) THEN
            CALL cotra3(bkpt,cartk,bmat)
            WRITE (85,FMT=8020) cartk(1),cartk(2),cartk(3),nbands,wk
         END IF
      END IF
 8020 FORMAT (1x,3e20.12,i6,e20.12)

      DO iband = 1,nbands
         IF (slice) THEN
            WRITE (6,FMT=8030) iband,eig(iband)
            WRITE (16,FMT=8030) iband,eig(iband)
 8030       FORMAT (' cdnval: slice for i=',i4,'  and energy',1e12.4)
         END IF

         qvacmt = 0.0
         qvact = 0.0
         IF (film) THEN
            DO ivac = 1,nvac
               qvact = qvact + qvac(iband,ivac,ikpt,jspin)
            END DO
            IF (invs .OR. zrfs) qvact = 2.0*qvact
            iqvacpc = nint(qvact*100.0)
            qvacmt = qvact
         END IF
         qalmax = 0.0
         lqmax = 0
         itypqmax = 0
         DO ityp = 1,ntyp
            DO l = 0,3
               iqalpc(l,ityp) = nint(qal(l,ityp,iband)*100.0)
               qvacmt = qvacmt + qal(l,ityp,iband)*neq(ityp)
               IF (qalmax.LT.qal(l,ityp,iband)) THEN
                  qalmax = qal(l,ityp,iband)
                  lqmax = l
                  itypqmax = ityp
               END IF
            END DO
         END DO
         qishlp = 1.0 - qvacmt
         IF (l_noco) qishlp = qis(iband,ikpt,jspin)
         iqispc = nint(qishlp*100.0)
         IF (film) THEN
            WRITE (6,FMT=8040) eig(iband),chstat(lqmax),itypqmax,
     +        iqispc,iqvacpc, ((iqalpc(l,ityp),l=0,3),ityp=1,ntyp)
            WRITE (16,FMT=8040) eig(iband),chstat(lqmax),itypqmax,
     +        iqispc,iqvacpc, ((iqalpc(l,ityp),l=0,3),ityp=1,ntyp)
 8040       FORMAT (f10.4,2x,a1,i2,2x,2i3, (t26,6 (4i3,1x)))
            IF (dos .AND. (.NOT.semic)) THEN
               IF (ndir.NE.0) THEN
                  WRITE (85,FMT=8050) eig(iband),ksym(iband),
     +              jsym(iband),qvact, ((qal(l,ityp,iband),l=0,3),
     +              ityp=1,ntyp)
 8050             FORMAT (f12.5,2i2,f12.5,/, (4f12.5))
               ELSE
                  WRITE (85,FMT=8060) eig(iband),
     +              ((qal(l,ityp,iband),l=0,3),ityp=1,ntyp),qvact
 8060             FORMAT (10f12.5)
               END IF
            END IF
C     ***************** for vacdos shz Jan.96
            IF (vacdos) THEN
                IF (.NOT.starcoeff) THEN
                   WRITE (86,FMT=8070) eig(iband),
     +               ((qvlay(iband,ilay,ivac),ilay=1,layers),
     +                               ivac=1,nvac)
               ELSE
                   WRITE (86,FMT=8070) eig(iband),
     +                 ((qvlay(iband,ilay,ivac),
     +                 (REAL(qstars(istar,iband,ilay,ivac)),
     +                 istar=1,nstars-1),ilay=1,layers),ivac=1,nvac)
                END IF
 8070          FORMAT (f10.4,2x,20(e16.8,1x))

            END IF
C     **************************************
         ELSE
            WRITE (6,FMT=8080) eig(iband),chstat(lqmax),itypqmax,
     +        iqispc, ((iqalpc(l,ityp),l=0,3),ityp=1,ntyp)
            WRITE (16,FMT=8080) eig(iband),chstat(lqmax),itypqmax,
     +        iqispc, ((iqalpc(l,ityp),l=0,3),ityp=1,ntyp)
 8080       FORMAT (f10.4,2x,a1,i2,2x,i3, (t26,6 (4i3,1x)))
            IF (dos .AND. (.NOT.semic)) THEN
               IF (ndir.NE.0) THEN
                  WRITE (85,FMT=8050) eig(iband),ksym(iband),
     +              jsym(iband),0.0, ((qal(l,ityp,iband),l=0,3),ityp=1,
     +              ntyp)
               ELSE
                  WRITE (85,FMT=8060) eig(iband),
     +              ((qal(l,ityp,iband),l=0,3),ityp=1,ntyp),0.0
               END IF
            END IF
         END IF
      END DO
      END


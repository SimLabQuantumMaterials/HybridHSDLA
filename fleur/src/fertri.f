      MODULE m_fertri
!
!     calculates fermi energy and weights using triangular method
!
      CONTAINS
      SUBROUTINE fertri(
     >                  neigd,nkptd,jspd,
     >                  ne,nkpt,jspins,zc,eig,bk,sfac,
     X                  ef,
     <                  seigv,w)

      USE m_triang
      USE m_maketetra
      USE m_tetraef
      USE m_dosef
      USE m_dosint
      USE m_doswt
!     USE m_bzints

      IMPLICIT NONE
!     ..
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN)    :: neigd,nkptd,jspd,nkpt,jspins
      REAL,    INTENT (IN)    :: zc,sfac
      REAL,    INTENT (OUT)   :: seigv
      REAL,    INTENT (INOUT) :: ef
!     ..
!     .. Array Arguments ..
      INTEGER, INTENT (IN)    :: ne(nkptd,jspd)
      REAL,    INTENT (IN)    :: bk(3,nkptd)
      REAL,    INTENT (OUT)   :: w(neigd,nkptd,jspd)
      REAL,    INTENT (INOUT) :: eig(neigd,nkptd,jspd)
!     ..
!     .. Local Scalars ..
      REAL chmom,ct,de,del,dez,ei,emax,emin,s,s1,workf
      REAL lb,ub,e_set
      LOGICAL film
      INTEGER i,ic,j,jsp,k,neig
      INTEGER ntria,ntetra      ! number of triangles & tetrahedrons
      REAL as                   ! total area covered by triangles
!     ..
!     .. Local Arrays ..
      INTEGER itria(3,2*nkptd)  ! index of k-points that are corner points of a triangle
      REAL    atr(2*nkptd)      ! area of a triangle
      INTEGER itetra(4,6*nkptd) ! ditto for tetrahedrons
      REAL    voltet(6*nkpt)
      INTEGER nemax(2)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC real,max0
!     ..
!     .. Data statements ..
      DATA de/5.0e-3/

      WRITE (16,FMT=8000)
      WRITE (6,FMT=8000)
 8000 FORMAT (/,/,10x,'linear triangular method')
c
      film = .true.
      CALL triang(
     >            bk,nkpt,
     <            itria,ntria,atr,as,film)
c
c--->   clear w and set eig=-9999.9
      e_set = -9999.9
      IF (.NOT.film) e_set = 1.0e10
      DO jsp = 1,jspins
         nemax(jsp) = 0.0
         DO k = 1,nkpt
            nemax(jsp) = max0(nemax(jsp),ne(k,jsp))
            DO i = 1,ne(k,jsp)
               w(i,k,jsp) = 0.
            ENDDO
            DO i = ne(k,jsp)+1,neigd
               w(i,k,jsp) = 0.
               eig(i,k,jsp) = e_set
            ENDDO
         ENDDO
      ENDDO
c
!      sfac = 2.0/real(jspins)
c
c--->   write results of triang
c
      IF (.not.film) THEN
         WRITE (6,*)  'reading tetrahedrons from file kpts'
         WRITE (16,*) 'reading tetrahedrons from file kpts'
         OPEN (41,file='kpts',FORM='formatted',STATUS='old')
         DO i = 1, nkpt+1
           READ (41,*)
         ENDDO
         READ (41,'(i5)',ERR=66,END=66) ntetra
         IF (ntetra.GT.6*nkpt) STOP 'fertri: ntetra > 6 nkpt'
         READ (41,'(4(4i6,4x))') ((itetra(i,j),i=1,4),j=1,ntetra)
         READ (41,'(4f20.13)') (voltet(j),j=1,ntetra)
         voltet(1:ntetra) = voltet(1:ntetra) / ntetra
         GOTO 67
 66      CONTINUE                       ! no tetrahedron-information of file
         CALL make_tetra(
     >                   nkptd,nkpt,bk,ntria,itria,atr,
     <                   ntetra,itetra,voltet)
 
 67      CONTINUE                       ! tetrahedron-information read or created
         CLOSE(41)
         lb = MINVAL(eig(:,:,:)) - 0.01
         ub = ef + 0.2
         CALL tetra_ef(
     >                 neigd,nkptd,jspd,jspins,nkpt,
     >                 lb,ub,eig,zc,sfac,
     >                 ntetra,itetra,voltet,
     <                 ef,w)
!         CALL bzints(
!     >               2,6,neigd,nkpt,ntetra,
!     >               ev,wt,itetra,voltet,lo,ub,ned,
!     >               efermi,
!     <               ferwe,g(1,1),g(1,2),sumwei,sume)
      ELSE

        DO i = 1,ntria
           atr(i) = atr(i)/as
        ENDDO
        WRITE (16,FMT=8010) ntria,as
        WRITE (6,FMT=8010) ntria,as
        DO i = 1,ntria
           WRITE (6,FMT=8020) i, (itria(j,i),j=1,3),atr(i)
           WRITE (16,FMT=8020) i, (itria(j,i),j=1,3),atr(i)
        ENDDO 
 8010   FORMAT (/,10x,'triangular decomposition of brillouin zone:',/,
     +          10x,'number of triangles=',i3,/,10x,
     +          'total area of triangles=',f12.6,/,10x,
     +          'no.,corners and (normalized) area of each triangle:',/)
 8020   FORMAT (10x,i3,3x,3i3,f14.6)
        WRITE (16,FMT=*) 'ef_hist=',ef
        WRITE (6,FMT=*) 'ef_hist=',ef
        ei = ef
cjr     emin = -9999.9
        emin = +9999.9
        emax = -emin
        ic = 1
   90   IF (ic.GT.100) GO TO 230
        ic = ic + 1
c
c     results from triang are included here
c
        CALL dosint(
     >              neigd,nkptd,jspd,2*nkptd,
     >              ei,nemax,jspins,sfac,ntria,itria,atr,eig,
     <              ct)
c
        WRITE (6,FMT=*) 'ct=',ct

        IF (ct.LT.zc) THEN            ! ei < ef
          emin = ei
          ei = ei + de
          IF (emin.GT.emax) GO TO 90
        ELSEIF (ct.GT.zc) THEN        ! ei > ef
          emax = ei
          ei = ei - de
          IF (emin.GT.emax) GO TO 90
        ENDIF
        IF (ct.NE.zc) THEN
          WRITE (16,FMT=*) '2nd dosint'
c--->     refine ef to a value of 5 mry * (2**-20)
          iterate : DO i = 1, 40
            ei = 0.5* (emin+emax)
c
            CALL dosint(
     >               neigd,nkptd,jspd,2*nkptd,
     >               ei,nemax,jspins,sfac,ntria,itria,atr,eig,
     <               ct)
c
            WRITE (6,FMT=*) 'i=',i,', ct=',ct
            IF ( ct == zc ) THEN
              EXIT iterate
            ELSEIF ( ct > zc ) THEN
              emax = ei
            ELSE
              emin = ei
            ENDIF
          ENDDO iterate
        ENDIF
        ef = ei
        del = emax - emin
        dez = zc - ct
        workf = -13.6058*2*ef
        WRITE (16,FMT=8030) ef,workf,del,dez
        WRITE (6,FMT=8030) ef,workf,del,dez
 8030   FORMAT(/,10x,'fermi energy=',f10.5,' har',/,10x,'work function='
     +         ,f10.5,' ev',/,10x,'uncertainity in energy and weights=',
     +         2e16.6)
c
c--->   obtain dos at ef
c
        CALL dosef(
     >             neigd,nkptd,jspd,2*nkptd,
     >             ei,nemax,jspins,sfac,ntria,itria,atr,eig)
c
c--->   obtain weights needed for integration
c
        CALL doswt(
     >             neigd,nkptd,jspd,2*nkptd,
     >             ei,nemax,jspins,ntria,itria,atr,eig,
     <             w)

      ENDIF ! .NOT.film
c
c--->   write weights
c
c      DO 190 jsp = 1,jspins
c         neig = nemax(jsp)
c         DO 180 i = 1,neig
c            DO 170 k = 1,nkpt
c               WRITE (6,FMT=*) 'w(',i,',',k,',',jsp,')=',w(i,k,jsp)
c  170       CONTINUE
c  180    CONTINUE
c  190 CONTINUE
c
c--->   obtain sum of weights and valence eigenvalues
c
      s1 = 0.
      seigv = 0.
      DO 220 jsp = 1,jspins
         s = 0.
         neig = nemax(jsp)
         DO 210 i = 1,neig
            DO 200 k = 1,nkpt
               s = s + w(i,k,jsp)
               seigv = seigv + w(i,k,jsp)*eig(i,k,jsp)
  200       CONTINUE
  210    CONTINUE
         s1 = s1 + s
  220 CONTINUE
      seigv = sfac*seigv
      chmom = s1 - jspins*s
      WRITE (16,FMT=8040) seigv,s1,chmom
      WRITE (6,FMT=8040) seigv,s1,chmom
 8040 FORMAT (/,10x,'sum of valence eigenvalues=',f20.6,5x,
     +       'sum of weights=',f10.6,/,10x,'moment=',f12.6)
      RETURN
c
  230 WRITE (16,FMT=8050) ei,ef,emin,emax,ct,zc
      WRITE (6,FMT=8050) ei,ef,emin,emax,ct,zc
 8050 FORMAT (/,/,10x,'error fertri: initial guess of ef off by 25 mry',
     +       ' ei,ef,emin,emax,ct,zc',/,10x,6e16.7,/,10x,
     +       'check number of bands')
      STOP 'fertri'
      END SUBROUTINE fertri
      END MODULE m_fertri

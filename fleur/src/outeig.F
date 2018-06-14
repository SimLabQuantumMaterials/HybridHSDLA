      MODULE m_outeig
c*********************************************************************
c     outputs the eigenvalues and eigenvectors to unit66.
c     if lpr.gt.0, then also list eigenvectors on output file;
c     otherwise only eigenvalues. form66=.true. gives a formatted
c     eigenvector file.
c             m. weinert
c*********************************************************************
      CONTAINS
      SUBROUTINE outeig(
     >     nrec,nvd,neigd,lmaxd,ntypd,jspd,irank,isize,
     >     ellow,wtkpt,lpr,film,l_noco,l_ss,l_J,form66,eonly,
     >     jspin,nw,el,evac,bkpt,nv,k1,k2,k3,z,eig,ne,nblw,nv2,
     >     nbasfcn,nmat,ello,nlod,nlotot,kveclo,nkpt_l,nkpt,
     <     eig_l)

#ifdef CPP_MPI
      USE m_mpi_col_eig
#endif
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nvd,neigd,lmaxd,ntypd,jspd,irank
      INTEGER, INTENT (IN) :: lpr,ne,nblw,jspin,nw,nrec,isize
      INTEGER, INTENT (IN) :: nbasfcn,nmat,nlod,nlotot,nkpt_l,nkpt
      REAL,    INTENT (IN) :: ellow,wtkpt
      LOGICAL, INTENT (IN) :: eonly,film,form66,l_noco,l_ss,l_J
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      INTEGER, INTENT (IN) :: nv(jspd),nv2(jspd),kveclo(nlotot)
      REAL,    INTENT (IN) :: bkpt(3)
      REAL,    INTENT (IN) :: eig(neigd)
      REAL,    INTENT (IN) :: el(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (IN) :: ello(nlod,ntypd,jspd)
      REAL,    INTENT (INOUT) :: evac(2,jspd)
#ifdef CPP_INVERSION
      REAL,    INTENT (IN) :: z(nbasfcn,neigd)
#else
      COMPLEX, INTENT (INOUT) :: z(nbasfcn,neigd)
#endif
      REAL,   INTENT(OUT):: eig_l(neigd+5,nkpt_l)
C     ..
C     .. Local Scalars ..
      INTEGER i,j,k,m,n4,nbg,nend,l,n,nfcn,nrec1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC min
C     ..
c
c--->    hsint, hssphn and hsvac set up the lower triangle of the
c--->    hamiltonian and overlap matrix row-wise, but franza expects
c--->    column-wise input, therefore it is used as the upper triangle.
c--->    Now, the complex conjugates of the eigenvectors have to be used.
c
#if ( !defined(CPP_INVERSION) )
        nfcn = nv(jspin)
        IF (l_noco) nfcn = 2*nv(jspin)
        IF (l_ss)   nfcn = nv(1) + nv(2)
        DO i=1,ne
          DO k=1,nmat
             z(k,i) = conjg(z(k,i))
          ENDDO
        ENDDO
#endif
c
      IF ( (irank.EQ.0).AND.(.not.l_J) ) THEN
        IF (film) THEN
           IF (l_ss) THEN
              WRITE (6,FMT=8005) (bkpt(i),i=1,3),nv,nv2,ne
              WRITE (16,FMT=8005) (bkpt(i),i=1,3),nv,nv2,ne
           ELSE
              WRITE (6,FMT=8000) (bkpt(i),i=1,3),nv(jspin),nv2(jspin),ne
              WRITE (16,FMT=8000)(bkpt(i),i=1,3),nv(jspin),nv2(jspin),ne
           ENDIF
        ELSE
           IF (l_ss) THEN
              WRITE (6,FMT=8015) (bkpt(i),i=1,3),nv,ne
              WRITE (16,FMT=8015) (bkpt(i),i=1,3),nv,ne
           ELSE
              WRITE (6,FMT=8010) (bkpt(i),i=1,3),nv(jspin),ne
              WRITE (16,FMT=8010) (bkpt(i),i=1,3),nv(jspin),ne
           ENDIF
        END IF
 8000 FORMAT (1x,/,/,' k=(',3f12.6,'):',i10,' 3-d basis functions and',
     +       i6,' 2-d basis functions',/,' the',i4,' eigenvalues are:')
 8005 FORMAT (1x,/,/,' k=(',3f12.6,'):',i10,',',i10,
     +       ' 3-d basis functions and',i6,',',i6,'
     +       2-d basis functions',/,' the',i4,' eigenvalues are:')
 8010 FORMAT (1x,/,/,' k=(',3f12.6,'):',i10,' 3-d basis functions',/,
     +       ' the',i4,' eigenvalues are:')
 8015 FORMAT (1x,/,/,' k=(',3f12.6,'):',i10,',',i10,
     +       ' 3-d basis functions',/,' the',i4,' eigenvalues are:')
        WRITE (6,FMT=8020) (eig(i),i=1,ne)
        WRITE (16,FMT=8020) (eig(i),i=1,ne)
 8020 FORMAT (5x,5f12.6)

        IF (nblw.GT.-9999) THEN
          WRITE (6,FMT=8030) nblw,ellow
          WRITE (16,FMT=8030) nblw,ellow
        ENDIF
 8030 FORMAT (5x,60 ('*'),/,10x,i3,' eigenvalue(s) below the energy',
     +       f12.6,/,5x,60 ('*'),/)
        IF ((lpr.GT.0) .AND. (.NOT.eonly)) THEN
           n4 = ne/4
           IF (4*n4.LT.ne) n4 = n4 + 1
           WRITE (6,FMT=8040)
           nend = 0
           DO 10 m = 1,n4
              nbg = nend + 1
              nend = min(nend+4,ne)
              WRITE (6,FMT=8050) (eig(i),i=nbg,nend)
              DO 1 k=1,nv(jspin)
                 WRITE (6,8060) k1(k,jspin),k2(k,jspin),k3(k,jspin),
     +                          (z(k,j),j=nbg,nend)
    1         CONTINUE
 10        CONTINUE
        END IF
      END IF
 8040 FORMAT (1x,/,' eigenvectors:')
 8050 FORMAT (1x,/,' eigenvalues   ',4 (3x,f12.6))
#ifdef CPP_INVERSION
 8060 FORMAT (3i5,5x,4 (3x,f11.6))
#else
 8060 FORMAT (3i5,5x,4 (3x,'(',f11.6,',',f11.6,')'))
#endif
c--->    write to unit66 (eigenvector file)
c
      IF (form66) THEN
c--->    formatted file
         IF (irank.EQ.0) THEN 
           IF (.NOT.eonly) THEN
             WRITE (66,FMT=8110) ((el(l,n,jspin),l=0,lmaxd),n=1,ntypd)
             IF (film) WRITE (66,FMT=8110) (evac(i,jspin),i=1,2)
           END IF
           WRITE (66,FMT=8070) bkpt,wtkpt,ne,nv
           WRITE (66,FMT=8080) (eig(i),i=1,ne)
           IF (.NOT.eonly) THEN
             WRITE (66,FMT=8090) (k1(k,jspin),k2(k,jspin),k3(k,jspin)
     +                             ,k=1,nv(jspin))
             WRITE (66,FMT=8090) (kveclo(k),k=1,nlotot)
            DO 20 j = 1,ne
               WRITE (66,FMT=8100) (z(k,j),k=1,nv(jspin))
   20       CONTINUE
           END IF
         ENDIF
 8110    FORMAT (1p,4e20.13)
 8070    FORMAT (3f10.6,f12.8,4i5)
 8080    FORMAT (6f15.6)
 8090    FORMAT (12i6)
 8100    FORMAT (1p,4e20.13)
#ifdef  CPP_MPI
         CALL mpi_col_eig(
     >                    irank,isize,jspd,neigd,lmaxd,ntypd,nlod,
     >                    l_ss,l_noco,l_J,jspin,form66,eonly,film,nkpt,
     X                    nrec,el,evac,ello,bkpt,wtkpt,ne,nv,
     X                    nmat,eig)
#endif

      ELSE ! => form66=F

#if ( defined(CPP_MPI) && !defined(CPP_T90) )
         IF (irank.ne.0) THEN
            nrec1 = (nrec - irank - 1) / isize + 1
            IF ((jspin == 2).AND.(irank<mod(nkpt,isize)) ) THEN
              nrec1 = nrec1 + 1
            ENDIF
         ELSE
            nrec1 = nrec
         ENDIF
#else
         nrec1 = nrec
#endif
c--->    standard unformatted file
         IF (l_J) THEN

           nrec1 = (nrec - irank - 1) / isize + 1
           eig_l(1:3,nrec1)    = bkpt(1:3)
           eig_l(4,nrec1)      = wtkpt
           eig_l(5,nrec1)      = REAL(ne)
           eig_l(6:ne+5,nrec1) = eig(1:ne)
           eig_l(ne+6:neigd+5,nrec1) = 0.0
!           IF (irank == 0) THEN
!             nrec1 = nrec
!             WRITE (66,rec=nrec1) bkpt,wtkpt,ne,eig
!           ENDIF

         ELSEIF (.NOT.eonly) THEN

           !write(*,*) 'process',irank,' wrote ',nrec1
           IF (l_ss) THEN
              WRITE (66,rec=nrec1)
     +             el,evac,ello,bkpt,wtkpt,ne,nv,nmat,
     +             eig,k1,k2,k3,kveclo,z
           ELSEIF (l_noco) THEN
              WRITE (66,rec=nrec1)
     +             el,evac,ello,bkpt,wtkpt,ne,nv(jspin),nmat,
     +             eig,k1(:,jspin),k2(:,jspin),
     +                 k3(:,jspin),kveclo, z
           ELSE
              WRITE (66,rec=nrec1) el(:,:,jspin),evac(:,jspin),
     +             ello(:,:,jspin),bkpt,wtkpt,ne,
     +             nv(jspin),nmat,eig,k1(:,jspin),k2(:,jspin),
     +                 k3(:,jspin),kveclo, z
           ENDIF
         ELSE                 ! => eonly=T
           IF (l_ss) THEN
              WRITE (66,rec=nrec)
     +             el,evac,ello,bkpt,wtkpt,ne,nv,nmat,
     +             eig,k1,k2,k3,kveclo
           ELSEIF (l_noco) THEN
              WRITE (66,rec=nrec)
     +             el,evac,ello,bkpt,wtkpt,ne,nv(jspin),nmat,
     +             eig,(k1(k,jspin),k=1,nvd),(k2(k,jspin),k=1,nvd),
     +             (k3(k,jspin),k=1,nvd),kveclo
           ELSE
              WRITE (66,rec=nrec) el(:,:,jspin),evac(:,jspin),
     +             ello(:,:,jspin),bkpt,wtkpt,ne,
     +             nv(jspin),nmat,eig,(k1(k,jspin),k=1,nvd),
     +             (k2(k,jspin),k=1,nvd),(k3(k,jspin),k=1,nvd),kveclo
           ENDIF
         END IF
#if ( defined(CPP_MPI) && !defined(CPP_T90) )
         IF (.not.l_J) THEN
           CALL mpi_col_eig(
     >                    irank,isize,jspd,neigd,lmaxd,ntypd,nlod,
     >                    l_ss,l_noco,l_J,jspin,form66,eonly,film,nkpt,
     X                    nrec,el,evac,ello,bkpt,wtkpt,ne,nv,
     X                    nmat,eig)
         ENDIF
#endif
      END IF

      END SUBROUTINE outeig
      END MODULE m_outeig

      MODULE m_rwnoco
c---------------------------------------------------------------------
c     read or write the nocoinp-file
c---------------------------------------------------------------------
      CONTAINS
      SUBROUTINE rw_noco(
     >                   ch_rw,noinpfile,ntype,l_J,l_soc,
     X                   alpha,beta,l_relax,b_con,
     X                   l_ss,l_mperp,l_constr,mix_b,qss,sso_opt,
     <                   thetaJ,l_magn,nmagn,M,mtypes,
     <                   magtype,nmagtype,nsh,l_disp)

      IMPLICIT NONE
C ..
C ..  Arguments ..
      CHARACTER*1,INTENT (IN) :: ch_rw
      INTEGER, INTENT    (IN) :: noinpfile,ntype
      LOGICAL, INTENT (INOUT) :: l_ss,l_mperp,l_constr
      REAL,    INTENT (INOUT) :: mix_b

      REAL,    INTENT (INOUT) :: alpha(ntype),beta(ntype),qss(3)
      REAL,    INTENT (INOUT) :: b_con(2,ntype)
      LOGICAL, INTENT (INOUT) :: l_relax(ntype)
      LOGICAL, INTENT (INOUT) :: sso_opt(2) 

      REAL,    INTENT (OUT)   :: thetaJ,M(ntype)
      INTEGER, INTENT (OUT)   :: nmagn,nsh,mtypes
      INTEGER, INTENT (OUT)   :: magtype(ntype),nmagtype(ntype)
      LOGICAL, INTENT (OUT)   :: l_magn(ntype),l_disp
      LOGICAL, INTENT (IN)    :: l_J,l_soc 
C ..
      INTEGER :: itype, j, fileend
      REAL    :: pi, qsc(3)
      CHARACTER(len=8) :: inpchar 

c---------------------------------------------------------------------
      IF (ch_rw.eq.'R') THEN
c---------------------------------------------------------------------

      pi= 4.*ATAN(1.)

      IF (l_J) nmagn=0
         magtype(:)=0
         nmagtype(:)=0
         mtypes=1
      DO itype = 1,ntype
         IF (l_J) THEN
         READ (noinpfile,'(22x,l1,8x,l1,3x,f8.5,9x,i4)') l_relax(itype),
     &        l_magn(itype),M(itype),magtype(itype)
         WRITE (6,8010)   itype,l_relax(itype),l_magn(itype),
     &                      M(itype),magtype(itype)

         ELSE
           READ (noinpfile,'(22x,l1)') l_relax(itype)
           l_magn(itype) = .false. ; M(itype) = 0.0
         ENDIF

         inpchar(1:2)= 'XX'
         READ (noinpfile,fmt='(19x,a2)') inpchar(1:2)
         BACKSPACE (noinpfile)
         IF (inpchar(1:2)/='pi') THEN
           READ (noinpfile,fmt='(7x,f14.10,11x,f14.10)') 
     &      alpha(itype),b_con(1,itype)
         ELSE
           READ (noinpfile,fmt='(7x,f12.8,13x,f14.10)') 
     &      alpha(itype),b_con(1,itype)
            alpha(itype)= alpha(itype)*pi
         ENDIF 
         inpchar(1:2)= 'XX'
         READ (noinpfile,fmt='(19x,a2)') inpchar(1:2)
         BACKSPACE (noinpfile)
         IF (inpchar(1:2)/='pi') THEN
           READ (noinpfile,fmt='(7x,f14.10,11x,f14.10)') 
     &      beta(itype),b_con(2,itype)
         ELSE
           READ (noinpfile,fmt='(7x,f12.8,13x,f14.10)') 
     &      beta(itype),b_con(2,itype)
            beta(itype)= beta(itype)*pi
         ENDIF 

         READ (noinpfile,*)
         IF (l_J) THEN
           IF (l_magn(itype)) THEN
             nmagn = nmagn+1
             IF(magtype(itype).eq.0) magtype(itype)=nmagn
             nmagtype(magtype(itype))=nmagtype(magtype(itype))+1
             IF((nmagn.GE.2).and.(nmagtype(magtype(itype)).eq.1))
     &           mtypes=mtypes+1
           ENDIF
         ELSE 
           WRITE (6,8010)        itype,l_relax(itype)
           WRITE (6,8020)        alpha(itype),b_con(1,itype)
           WRITE (6,8025)         beta(itype),b_con(2,itype)
           WRITE (6,*)
        ENDIF
      ENDDO
           WRITE (6,8026) nmagn,mtypes
 8010 FORMAT ('atom-type',i4,',l_relax=',l1,',l_magn=',l1,
     &',M=',f8.5,',magtype=',i4)
 8020 FORMAT ('alpha =',f14.10,',b_cons_x =',f14.10)
 8025 FORMAT ('beta  =',f14.10,',b_cons_y =',f14.10)
 8026 FORMAT ('Total number of magnetic atoms',i4,',magnetic types'i4)

      READ (noinpfile,*)
      IF (l_J) THEN
         READ (noinpfile,8035) l_ss,l_mperp,l_constr,l_disp
      ELSE
         READ (noinpfile,8036) l_ss,l_mperp,l_constr
      ENDIF
      BACKSPACE (noinpfile)
      inpchar= 'XXXXXXXX' 
      READ (noinpfile,fmt='(37x,a8)',ERR=200) inpchar(1:8)
      IF ( (inpchar(1:8)=='sso_opt=') .or. (l_ss .and. l_soc) ) THEN
        BACKSPACE (noinpfile)
        READ (noinpfile,fmt='(45x,2l1)') sso_opt(1),sso_opt(2)
      ELSE
        sso_opt(1)= .false. 
        sso_opt(2)= .false.
      ENDIF 
 200  IF (l_J) THEN
         READ (noinpfile,8045) mix_b,thetaJ,nsh
      ELSE
         READ (noinpfile,8046) mix_b
         l_disp = .false. ; thetaJ = 0.0 ; nsh = 0
      ENDIF

c--- J constants
      IF (l_J) THEN 
        l_ss = .true.
        IF (l_disp) THEN
          WRITE(6,*) 'Calculating magnon spectrum'
        ELSE
         WRITE(6,*)'This is a calculation of the coupling constants Jij'
          WRITE(6,*)'The cone angle used for this calculation is'
          WRITE(6,*)'thetaJ=',thetaJ
          WRITE(6,*)'The interactions are calculated for the following'
          WRITE(6,*)'atom types (using the given magnetic moments):'
          DO itype = 1,ntype
            IF (l_magn(itype)) THEN
               WRITE(6,*) 'atom type',itype,',M=',M(itype),',one of ',
     &nmagtype(magtype(itype)),' atoms of magnetic type ',magtype(itype)
            ENDIF
          ENDDO 
        ENDIF
        RETURN
      ENDIF
      WRITE (6,fmt='(5(A,l1),l1)') 
     & 'l_ss=',l_ss,',l_mperp=',l_mperp,',l_constr=',l_constr,
     & ',l_disp=',l_disp,',sso_opt=',sso_opt(1),sso_opt(2)
      WRITE (6,8040) mix_b,thetaJ
 8030 FORMAT ('l_ss=',l1,',l_mperp=',l1,',l_constr=',l1,',l_disp=',l1)
 8035 FORMAT (5x,l1,9x,l1,10x,l1,8x,l1)
 8036 FORMAT (5x,l1,9x,l1,10x,l1)
 8040 FORMAT ('mix_b=',f6.3,',thetaJ=',f14.10,',nsh=',i4)
 8045 FORMAT (6x,f6.3,8x,f14.10,5x,i4)
 8046 FORMAT (6x,f6.3)

      IF (l_ss) THEN
         READ (noinpfile,fmt='(5x,3(f14.10,1x))') 
     &    qss(1),qss(2),qss(3)
         inpchar(1:3)= 'XXX'
         READ (noinpfile,fmt='(a4)',iostat=fileend) inpchar(1:4)
         IF (fileend==0) THEN ; IF (inpchar(1:4)=='qsc=') THEN
           BACKSPACE (noinpfile)
           READ (noinpfile,fmt='(5x,3(f14.10,1x))') 
     &      qsc(1),qsc(2),qsc(3)
           DO j= 1,3
             IF ( ABS(qsc(j)) < 1.e-6 ) THEN
               WRITE (6,fmt='(A,i1,A,x,f14.10)')
     &          'Error reading nocoinp: qsc(',j,') =',qsc(j)
               STOP 'Stop in rw_noco'
             ENDIF
             qss(j)= qss(j)/qsc(j)
           ENDDO
         ENDIF ; ENDIF 
         WRITE (6,*) 'This is a Spin-Spiral (SS) calculation. The'
         WRITE (6,*) 'q-vector of the Spin-Spiral is:'
         WRITE (6,8060) qss(1),qss(2),qss(3)
 8060    FORMAT ('qss=(',f14.10,',',f14.10,',',f14.10,')')
      ENDIF


c---------------------------------------------------------------------
      ELSEIF (ch_rw.eq.'W') THEN
c---------------------------------------------------------------------
 
      DO itype = 1,ntype
         IF ( l_relax(itype) ) b_con(:,itype) = 0.0
         WRITE (noinpfile,8010) itype,l_relax(itype),l_magn(itype),
     &          M(itype)
         WRITE (noinpfile,8020) alpha(itype),b_con(1,itype)
         WRITE (noinpfile,8025)  beta(itype),b_con(2,itype)
         WRITE (noinpfile,*)
      ENDDO
 
      WRITE (noinpfile,*) '-- logical parameters --'
      WRITE (noinpfile,fmt='(5(A,l1),l1)') 
     & 'l_ss=',l_ss,',l_mperp=',l_mperp,',l_constr=',l_constr,
     & ',l_disp=',l_disp,',sso_opt=',sso_opt(1),sso_opt(2)
      WRITE (noinpfile,8040) mix_b,thetaJ,nsh

      IF (l_ss) THEN
       WRITE (noinpfile,8060) qss(1),qss(2),qss(3)
       WRITE (6,8060) qss(1),qss(2),qss(3)
      ENDIF

c---------------------------------------------------------------------
      ELSE
        STOP "rw_noco: choose 'R' for read or 'W' to write !"
      ENDIF

      END SUBROUTINE rw_noco
      END MODULE m_rwnoco

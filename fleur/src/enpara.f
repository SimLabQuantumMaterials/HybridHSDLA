      MODULE m_enpara
c     *************************************************************
c     Module containing three subroutines
c     r_enpara: read enpara file
c     w_enpara: write enpara file
c     mix_enpara: calculate new energy parameters
c     *************************************************************
      CONTAINS
      SUBROUTINE w_enpara(   
     >                    lmaxd,nlod,ntype,nw,jspin,film,nlo,
     >                    skiplo,ello0,el0,evac0,lchange,
     >                    llochg,lchg_v,enmix,id)
c
c write enpara-file
c
      IMPLICIT NONE

      INTEGER, INTENT (IN) :: lmaxd,nlod,ntype,nw,jspin,id
      LOGICAL, INTENT (IN) :: lchg_v,film

      INTEGER, INTENT (IN) :: skiplo(ntype),nlo(ntype)
      REAL,    INTENT (IN) :: ello0(nlod,ntype),el0(0:lmaxd,ntype)
      REAL,    INTENT (IN) :: evac0(2),enmix
      LOGICAL, INTENT (IN) :: lchange(0:lmaxd,ntype),llochg(nlod,ntype)

      INTEGER n,l,lo


      WRITE (6,FMT=8030) nw
      WRITE (id,FMT=8030) nw

      WRITE (40,FMT=8035) nw,jspin,enmix
 8030 FORMAT (/,5x,'energy parameters for window',i2,/,t6,'atom',t15,
     +       's',t24,'p',t33,'d',t42,'f')
 8035 FORMAT (5x,'energy parameters for window',i2,' spin ',i1,
     +     ' mix=',f10.6,/,t6,'atom',t15,
     +       's',t24,'p',t33,'d',t42,'f')
      DO n = 1,ntype
         WRITE (6,FMT=8040)  n, (el0(l,n),l=0,3),
     +                          (lchange(l,n),l=0,3),skiplo(n)
         WRITE (id,FMT=8040) n, (el0(l,n),l=0,3),
     +                          (lchange(l,n),l=0,3),skiplo(n)
         WRITE (40,FMT=8040) n, (el0(l,n),l=0,3),
     +                          (lchange(l,n),l=0,3),skiplo(n)
c--->    energy parameters for the local orbitals
         IF (nlo(n).GE.1) THEN
            WRITE (6,FMT=8039) (ello0(lo,n),lo=1,nlo(n))
            WRITE (6,FMT=8038) (llochg(lo,n),lo=1,nlo(n))
            WRITE (id,FMT=8039) (ello0(lo,n),lo=1,nlo(n))
            WRITE (id,FMT=8038) (llochg(lo,n),lo=1,nlo(n))
            WRITE (40,FMT=8039) (ello0(lo,n),lo=1,nlo(n))
            WRITE (40,FMT=8038) (llochg(lo,n),lo=1,nlo(n))
         END IF

      ENDDO
 8038 FORMAT (' --> change   ',30(l1,8x))
 8039 FORMAT (' --> lo ',10f9.5)
 8040 FORMAT (' -->',i2,2x,4f9.5,' change: ',4l1,' skiplo: ',i3)

      IF (film) THEN
         WRITE (40,FMT=8050) evac0(1),lchg_v,evac0(2)
         WRITE (6,FMT=8050)  evac0(1),lchg_v,evac0(2)
         WRITE (id,FMT=8050) evac0(1),lchg_v,evac0(2)
 8050    FORMAT ('  vacuum parameter=',f9.5,' change: ',l1,
     +           ' second vacuum=',f9.5)
      ENDIF

      RETURN
      END SUBROUTINE w_enpara
!
!------------------------------------------------------------------
      SUBROUTINE r_enpara(
     >                    lmaxd,nlod,ntype,film,jsp,
     >                    nw,nlo,lmax,neq,
     <                    skiplo,ello0,el0,evac0,lchange,llochg,
     <                    lchg_v,enmix)
!------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT (IN)  :: lmaxd,nlod,ntype,nw,jsp
      LOGICAL, INTENT (IN)  :: film
      REAL,    INTENT (OUT) :: enmix
      LOGICAL, INTENT (OUT) :: lchg_v

      INTEGER, INTENT (IN)  :: nlo(ntype),lmax(ntype),neq(ntype)
      INTEGER, INTENT (OUT) :: skiplo(ntype)
      REAL,    INTENT (OUT) :: ello0(nlod,ntype)
      REAL,    INTENT (OUT) :: el0(0:lmaxd,ntype)
      REAL,    INTENT (OUT) :: evac0(2)
      LOGICAL, INTENT (OUT) :: lchange(0:lmaxd,ntype)
      LOGICAL, INTENT (OUT) :: llochg(nlod,ntype)

      INTEGER n,l,lo,skip_t
      lchange=.false.
 
!-->  first line contains mixing parameter!

      enmix = 0.0
      READ (40,FMT='(48x,f10.6)',END=100) enmix
 100  READ (40,*)                       ! skip next line
      IF (enmix.EQ.0.0) enmix = 1.0
      WRITE (6,FMT=8001) jsp,nw
      WRITE (6,FMT=8000)
      skip_t = 0
      DO n = 1,ntype
         READ (40,FMT=8040,END=200) (el0(l,n),l=0,3),
     +                          (lchange(l,n),l=0,3),skiplo(n)    
         WRITE (6,FMT=8140) n,(el0(l,n),l=0,3),
     +                          (lchange(l,n),l=0,3),skiplo(n)    
!
!--->    energy parameters for the local orbitals
!
         IF (nlo(n).GE.1) THEN
             skip_t = skip_t + skiplo(n) * neq(n)
             IF (nw.GT.1) THEN
                 WRITE (6,*) 'multiple windows and local orbitals',
     +                       'cannot be used at the same time'
                 STOP 'inpeig: nlo > 0 and nwd > 1'
             END IF
             READ (40,FMT=8039,END=200)  (ello0(lo,n),lo=1,nlo(n))
             READ (40,FMT=8038,END=200) (llochg(lo,n),lo=1,nlo(n))
             WRITE (6,FMT=8139)          (ello0(lo,n),lo=1,nlo(n))
             WRITE (6,FMT=8138)         (llochg(lo,n),lo=1,nlo(n))
         ELSEIF (skiplo(n).GT.0) THEN
             WRITE (6,*) "for atom",n," no LO's were specified"
             WRITE (6,*) 'but skiplo was set to',skiplo 
             STOP "inpeig: no LO's but skiplo > 0 !"
         END IF
!
!--->    set the energy parameters with l>3 to the value of l=3
!
         DO  l = 4,lmax(n)
             el0(l,n) = el0(3,n)
         ENDDO
      ENDDO   ! ntype
 
      IF (film) THEN
         lchg_v = .true.
         READ (40,FMT=8050,END=200) evac0(1),lchg_v,evac0(2)
         WRITE (6,FMT=8150)         evac0(1),lchg_v,evac0(2)
      ENDIF
      IF (nlod.GE.1) THEN               
         WRITE (6,FMT=8090) jsp,skip_t
         WRITE (6,FMT=8091) 
      END IF

! input formats

 8038 FORMAT (14x,30(l1,8x))
 8039 FORMAT (8x,30f9.5)
 8040 FORMAT (8x,4f9.5,9x,4l1,9x,i3)
 8050 FORMAT (19x,f9.5,9x,l1,15x,f9.5)

! output formats

 8138 FORMAT (' --> change   ',30(l1,8x))
 8139 FORMAT (' --> lo ',10f9.5)
 8140 FORMAT (' -->',i2,2x,4f9.5,' change: ',4l1,' skiplo: ',i3)
 8150 FORMAT ('  vacuum parameter=',f9.5,' change: ',l1,
     +           ' second vacuum=',f9.5)
 8001 FORMAT ('READING enpara for spin: ',i1,' window:',i2)
 8000 FORMAT (/,' energy parameters:',/,t10,'s',t20,
     +        'p',t30,'d',t37,'higher l - - -')
 8090 FORMAT ('Spin: ',i1,' -- ',i3,'eigenvalues')
 8091 FORMAT ('will be skipped for energyparameter computation')

      RETURN

 200  WRITE (6,*) 'the end of the file enpara has been reached while'
      WRITE (6,*) 'reading the energy-parameters.'
      WRITE (6,*) 'possible reason: energy parameters have not been'
      WRITE (6,*) 'specified for all atom types.'
      WRITE (6,FMT='(a,i4)')
     +     'the actual number of atom-types is: ntype=',ntype
      STOP 'end of file enpara reached while reading'

      END SUBROUTINE r_enpara
!
!------------------------------------------------------------------
      SUBROUTINE mix_enpara(
     >                      ntype,nlod,lmaxd,nmzd,jmtd,
     >                      lepr,film,nvac,enmix,
     >                      nlo,llo,jri,dx,rmt,vr,vz,pvac,svac,lchg_v,
     >                      ener,sqal,lchange,enerlo,sqlo,llochg,l_dulo,
     X                      ello0,el0,evac0)
!------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nmzd,nlod,lmaxd,jmtd
      INTEGER, INTENT(IN) :: lepr,nvac,ntype
      LOGICAL, INTENT(IN) :: film

      INTEGER, INTENT(IN) :: jri(ntype),nlo(ntype)
      INTEGER, INTENT(IN) :: llo(nlod,ntype)
      REAL,    INTENT(IN) :: enmix
      REAL,    INTENT(IN) :: dx(ntype),rmt(ntype),vr(jmtd,ntype)
      REAL,    INTENT(IN) :: ener(0:3,ntype),sqal(0:3,ntype)
      REAL,    INTENT(IN) :: enerlo(nlod,ntype),sqlo(nlod,ntype)
      REAL,    INTENT(IN) :: pvac(2),svac(2),vz(nmzd,2)
      LOGICAL, INTENT(IN) :: lchange(0:lmaxd,ntype),lchg_v
      LOGICAL, INTENT(IN) :: llochg(nlod,ntype),l_dulo(nlod,ntype)
      REAL,    INTENT(INOUT) :: ello0(nlod,ntype),el0(0:lmaxd,ntype)
      REAL,    INTENT(INOUT) :: evac0(2)

      INTEGER ityp,j,l,lo
      REAl    vbar,maxdist
      INTEGER same(nlod)

      INTRINSIC exp,log

      maxdist=0.0
      DO ityp = 1,ntype
c        look for LO's energy parameters equal to the LAPW (and previous LO) ones
         same = 0
         DO lo = 1,nlo(ityp)
           IF(el0(llo(lo,ityp),ityp).eq.ello0(lo,ityp)) same(lo)=-1
           DO l = 1,lo-1
             IF(llo(l,ityp).ne.llo(lo,ityp)) cycle
             IF(ello0(l,ityp).eq.ello0(lo,ityp).and.same(lo).eq.0)
     &         same(lo)=l
           ENDDO
         ENDDO
c
c--->   change energy parameters
c
         IF ( lepr.EQ.1) THEN
            j = jri(ityp) - (log(4.0)/dx(ityp)+1.51)
            vbar = vr(j,ityp)/( rmt(ityp)*exp(dx(ityp)*(j-jri(ityp))) )
         ELSE
            vbar = 0.0
         END IF
 
         DO l = 0,3
            IF ( lchange(l,ityp) ) THEN
               write(6,*) 'Type:',ityp,' l:',l
               write(6,FMT=777) el0(l,ityp),
     +              (ener(l,ityp)/sqal(l,ityp) - vbar),
     +              abs(el0(l,ityp)-(ener(l,ityp)/sqal(l,ityp) - vbar))
               maxdist=max(maxdist,
     +              abs(el0(l,ityp)-(ener(l,ityp)/sqal(l,ityp) - vbar)))
               el0(l,ityp) =(1.0-enmix)*el0(l,ityp) + 
     +              enmix*(ener(l,ityp)/sqal(l,ityp) - vbar)
           ENDIF
         ENDDO
         DO l = 4, lmaxd
            IF ( lchange(3,ityp) ) THEN
              el0(l,ityp) = el0(3,ityp)
            ENDIF
         ENDDO
c
c--->    determine and change local orbital energy parameters
c
         DO lo = 1,nlo(ityp)
            IF (l_dulo(lo,ityp)) THEN
               ello0(lo,ityp) =el0(llo(lo,ityp),ityp)
            ELSEIF (llochg(lo,ityp) ) THEN
               IF(same(lo).eq.-1) THEN
                 ello0(lo,ityp) = el0(llo(lo,ityp),ityp)
                 cycle
               ELSE IF(same(lo).gt.0) THEN
                 ello0(lo,ityp) = ello0(same(lo),ityp)
                 cycle
               ENDIF 
               write(6,*) 'Type:',ityp,' lo:',lo
               write(6,FMT=777) ello0(lo,ityp),
     +           (enerlo(lo,ityp)/sqlo(lo,ityp) - vbar),
     +          abs(ello0(lo,ityp)-(enerlo(lo,ityp)/sqlo(lo,ityp)-vbar))
               maxdist=max(maxdist,
     +         abs(ello0(lo,ityp)-(enerlo(lo,ityp)/sqlo(lo,ityp)-vbar)))
               ello0(lo,ityp) =(1.0-enmix)*ello0(lo,ityp)+
     +              enmix*(enerlo(lo,ityp)/sqlo(lo,ityp) - vbar)
            END IF
         END DO

      END DO

      IF (film) THEN

         IF (lepr.eq.1) THEN
           vbar = vz(1,1)
         ELSE
           vbar = 0.0
         ENDIF

         IF (lchg_v ) THEN
            write(6,*) 'Vacuum:'
            write(6,FMT=777) evac0(1),(pvac(1)/svac(1) - vbar),
     +              abs(evac0(1)-(pvac(1)/svac(1) - vbar))
            maxdist=max(maxdist,abs(evac0(1)-(pvac(1)/svac(1) - vbar)))
            evac0(1) =(1.0-enmix)*evac0(1)+
     +              enmix*(pvac(1)/svac(1) - vbar)
            IF (nvac.EQ.2) THEN
               IF (lepr.eq.1) vbar = vz(1,nvac)
               evac0(2) = (1.0-enmix)*evac0(2)+
     +              enmix*(pvac(2)/svac(2) - vbar)
            ELSE
               evac0(2) = evac0(1)
            ENDIF
         ENDIF

      ENDIF
      WRITE(6,'(a36,f12.6)') 'Max. missmatch of energy parameters:',
     +                                                       maxdist
      RETURN
 777  FORMAT('Old:',f8.5,' new:',f8.5,' diff:',f8.5)

      END SUBROUTINE mix_enpara
      END MODULE m_enpara

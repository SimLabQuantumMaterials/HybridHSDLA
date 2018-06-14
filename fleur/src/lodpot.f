      MODULE m_lodpot
      CONTAINS
      SUBROUTINE lodpot(
     >                  irank,jspd,lmaxd,jmtd,nlhd,ntypd,nwdd,nmzd,
     >                  jspins,ntype,nwd,film,nvac,nmz,lepr,
     >                  jri,dx,rmt,rmsh,lmax,vr,vz,llo,zatom,
     >                  el0,evac0,ello0,nlo,nlod,l_dulo,ellow,elup,
     <                  el,evac,ello,bound_lo,bound_up)
c*********************************************************************
c
c set el and evac from el0 & evac0 (depending on lepr)
c
c*********************************************************************
      USE m_constants, ONLY : c_light
      USE m_radsra
      USE m_differ
      IMPLICIT NONE
c
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: irank
      INTEGER, INTENT (IN) :: jmtd,nlhd,ntypd,nwdd,nmzd,jspd,lmaxd
      INTEGER, INTENT (IN) :: jspins,ntype,nwd,nvac,nmz,lepr,nlod
      LOGICAL, INTENT (IN) :: film
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: jri(ntypd),lmax(ntypd),nlo(ntypd)
      INTEGER, INTENT (IN) :: llo(nlod,ntypd)
      LOGICAL, INTENT (IN) :: l_dulo(nlod,ntypd)
      REAL,    INTENT (IN) :: dx(ntypd),rmt(ntypd)
      REAL,    INTENT (IN) :: vr(jmtd,0:nlhd,ntypd,jspd),vz(nmzd,2,4)
      REAL,    INTENT (IN) :: el0(0:lmaxd,ntypd,jspd,nwdd)
      REAL,    INTENT (IN) :: evac0(2,jspd,nwdd)
      REAL,    INTENT (IN) :: ello0(nlod,ntypd,jspd)
      REAL,    INTENT (IN) :: elup(nwdd),ellow(nwdd)
      REAL,    INTENT (OUT):: bound_lo(nwdd),bound_up(nwdd)
      REAL,    INTENT (OUT):: el(0:lmaxd,ntypd,jspd,nwdd) 
      REAL,    INTENT (OUT):: evac(2,jspd,nwdd),ello(nlod,ntypd,jspd)
      REAL,    INTENT (IN) :: rmsh(jmtd,ntypd),zatom(ntypd)

C     ..
C     .. Local Scalars ..
      INTEGER jsp,n,nw,ivac,j,l,ilo
      INTEGER nodeu,node,ierr,msh
      REAL vbar,vz0,rj,e_up,e_lo
      REAL us,dus,c,e,d,rn,fl,fn,fj,t2,rr,t1,test
      LOGICAL start
C     ..
C     .. Local Arrays ..
      INTEGER nqn(0:3),nqn_lo(nlod)
      REAL, ALLOCATABLE :: f(:,:),vrd(:)
      LOGICAL l_done(ntypd,jspd,nwdd),lo_done(nlod,ntypd,jspd)
      CHARACTER(len=1) :: ch(0:3)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC exp,log,min,max
C     ..
c
c initialise also if film=.false. to use as marker for ev-parallelisation
c
      el(:,:,:,:) = 0.0 ; evac(:,:,:) = 0.0 ; ello(:,:,:) = 0.0

!test_start
      IF ( lepr == 0 ) THEN ! not for floating energy parameters
      c = c_light(1.0)
      ch(0:3) = (/'s','p','d','f'/)
      DO nw = 1,nwd
        DO jsp = 1,jspins
          DO n = 1, ntype
! check what to do ( for 'normal' energy parameters )
            nqn(0:3) = nint( el0(0:3,n,jsp,nwd) )
            test = 0.0
            DO l = 0,3
              test = test + abs ( el0(l,n,jsp,nwd) - nqn(l) )
            ENDDO

            nqn_lo( 1:nlo(n)) = nint( ello0(1:nlo(n),n,jsp) )
            DO ilo = 1, nlo(n)
              test = test + abs( ello0(ilo,n,jsp) - nqn_lo(ilo) )
            ENDDO
            
            IF ( (abs(test) < 0.00001 ) .AND.
     +         ( maxval(nqn) < 9) .AND.
     +         ( minval(nqn) > 0) ) THEN ! determine the energy parameters
            l_done(n,jsp,nw) = .true.

            d = exp(dx(n))
! set up core-mesh
            rn = rmt(n)
            msh = jri(n)
            DO WHILE (rn < rmt(n) + 20.0)
               msh = msh + 1
               rn = rn*d
            ENDDO
            rn = rmsh(1,n)*( d**(msh-1) )
            ALLOCATE ( f(msh,2),vrd(msh) )

! extend core potential (linear with slope t1 / a.u.)

            DO j = 1, jri(n)
              vrd(j) = vr(j,0,n,jsp)
            ENDDO      
            t1=0.125
            t2 = vrd(jri(n))/rmt(n) - rmt(n)*t1
            rr = rmt(n)
            DO j = jri(n) + 1, msh
               rr = d*rr
               vrd(j) = rr*( t2 + rr*t1 )
            ENDDO

! search for branches

            DO l = 0, 3
              node = max(nqn(l) - (l+1),0)
              e = 0.0
! determine upper edge
              nodeu = -1 ; start = .true.
              DO WHILE ( nodeu <= node ) 
                CALL radsra(
     >                      e,l,vr(1,0,n,jsp),rmsh(1,n),
     >                      dx(n),jri(n),jmtd,c,
     <                      us,dus,nodeu,f(1,1),f(1,2))
                IF  ( ( nodeu > node ) .AND. start ) THEN
                  e = e - 1.0
                  nodeu = -1
                ELSE
                  e = e + 0.01
                  start = .false.
                ENDIF
              ENDDO 
              e_up = e
              IF (node /= 0) THEN
! determine lower edge
                nodeu = node + 1
                DO WHILE ( nodeu >= node ) 
                  CALL radsra(
     >                        e,l,vr(1,0,n,jsp),rmsh(1,n),
     >                        dx(n),jri(n),jmtd,c,
     <                        us,dus,nodeu,f(1,1),f(1,2))
                  e = e - 0.01
                ENDDO 
                e_lo = e
              ELSE
                e_lo = -9.99 
              ENDIF

! determine notches
              test = 99.0
              e = (e_up+e_lo)/2
              DO WHILE  (test < -l-1 ) 
                e = e - 0.05
                 IF (e < e_lo) THEN
                  e = (e_up+e_lo)/2
                  exit
                ENDIF
                CALL radsra(
     >                      e,l,vr(1,0,n,jsp),rmsh(1,n),
     >                      dx(n),jri(n),jmtd,c,
     <                      us,dus,nodeu,f(1,1),f(1,2))
                test = dus/us
              ENDDO
              DO WHILE (test > -l-1 )
                e = e + 0.005
                IF (e > e_up) THEN
                  e = (e_up+e_lo)/2
                  exit
                ENDIF
                CALL radsra(
     >                      e,l,vr(1,0,n,jsp),rmsh(1,n),
     >                      dx(n),jri(n),jmtd,c,
     <                      us,dus,nodeu,f(1,1),f(1,2))
                test = dus/us
              ENDDO
              IF (irank == 0) THEN
              WRITE (6,'(a5,i3,i2,a1,a12,f7.2,a4,f7.2,a5)') "Atom ",n,
     +              nqn(l),ch(l)," branch, D = ",test," at ",e," htr."
              ENDIF
! calculate core

!              e = (e_up+e_lo)/2
              fn = real(nqn(l)) ; fl = real(l) ; fj = fl + 0.5
              CALL differ(
     >                  fn,fl,fj,c,zatom(n),dx(n),rmsh(1,n),
     >                  rn,d,msh,vrd,
     X                  e,
     <                  f(1,1),f(1,2),ierr)
              el(l,n,jsp,nwd) = e
              IF (irank  == 0) THEN
                 WRITE(6,'(a31,f6.2,a3,f6.2,a13,f8.4)') 
     $               '            branch reaches from',e_lo,' to',e_up
     $               ,' htr. ; e_l =',e
              ENDIF
            ENDDO
            el(4:lmax(n),n,jsp,nwd) = el(3,n,jsp,nwd)

! Now for the lo's
!
            DO ilo = 1, nlo(n)
              l = llo(ilo,n)
! search for branches
              node = nqn_lo(ilo) - (l+1)
              e = 0.0
! determine upper edge
              nodeu = -1 ; start = .true.
              DO WHILE ( nodeu <= node )
                CALL radsra(
     >                      e,l,vr(1,0,n,jsp),rmsh(1,n),
     >                      dx(n),jri(n),jmtd,c,
     <                      us,dus,nodeu,f(1,1),f(1,2))
                IF  ( ( nodeu > node ) .AND. start ) THEN
                  e = e - 1.0
                  nodeu = -1
                ELSE
                  e = e + 0.01
                  start = .false.
                ENDIF
              ENDDO
              e_up = e
              e = 0.0
              IF (node /= 0) THEN
! determine lower edge
                nodeu = node + 1
                DO WHILE ( nodeu >= node )
                  CALL radsra(
     >                        e,l,vr(1,0,n,jsp),rmsh(1,n),
     >                        dx(n),jri(n),jmtd,c,
     <                        us,dus,nodeu,f(1,1),f(1,2))
                  e = e - 0.01
                ENDDO
                e_lo = e
              ELSE
                e_lo = -9.99
              ENDIF

! calculate core

              e = (e_up+e_lo)/2
              fn = real(nqn_lo(ilo)) ; fl = real(l) ; fj = fl + 0.5
              CALL differ(
     >                  fn,fl,fj,c,zatom(n),dx(n),rmsh(1,n),
     >                  rn,d,msh,vrd,
     X                  e,
     <                  f(1,1),f(1,2),ierr)
              ello(ilo,n,jsp) = e
              IF (irank == 0) THEN
                 WRITE(6,'(a4,i3,i2,a1,a12,f6.2,a3,f6.2,a13,f8.4)')
     $                'Atom',n,nqn_lo(ilo),ch(l),' branch from',e_lo
     $                ,' to',e_up,' htr. ; e_l =',e
              ENDIF
            ENDDO

            DEALLOCATE ( f,vrd )
            ELSE ! set eparas below in 'nomal' way
              l_done(n,jsp,nw) = .false.
            ENDIF 

          ENDDO ! n
        ENDDO   ! jsp
      ENDDO     ! nwd
      ELSE
       l_done(:,:,:) = .false.
      ENDIF ! lepr == 0
!test_end
c
      IF ((lepr.eq.1).AND.(irank.EQ.0)) THEN
         WRITE ( 6,'(//,'' Reference energies for energy parameters'')')
         WRITE ( 6,'('' ----------------------------------------'')')
         WRITE (16,'(//,'' Reference energies for energy parameters'')')
         WRITE (16,'('' ----------------------------------------'')')
      ENDIF
c
      spins: DO jsp = 1,jspins
         types: DO n = 1,ntype
         
c
c--->    determine energy parameters if lepr=1. the reference energy
c--->    is the value of the l=0 potential at approximately rmt/4.
c
            IF (lepr.EQ.1) THEN
               j = jri(n) - (log(4.0)/dx(n)+1.51)
               rj = rmt(n)*exp(dx(n)* (j-jri(n)))
               vbar = vr(j,0,n,jsp)/rj
               IF (irank.EQ.0) THEN
               WRITE ( 6,'('' spin'',i2,'', atom type'',i3,'' ='',
     &               f12.6,''   r='',f8.5)') jsp,n,vbar,rj
               WRITE (16,'('' spin'',i2,'', atom type'',i3,'' ='',
     &               f12.6,''   r='',f8.5)') jsp,n,vbar,rj
               ENDIF
            ELSE
               vbar = 0.0
            END IF
            DO nw = 1,nwd
            IF ( .not.l_done(n,jsp,1) ) THEN
               DO l = 0,lmax(n)
                  el(l,n,jsp,nw) = vbar + el0(l,n,jsp,nw)
               ENDDO
            ENDIF
            ENDDO
            IF ( .not.l_done(n,jsp,nwd) ) THEN
            IF (nlo(n).GE.1) THEN
              DO ilo = 1,nlo(n)
                 ello(ilo,n,jsp) = vbar + ello0(ilo,n,jsp)
c+apw+lo
                 IF (l_dulo(ilo,n)) THEN
                     ello(ilo,n,jsp) = el(llo(ilo,n),n,jsp,1)
                 ENDIF
c-apw+lo
              END DO
            ENDIF ! not done
            ENDIF


         ENDDO types

         IF (film) THEN
c
c--->    vacuum energy parameters: for lepr=1, relative to potential
c--->    at vacuum-interstitial interface (better for electric field)
c
            DO ivac = 1,nvac
               vz0 = 0.0
               IF (lepr.eq.1) THEN
                  vz0 = vz(1,ivac,jsp)
                  IF (irank.EQ.0) THEN
                    WRITE ( 6,'('' spin'',i2,'', vacuum   '',i3,'' ='',
     &                  f12.6)') jsp,ivac,vz0
                    WRITE (16,'('' spin'',i2,'', vacuum   '',i3,'' ='',
     &                  f12.6)') jsp,ivac,vz0
                  ENDIF
               ENDIF
               DO nw = 1,nwd
                  evac(ivac,jsp,nw) = evac0(ivac,jsp,nw) + vz0
               ENDDO
            ENDDO
            IF (nvac.EQ.1) THEN
               DO nw = 1,nwd
                  evac(2,jsp,nw) = evac(1,jsp,nw)
               ENDDO
            END IF
         END IF
      ENDDO spins
c
c--->    determine lower and upper bounds for energy window
c
      IF (lepr.eq.0) THEN
         DO nw = 1,nwd
            bound_up(nw) = elup(nw)
            bound_lo(nw) = ellow(nw)
         ENDDO
      ELSE
         DO nw =1,nwd
            e_up = -99999.9
            e_lo =  99999.9
            DO jsp=1,jspins
               DO n=1,ntype
                  DO l = 0,lmax(n)
                     e_up = max( e_up , el(l,n,jsp,nw) )
                     e_lo = min( e_lo , el(l,n,jsp,nw) )
                  ENDDO
                  DO ilo=1,nlo(n)
                     e_up = max( e_up , ello(ilo,n,jsp ) )
                     e_lo = min( e_lo , ello(ilo,n,jsp ) )
                  ENDDO
               ENDDO
            ENDDO
            bound_up(nw) = e_up + elup(nw)
            bound_lo(nw) = e_lo + ellow(nw)
         ENDDO
      ENDIF

      END SUBROUTINE lodpot
      END MODULE m_lodpot

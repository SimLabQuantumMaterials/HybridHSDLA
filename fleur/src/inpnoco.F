      MODULE m_inpnoco
c**********************************************************************
c     This subroutine reads the  Euler angles of the  magnetic field 
c     directions and other noncollinear parameters from the file nocoinp.
c
c                                                 Philipp Kurz 98/01/28
c******** ABBREVIATIONS ***********************************************
c     alpha,beta:Euler angles of the local magnetic field direction of
c                each atom(-type). 
c**********************************************************************
      CONTAINS
      SUBROUTINE inpnoco(
     > natd,noinpfile,itmax,isec1,secvar,gauss,tria,
     > ctail,score,ntype,neq,taual,nstm,starcoeff,l_f,l_J,l_soc,
     < l_ss,l_mperp,l_constr,mix_b,alpha,beta,b_con,
     < l_relax,qss,thetaJ,l_magn,nmagn,M,mtypes,magtype,
     < nmagtype,nsh,l_disp,sso_opt)

      USE m_constants, ONLY : pimach
      USE m_rwnoco
      IMPLICIT NONE

C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: natd,noinpfile,itmax,isec1,ntype,nstm
      LOGICAL, INTENT (IN) :: secvar,gauss,tria,ctail,score
      LOGICAL, INTENT (IN) :: starcoeff,l_f
      LOGICAL, INTENT (OUT):: l_ss,l_mperp,l_constr
      REAL,    INTENT (OUT):: mix_b
c--- J constants
      INTEGER, INTENT (OUT) :: nmagn,nsh,mtypes
      REAL,    INTENT (OUT) :: thetaJ
      LOGICAL, INTENT (OUT) :: l_disp
      LOGICAL, INTENT (IN)  :: l_J,l_soc
c--- J constants

C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntype)
      REAL,    INTENT (IN) :: taual(3,natd)
      REAL,    INTENT (OUT):: alpha(ntype),beta(ntype),qss(3)
      REAL,    INTENT (OUT):: b_con(2,ntype)
      LOGICAL, INTENT (OUT):: l_relax(ntype),sso_opt(2) 
c--- J constants
      LOGICAL, INTENT (OUT) :: l_magn(ntype)
      REAL,    INTENT (OUT) :: M(ntype)
      INTEGER, INTENT (OUT) :: magtype(ntype),nmagtype(ntype)
c--- J constants

C     ..
C     .. Local Scalars ..
      INTEGER itype,iatom
      REAL tpi,phi
      LOGICAL l_relax_any

c---> define constant 2*Pi
      tpi = 2 * pimach()

c---> make sure that the complex program has been compiled
#ifdef CPP_INVERSION
         WRITE (6,*) 'Non-collinear calculations can only be done with'
         WRITE (6,*) 'the complex program. Please compile without the'
         WRITE (6,*) 'option "CPP_INVERSION".'
         STOP 'for l_noco=T: recompile without CPP_INVERSION!'
#endif

c---> make sure Wu-diagonalization is switched off
      IF (isec1 .LE. itmax) THEN
         WRITE (6,*) 'This non-collinear version of the flapw program'
         WRITE (6,*) 'cannot be used with the Wu-diagonalization!!'
         WRITE (6,*) 'itmax = ',itmax,'isec1 = ',isec1
         STOP 'Wu-diagonalization cannot be used!!!'
      ENDIF

c---> make sure second variation is switched off
      IF (secvar) THEN
         WRITE (6,*) 'This non-collinear version of the flapw program'
         WRITE (6,*) 'cannot be used with the second variation!!'
         STOP 'Second variation cannot be used!!!'
      ENDIF

c---> make sure histogram method is used
      IF (gauss) THEN
         WRITE (6,*) 'This non-collinear version of the flapw program'
         WRITE (6,*) 'cannot be used with the Gaussian smearing for '
         WRITE (6,*) 'the Brillouin zone integration!!'
         WRITE (6,*) 'Please use the histogram method.'
         STOP 'Only histogram Brillouin zone integration can be used!!!'
      ENDIF

c---> make sure force is switched off
      IF (l_f) THEN
         WRITE (6,*) 'This non-collinear version of the flapw program'
         WRITE (6,*) 'does not support force calculations.'
         STOP 'force calculations not supported!!!'
      ENDIF

c---> make sure nstm equals zero
      IF (nstm.NE.0) THEN
         WRITE (6,*) 'This non-collinear version of the flapw program'
         WRITE (6,*) 'does not support STM calculations(nstm .NE. 0).'
         STOP 'nstm .NE. 0 not supported!!!'
      ENDIF

c---> make sure starcoeff is switched off
c      IF (starcoeff) THEN
c         WRITE (6,*) 'This non-collinear version of the flapw program'
c         WRITE (6,*) 'does not support starcoefficients output.'
c         STOP 'starcoefficients output (for STM) cannot be generated!!!'
c      ENDIF

c---> make sure coretails are switched off
      IF (ctail) THEN
         WRITE (6,*) 'This non-collinear version of the flapw program'
         WRITE (6,*) 'cannot be used with the coretail option!! '
         STOP 'Coretail option cannot be used!!!'
      ENDIF

c---> make sure score is false
      IF (score) THEN
         WRITE (6,*) 'This non-collinear version of the flapw program'
         WRITE (6,*) 'cannot be used with the score option!! '
         STOP 'score must be false!!!'
      ENDIF

      OPEN (noinpfile,file='nocoinp',form='formatted',status='old')

      WRITE (6,*)'This is a non-collinear calculation. The magnetic'
      WRITE (6,*)'moments of the atoms have a fixed direction.'
      WRITE (6,*)'The Euler-angles alpha and beta of this direction'
      WRITE (6,*)'are equal to the polar angles of the magnetic'
      WRITE (6,*)'moment vector phi and theta respectively.'
      WRITE (6,*)

      CALL rw_noco(
     >             'R',noinpfile,ntype,l_J,l_soc,
     <              alpha,beta,l_relax,b_con,
     <              l_ss,l_mperp,l_constr,mix_b,qss,sso_opt,
     <              thetaJ,l_magn,nmagn,M,mtypes,
     <              magtype,nmagtype,nsh,l_disp)

c---> make sure that moments are not relaxed and constrained
      l_relax_any = .false.
      DO itype = 1,ntype
         l_relax_any = l_relax_any.OR.l_relax(itype)
      ENDDO
      IF (l_relax_any.AND.l_constr) THEN
         WRITE (6,*)'The relaxation of the moment is switched on for at'
         WRITE (6,*)'least one atom. At the same time the constrained'
         WRITE (6,*)'moment option has been switched on!!!'
!         STOP 'relaxation of moments and constraint are switched on!!'
      ENDIF
c---> make sure that perp. component of mag. is calculated if needed
      IF ( (l_relax_any .or. l_constr) .and. (.not. l_mperp) ) THEN
        WRITE (6,*)'The relaxation of the moment is switched on for at'
        WRITE (6,*)'least one atom or the constrained moment option is'
        WRITE (6,*)'switched on. In either case, you need to set'
        WRITE (6,*)'l_mperp=T !!'
        STOP 'Stop: Set l_mperp=T to relax or constrain the moments!!'
      ENDIF
c---> make sure l_constr is switched off in the case of spin spirals
      IF (l_constr .and. l_ss) THEN
        WRITE (6,*)'The constraint moment option is not implemeted'
        WRITE (6,*)'for spin spirals.'
        STOP 'Stop: constraint not implemented for spin spirals!!' 
      ENDIF
         
      IF (l_J) RETURN
      IF (l_ss) THEN
c
c--->    the angle beta is relative to the spiral in a spin-spiral
c--->    calculation, i.e. if beta = 0 for all atoms in the unit cell
c--->    that means that the moments are "in line" with the spin-spiral
c--->    (beta = qss * taual). note: this means that only atoms within
c--->    a plane perpendicular to qss can be equivalent!
         iatom = 1
         DO itype = 1,ntype
            phi = tpi*(  qss(1)*taual(1,iatom)
     +                 + qss(2)*taual(2,iatom)
     +                 + qss(3)*taual(3,iatom) )
            alpha(itype) = alpha(itype) + phi
            iatom = iatom + neq(itype)
         ENDDO
      ENDIF

      WRITE (6,*)

      CLOSE (noinpfile)

      END SUBROUTINE inpnoco
      END MODULE m_inpnoco

      MODULE m_forcea21
      CONTAINS
      SUBROUTINE force_a21(
     >                     lmaxd,ntypd,neigd,nobd,natd,nop,
     >                     nlod,llod,loplod,lmd,
     >                     jspd,mrot,ngopr,lmax,invarop,invarind,
     >                     ntype,neq,invtab,multab,amat,bmat,
     >                     we,jsp,nw,epar,ne,eig,ddn,l_geo,lda_u,
     >                     nlo,llo,nlol,lo1l,dulon,uulon,uloulopn,
     >                     acof,bcof,ccof,aveccof,bveccof,cveccof,
     >                     us,uds,dus,duds,ulos,dulos,rmt,odi,ods,
     X                     force,f_a21,f_b4)

CHANGE FOR 2 WINDOW FORCES
c ************************************************************
c Pulay 2nd and 3rd (A17+A20) term force contribution a la Rici
c combined
c NOTE: we do NOT include anymore  the i**l factors
c in the alm,blm coming from to_pulay. Therefore, we can
c use matrixelements from file 28,38 DIRECTLY
c note: present version only yields forces for
c highest energy window (=valence states)
c if also semicore forces are wanted the tmas and tmat files
c have to be saved, indexed and properly used here in force_a21
C 22/june/97: probably we found symmetrization error replacing
C now S^-1 by S (IS instead of isinv)
c ************************************************************
C
C Force contribution B4 added
C Madsen, Blaha, Schwarz, Sjostedt, Nordstrom
C GMadsen FZJ 20/3-01
c
      USE m_cotra, ONLY : cotra0,cotra1
      USE m_od_types, ONLY : od_inp, od_sym
      USE m_forcea21lo
      USE m_forcea21U
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: lmaxd,ntypd,neigd,nobd,natd,nop
      INTEGER, INTENT (IN) :: nlod,llod,loplod,lmd
      INTEGER, INTENT (IN) :: jspd,ne,ntype,jsp,nw
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntypd),lda_u(ntypd)
      INTEGER, INTENT (IN) :: mrot(3,3,nop),ngopr(natd),lmax(ntypd)
      INTEGER, INTENT (IN) :: invarop(natd,nop),invarind(natd)
      INTEGER, INTENT (IN) :: multab(nop,nop),invtab(nop)
      INTEGER, INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd)
      INTEGER, INTENT (IN) :: nlol(0:llod,ntypd),lo1l(0:llod,ntypd)
      REAL,    INTENT (IN) :: dulon(nlod,ntypd),uulon(nlod,ntypd)
      REAL,    INTENT (IN) :: uloulopn(nlod,nlod,ntypd)
      REAL,    INTENT (IN) :: us(0:lmaxd,ntypd)
      REAL,    INTENT (IN) :: uds(0:lmaxd,ntypd)
      REAL,    INTENT (IN) :: dus(0:lmaxd,ntypd)
      REAL,    INTENT (IN) :: duds(0:lmaxd,ntypd)
      REAL,    INTENT (IN) :: ulos(nlod,ntypd),dulos(nlod,ntypd)
      REAL,    INTENT (IN) :: rmt(ntypd)
      REAL,    INTENT (IN) :: we(nobd),epar(0:lmaxd,ntypd)
      REAL,    INTENT (IN) :: eig(neigd),amat(3,3),bmat(3,3)
      REAL,    INTENT (IN) :: ddn(0:lmaxd,ntypd)
      REAL,    INTENT (INOUT) :: force(3,ntypd,jspd)
      COMPLEX, INTENT (INOUT) :: f_a21(3,ntypd),f_b4(3,ntypd)
      COMPLEX, INTENT (IN) ::      acof(nobd,0:lmaxd*(lmaxd+2),natd)
      COMPLEX, INTENT (IN) ::      bcof(nobd,0:lmaxd*(lmaxd+2),natd)
      COMPLEX, INTENT (IN) ::      ccof(-llod:llod,nobd,nlod,natd)
      COMPLEX, INTENT (IN) :: aveccof(3,nobd,0:lmaxd*(lmaxd+2),natd)
      COMPLEX, INTENT (IN) :: bveccof(3,nobd,0:lmaxd*(lmaxd+2),natd)
      COMPLEX, INTENT (IN) :: cveccof(3,-llod:llod,nobd,nlod,natd)
      LOGICAL, INTENT (IN) :: l_geo(ntypd)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
C     ..
C     .. Local Scalars ..
      INTEGER, PARAMETER :: lmaxb=3
      COMPLEX czero,dtd,dtu,utd,utu
      INTEGER lmplmd,tfile,lo,m
      INTEGER i,ie,im,in,l1,l2,ll1,ll2,lm1,lm2,m1,m2,n,natom,
     +        natrun,nb0,nrec,is,isinv,j,irinv,it
      REAL    zero
C     ..
C     .. Local Arrays ..
      COMPLEX, ALLOCATABLE :: tdd(:),tuu(:),tdu(:),tdulo(:,:,:)
      COMPLEX, ALLOCATABLE :: tuulo(:,:,:),tuloulo(:,:,:),tud(:)
      COMPLEX, ALLOCATABLE :: v_mmp(:,:)
      INTEGER, ALLOCATABLE :: ind(:,:)
      REAL,    ALLOCATABLE :: a21(:,:),b4(:,:)
      COMPLEX forc_a21(3),forc_b4(3)
      REAL starsum(3),starsum2(3),gvint(3),gvint2(3)
      REAL vec(3),vec2(3),vecsum(3),vecsum2(3)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC conjg,aimag
C     ..
C     .. Data statements ..
      DATA czero/ (0.000,0.000)/,zero/0.0/
C     ..
      lmplmd = (lmd* (lmd+3))/2
      ALLOCATE ( tdd(0:lmplmd),tuu(0:lmplmd),
     +           tdu(0:lmplmd),tud(0:lmplmd),
     +           tuulo(0:lmd,-llod:llod,nlod),
     +           tdulo(0:lmd,-llod:llod,nlod),
     +           tuloulo(-llod:llod,-llod:llod,loplod),
     +           v_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb),
     +           a21(3,natd),b4(3,natd),ind(0:lmd,0:lmd) )
c
CHANGE FOR 2 WINDOW FORCES
c     nb0=0
      nb0 = (nw-1)*ntype
c
      natom = 1
      DO 10 n = 1,ntype
         IF (l_geo(n)) THEN
         DO i = 1,3
            forc_a21(i) = czero
            forc_b4(i) = czero
         END DO
         nrec = nb0 + n

         IF (jsp.EQ.1) tfile = 28
         IF (jsp.EQ.2) tfile = 38
         IF (nlo(n).GE.1) THEN
           IF (lda_u(n).GE.0) THEN
           READ (tfile,rec=nrec) tuu(:),tdd(:),tud(:),tdu(:),ind(:,:),
     +                        tuulo(:,:,1:nlo(n)),tdulo(:,:,1:nlo(n)),
     +                   tuloulo(:,:,1:nlo(n)*(nlo(n)+1)/2),v_mmp(:,:)
           ELSE
           READ (tfile,rec=nrec) tuu(:),tdd(:),tud(:),tdu(:),ind(:,:),
     +                        tuulo(:,:,1:nlo(n)),tdulo(:,:,1:nlo(n)),
     +                              tuloulo(:,:,1:nlo(n)*(nlo(n)+1)/2)
           ENDIF
         ELSE
           IF (lda_u(n).GE.0) THEN
            READ (tfile,rec=nrec) tuu(:),tdd(:),tud(:),tdu(:),
     +                            ind(:,:),v_mmp(:,:)
           ELSE
            READ (tfile,rec=nrec) tuu(:),tdd(:),tud(:),tdu(:),ind(:,:)
           ENDIF
         ENDIF
         DO natrun = natom,natom + neq(n) - 1
           DO i = 1,3
              a21(i,natrun) = zero
              b4(i,natrun) = zero
           END DO
         END DO
c
         DO ie = 1,ne
c
c
            DO l1 = 0,lmax(n)
               ll1 = l1* (l1+1)
               DO m1 = -l1,l1
                  lm1 = ll1 + m1
                  DO l2 = 0,lmax(n)
c
                     ll2 = l2* (l2+1)
                     DO m2 = -l2,l2
                        lm2 = ll2 + m2
                        DO natrun = natom,natom + neq(n) - 1
                           in = ind(lm1,lm2)
                           IF (in.NE.-9999) THEN
                              IF (in.GE.0) THEN
c
C ATTENTION: the matrix elements tuu,tdu,tud,tdd
C as calculated in tlmplm are the COMPLEX CONJUGATE
C of the non-spherical matrix elements because in the
C matrix building routine hssphn (or similar routines)
C the COMPLEX CONJUGATE of alm,blm is calculated (to
C save complex operations presumably)
C Her, A20 is formulated in the usual way therefore
C we have to take the COMPLEX CONJUGATE versions
C of tuu,tdu,tud,tdd as compared to hssphn!
c
                                 utu = tuu(in)
                                 dtu = tdu(in)
                                 utd = tud(in)
                                 dtd = tdd(in)
                              ELSE
                                 im = -in
                                 utu = conjg(tuu(im))
                                 dtd = conjg(tdd(im))
                                 utd = conjg(tdu(im))
                                 dtu = conjg(tud(im))
                              END IF
                              DO i = 1,3
                                 a21(i,natrun) = a21(i,natrun) + 2.0*
     +                               aimag( conjg(acof(ie,lm1,natrun))
     +                                 *utu*aveccof(i,ie,lm2,natrun)
     +                                    +conjg(acof(ie,lm1,natrun))
     +                                 *utd*bveccof(i,ie,lm2,natrun)
     +                                    +conjg(bcof(ie,lm1,natrun))
     +                                 *dtu*aveccof(i,ie,lm2,natrun)
     +                                    +conjg(bcof(ie,lm1,natrun))
     +                                 *dtd*bveccof(i,ie,lm2,natrun)
     +                                           )*we(ie)/neq(n)
c   END i loop
                              END DO
                           END IF
c   END natrun
                        END DO
c
c   END m2 loop
                     END DO
c   END l2 loop
                  END DO
c+gu 20.11.97
                  utu = epar(l1,n)-eig(ie)
                  utd = 0.5
                  dtu = 0.5
                  dtd = utu*ddn(l1,n)
                  DO i = 1,3
                     DO natrun = natom,natom + neq(n) - 1
                        a21(i,natrun) = a21(i,natrun) + 2.0*
     +                                 aimag(conjg(acof(ie,lm1,natrun))
     +                                  *utu*aveccof(i,ie,lm1,natrun)
     +                                     +conjg(acof(ie,lm1,natrun))
     +                                  *utd*bveccof(i,ie,lm1,natrun)
     +                                     +conjg(bcof(ie,lm1,natrun))
     +                                  *dtu*aveccof(i,ie,lm1,natrun)
     +                                     +conjg(bcof(ie,lm1,natrun))
     +                                  *dtd*bveccof(i,ie,lm1,natrun)
     +                                           )*we(ie) /neq(n)
                     END DO
c
c-gu
c END  i loop
                  END DO
c   END m1 loop
               END DO
c   END l1 loop
            END DO
c   END ie loop
         END DO
c
c--->    add the local orbital and U contribution to a21
c
         CALL force_a21_lo(
     >                     nobd,natd,ntypd,neigd,nlod,llod,lmd,
     >                     loplod,n,natom,neq,lmax,we,eig,ne,
     >                     nlo,llo,nlol,lo1l,ind,
     >                     acof,bcof,ccof,aveccof,bveccof,cveccof,
     >                     tuulo,tdulo,tuloulo,uulon,dulon,uloulopn,
     X                     a21)

         CALL force_a21_U(
     >                     nobd,natd,nlod,llod,lmd,lmaxb,
     >                     natom,neq(n),lmax(n),we,ne,
     >                     nlo(n),llo(1,n),lda_u(n),ddn(0,n),v_mmp,
     >                     acof,bcof,ccof,aveccof,bveccof,cveccof,
     >                     uulon(1,n),dulon(1,n),
     X                     a21)

#ifdef CPP_APW
c -> B4 force
         DO ie = 1,ne
            DO l1 = 0,lmax(n)
               ll1 = l1* (l1+1)
               DO m1 = -l1,l1
                  lm1 = ll1 + m1
                  DO i = 1,3
                     DO natrun = natom,natom + neq(n) - 1
                        b4(i,natrun) = b4(i,natrun) + 0.5 *
     $                            we(ie)/neq(n)*rmt(n)**2*aimag(
     $                    conjg(acof(ie,lm1,natrun)*us(l1,n)
     $                         +bcof(ie,lm1,natrun)*uds(l1,n))*
     $                    (aveccof(i,ie,lm1,natrun)*dus(l1,n)
     $                    +bveccof(i,ie,lm1,natrun)*duds(l1,n) )
     $              -conjg(aveccof(i,ie,lm1,natrun)*us(l1,n)
     $                    +bveccof(i,ie,lm1,natrun)*uds(l1,n) )*
     $                         (acof(ie,lm1,natrun)*dus(l1,n)
     $                         +bcof(ie,lm1,natrun)*duds(l1,n)) )
                     END DO
                  END DO
               END DO
            END DO
            DO lo = 1,nlo(n)
               l1 = llo(lo,n)
               DO m = -l1,l1
                  lm1 = l1* (l1+1) + m
                  do i=1,3
                    DO natrun = natom,natom + neq(n) - 1
                       b4(i,natrun) = b4(i,natrun) + 0.5 *
     $                        we(ie)/neq(n)*rmt(n)**2*aimag(
     $       conjg( acof(ie,lm1,natrun)* us(l1,n)
     $            + bcof(ie,lm1,natrun)* uds(l1,n) ) *
     $          cveccof(i,m,ie,lo,natrun)*dulos(lo,n)
     $        + conjg(ccof(m,ie,lo,natrun)*ulos(lo,n)) *
     $            ( aveccof(i,ie,lm1,natrun)* dus(l1,n)
     $            + bveccof(i,ie,lm1,natrun)* duds(l1,n)
     $            + cveccof(i,m,ie,lo,natrun)*dulos(lo,n) )  
     $     - (conjg( aveccof(i,ie,lm1,natrun) *us(l1,n)
     $             + bveccof(i,ie,lm1,natrun) *uds(l1,n) ) *
     $                  ccof(m,ie,lo,natrun)  *dulos(lo,n)
     $        + conjg(cveccof(i,m,ie,lo,natrun)*ulos(lo,n)) *
     $                 ( acof(ie,lm1,natrun)*dus(l1,n)
     $                 + bcof(ie,lm1,natrun)*duds(l1,n)
     $                 + ccof(m,ie,lo,natrun)*dulos(lo,n) ) ) )  
                     END DO
                  ENDDO
              ENDDO
            ENDDO
         END DO
#endif   
c
         DO natrun = natom,natom + neq(n) - 1
c
c  to complete summation over stars of k now sum
c  over all operations which leave (k+G)*R(natrun)*taual(natrun)
c  invariant. We sum over ALL these operations and not only
c  the ones needed for the actual star of k. Should be
c  ok if we divide properly by the number of operations
c  First, we find operation S where RS=T. T -like R- leaves
c  the above scalar product invariant (if S=1 then R=T).
c  R is the operation which generates position of equivalent atom
c  out of position of representative
c  S=R^(-1) T
c  number of ops which leave (k+G)*op*taual invariant: invarind
c  index of inverse operation of R: irinv
c  index of operation T: invarop
c  now, we calculate index of operation S: is
c
c  note, that vector in expression A17,A20 + A21 is a
c  reciprocal lattice vector! other transformation rules
c
c  transform recip vector g-g' into internal coordinates
            DO i = 1,3
               vec(i) = a21(i,natrun)
               vec2(i) = b4(i,natrun)
            END DO
            CALL cotra1(vec,gvint,bmat)
            CALL cotra1(vec2,gvint2,bmat)
            DO i = 1,3
               vecsum(i) = zero
               vecsum2(i) = zero
            END DO
!-gb2002
!            irinv = invtab(ngopr(natrun))
!            DO it = 1,invarind(natrun)
!               is = multab(irinv,invarop(natrun,it))
!c  note, actually we need the inverse of S but -in principle
!c  because {S} is a group and we sum over all S- S should also
!c  work; to be lucid we take the inverse:
!                isinv = invtab(is)
!!               isinv = is
! Rotation is alreadt done in to_pulay, here we work only in the
! coordinate system of the representative atom (natom):
!!        
            DO it = 1,invarind(natom)
               is =invarop(natom,it)
               isinv = invtab(is)
               IF (odi%d1) isinv = ods%ngopr(natom)
!-gb 2002
c  now we have the wanted index of operation with which we have
c  to rotate gv. Note gv is given in cart. coordinates but
c  mrot acts on internal ones
               DO i = 1,3
                  vec(i) = zero
                  vec2(i) = zero
                  DO j = 1,3
                     IF (.NOT.odi%d1) THEN
                        vec(i) = vec(i) + mrot(i,j,isinv)*gvint(j)
                        vec2(i) = vec2(i) + mrot(i,j,isinv)*gvint2(j)
                     ELSE
                       vec(i) = vec(i) + ods%mrot(i,j,isinv)*gvint(j)
                       vec2(i) = vec2(i) + ods%mrot(i,j,isinv)*gvint2(j)
                     END IF
                  END DO
               END DO
               DO i = 1,3
                  vecsum(i) = vecsum(i) + vec(i)
                  vecsum2(i) = vecsum2(i) + vec2(i)
               END DO
c   end operator loop
            END DO
c
c   transform from internal to cart. coordinates
            CALL cotra0(vecsum,starsum,amat)
            CALL cotra0(vecsum2,starsum2,amat)
            DO i = 1,3
               forc_a21(i) = forc_a21(i) + starsum(i)/invarind(natrun)
               forc_b4(i) = forc_b4(i) + starsum2(i)/invarind(natrun)
            END DO
c
c  natrun loop end
         END DO
c
c     sum to existing forces
c
C  NOTE: force() IS REAL AND THEREFORE TAKES ONLY THE
C  REAL PART OF forc_a21(). IN GENERAL, FORCE MUST BE
C  REAL AFTER k-STAR SUMMATION. NOW, WE PUT THE PROPER
C  OPERATIONS INTO REAL SPACE. PROBLEM: WHAT HAPPENS
C  IF IN REAL SPACE THERE IS NO INVERSION ANY MORE?
C  BUT WE HAVE INVERSION IN k-SPACE DUE TO TIME REVERSAL
C  SYMMETRY, E(k)=E(-k)
C  WE ARGUE THAT k-SPACE INVERSION IS AUTOMATICALLY TAKEN
C  INTO ACCOUNT IF FORCE = (1/2)(forc_a21+conjg(forc_a21))
C  BECAUSE TIME REVERSAL SYMMETRY MEANS THAT conjg(PSI)
C  IS ALSO A SOLUTION OF SCHR. EQU. IF PSI IS ONE.
         DO i = 1,3
            force(i,n,jsp) = force(i,n,jsp) + real(forc_a21(i) + 
     +                                              forc_b4(i))
            f_a21(i,n)     = f_a21(i,n)     + forc_a21(i)
            f_b4(i,n)      = f_b4(i,n)      + forc_b4(i)
         END DO
c
c     write result moved to force_a8
c
!         write(*,*) a21(:,n) 
         ENDIF                                            !  IF (l_geo(n)) ...
         natom = natom + neq(n)
   10 CONTINUE
c
      DEALLOCATE ( tdd,tuu,tdu,tud,tuulo,tdulo,tuloulo,ind,a21,b4 )
    
      END SUBROUTINE force_a21
      END MODULE m_forcea21

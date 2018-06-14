      MODULE m_qpwtonmt
c***************************************************************
c     This subroutine calculates the lattice harmonic expansions
c     for the plane wave density inside the spheres.      
c
c             Stefan Bl"ugel  , IFF, Nov. 1997
c***************************************************************
      CONTAINS
      SUBROUTINE qpw_to_nmt(
     >                      memd,nlhd,ntypsd,jmtd,ntypd,n3d,jspd,lmaxd,
     >                      lmax,ntypsy,jri,dx,rmsh,ntype,nop,natd,
     >                      symor,bmat,tau,taual,neq,kv3,mrot,invtab,
     >                      nmem,nlh,mlh,llh,clnu,odi,ods,irank,isize,
     >                      jspin,l_cutoff,ng3,nstr,sk3,fpi,sfp,qpwc,
     X                      rho)
c
      USE m_phasy1
      USE m_sphbes
      USE m_cputime
      USE m_outtime
      USE m_od_types, ONLY : od_inp, od_sym
      USE m_od_phasy

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: irank,isize
      INTEGER, INTENT (IN) :: memd,nlhd,ntypsd,jmtd,ntypd,n3d,jspd,lmaxd
      INTEGER, INTENT (IN) :: jspin,l_cutoff,ng3,ntype,nop,natd
      REAL,    INTENT (IN) :: sfp,fpi
      LOGICAL, INTENT (IN) :: symor
C     ..
C     .. Array Arguments ..
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd),qpwc(n3d)
      INTEGER, INTENT (IN) :: lmax(ntypd),ntypsy(natd),jri(ntypd)
      INTEGER, INTENT (IN) :: nmem(0:nlhd,ntypsd),nlh(ntypsd),nstr(n3d)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),llh(0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: neq(ntypd),kv3(3,n3d),mrot(3,3,nop)
      INTEGER, INTENT (IN) :: invtab(nop)
      REAL,    INTENT (IN) :: bmat(3,3),tau(3,nop),taual(3,natd)
      REAL,    INTENT (IN) :: dx(ntypd),rmsh(jmtd,ntypd),sk3(n3d)
      REAL,    INTENT (INOUT) :: rho(jmtd,0:nlhd,ntypd,jspd)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
C     ..
C     .. Local Scalars ..
      REAL,    PARAMETER :: zero = 0.0
      COMPLEX, PARAMETER :: czero = (0.0,0.0)
      INTEGER in,j,jl,j0,jm,k,l,lh,m,n,nd,nrm,n1,n2,na,lm
      REAL    d0,gr,r0,rr,time1,time2
      COMPLEX cp,sm,cprr2
C     ..
C     .. Local Arrays ..
      COMPLEX pylm( (lmaxd+1)**2, ntypd) !,bsl_c(jmtd,0:lmaxd)
      REAL    bsl(0:lmaxd),bsl_r(jmtd,0:lmaxd),bsl_i(jmtd,0:lmaxd)
      INTEGER mr(ntypd),lmx(ntypd),lmxx(ntypd),ntypsy_o(ntypd)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC conjg,max,min
c
c----> cut-off l-expansion of non-spherical charge contribution
c      from coretails of neighboring atom for l> l_cutoff
c
      na = 1
      DO n = 1,ntype
         lmx(n) = min( lmax(n) , l_cutoff )
         ntypsy_o(n) = ntypsy(na)
         na = na + neq(n)
      END DO
c
c----> identify atoms with the same radial mesh
c
      j0 = 0
      r0 = zero
      d0 = zero
      nrm= 0
      DO n = 1 , ntype
         IF (.not.(jri(n).eq.j0 .and. rmsh(1,n).eq.r0 
     +                          .and. dx(n).eq.d0)) THEN   
             j0 = jri(n)
             r0 = rmsh(1,n)
             d0 = dx(n)
             nrm= nrm + 1
             lmxx(nrm) = lmx(n)
         END IF
         mr(nrm)=n
         lmxx(nrm) = max( lmx(n) , lmx(nrm) )
      END DO
c
c=====> Loop over the g-vectors
c
c ----> g=0 component

      IF (irank == 0) THEN
         CALL cpu_time(time1)
         cp = qpwc(1)*nstr(1)
         DO  n = 1 , ntype
           DO j = 1,jri(n)
             rr = rmsh(j,n)*rmsh(j,n)
             rho(j,0,n,jspin) = rho(j,0,n,jspin) + rr*sfp*cp
           ENDDO
         ENDDO
      ELSE
         rho(:,:,:,jspin) = 0.0
      ENDIF

c ----> g.ne.0 components
C
C     g=|=0 terms : \sum_{g =|= 0} 4 \pi i^l \rho_int(g) r_i^{2} \times
C                    j_{l} (gr_i) \exp{ig\xi_i} Y^*_{lm} (g)
c
c ---->     calculate structure constant for each atom
c     (4pi*i**l/nop(3)*sum(R){exp(iRG(taual-taur)*conjg(ylm(RG)) 
c
      DO 100 k = irank+2,ng3,isize
            cp = qpwc(k)*nstr(k)
            IF (.NOT.odi%d1) THEN
            CALL phasy1(
     >                  ntypd,n3d,natd,nop,lmaxd,ntype,neq,lmax,
     >                  fpi,taual,bmat,kv3,tau,mrot,symor,k,invtab,
     <                  pylm)
            ELSE
c-odim
               CALL od_phasy(
     >              ntypd,n3d,natd,lmaxd,ntype,neq,lmax,
     >              fpi,taual,bmat,kv3,k,odi,ods,
     <              pylm)
c+odim
            END IF
c
            n1 = 1
            DO 200 in = 1 , nrm
               n2 = mr(in)
               DO j = 1,jri(n1)
                  cprr2 = cp*rmsh(j,n1)*rmsh(j,n1)
                  gr = sk3(k)*rmsh(j,n1)
                  CALL sphbes(lmxx(in),gr,bsl)
                  DO jl=0,lmxx(in)
                    bsl_r(j,jl) = bsl(jl) *  real(cprr2)
                    bsl_i(j,jl) = bsl(jl) * aimag(cprr2)
                  END DO
               END DO
               DO n = n1,n2 
                  nd = ntypsy_o(n)
                  DO lh = 0,nlh(nd)
                     l = llh(lh,nd)
                     IF ( l.le.lmx(n) ) THEN
                        sm = czero
                        DO jm = 1,nmem(lh,nd)
                           m  = mlh(jm,lh,nd)
                           lm = l*(l+1) + m + 1
                           sm = sm + conjg(clnu(jm,lh,nd))
     *                              *pylm(lm,n)
                        END DO
                        DO j = 1,jri(n1)
                           rho(j,lh,n,jspin) = rho(j,lh,n,jspin)
     +                            + bsl_r(j,l) *  real(sm)
                        END DO
                        DO j = 1,jri(n1)
                           rho(j,lh,n,jspin) = rho(j,lh,n,jspin)
     +                            - bsl_i(j,l) * aimag(sm)
                        END DO
                     END IF
                  END DO
               END DO
               n1 = n2 + 1 
  200       CONTINUE

  100 CONTINUE
c
      IF (irank == 0) THEN
        CALL cpu_time(time2)
        CALL outtime('qpw_to_nmt: ',time2-time1)
      ENDIF
c
      END SUBROUTINE qpw_to_nmt
      END MODULE m_qpwtonmt

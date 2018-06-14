      MODULE m_stern
c     **************************************************************
c     returns star of recipocal space vector g
c     called by force_a8 - APW+LO package
c     *************************************************************
      CONTAINS
      SUBROUTINE stern(
     >                 nop,tau,g,mrot,bbmat,bmat,
     <                 nst,stg,taup,gl,rstg)

      USE m_constants, ONLY : pimach
      USE m_dotir, ONLY : dotirp
      IMPLICIT NONE
C     ..
C     .. Arguments
      INTEGER, INTENT (IN)  :: nop,g(3),mrot(3,3,nop)
      REAL,    INTENT (IN)  :: tau(3,nop),bbmat(3,3),bmat(3,3)

      INTEGER, INTENT (OUT) :: nst,stg(3,nop)
      REAL,    INTENT (OUT) :: gl,rstg(3,nop)
      COMPLEX, INTENT (OUT) :: taup(nop)
C     ..
C     .. Local Variables
      INTEGER               :: i,m,j,k,l,ind(nop)
      REAL                  :: tpi,tk,s,rg(3)

      tpi = 2.0 * pimach()
      ind(1:nop)  = 0
      taup(1:nop) = 0.0
      nst = 0                                                             

      DO i = 1,3
        rg(i) = real( g(i) )
      ENDDO
      gl = sqrt( dotirp(rg,rg,bbmat) )

      DO i = 1,nop
         tk=0.                                                          
         DO j=1,3                                                     
            tk=tk+tau(j,i)*g(j)*tpi                                     
            k=0                                                         
            DO l=1,3
               k=mrot(l,j,i)*g(l)+k                                       
            ENDDO
            stg(j,i)=k                                                  
         ENDDO
         IF (nst.NE.0) THEN                                              
            DO m = 1,nst                                                   
               DO j = 1,3                                                  
                  IF (stg(j,m).NE.stg(j,i)) GOTO 4                          
               ENDDO
               ind(m)=ind(m)+1
               taup(m)=taup(m) + cmplx(cos(tk),sin(tk))
               GOTO 1
 4             CONTINUE
            ENDDO
         ENDIF 
         nst=nst+1
         DO j=1,3
            stg(j,nst)=stg(j,i)
         ENDDO
         DO j = 1,3
            s = 0
            DO l = 1,3
               s = s + stg(l,nst)*bmat(l,j)
            ENDDO
            rstg(j,nst) = s
         ENDDO

         ind(nst)  = 1
         taup(nst) = cmplx(cos(tk),sin(tk))
 1       CONTINUE
      ENDDO                                                        

      DO i=1,nst                                                     
         taup(i)=taup(i)/ind(i)                                            
      ENDDO

      RETURN
      END SUBROUTINE stern
      END MODULE m_stern

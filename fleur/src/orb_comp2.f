      MODULE m_orbcomp
      CONTAINS
      SUBROUTINE orb_comp(
     >                     nobd,lmd,lmaxd,natd,ntypd,nlod,llod,
     >                     ntype,ne,neq,ddn,uulon,dulon,uloulopn,
     >                     acof,bcof,ccof,nlo,llo,neigd,
     <                     orbcomp,qmtp)
c***********************************************************************
c     Calculates an orbital composition of eigen states
c     
c                                   Yury  Koroteev  2003-12-24 
c***********************************************************************
C                     ABBREVIATIONS
C          dimentions
C nobd                  : in, number of considered bands   
C lmd                   : in, (lmaxd + 1)**2
C natd                  : in, number of atoms in a film
C lmaxd                 : in, max of l    
C ntypd                 : in, number of mt-sphere types
C nlod                  : in, number of local orbitals in mt-sphere types
C llod                  : in, l max for local orbitals in mt-sphere types
C ---------------------------------------------------------------------- 
C neq(ntypd)            : in, number of mt-spheres of the same type
C acof(nobd,0:lmd,natd) : in, a,b  coefficients of linearized 
C bcof(nobd,0:lmd,natd) : in, mt-wavefunctions for each band and atom
C ccof(-llod:llod,nobd, :
C     :      nobd,natd) : in, c coefficients for local orbitals
C ddn(16,ntypd)         : in,   
C uulon(16,ntypd)       : in,   
C dulon(16,ntypd)       : in,  
C uloulopn(16,ntypd)    : in,   
C nlo(ntypd)            : in, 
C llo(nlod,ntypd)       : in,
c-----------------------------------------------------------------------
C orbcomp(nobd,16,natd) : out, an orbital composition of  states
C qmtp(nobd,natd)       : out, the portion of the state in mt-sphere
C-----------------------------------------------------------------------
c
	IMPLICIT NONE
c	..
c	..Scalar Argument
	INTEGER, INTENT  (IN) :: nobd,lmd,lmaxd,natd,ntypd
	INTEGER, INTENT  (IN) :: ntype,ne,nlod,llod,neigd
c	..
c	..Array Arguments
	INTEGER, INTENT  (IN) :: neq(ntypd),nlo(ntypd),llo(nlod,ntypd)
	REAL,    INTENT  (IN) :: ddn(0:lmaxd,ntypd)
	REAL,    INTENT  (IN) :: uulon(nlod,ntypd)
	REAL,    INTENT  (IN) :: dulon(nlod,ntypd)
	REAL,    INTENT  (IN) :: uloulopn(nlod,nlod,ntypd)
        COMPLEX, INTENT  (IN) :: acof(nobd,0:lmd,natd)
        COMPLEX, INTENT  (IN) :: bcof(nobd,0:lmd,natd)
        COMPLEX, INTENT  (IN) :: ccof(-llod:llod,nobd,nlod,natd)
	REAL,    INTENT (OUT) :: orbcomp(neigd,23,natd) 
	REAL,    INTENT (OUT) :: qmtp(neigd,natd) 
c	..
c	..Local Scalars 
	INTEGER  n,mt,ityp,imt,lm,lo
	INTEGER  l,lme,nate,lmaxe,jspe,jspd,nobc,nei
	REAL     c3,c5,g,h,sum,cf
	REAL     ddn0,ddn1,ddn2,ddn3,ddn12,ddn22,ddn32
	COMPLEX  ca00,ca01,ca02,ca03,ca04,ca05,ca06,ca07,ca08,ca09
	COMPLEX  ca10,ca11,ca12,ca13,ca14,ca15,ca16,ca17,ca18,ca19
	COMPLEX  ca20,ca21,ca22
	COMPLEX  cb00,cb01,cb02,cb03,cb04,cb05,cb06,cb07,cb08,cb09
	COMPLEX  cb10,cb11,cb12,cb13,cb14,cb15,cb16,cb17,cb18,cb19
	COMPLEX  cb20,cb21,cb22
	COMPLEX  cc00,cc01,cc02,cc03,cc04,cc05,cc06,cc07,cc08,cc09
	COMPLEX  cc10,cc11,cc12,cc13,cc14,cc15,cc16,cc17,cc18,cc19
	COMPLEX  cc20,cc21,cc22
	COMPLEX  ck00,ck01,ck02,ck03,ck04,ck05,ck06,ck07,ck08,ck09
	COMPLEX  ck10,ck11,ck12,ck13,ck14,ck15,ck16,ck17,ck18,ck19
	COMPLEX  ck20,ck21,ck22
c	..
c	..Local Arrays 
	REAL     comp(23)
c	..
c	Intrinsic Function
	Intrinsic sqrt,conjg
c
	DATA h/0.50/ g/0.0625/
c**************************************************** 
	c5 = sqrt(5.0)
	c3 = sqrt(3.0)
c
	mt=0
	DO 10  ityp = 1,ntype
	   ddn0 = ddn(0,ityp)
	   ddn1 = ddn(1,ityp) 	 
	   ddn2 = ddn(2,ityp) 	 
	   ddn3 = ddn(3,ityp) 
	   DO 20 imt=1,neq(ityp)
	      mt=mt+1
              DO 30 n=1,ne
c
c acof
c   s-states
	         ca00 = acof(n,0,mt)
c   p-states
	         ca01 = acof(n,1,mt) - acof(n,3,mt)
	         ca02 = acof(n,1,mt) + acof(n,3,mt)
	         ca03 = acof(n,2,mt)
c   d-states
	         ca04 = acof(n,4,mt) - acof(n,8,mt)
	         ca05 = acof(n,5,mt) + acof(n,7,mt)
	         ca06 = acof(n,5,mt) - acof(n,7,mt)
	         ca07 = acof(n,4,mt) + acof(n,8,mt)
	         ca08 = acof(n,6,mt)
c
c   f-states: a cubic set (cub) 
c 
	         ca09 = ( acof(n,9,mt)  - acof(n,15,mt) )*c5 -
     -                  ( acof(n,11,mt) - acof(n,13,mt) )*c3
	         ca10 = ( acof(n,9,mt)  + acof(n,15,mt) )*c5 +
     +                  ( acof(n,11,mt) + acof(n,13,mt) )*c3 
	         ca11 =   acof(n,12,mt)
	         ca12 = ( acof(n,9,mt)  + acof(n,15,mt) )*c3 -
     -                  ( acof(n,11,mt) + acof(n,13,mt) )*c5 
	         ca13 =   acof(n,10,mt) + acof(n,14,mt)
	         ca14 = ( acof(n,9,mt)  - acof(n,15,mt) )*c3 +
     +                  ( acof(n,11,mt) - acof(n,13,mt) )*c5 
	         ca15 =   acof(n,10,mt) - acof(n,14,mt) 
c
c   f-states:	a low symmetry set (lss)
c
	         ca16 =  acof(n,11,mt) - acof(n,13,mt)
	         ca17 =  acof(n,11,mt) + acof(n,13,mt)
	         ca18 =  acof(n,12,mt)
	         ca19 =  acof(n,10,mt) - acof(n,14,mt)
 	         ca20 =  acof(n,10,mt) + acof(n,14,mt)
	         ca21 =  acof(n,9,mt)  - acof(n,15,mt)
	         ca22 =  acof(n,9,mt)  + acof(n,15,mt)
c
c bcof
c   s-states
	         cb00 =  bcof(n,0,mt)
c   p-states
	         cb01 =  bcof(n,1,mt) - bcof(n,3,mt) 
	         cb02 =  bcof(n,1,mt) + bcof(n,3,mt) 
	         cb03 =  bcof(n,2,mt)
c   d-states
	         cb04 =  bcof(n,4,mt) - bcof(n,8,mt) 
	         cb05 =  bcof(n,5,mt) + bcof(n,7,mt) 
	         cb06 =  bcof(n,5,mt) - bcof(n,7,mt) 
	         cb07 =  bcof(n,4,mt) + bcof(n,8,mt) 
	         cb08 =  bcof(n,6,mt)
c
c   f-states: a cubic set (cub)
c
	         cb09 = ( bcof(n,9,mt)  - bcof(n,15,mt) )*c5 -
     -                  ( bcof(n,11,mt) - bcof(n,13,mt) )*c3
	         cb10 = ( bcof(n,9,mt)  + bcof(n,15,mt) )*c5 +
     +                  ( bcof(n,11,mt) + bcof(n,13,mt) )*c3 
	         cb11 =   bcof(n,12,mt)
	         cb12 = ( bcof(n,9,mt)  + bcof(n,15,mt) )*c3 -
     -                  ( bcof(n,11,mt) + bcof(n,13,mt) )*c5 
	         cb13 =   bcof(n,10,mt) + bcof(n,14,mt)
	         cb14 = ( bcof(n,9,mt)  - bcof(n,15,mt) )*c3 +
     +                  ( bcof(n,11,mt) - bcof(n,13,mt) )*c5
	         cb15 =   bcof(n,10,mt) - bcof(n,14,mt) 
c
c   f-states:	a low symmetry set (lss)
c
	         cb16 =  bcof(n,11,mt) - bcof(n,13,mt)
	         cb17 =  bcof(n,11,mt) + bcof(n,13,mt)
	         cb18 =  bcof(n,12,mt)
	         cb19 =  bcof(n,10,mt) - bcof(n,14,mt)
	         cb20 =  bcof(n,10,mt) + bcof(n,14,mt)
	         cb21 =  bcof(n,9,mt)  - bcof(n,15,mt)
	         cb22 =  bcof(n,9,mt)  + bcof(n,15,mt)
c------------------------------------------------------------------
c  s
	 comp(1)  =   ca00*conjg(ca00) + cb00*conjg(cb00)*ddn0 
c  p
	 comp(2)  = ( ca01*conjg(ca01) + cb01*conjg(cb01)*ddn1 )*h
	 comp(3)  = ( ca02*conjg(ca02) + cb02*conjg(cb02)*ddn1 )*h
	 comp(4)  =   ca03*conjg(ca03) + cb03*conjg(cb03)*ddn1 
c  d
	 comp(5)  = ( ca04*conjg(ca04) + cb04*conjg(cb04)*ddn2 )*h
	 comp(6)  = ( ca05*conjg(ca05) + cb05*conjg(cb05)*ddn2 )*h
	 comp(7)  = ( ca06*conjg(ca06) + cb06*conjg(cb06)*ddn2 )*h
	 comp(8)  = ( ca07*conjg(ca07) + cb07*conjg(cb07)*ddn2 )*h
	 comp(9)  =   ca08*conjg(ca08) + cb08*conjg(cb08)*ddn2 
c  f: a cubic set
	 comp(10) = ( ca09*conjg(ca09) + cb09*conjg(cb09)*ddn3 )*g       
         comp(11) = ( ca10*conjg(ca10) + cb10*conjg(cb10)*ddn3 )*g
	 comp(12) =   ca11*conjg(ca11) + cb11*conjg(cb11)*ddn3 
	 comp(13) = ( ca12*conjg(ca12) + cb12*conjg(cb12)*ddn3 )*g
	 comp(14) = ( ca13*conjg(ca13) + cb13*conjg(cb13)*ddn3 )*h
	 comp(15) = ( ca14*conjg(ca14) + cb14*conjg(cb14)*ddn3 )*g
	 comp(16) = ( ca15*conjg(ca15) + cb15*conjg(cb15)*ddn3 )*h
c  f: a low symmetry set
	 comp(17) = ( ca16*conjg(ca16) + cb16*conjg(cb16)*ddn3 )*h    
         comp(18) = ( ca17*conjg(ca17) + cb17*conjg(cb17)*ddn3 )*h
	 comp(19) =   ca18*conjg(ca18) + cb18*conjg(cb18)*ddn3 
	 comp(20) = ( ca19*conjg(ca19) + cb19*conjg(cb19)*ddn3 )*h
	 comp(21) = ( ca20*conjg(ca20) + cb20*conjg(cb20)*ddn3 )*h
	 comp(22) = ( ca21*conjg(ca21) + cb21*conjg(cb21)*ddn3 )*h
	 comp(23) = ( ca22*conjg(ca22) + cb22*conjg(cb22)*ddn3 )*h
c--------------------------------------------------------------------
c ccof   ( contributions from local orbitals )
c
	 DO 60 lo = 1,nlo(ityp)
	       l = llo(lo,ityp)
c lo-s
	       IF ( l.EQ.0 )  THEN
	           cc00 = ccof(0,n,lo,mt)
	           ck00 = conjg(cc00)
c
	           comp(1)  =  comp(1)  +
     +           ( ca00*ck00 + cc00*conjg(ca00) )*uulon(lo,ityp) +
     +           ( cb00*ck00 + cc00*conjg(cb00) )*dulon(lo,ityp) +
     +             cc00*ck00*uloulopn(lo,lo,ityp) 
	           GOTO 60
	       ENDIF
c lo-p
	       IF ( l.EQ.1 )  THEN
	           cc01 = ccof(-1,n,lo,mt) - ccof(1,n,lo,mt)
	           cc02 = ccof(-1,n,lo,mt) + ccof(1,n,lo,mt)
	           cc03 = ccof( 0,n,lo,mt)
c
	           ck01 = conjg(cc01)
		   ck02 = conjg(cc02)
		   ck03 = conjg(cc03)
c
		   comp(2) = comp(2)  +
     +             (( ca01*ck01 + cc01*conjg(ca01) )*uulon(lo,ityp) +
     +              ( cb01*ck01 + cc01*conjg(cb01) )*dulon(lo,ityp) +
     +                cc01*ck01*uloulopn(lo,lo,ityp) )*h 	
	           comp(3) = comp(3)  +
     +		   (( ca02*ck02 + cc02*conjg(ca02) )*uulon(lo,ityp) +
     +              ( cb02*ck02 + cc02*conjg(cb02) )*dulon(lo,ityp) +
     +                cc02*ck02*uloulopn(lo,lo,ityp) )*h 
	           comp(4) = comp(4)  +
     +		    ( ca03*ck03 + cc03*conjg(ca03) )*uulon(lo,ityp) +
     +              ( cb03*ck03 + cc03*conjg(cb03) )*dulon(lo,ityp) +
     +                cc03*ck03*uloulopn(lo,lo,ityp) 
	           GOTO 60
	       ENDIF
c lo-d
	       IF ( l.EQ.2 )  THEN
	           cc04 = ccof(-2,n,lo,mt) - ccof(2,n,lo,mt)
	           cc05 = ccof(-1,n,lo,mt) + ccof(1,n,lo,mt)
	           cc06 = ccof(-1,n,lo,mt) - ccof(1,n,lo,mt)
	           cc07 = ccof(-2,n,lo,mt) + ccof(2,n,lo,mt)
	           cc08 = ccof( 0,n,lo,mt)
c
		   ck04 = conjg(cc04)
		   ck05 = conjg(cc05)
		   ck06 = conjg(cc06)
                   ck07 = conjg(cc07)
                   ck08 = conjg(cc08)
c
		   comp(5) = comp(5)  +
     +             (( ca04*ck04 + cc04*conjg(ca04) )*uulon(lo,ityp) +
     +              ( cb04*ck04 + cc04*conjg(cb04) )*dulon(lo,ityp) +
     +                cc04*ck04*uloulopn(lo,lo,ityp) )*h 	
	           comp(6) = comp(6)  +
     +		   (( ca05*ck05 + cc05*conjg(ca05) )*uulon(lo,ityp) +
     +              ( cb05*ck05 + cc05*conjg(cb05) )*dulon(lo,ityp) +
     +                cc05*ck05*uloulopn(lo,lo,ityp) )*h 
	           comp(7) = comp(7)  +
     +		   (( ca06*ck06 + cc06*conjg(ca06) )*uulon(lo,ityp) +
     +              ( cb06*ck06 + cc06*conjg(cb06) )*dulon(lo,ityp) +
     +                cc06*ck06*uloulopn(lo,lo,ityp) )*h
     		   comp(8) = comp(8)  +
     +             (( ca07*ck07 + cc07*conjg(ca07) )*uulon(lo,ityp) +
     +              ( cb07*ck07 + cc07*conjg(cb07) )*dulon(lo,ityp) +
     +                cc07*ck07*uloulopn(lo,lo,ityp) )*h 	
	           comp(9) = comp(9)  +
     +		    ( ca08*ck08 + cc08*conjg(ca08) )*uulon(lo,ityp) +
     +              ( cb08*ck08 + cc08*conjg(cb08) )*dulon(lo,ityp) +
     +                cc08*ck08*uloulopn(lo,lo,ityp)  
		   GOTO 60				
	       ENDIF
c lo-f
	       IF ( l.EQ.3 )  THEN
c
c  a cubic set (cub)
c
	           cc09 = ( ccof(-3,n,lo,mt) - ccof(3,n,lo,mt) )*c5 -
     -                    ( ccof(-1,n,lo,mt) - ccof(1,n,lo,mt) )*c3 
	           cc10 = ( ccof(-3,n,lo,mt) + ccof(3,n,lo,mt) )*c5 +
     +                    ( ccof(-1,n,lo,mt) + ccof(1,n,lo,mt) )*c3 
	           cc11 =   ccof( 0,n,lo,mt)
	           cc12 = ( ccof(-3,n,lo,mt) + ccof(3,n,lo,mt) )*c3 -
     -                    ( ccof(-1,n,lo,mt) + ccof(1,n,lo,mt) )*c5 
	           cc13 =   ccof(-2,n,lo,mt) + ccof(2,n,lo,mt) 
	           cc14 = ( ccof(-3,n,lo,mt) - ccof(3,n,lo,mt) )*c3 +
     +                    ( ccof(-1,n,lo,mt) - ccof(1,n,lo,mt) )*c5
	           cc15 =   ccof(-2,n,lo,mt) - ccof(2,n,lo,mt)
c
                   ck09 = conjg(cc09)
		   ck10 = conjg(cc10)
		   ck11 = conjg(cc11)
                   ck12 = conjg(cc12)
		   ck13 = conjg(cc13)
		   ck14 = conjg(cc14)
		   ck15 = conjg(cc15)
c
		   comp(10) = comp(10)  +
     +             (( ca09*ck09 + cc09*conjg(ca09) )*uulon(lo,ityp) +
     +              ( cb09*ck09 + cc09*conjg(cb09) )*dulon(lo,ityp) +
     +                cc09*ck09*uloulopn(lo,lo,ityp) )*g 	
	           comp(11) = comp(11)  +
     +		   (( ca10*ck10 + cc10*conjg(ca10) )*uulon(lo,ityp) +
     +              ( cb10*ck10 + cc10*conjg(cb10) )*dulon(lo,ityp) +
     +                cc10*ck10*uloulopn(lo,lo,ityp) )*g 
	           comp(12) = comp(12)  +
     +		    ( ca11*ck11 + cc11*conjg(ca11) )*uulon(lo,ityp) +
     +              ( cb11*ck11 + cc11*conjg(cb11) )*dulon(lo,ityp) +
     +                cc11*ck11*uloulopn(lo,lo,ityp) 
	           comp(13) = comp(13)  +
     +             (( ca12*ck12 + cc12*conjg(ca12) )*uulon(lo,ityp) +
     +              ( cb12*ck12 + cc12*conjg(cb12) )*dulon(lo,ityp) +
     +                cc12*ck12*uloulopn(lo,lo,ityp) )*g
	           comp(14) = comp(14)  +
     +		   (( ca13*ck13 + cc13*conjg(ca13) )*uulon(lo,ityp) +
     +              ( cb13*ck13 + cc13*conjg(cb13) )*dulon(lo,ityp) +
     +                cc13*ck13*uloulopn(lo,lo,ityp) )*h 
	           comp(15) = comp(15)  +
     +		   (( ca14*ck14 + cc14*conjg(ca14) )*uulon(lo,ityp) +
     +              ( cb14*ck14 + cc14*conjg(cb14) )*dulon(lo,ityp) +
     +                cc14*ck14*uloulopn(lo,lo,ityp) )*g
     		   comp(16) = comp(16)  +
     +             (( ca15*ck15 + cc15*conjg(ca15) )*uulon(lo,ityp) +
     +              ( cb15*ck15 + cc15*conjg(cb15) )*dulon(lo,ityp) +
     +                cc15*ck15*uloulopn(lo,lo,ityp) )*h
c
c  a low symmetry set (lss)
c
	           cc16 = ccof(-1,n,lo,mt) - ccof(1,n,lo,mt)
	           cc17 = ccof(-1,n,lo,mt) + ccof(1,n,lo,mt)
	           cc18 = ccof( 0,n,lo,mt)
	           cc19 = ccof(-2,n,lo,mt) - ccof(2,n,lo,mt)
	           cc20 = ccof(-2,n,lo,mt) + ccof(2,n,lo,mt)
	           cc21 = ccof(-3,n,lo,mt) - ccof(3,n,lo,mt)
	           cc22 = ccof(-3,n,lo,mt) + ccof(3,n,lo,mt)
c
                   ck16 = conjg(cc16)
                   ck17 = conjg(cc17)
		   ck18 = conjg(cc18)
		   ck19 = conjg(cc19)
		   ck20 = conjg(cc20)
		   ck21 = conjg(cc21)
                   ck22 = conjg(cc22)
c
	           comp(17) = comp(17)  +
     +		   (( ca16*ck16 + cc16*conjg(ca16) )*uulon(lo,ityp) +
     +              ( cb16*ck16 + cc16*conjg(cb16) )*dulon(lo,ityp) +
     +                cc16*ck16*uloulopn(lo,lo,ityp) )*h 
	           comp(18) = comp(18)  +
     +             (( ca17*ck17 + cc17*conjg(ca17) )*uulon(lo,ityp) +
     +              ( cb17*ck17 + cc17*conjg(cb17) )*dulon(lo,ityp) +
     +                cc17*ck17*uloulopn(lo,lo,ityp) )*h
	           comp(19) = comp(19)  +
     +		    ( ca18*ck18 + cc18*conjg(ca18) )*uulon(lo,ityp) +
     +              ( cb18*ck18 + cc18*conjg(cb18) )*dulon(lo,ityp) +
     +                cc18*ck18*uloulopn(lo,lo,ityp)
	           comp(20) = comp(20)  +
     +		   (( ca19*ck19 + cc19*conjg(ca19) )*uulon(lo,ityp) +
     +              ( cb19*ck19 + cc19*conjg(cb19) )*dulon(lo,ityp) +
     +                cc19*ck19*uloulopn(lo,lo,ityp) )*h
     		   comp(21) = comp(21)  +
     +             (( ca20*ck20 + cc20*conjg(ca20) )*uulon(lo,ityp) +
     +              ( cb20*ck20 + cc20*conjg(cb20) )*dulon(lo,ityp) +
     +                cc20*ck20*uloulopn(lo,lo,ityp) )*h 
	           comp(22) = comp(22)  +
     +		   (( ca21*ck21 + cc21*conjg(ca21) )*uulon(lo,ityp) +
     +              ( cb21*ck21 + cc21*conjg(cb21) )*dulon(lo,ityp) +
     +                cc21*ck21*uloulopn(lo,lo,ityp) )*h
	           comp(23) = comp(23)  +
     +		   (( ca22*ck22 + cc22*conjg(ca22) )*uulon(lo,ityp) +
     +              ( cb22*ck22 + cc22*conjg(cb22) )*dulon(lo,ityp) +
     +                cc22*ck22*uloulopn(lo,lo,ityp) )*h
	       ENDIF
  60	 CONTINUE
c-------------------------------------------------------------------
c    calculate an orbital cnomposition in percets
c
	         sum = 0.0
	         DO 40  lm=1,16
	            sum = sum + comp(lm)
 40	         CONTINUE
	         cf = 100.0/sum 
	         qmtp(n,mt) = sum*100.0
	         DO 50  lm=1,23
	            orbcomp(n,lm,mt) = comp(lm)*cf
 50	         CONTINUE
c----------------------------------------------------
 30	      CONTINUE ! bands (n)
 20	   CONTINUE    ! atoms (imt) -> mt (=nat)
 10	CONTINUE       ! types (ityp)
c
      END SUBROUTINE orb_comp
      END MODULE m_orbcomp

      subroutine rhzgrd(ied,dx,ro,ndvgrd,drr,ddrr)
c.....-----------------------------------------------------------------
c     input: equal interval (dx) mesh rv.
cc           ro: for example charge

c     evaluates d(ro)/dr,d{d(ro)/dr}/dr.
c     drr=d(ro)/dr, ddrr=d(drr)/dr.
c     by ndvgrd points numerical derivative.
c     coded by t.asada. june,1995.
c.....-----------------------------------------------------------------
      implicit none
      integer ied,ndvgrd
c.....-----------------------------------------------------------------
c     .. array arguments ..
      real ddrr(ied),drr(ied),ro(ied)
c     ..
c     .. intrinsic functions ..
      intrinsic real
c     ..
c     .. scalar arguments ..
      real dx
c     ..
c     .. scalars in common ..
c     ..
c     .. local scalars ..
      real d,drx0,drx1,drx2,drx3,drxx0,drxx1,drxx2,drxx3,f0,f1,f2,f3,f4,
     +     f5,g1,g2,g3,g4,g5
      integer i,i1,i2,i3,i4,i5,i6,ist,j,nred
c     ..
c     .. statement functions ..
      real f131,f132,f133,f141,f142,f143,f144,f151,f152,f153,f154,f155,
     +     f161,f162,f163,f164,f165,f166,f231,f232,f233,f241,f242,f243,
     +     f244,f251,f252,f253,f254,f255,f261,f262,f263,f264,f265,f266
c.....-----------------------------------------------------------------
c     .. statement function definitions ..

c.....three point formula for the 1st deriv.
      f131(f0,f1,f2,d) = (-3*f0+4*f1-f2)/ (2*d)
      f132(g1,f0,f1,d) = (-1*g1-0*f0+f1)/ (2*d)
      f133(g2,g1,f0,d) = (g2-4*g1+3*f0)/ (2*d)

c.....four point formula for the 1st deriv.
      f141(f0,f1,f2,f3,d) = (-11*f0+18*f1-9*f2+2*f3)/ (6*d)
      f142(g1,f0,f1,f2,d) = (-2*g1-3*f0+6*f1-f2)/ (6*d)
      f143(g2,g1,f0,f1,d) = (g2-6*g1+3*f0+2*f1)/ (6*d)
      f144(g3,g2,g1,f0,d) = (-2*g3+9*g2-18*g1+11*f0)/ (6*d)

c.....five point formula for the 1st deriv.
      f151(f0,f1,f2,f3,f4,d) = (-50*f0+96*f1-72*f2+32*f3-6*f4)/ (24*d)
      f152(g1,f0,f1,f2,f3,d) = (-6*g1-20*f0+36*f1-12*f2+2*f3)/ (24*d)
      f153(g2,g1,f0,f1,f2,d) = (2*g2-16*g1-0*f0+16*f1-2*f2)/ (24*d)
      f154(g3,g2,g1,f0,f1,d) = (-2*g3+12*g2-36*g1+20*f0+6*f1)/ (24*d)
      f155(g4,g3,g2,g1,f0,d) = (6*g4-32*g3+72*g2-96*g1+50*f0)/ (24*d)

c.....six point formula for the 1st deriv.
      f161(f0,f1,f2,f3,f4,f5,d) = (-274*f0+600*f1-600*f2+400*f3-150*f4+
     +                            24*f5)/ (120*d)
      f162(g1,f0,f1,f2,f3,f4,d) = (-24*g1-130*f0+240*f1-120*f2+40*f3-
     +                            6*f4)/ (120*d)
      f163(g2,g1,f0,f1,f2,f3,d) = (6*g2-60*g1-40*f0+120*f1-30*f2+4*f3)/
     +                            (120*d)
      f164(g3,g2,g1,f0,f1,f2,d) = (-4*g3+30*g2-120*g1+40*f0+60*f1-6*f2)/
     +                            (120*d)
      f165(g4,g3,g2,g1,f0,f1,d) = (6*g4-40*g3+120*g2-240*g1+130*f0+
     +                            24*f1)/ (120*d)
      f166(g5,g4,g3,g2,g1,f0,d) = (-24*g5+150*g4-400*g3+600*g2-600*g1+
     +                            274*f0)/ (120*d)

c.....three point formula for the 2nd deriv.
      f231(f0,f1,f2,d) = (f0-2*f1+f2)/ (d*d)
      f232(g1,f0,f1,d) = (g1-2*f0+f1)/ (d*d)
      f233(g2,g1,f0,d) = (g2-2*g1+f0)/ (d*d)

c.....four point formula for the 2nd deriv.
      f241(f0,f1,f2,f3,d) = (6*f0-15*f1+12*f2-3*f3)/ (3*d*d)
      f242(g1,f0,f1,f2,d) = (3*g1-6*f0+3*f1+0*f2)/ (3*d*d)
      f243(g2,g1,f0,f1,d) = (0*g2+3*g1-6*f0+3*f1)/ (3*d*d)
      f244(g3,g2,g1,f0,d) = (-3*g3+2*g2+15*g1+6*f0)/ (3*d*d)
c.....five point formula for the 2nd deriv.
      f251(f0,f1,f2,f3,f4,d) = (35*f0-104*f1+114*f2-56*f3+11*f4)/
     +                         (12*d*d)
      f252(g1,f0,f1,f2,f3,d) = (11*g1-20*f0+6*f1+4*f2-f3)/ (12*d*d)
      f253(g2,g1,f0,f1,f2,d) = (-g2+16*g1-30*f0+16*f1-f2)/ (12*d*d)
      f254(g3,g2,g1,f0,f1,d) = (-g3+4*g2+6*g1-20*f0+11*f1)/ (12*d*d)
      f255(g4,g3,g2,g1,f0,d) = (11*g4-56*g3+114*g2-104*g1+35*f0)/
     +                         (12*d*d)

c.....six point formula for the 2nd deriv.
      f261(f0,f1,f2,f3,f4,f5,d) = (225*f0-770*f1+1070*f2-780*f3+305*f4-
     +                            50*f5)/ (60*d*d)
      f262(g1,f0,f1,f2,f3,f4,d) = (50*g1-75*f0-20*f1+70*f2-30*f3+5*f4)/
     +                            (60*d*d)
      f263(g2,g1,f0,f1,f2,f3,d) = (-5*g2+80*g1-150*f0+80*f1-5*f2+0*f3)/
     +                            (60*d*d)
      f264(g3,g2,g1,f0,f1,f2,d) = (0*g3-5*g2+80*g1-150*f0+80*f1-5*f2)/
     +                            (60*d*d)
      f265(g4,g3,g2,g1,f0,f1,d) = (5*g4-30*g3+70*g2-20*g1-75*f0+50*f1)/
     +                            (60*d*d)
      f266(g5,g4,g3,g2,g1,f0,d) = (-50*g5+305*g4-780*g3+1070*g2-770*g1+
     +                            225*f0)/ (60*d*d)
c.....-----------------------------------------------------------------
c     ..
      ist = 1
c     ..
      if(ied-ist.lt.3) then
          write(16,fmt='(/'' ied too small. ied='',i4)') ied
          STOP 'ied-ist.lt.3'
      endif

      if(ndvgrd.lt.3 .or. ndvgrd.gt.6) then
          write(16,fmt=126) ndvgrd
  126     format (/,' ndvgrd should be ge.4 .or. le.6. ndvgrd=',i3)
          STOP 'ndvgrd.lt.3 .or. ndvgrd.gt.6'
      endif
c.....

      DO i = ist,ied
          drr(i) = 0
          ddrr(i) = 0
      ENDDO

c.....
      i1 = ist
      i2 = ist + 1
      i3 = ist + 2
      i4 = ist + 3
      i5 = ist + 4
      i6 = ist + 5

c.....ro: total(core+val)(up+down) charge density.
c.....drr:d(ro)/dr, ddrr=d(d(ro)/dr)/dr
c.....
      if(ndvgrd.eq.3) then

          drx1 = f131(ro(i1),ro(i2),ro(i3),dx)
          drxx1 = f231(ro(i1),ro(i2),ro(i3),dx)

      elseif(ndvgrd.eq.4) then

          drx1 = f141(ro(i1),ro(i2),ro(i3),ro(i4),dx)
          drxx1 = f241(ro(i1),ro(i2),ro(i3),ro(i4),dx)
          drx2 = f142(ro(i1),ro(i2),ro(i3),ro(i4),dx)
          drxx2 = f242(ro(i1),ro(i2),ro(i3),ro(i4),dx)

      elseif(ndvgrd.eq.5) then

          drx1 = f151(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),dx)
          drxx1 = f251(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),dx)
          drx2 = f152(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),dx)
          drxx2 = f252(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),dx)

      elseif(ndvgrd.eq.6) then

          drx1 = f161(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)
          drxx1 = f261(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)
          drx2 = f162(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)
          drxx2 = f262(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)
          drx3 = f163(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)
          drxx3 = f263(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)

      endif

      drr(i1) = drx1
      ddrr(i1) = drxx1

      if(ndvgrd.gt.3) then

          drr(i2) = drx2
          ddrr(i2) = drxx2

          if(ndvgrd.eq.6) then
              drr(i3) = drx3
              ddrr(i3) = drxx3
          endif

      endif

      nred = real(ndvgrd)/2 + .1

      IF (ied-nred.LE.ist) THEN
          WRITE(16,fmt='(/'' ied-nred.lt.ist. ied,nred,ist='',3i4)')
     +      ied,nred,ist
          STOP 'ied-nred.le.ist'
      ENDIF


      IF (ndvgrd.EQ.3) THEN

         DO j = nred + ist,ied - nred
             drr(j) = f132(ro(j-1),ro(j),ro(j+1),dx)
            ddrr(j) = f232(ro(j-1),ro(j),ro(j+1),dx)
         ENDDO

      ELSEIF (ndvgrd.EQ.4) THEN

         DO j = nred + ist,ied - nred
             drr(j) = f142(ro(j-1),ro(j),ro(j+1),ro(j+2),dx)
            ddrr(j) = f242(ro(j-1),ro(j),ro(j+1),ro(j+2),dx)
         ENDDO

      ELSEIF (ndvgrd.EQ.5) THEN

         DO j = nred + ist,ied - nred
            drr(j) = f153(ro(j-2),ro(j-1),ro(j),ro(j+1),ro(j+2),dx)
           ddrr(j) = f253(ro(j-2),ro(j-1),ro(j),ro(j+1),ro(j+2),dx)
         ENDDO

      ELSEIF (ndvgrd.EQ.6) THEN

         DO j = nred + ist,ied - nred
            drr(j) = f164(ro(j-3),ro(j-2),ro(j-1),ro(j),ro(j+1),
     +                 ro(j+2),dx)
           ddrr(j) = f264(ro(j-3),ro(j-2),ro(j-1),ro(j),ro(j+1),
     +                 ro(j+2),dx)
         ENDDO

      ENDIF

c.....
      if(ndvgrd.eq.3) then

          drx0 = f133(ro(ied-2),ro(ied-1),ro(ied),dx)
          drxx0 = f233(ro(ied-2),ro(ied-1),ro(ied),dx)

      elseif(ndvgrd.eq.4) then

          drx1 = f143(ro(ied-3),ro(ied-2),ro(ied-1),ro(ied),dx)
          drxx1 = f243(ro(ied-3),ro(ied-2),ro(ied-1),ro(ied),dx)
          drx0 = f144(ro(ied-3),ro(ied-2),ro(ied-1),ro(ied),dx)
          drxx0 = f244(ro(ied-3),ro(ied-2),ro(ied-1),ro(ied),dx)

      elseif(ndvgrd.eq.5) then

          drx1 = f154(ro(ied-4),ro(ied-3),ro(ied-2),ro(ied-1),ro(ied),
     +           dx)
          drxx1 = f254(ro(ied-4),ro(ied-3),ro(ied-2),ro(ied-1),ro(ied),
     +            dx)
          drx0 = f155(ro(ied-4),ro(ied-3),ro(ied-2),ro(ied-1),ro(ied),
     +           dx)
          drxx0 = f255(ro(ied-4),ro(ied-3),ro(ied-2),ro(ied-1),ro(ied),
     +            dx)

      elseif(ndvgrd.eq.6) then

          drx2 = f164(ro(ied-5),ro(ied-4),ro(ied-3),ro(ied-2),ro(ied-1),
     +           ro(ied),dx)
          drxx2 = f264(ro(ied-5),ro(ied-4),ro(ied-3),ro(ied-2),
     +            ro(ied-1),ro(ied),dx)

          drx1 = f165(ro(ied-5),ro(ied-4),ro(ied-3),ro(ied-2),ro(ied-1),
     +           ro(ied),dx)
          drxx1 = f265(ro(ied-5),ro(ied-4),ro(ied-3),ro(ied-2),
     +            ro(ied-1),ro(ied),dx)

          drx0 = f166(ro(ied-5),ro(ied-4),ro(ied-3),ro(ied-2),ro(ied-1),
     +           ro(ied),dx)
          drxx0 = f266(ro(ied-5),ro(ied-4),ro(ied-3),ro(ied-2),
     +            ro(ied-1),ro(ied),dx)

      endif

      if(ndvgrd.gt.3) then

          if(ndvgrd.eq.6) then
              drr(ied-2) = drx2
              ddrr(ied-2) = drxx2
          endif

          drr(ied-1) = drx1
          ddrr(ied-1) = drxx1

      endif

      drr(ied) = drx0
      ddrr(ied) = drxx0

      return
      end

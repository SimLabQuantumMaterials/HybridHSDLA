      MODULE m_outtime
      CONTAINS
      SUBROUTINE outtime(name,time)

c**********************************************************************
c     every time this subroutine is called it produces one line of
c     output on the file with the unit-number 2, which consitst
c     of up to 59 characters (from the variable name) and the time.
c     here the time must be in seconds. in the output the time is
c     given in seconds, but also in hours minutes and seconds.
c                                              p.kurz   8.2.96
c**********************************************************************

      IMPLICIT NONE

C     .. Scalar Arguments ..
      REAL time
      CHARACTER*(*) name
C     ..
C     .. Local Scalars ..
      REAL rest,seconds
      INTEGER ihours,iminutes
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC real,int
C     ..

c     calculate time in hours, minutes and seconds

      rest = time

      ihours = int(rest/3600.0)
      rest = rest - real(ihours)*3600

      iminutes = int(rest/60.0)
      seconds = rest - real(iminutes)*60

c     output of the results

      WRITE (2,FMT=8000) name,time,ihours,iminutes,seconds
 8000 FORMAT (a,t60,f9.2,' sec = ',i3,' h ',i2,' min ',f5.2,' sec')

      END SUBROUTINE outtime
      END MODULE m_outtime

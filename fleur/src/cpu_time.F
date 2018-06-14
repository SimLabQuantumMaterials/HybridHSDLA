      MODULE m_cputime
      CONTAINS
      SUBROUTINE cpu_time(time)

c**********************************************************************
c     this soubroutine determines the CPU-time in seconds the 
c     program needed since it has been started. the result is 
c     stored in the variable time.
c     on IBM-AIX using xl-fortran this is done with the mclock
c     command, which has an accuracy of 1/100 of a second.
c     on other machines this command probably needs to be replaced.
c                                             p.kurz  8.2.96
c**********************************************************************

      IMPLICIT NONE
C     .. Scalar Arguments ..
      REAL time
C     ..
C     .. External Functions ..
#if (defined(CPP_AIX)||defined(CPP_APC))
c---> IBM AIX or ATHLON PC
c      INTEGER mclock
c      EXTERNAL mclock
#elif defined (CPP_T90)
c---> Cray T90, J90
      REAL second
      EXTERNAL second
#elif defined (CPP_VPP)
c---> vpp300
      REAL Y1
      EXTERNAL clock
#else
c---> DEC ALPHA
      REAL*4 etime
      REAL*4 tt(2)
      EXTERNAL etime
#endif
C     ..

#ifdef CPP_AIX
c---> IBM RS6000
      time = mclock()*0.01
#elif defined (CPP_APC)
      time = mclock()*0.001
#elif defined (CPP_T90)
c---> Cray T90, J90
      time = second()
#elif defined (CPP_VPP)
c---> vpp300
      CALL CLOCK(Y1,2,2)
      time = Y1*0.000001
#else
c---> DEC ALPHA
      time = etime(tt)
#endif

      END SUBROUTINE cpu_time
      END MODULE m_cputime

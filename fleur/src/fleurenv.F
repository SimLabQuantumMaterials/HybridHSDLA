      MODULE m_fleurenv
!-----------------------------------------------
!    fleur_err: print various error logs and stop fleur
!                 Daniel Wortmann, (05-02-25)
!-----------------------------------------------
      CONTAINS
      !<--S:fleur_err
#if defined(__INTEL_COMPILER)
      !<-- IFC-version
      SUBROUTINE fleur_err(message)
!-----------------------------------------------
!     Version for LINUX compiled with IFC
!             (last modified: 05-02-25) D. Wortmann
!-----------------------------------------------
      USE ifcore
      IMPLICIT NONE
      !<-- Arguments
      CHARACTER*(*)        ::message

      !>
      WRITE(*,*) message
      CALL tracebackqq()
      STOP "FLEUR-stopped"
      END SUBROUTINE
      !>
#elif (defined(CPP_AIX)&&defined(CPP_MPI))
      !<-- Parallel-AIX-Version
      SUBROUTINE fleur_err(message)
!-----------------------------------------------
!     Version for IBM sp4
!             (last modified: 05-02-25) D. Wortmann
!-----------------------------------------------
      include 'mpif.h'
      IMPLICIT NONE
      !<-- Arguments
      CHARACTER*(*)        ::message
      !>
      INTEGER :: ierr
      WRITE(*,*) message
!      CALL xl__trbk  !generate stacktrace
      CALL MPI_ABORT(MPI_COMM_WORLD,1,ierr)
      STOP "FLEUR-STOPPED"
      END SUBROUTINE
      !>
#else
      !<-- Default
      SUBROUTINE fleur_err(message)
!-----------------------------------------------
!     Version for LINUX compiled with IFC
!             (last modified: 05-02-25) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<-- Arguments
      CHARACTER*(*)        ::message
      !>
      WRITE(*,*) message
      STOP "FLEUR-STOPPED"
      END SUBROUTINE
      !>
#endif
      !>
      END MODULE m_fleurenv

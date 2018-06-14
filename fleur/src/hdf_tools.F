      MODULE m_hdf_tools
      INTEGER,PARAMETER::rkind=SELECTED_REAL_KIND(12)
      INTEGER,PARAMETER::ckind=SELECTED_REAL_KIND(12)
c*****************************************************************
c DESC:Some subroutines that wrap around simple hdf-calls:
!    hdf_init()
!    hdf_close()
!    io_[read/write]_real[0-6](did,start,count,data,transprop)
!    io_[read/write]_integer[0-6](did,start,count,data,transprop)
!    io_[read/write]_complex[0-5](did,start,count,data,trans)
!    cleartransprop(trans)
c                          Daniel Wortmann, Thu Nov  7 16:12:40 2002
c*****************************************************************

c*****************************************************************
c     These interfaces define easy to USE generic procedures
c
c*****************************************************************

      INTERFACE io_write_att
      MODULE PROCEDURE io_write_attreal0,io_write_attreal1
     +     ,io_write_attreal2,io_write_attreal3
      MODULE PROCEDURE io_write_attint0,io_write_attint1
     +     ,io_write_attint2,io_write_attint3
      MODULE PROCEDURE io_write_attlog0
      MODULE PROCEDURE io_write_attchar0
      END INTERFACE
      INTERFACE io_read_att
      MODULE PROCEDURE io_read_attreal0,io_read_attreal1
     +     ,io_read_attreal2,io_read_attreal3
      MODULE PROCEDURE io_read_attint0,io_read_attint1
     +     ,io_read_attint2,io_read_attint3
      MODULE PROCEDURE io_read_attlog0
      MODULE PROCEDURE io_read_attchar0
      END INTERFACE
      INTERFACE io_read
      MODULE PROCEDURE io_read_real0,io_read_real1,io_read_real2
     +     ,io_read_real3,io_read_real4,io_read_real5,io_read_real6
      MODULE PROCEDURE io_read_integer0,io_read_integer1
     +     ,io_read_integer2,io_read_integer3,io_read_integer4
     +     ,io_read_integer5,io_read_integer6
      MODULE PROCEDURE io_read_complex0,io_read_complex1
     +     ,io_read_complex2,io_read_complex3,io_read_complex4
     +     ,io_read_complex5
      END INTERFACE
      INTERFACE io_write
      MODULE PROCEDURE io_write_real0,io_write_real1,io_write_real2
     +     ,io_write_real3,io_write_real4,io_write_real5,io_write_real6
      MODULE PROCEDURE io_write_integer0,io_write_integer1
     +     ,io_write_integer2,io_write_integer3,io_write_integer4
     +     ,io_write_integer5,io_write_integer6
      MODULE PROCEDURE io_write_complex0,io_write_complex1
     +     ,io_write_complex2,io_write_complex3,io_write_complex4
     +     ,io_write_complex5
      END INTERFACE



      CONTAINS 

      SUBROUTINE hdf_init()
c*****************************************************************
c DESC:Opens library    No longer needed?!
c*****************************************************************
      USE hdf5
      IMPLICIT NONE
      INTEGER::hdferr
      CALL h5open_f(hdferr)
      !Turn automatic error checking on!
      CALL h5eset_auto_f(1,hdferr)
      !CALL h5init_types(hdferr)
#ifndef CPP_AIX
!      CALL checklib()
#endif
      END SUBROUTINE
      SUBROUTINE hdf_close()
c*****************************************************************
c DESC:Closes library   No longer needed?!
c*****************************************************************
      USE hdf5
      IMPLICIT NONE
      INTEGER::hdferr
      !CALL h5close_types(hdferr)
      CALL h5close_f(hdferr)
      END SUBROUTINE

      !<-- S: io_groupexists(gid,name)RESULT(exists)
      FUNCTION io_dataexists(gid,name)RESULT(exist)
c******************************************
c     checks if dataset called 'name' exists at position gid
c                          D. Wortmann
c******************************************
      USE hdf5
      IMPLICIT NONE
      !<--Arguments
      INTEGER(HID_T),INTENT(IN)  ::gid
      CHARACTER,INTENT(IN)       ::name*(*)
      LOGICAL                    ::exist
      !>
      !<--Locals
      INTEGER(HID_t)::testid
      INTEGER       ::hdferr
      !>
      CALL h5eclear_f(hdferr) 
      CALL h5eset_auto_f(0, hdferr) !No automatic error checking!
      CALL h5dopen_f(gid,name,testid,hdferr)
      exist=(hdferr.EQ.0)
      CALL h5dclose_f(testid,hdferr)
      CALL h5eclear_f(hdferr) 
      CALL h5eset_auto_f(1, hdferr) !Resume automatic error checking!
      END FUNCTION
      !> 

      !<-- S: io_groupexists(gid,name)RESULT(exists)
      FUNCTION io_attexists(gid,name)RESULT(exist)
c******************************************
c     checks if attribute called 'name' exists at position gid
c                          D. Wortmann
c******************************************
      USE hdf5
      IMPLICIT NONE
      !<--Arguments
      INTEGER(HID_T),INTENT(IN)  ::gid
      CHARACTER,INTENT(IN)       ::name*(*)
      LOGICAL                    ::exist
      !>
      !<--Locals
      INTEGER(HID_t)::testid
      INTEGER       ::hdferr
      !>
      CALL h5eclear_f(hdferr) 
      CALL h5eset_auto_f(0, hdferr) !No automatic error checking!
      CALL h5aopen_name_f(gid,name,testid,hdferr)
      exist=(hdferr.EQ.0)
      CALL h5aclose_f(testid,hdferr)
      CALL h5eclear_f(hdferr) 
      CALL h5eset_auto_f(1, hdferr) !Resume automatic error checking!
      END FUNCTION
      !> 

      !<-- S: io_groupexists(gid,name)RESULT(exists)
      FUNCTION io_groupexists(gid,name)RESULT(exist)
c******************************************
c     checks if group called 'name' exists at position gid
c                          D. Wortmann
c******************************************
      USE hdf5
      IMPLICIT NONE
      !<--Arguments
      INTEGER(HID_T),INTENT(IN)  ::gid
      CHARACTER,INTENT(IN)       ::name*(*)
      LOGICAL                    ::exist
      !>
      !<--Locals
      INTEGER(HID_t)::testid
      INTEGER       ::hdferr
      !>
      CALL h5eclear_f(hdferr) 
      CALL h5eset_auto_f(0, hdferr) !No automatic error checking!
      CALL h5gopen_f(gid,name,testid,hdferr)
      exist=(hdferr.EQ.0)
      CALL h5gclose_f(testid,hdferr)
      CALL h5eclear_f(hdferr) 
      CALL h5eset_auto_f(1, hdferr) !Resume automatic error checking!
      END FUNCTION
      !> 

c*****************************************************************
c
c     The following subroutines create a var 
c     
c
c************************************************************
      SUBROUTINE io_createvar(did,name,TYPE,dims,vid)
c*****************************************************************
      USE hdf5
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      CHARACTER,INTENT(IN)       ::name*(*)
      INTEGER(HID_T),INTENT(IN)  ::TYPE
      INTEGER,INTENT(IN)         ::dims(:)
      INTEGER(HID_T),INTENT(OUT) ::vid

      !locals
      INTEGER(HSIZE_T)::DIM(SIZE(dims))
      INTEGER(HID_t)::spaceid
      INTEGER       ::hdferr
      dim=dims
      CALL h5screate_simple_f(SIZE(dim),dim ,spaceid,hdferr)
      CALL h5dcreate_f(did,name,TYPE,spaceid, 
     +                 vid, hdferr)
      CALL h5sclose_f(spaceid,hdferr)
      CALL io_check('io_createvar:'//name,hdferr)
      END SUBROUTINE


c*****************************************************************
c
c     The following subroutines READ or WRITE a CHARACTER*-attributes 
c     
c
c*****************************************************************
      SUBROUTINE io_write_attchar0(did,name,DATA)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      CHARACTER,INTENT(IN)       ::name*(*)
      CHARACTER,INTENT(IN)         ::DATA*(*)
      !locals
      INTEGER(HSIZE_T)::dims(7)
      INTEGER(HID_t)::atid,sid
      INTEGER       ::hdferr
      dims=(/LEN(DATA),0,0,0,0,0,0/)
      CALL h5screate_simple_f(1,dims,sid,hdferr)
      CALL h5acreate_f(did, name,H5T_NATIVE_CHARACTER,sid,atid,hdferr) 
      CALL h5awrite_f(atid,H5T_NATIVE_CHARACTER,DATA,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr)
      CALL h5sclose_f(sid,hdferr)
      CALL io_check('io_write_attchar0:'//name//'-'//DATA,hdferr)
      END SUBROUTINE
c*****************************************************************
      SUBROUTINE io_read_attchar0(did,name,DATA)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      CHARACTER,INTENT(IN)       ::name*(*)
      CHARACTER,INTENT(OUT)      ::DATA*(*)
 
      !locals
      INTEGER(HSIZE_T)::dims(7)
      INTEGER(HID_t)::atid
      INTEGER       ::hdferr,dum
      dims=(/LEN(DATA),0,0,0,0,0,0/)
      CALL h5aopen_name_f(did, name, atid, hdferr)
      CALL h5aread_f(atid,H5T_NATIVE_CHARACTER,DATA,dims, hdferr)
      CALL h5aclose_f(atid,hdferr)
      CALL io_check('io_write_attchar0:'//name//'-'//DATA,hdferr)
      END SUBROUTINE    

   
c*****************************************************************
c
c     The following subroutines READ or WRITE LOGICAL-attributes 
c     
c
c*****************************************************************
      SUBROUTINE io_write_attlog0(did,name,DATA)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      CHARACTER,INTENT(IN)       ::name*(*)
      LOGICAL,INTENT(IN)         ::DATA
      !locals
      INTEGER(HSIZE_T)::dims(7)
      INTEGER(HID_t)::atid,sid
      INTEGER       ::hdferr,dum
      dims=(/1,0,0,0,0,0,0/)
      IF (DATA) THEN
         dum=1
      ELSE
         dum=0
      ENDIF
      CALL h5screate_simple_f(1,dims,sid,hdferr)
      CALL h5acreate_f(did, name,H5T_NATIVE_INTEGER,sid,atid,hdferr) 
      CALL h5awrite_f(atid,H5T_NATIVE_INTEGER,dum,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr)
      CALL h5sclose_f(sid,hdferr)
      CALL io_check('io_write_attlog0'//name,hdferr)
      END SUBROUTINE
c*****************************************************************
      SUBROUTINE io_read_attlog0(did,name,DATA)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      CHARACTER,INTENT(IN)       ::name*(*)
      LOGICAL,INTENT(OUT)        ::DATA
 
      !locals
      INTEGER(HSIZE_T)::dims(7)
      INTEGER(HID_t)::atid
      INTEGER       ::hdferr,dum
      dims=(/1,0,0,0,0,0,0/)
      CALL h5aopen_name_f(did, name, atid, hdferr)
      CALL h5aread_f(atid,H5T_NATIVE_INTEGER,dum,dims, hdferr)
      DATA=(dum==1)
      CALL h5aclose_f(atid,hdferr)
      CALL io_check('io_read_attlog0'//name,hdferr)
      END SUBROUTINE    

   

c*****************************************************************
c
c     The following subroutines READ or WRITE REAL-attributes 
c     
c
c*****************************************************************
      SUBROUTINE io_write_attreal0(did,name,DATA)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      CHARACTER,INTENT(IN)       ::name*(*)
      REAL(rkind),INTENT(IN)    ::DATA
 
      CALL io_write_attreal1(did,name,(/DATA/))
      END SUBROUTINE
c*****************************************************************
      SUBROUTINE io_write_attreal1(did,name,DATA)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      CHARACTER,INTENT(IN)       ::name*(*)
      REAL(rkind),INTENT(IN)    ::DATA(:)
      !locals
      INTEGER(HSIZE_T)::dims(7)
      INTEGER(HID_t)::atid,sid
      INTEGER       ::hdferr
      dims=(/SIZE(DATA),0,0,0,0,0,0/)
      CALL h5screate_simple_f(1,dims,sid,hdferr)
      CALL h5acreate_f(did, name,H5T_NATIVE_DOUBLE,sid,atid,hdferr) 
      CALL h5awrite_f(atid,H5T_NATIVE_DOUBLE,DATA,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr)
      CALL h5sclose_f(sid,hdferr)
      CALL io_check('io_write_attreal1'//name,hdferr)
      END SUBROUTINE
c*****************************************************************
      SUBROUTINE io_write_attreal2(did,name,DATA)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      CHARACTER,INTENT(IN)       ::name*(*)
      REAL(rkind),INTENT(IN)    ::DATA(:,:)
      !locals
      INTEGER(HSIZE_T)::dims(7)
      INTEGER(HID_t)::atid,sid
      INTEGER       ::hdferr
      dims=(/SIZE(DATA,1),SIZE(DATA,2),0,0,0,0,0/)
      CALL h5screate_simple_f(2,dims,sid,hdferr)
      CALL h5acreate_f(did, name,H5T_NATIVE_DOUBLE,sid,atid,hdferr) 
      CALL h5awrite_f(atid,H5T_NATIVE_DOUBLE,DATA,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr)
      CALL h5sclose_f(sid,hdferr)
      CALL io_check('io_write_attreal2'//name,hdferr)
      END SUBROUTINE
c*****************************************************************
      SUBROUTINE io_write_attreal3(did,name,DATA)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      CHARACTER,INTENT(IN)       ::name*(*)
      REAL(rkind),INTENT(IN)    ::DATA(:,:,:)
      !locals
      INTEGER(HSIZE_T)::dims(7)
      INTEGER(HID_t)::atid,sid
      INTEGER       ::hdferr
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),0,0,0,0/)
      CALL h5screate_simple_f(3,dims,sid,hdferr)
      CALL h5acreate_f(did, name,H5T_NATIVE_DOUBLE,sid,atid,hdferr) 
      CALL h5awrite_f(atid,H5T_NATIVE_DOUBLE,DATA,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr)
      CALL h5sclose_f(sid,hdferr)
      CALL io_check('io_write_attreal3'//name,hdferr)
      END SUBROUTINE
c*****************************************************************
      SUBROUTINE io_read_attreal0(did,name,DATA)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      CHARACTER,INTENT(IN)       ::name*(*)
      REAL(rkind),INTENT(OUT)    ::DATA
      !locals
      INTEGER(HSIZE_T)::dims(7)
      INTEGER(HID_t)::atid
      INTEGER       ::hdferr
      REAL          ::buf(1)
      dims=(/1,0,0,0,0,0,0/)
      CALL h5aopen_name_f(did, name, atid, hdferr)
      CALL h5aread_f(atid,H5T_NATIVE_DOUBLE,buf,dims, hdferr)
      DATA=buf(1)
      CALL h5aclose_f(atid,hdferr)
      CALL io_check('io_read_attreal0'//name,hdferr)
      END SUBROUTINE
c*****************************************************************
      SUBROUTINE io_read_attreal1(did,name,DATA)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      CHARACTER,INTENT(IN)       ::name*(*)
      REAL(rkind),INTENT(OUT)    ::DATA(:)
 
      !locals
      INTEGER(HSIZE_T)::dims(7)
      INTEGER(HID_t)::atid
      INTEGER       ::hdferr
      dims=(/SIZE(DATA,1),0,0,0,0,0,0/)
      CALL h5aopen_name_f(did, name, atid, hdferr)
      CALL h5aread_f(atid,H5T_NATIVE_DOUBLE,DATA,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr)
      CALL io_check('io_read_attreal1'//name,hdferr)
      END SUBROUTINE     
c*****************************************************************
      SUBROUTINE io_read_attreal2(did,name,DATA)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      CHARACTER,INTENT(IN)       ::name*(*)
      REAL(rkind),INTENT(OUT)    ::DATA(:,:)

      !locals
      INTEGER(HSIZE_T)::dims(7)
      INTEGER(HID_t)::atid
      INTEGER       ::hdferr
      
      dims=(/SIZE(DATA,1),SIZE(DATA,2),0,0,0,0,0/)
      CALL h5aopen_name_f(did, name, atid, hdferr)
      CALL h5aread_f(atid,H5T_NATIVE_DOUBLE,DATA,dims,hdferr) 
      CALL h5aclose_f(atid,hdferr)
      CALL io_check('io_real_attreal2'//name,hdferr)
      END SUBROUTINE
c*****************************************************************
      SUBROUTINE io_read_attreal3(did,name,DATA)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      CHARACTER,INTENT(IN)       ::name*(*)
      REAL(rkind),INTENT(OUT)    ::DATA(:,:,:)

      !locals
      INTEGER(HSIZE_T)::dims(7)
      INTEGER(HID_t)::atid
      INTEGER       ::hdferr
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),0,0,0,0/)
      CALL h5aopen_name_f(did, name, atid, hdferr)
      CALL h5aread_f(atid,H5T_NATIVE_DOUBLE,DATA,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr)
      CALL io_check('io_read_attreal3'//name,hdferr)
      END SUBROUTINE
c*****************************************************************
c
c     The following subroutines READ or WRITE INTEGER-attributes 
c     
c
c*****************************************************************
      SUBROUTINE io_write_attint0(did,name,DATA)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      CHARACTER,INTENT(IN)       ::name*(*)
      INTEGER,INTENT(IN)    ::DATA
      !locals
      INTEGER(HSIZE_T)::dims(7)
      INTEGER(HID_t)::atid,sid
      INTEGER       ::hdferr
      dims=(/1,0,0,0,0,0,0/)
      CALL h5screate_simple_f(1,dims,sid,hdferr)
      CALL h5acreate_f(did, name,H5T_NATIVE_INTEGER,sid,atid,hdferr) 
      CALL h5awrite_f(atid,H5T_NATIVE_INTEGER,DATA,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr)
      CALL h5sclose_f(sid,hdferr)
      CALL io_check('io_write_attint0'//name,hdferr)
      END SUBROUTINE
c*****************************************************************
      SUBROUTINE io_write_attint1(did,name,DATA)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      CHARACTER,INTENT(IN)       ::name*(*)
      INTEGER,INTENT(IN)    ::DATA(:)
      !locals
      INTEGER(HSIZE_T)::dims(7)
      INTEGER(HID_t)::atid,sid
      INTEGER       ::hdferr
      dims=(/SIZE(DATA),0,0,0,0,0,0/)
      CALL h5screate_simple_f(1,dims,sid,hdferr)
      CALL h5acreate_f(did, name,H5T_NATIVE_INTEGER,sid,atid,hdferr) 
      CALL h5awrite_f(atid,H5T_NATIVE_INTEGER,DATA,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr)
      CALL h5sclose_f(sid,hdferr)
      CALL io_check('io_write_attint1'//name,hdferr)
      END SUBROUTINE
c*****************************************************************
      SUBROUTINE io_write_attint2(did,name,DATA)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      CHARACTER,INTENT(IN)       ::name*(*)
      INTEGER,INTENT(IN)    ::DATA(:,:)
      !locals
      INTEGER(HSIZE_T)::dims(7)
      INTEGER(HID_t)::atid,sid
      INTEGER       ::hdferr
      dims=(/SIZE(DATA,1),SIZE(DATA,2),0,0,0,0,0/)
      CALL h5screate_simple_f(2,dims,sid,hdferr)
      CALL h5acreate_f(did, name,H5T_NATIVE_INTEGER,sid,atid,hdferr) 
      CALL h5awrite_f(atid,H5T_NATIVE_INTEGER,DATA,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr)
      CALL h5sclose_f(sid,hdferr)
      CALL io_check('io_write_attint2'//name,hdferr)
      END SUBROUTINE
c*****************************************************************
      SUBROUTINE io_write_attint3(did,name,DATA)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      CHARACTER,INTENT(IN)       ::name*(*)
      INTEGER,INTENT(IN)    ::DATA(:,:,:)
      !locals
      INTEGER(HSIZE_T)::dims(7)
      INTEGER(HID_t)::atid,sid
      INTEGER       ::hdferr
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),0,0,0,0/)
      CALL h5screate_simple_f(3,dims,sid,hdferr)
      CALL h5acreate_f(did, name,H5T_NATIVE_INTEGER,sid,atid,hdferr) 
      CALL h5awrite_f(atid,H5T_NATIVE_INTEGER,DATA,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr)
      CALL h5sclose_f(sid,hdferr)
      CALL io_check('io_write_attint3'//name,hdferr)
      END SUBROUTINE
c*****************************************************************
      SUBROUTINE io_read_attint0(did,name,DATA)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      CHARACTER,INTENT(IN)       ::name*(*)
      INTEGER,INTENT(OUT)    ::DATA
      !locals
      INTEGER(HSIZE_T)::dims(7)
      INTEGER(HID_t)::atid
      INTEGER       ::hdferr
      dims=(/1,0,0,0,0,0,0/)
      CALL h5aopen_name_f(did, name, atid, hdferr)
      CALL h5aread_f(atid,H5T_NATIVE_INTEGER,DATA,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr)
      CALL io_check('io_read_attint0'//name,hdferr)
      END SUBROUTINE
c*****************************************************************
      SUBROUTINE io_read_attint1(did,name,DATA)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      CHARACTER,INTENT(IN)       ::name*(*)
      INTEGER,INTENT(OUT)    ::DATA(:)
 
      !locals
      INTEGER(HSIZE_T)::dims(7)
      INTEGER(HID_t)::atid
      INTEGER       ::hdferr
      dims=(/SIZE(DATA,1),0,0,0,0,0,0/)
      CALL h5aopen_name_f(did, name, atid, hdferr)
      CALL h5aread_f(atid,H5T_NATIVE_INTEGER,DATA,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr)
      CALL io_check('io_read_attint1'//name,hdferr)
      END SUBROUTINE     
c*****************************************************************
      SUBROUTINE io_read_attint2(did,name,DATA)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      CHARACTER,INTENT(IN)       ::name*(*)
      INTEGER,INTENT(OUT)    ::DATA(:,:)

      !locals
      INTEGER(HSIZE_T)::dims(7)
      INTEGER(HID_t)::atid
      INTEGER       ::hdferr
      
      dims=(/SIZE(DATA,1),SIZE(DATA,2),0,0,0,0,0/)
      CALL h5aopen_name_f(did, name, atid, hdferr)
      CALL h5aread_f(atid,H5T_NATIVE_INTEGER,DATA,dims,hdferr) 
      CALL h5aclose_f(atid,hdferr)
      CALL io_check('io_read_attint2'//name,hdferr)
      END SUBROUTINE
c*****************************************************************
      SUBROUTINE io_read_attint3(did,name,DATA)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      CHARACTER,INTENT(IN)       ::name*(*)
      INTEGER,INTENT(OUT) ::DATA(:,:,:)

      !locals
      INTEGER(HSIZE_T)::dims(7)
      INTEGER(HID_t)::atid
      INTEGER       ::hdferr
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),0,0,0,0/)
      CALL h5aopen_name_f(did, name, atid, hdferr)
      CALL h5aread_f(atid,H5T_NATIVE_INTEGER,DATA,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr)
      CALL io_check('io_read_attint3'//name,hdferr)
      END SUBROUTINE

c*****************************************************************
c
c     The following subroutines READ or WRITE int values 
c     from hdf-file
c
c*****************************************************************
      SUBROUTINE io_read_real0(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      REAL(rkind),INTENT(OUT)    ::DATA
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop
      !locals
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n
      IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims=(/1,0,0,0,0,0,0/)   !
      !check if size of count is ok!
      s=1
      trans=H5p_DEFAULT_F
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(1,dims,memspace,hdferr)
      CALL h5dread_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_read_real0!",hdferr)
      END SUBROUTINE

      SUBROUTINE io_read_real1(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  :: did
      INTEGER,       INTENT(IN)  :: start(:),COUNT(:)
      REAL(rkind),   INTENT(OUT) :: DATA(:)
      INTEGER(HID_T), INTENT(IN), OPTIONAL :: transprop
      !locals
      INTEGER(HSIZE_t) :: dims(7),foffset(SIZE(start))
      INTEGER(HSIZE_t) :: fcount(SIZE(count))      
      INTEGER(HID_t)   :: trans,fspace,memspace
      INTEGER          :: hdferr,s,n
      IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims=(/SIZE(DATA,1),0,0,0,0,0,0/)   !
      !check if size of count is ok!
      s=1
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(1,dims,memspace,hdferr)
      CALL h5dread_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_read_real1!",hdferr)

      END SUBROUTINE

      SUBROUTINE io_read_real2(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      IMPLICIT NONE

      INTEGER(HID_T), INTENT(IN)  :: did
      INTEGER,        INTENT(IN)  :: start(:),COUNT(:)
      REAL,           INTENT(OUT) :: DATA(:,:)
      INTEGER(HID_T), INTENT(IN), OPTIONAL :: transprop
! locals
      INTEGER(HSIZE_t) :: dims(7),foffset(SIZE(start))
      INTEGER(HSIZE_t) :: fcount(SIZE(count))      
      INTEGER(HID_t)   :: trans,fspace,memspace
      INTEGER          :: hdferr,s,n

       IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF

      foffset=start-1  !
      fcount=count     !
      dims=(/SIZE(DATA,1),SIZE(DATA,2),0,0,0,0,0/)   !
      !check if size of count is ok!
      IF (ANY(count<1)) RETURN !read nothing
      s=1
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(2,dims,memspace,hdferr)
      data = 0.0
!     write(*,*) dims,hdferr,memspace,fspace,trans
      CALL h5dread_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_read_real2!",hdferr)

      END SUBROUTINE

      SUBROUTINE io_read_real3(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      REAL(rkind),INTENT(OUT)    ::DATA(:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop
      !locals
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n

      IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),0,0,0,0/)   !
      !check if size of count is ok!
      s=1
      IF (ANY(count<1)) RETURN !read nothing
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(3,dims,memspace,hdferr)
      CALL h5dread_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_read_real3!",hdferr)

      END SUBROUTINE

      SUBROUTINE io_read_real4(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      IMPLICIT NONE

      INTEGER(HID_T),INTENT(IN)  :: did
      INTEGER,       INTENT(IN)  :: start(:),COUNT(:)
      REAL,          INTENT(OUT) :: DATA(:,:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL :: transprop
! locals
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n

      IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4),0,0,0/)
      !check if size of count is ok!
      IF (ANY(count<1)) RETURN !read nothing
      s=1
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(4,dims,memspace,hdferr)
      CALL h5dread_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_read_real4!",hdferr)
      END SUBROUTINE

      SUBROUTINE io_read_real5(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      REAL(rkind),INTENT(OUT)    ::DATA(:,:,:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop
      !locals
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n
       IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4) !
     $     ,SIZE(DATA,5),0,0/) 
      !check if size of count is ok!
      IF (ANY(count<1)) RETURN !read nothing
      s=1
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(5,dims,memspace,hdferr)
      CALL h5dread_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_read_real5!",hdferr)

      END SUBROUTINE

      SUBROUTINE io_read_real6(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      REAL(rkind),INTENT(OUT)    ::DATA(:,:,:,:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop
      !locals
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n
      IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims(:)=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4) !
     $     ,SIZE(DATA,5),SIZE(DATA,6),0/) 
      !check if size of count is ok!
      IF (ANY(count<1)) RETURN !read nothing
      s=1
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(6,dims,memspace,hdferr)
      CALL h5dread_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_read_real6!",hdferr)

      END SUBROUTINE

!----------------------------------------------------------------------

      SUBROUTINE io_write_real0(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE

!arguments
      INTEGER(HID_T), INTENT(IN)  :: did
      INTEGER,        INTENT(IN)  :: start(:),COUNT(:)
      REAL(rkind),    INTENT(IN)  :: DATA
      INTEGER(HID_T), INTENT(IN), OPTIONAL :: transprop
!locals
      INTEGER(HSIZE_t) :: dims(7)
      INTEGER(HSIZE_t) :: foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)   :: trans,fspace,memspace
      INTEGER          :: hdferr,s,n

      dims = (/1,0,0,0,0,0,0/)  ! write a single real 

      IF (.NOT.PRESENT(transprop)) THEN
         trans = gettransprop()
      ELSE
         trans = transprop
      ENDIF
      foffset = start-1
      fcount  = count
      IF (ANY(count<1)) RETURN ! write nothing
      s=1                      ! check if size of count is ok
      DO n = 1, SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.1) STOP 'Missmatch of sizes'
!do I/O
      CALL h5dget_space_f(
     >                    did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(
     >                           fspace,H5S_SELECT_SET_F,
     >                           foffset,fcount,hdferr)
      CALL h5screate_simple_f(
     >                        1,dims,memspace,hdferr)
      CALL h5dwrite_f(
     >                did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,
     >                memspace,fspace,trans)
      CALL h5sclose_f(
     >                memspace,hdferr)
      CALL h5sclose_f(
     >                fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_write_real0!",hdferr)
      END SUBROUTINE io_write_real0

!----------------------------------------------------------------------

      SUBROUTINE io_write_real1(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE

!arguments
      INTEGER(HID_T), INTENT(IN)  :: did
      INTEGER,        INTENT(IN)  :: start(:),COUNT(:)
      REAL(rkind),    INTENT(IN)  :: DATA(:)
      INTEGER(HID_T), INTENT(IN), OPTIONAL :: transprop
!locals
      INTEGER(HSIZE_t) :: dims(7)
      INTEGER(HSIZE_t) :: foffset(SIZE(start)),fcount(SIZE(count))
      INTEGER(HID_t)   :: trans,fspace,memspace
      INTEGER          :: hdferr,s,n

      dims = (/SIZE(DATA,1),0,0,0,0,0,0/)   ! write 1-dim array 

      IF (.NOT.PRESENT(transprop)) THEN
         trans = gettransprop()
      ELSE
         trans = transprop
      ENDIF
      foffset = start-1
      fcount  = count
      IF (ANY(count<1)) RETURN ! write nothing
      s=1                      ! check if size of count is ok
      DO n = 1, SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'

!do I/O
      CALL h5dget_space_f(
     >                    did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(
     >                           fspace,H5S_SELECT_SET_F,
     >                           foffset,fcount,hdferr)
      CALL h5screate_simple_f(
     >                        1,dims,memspace,hdferr)
      CALL h5dwrite_f(
     >                did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,
     >                memspace,fspace,trans)
      CALL h5sclose_f(
     >                memspace,hdferr)
      CALL h5sclose_f(
     >                fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_write_real1!",hdferr)

      END SUBROUTINE io_write_real1

!----------------------------------------------------------------------

      SUBROUTINE io_write_real2(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      REAL(rkind),INTENT(IN)    ::DATA(:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop
      !locals
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n
      IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims=(/SIZE(DATA,1),SIZE(DATA,2),0,0,0,0,0/)   !
      !check if size of count is ok!
      IF (ANY(count<1)) RETURN !read nothing
      s=1
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(2,dims,memspace,hdferr)
      CALL h5dwrite_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_write_real2!",hdferr)

      END SUBROUTINE

      SUBROUTINE io_write_real3(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      IMPLICIT NONE

      INTEGER(HID_T),INTENT(IN) :: did
      INTEGER,       INTENT(IN) :: start(:),COUNT(:)
      REAL,          INTENT(IN) :: DATA(:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL :: transprop
!  locals
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n
      IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),0,0,0,0/)   !
      !check if size of count is ok!
      IF (ANY(count<1)) RETURN !read nothing
      s=1
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(3,dims,memspace,hdferr)
      CALL h5dwrite_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_write_real3!",hdferr)
      END SUBROUTINE

      SUBROUTINE io_write_real4(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      REAL(rkind),INTENT(IN)    ::DATA(:,:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop
      !locals
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n
      IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4),0,0,0/)
      !check if size of count is ok!
      IF (ANY(count<1)) RETURN !read nothing
      s=1
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(4,dims,memspace,hdferr)
      CALL h5dwrite_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_write_real4!",hdferr)

      END SUBROUTINE

      SUBROUTINE io_write_real5(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      REAL(rkind),INTENT(IN)    ::DATA(:,:,:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop
      !locals
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n
      IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4) !
     $     ,SIZE(DATA,5),0,0/) 
      !check if size of count is ok!
      IF (ANY(count<1)) RETURN !read nothing
      s=1
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(5,dims,memspace,hdferr)
      CALL h5dwrite_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_write_real5!",hdferr)
      END SUBROUTINE

      SUBROUTINE io_write_real6(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      REAL(rkind),INTENT(IN)    ::DATA(:,:,:,:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop
      !locals
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n
      IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims(:)=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4) !
     $     ,SIZE(DATA,5),SIZE(DATA,6),0/) 
      !check if size of count is ok!
      IF (ANY(count<1)) RETURN !read nothing
      s=1
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(6,dims,memspace,hdferr)
      CALL h5dwrite_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_write_real6!",hdferr)

      END SUBROUTINE

c*****************************************************************
c
c     The following subroutines READ or WRITE INTEGER values 
c     from hdf-file
c
c*****************************************************************
      SUBROUTINE io_read_integer0(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      INTEGER,INTENT(OUT)    ::DATA
 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop
      !locals
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n
      IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims=(/1,0,0,0,0,0,0/)   !
      !check if size of count is ok!
      IF (ANY(count<1)) RETURN !read nothing
      s=1
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(1,dims,memspace,hdferr)
      CALL h5dread_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_read_integer0!",hdferr)
      END SUBROUTINE

      SUBROUTINE io_read_integer1(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      INTEGER,INTENT(OUT)    ::DATA(:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop
      !locals
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n
      IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims=(/SIZE(DATA,1),0,0,0,0,0,0/)   !
      !check if size of count is ok!
      IF (ANY(count<1)) RETURN !read nothing
      s=1
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(1,dims,memspace,hdferr)
      CALL h5dread_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_read_integer1!",hdferr)
      END SUBROUTINE

      SUBROUTINE io_read_integer2(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      INTEGER,INTENT(OUT)    ::DATA(:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop
      !locals
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n
      IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims=(/SIZE(DATA,1),SIZE(DATA,2),0,0,0,0,0/)   !
      !check if size of count is ok!
      IF (ANY(count<1)) RETURN !read nothing
      s=1
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(2,dims,memspace,hdferr)
      CALL h5dread_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_read_integer2!",hdferr)
      END SUBROUTINE

      SUBROUTINE io_read_integer3(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      INTEGER,INTENT(OUT)    ::DATA(:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop
      !locals
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n
      IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),0,0,0,0/)   !
      !check if size of count is ok!
      IF (ANY(count<1)) RETURN !read nothing
      s=1
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(3,dims,memspace,hdferr)
      CALL h5dread_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_read_integer3!",hdferr)
      END SUBROUTINE

      SUBROUTINE io_read_integer4(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      INTEGER,INTENT(OUT)    ::DATA(:,:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop
      !locals
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n
      IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4),0,0,0/)
      !check if size of count is ok!
      IF (ANY(count<1)) RETURN !read nothing
      s=1
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(4,dims,memspace,hdferr)
      CALL h5dread_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_read_integer4!",hdferr)

      END SUBROUTINE

      SUBROUTINE io_read_integer5(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      INTEGER,INTENT(OUT)    ::DATA(:,:,:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop
      !locals
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n
      IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4) !
     $     ,SIZE(DATA,5),0,0/) 
      !check if size of count is ok!
      IF (ANY(count<1)) RETURN !read nothing
      s=1
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(5,dims,memspace,hdferr)
      CALL h5dread_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_read_integer5!",hdferr)
      END SUBROUTINE

      SUBROUTINE io_read_integer6(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      INTEGER,INTENT(OUT)    ::DATA(:,:,:,:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop
      !locals
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n
      IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims(:)=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4) !
     $     ,SIZE(DATA,5),SIZE(DATA,6),0/) 
      !check if size of count is ok!
      IF (ANY(count<1)) RETURN !read nothing
      s=1
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(6,dims,memspace,hdferr)
      CALL h5dread_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_read_integer6!",hdferr)
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE io_write_integer0(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE

! arguments
      INTEGER(HID_T),INTENT(IN)  :: did                ! setid
      INTEGER,       INTENT(IN)  :: start(:),COUNT(:)
      INTEGER,       INTENT(IN)  :: DATA
      INTEGER(HID_T),INTENT(IN),OPTIONAL:: transprop

! locals
      INTEGER(HSIZE_T) :: dims(7)
      INTEGER(HSIZE_T) :: foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_T)   :: trans,fspace,memspace
      INTEGER          :: hdferr,s,n

      dims = (/1,0,0,0,0,0,0/) ! write a single integer

      IF (.NOT.PRESENT(transprop)) THEN
         trans = gettransprop()          ! set default
      ELSE                               ! transfer-properties
         trans = transprop
      ENDIF
      foffset = start-1        
      fcount  = count 
      IF (ANY(count<1)) RETURN ! write nothing
      s=1                      ! check if size of count is ok
      DO n = 1, SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.1) STOP 'Missmatch of sizes'

!do I/O
      CALL h5dget_space_f(
     >                    did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(
     >                           fspace,H5S_SELECT_SET_F,
     >                           foffset,fcount,hdferr)
      CALL h5screate_simple_f(
     >                        1,dims,memspace,hdferr)
      CALL h5dwrite_f(
     >                did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,
     >                memspace,fspace,trans)
      CALL h5sclose_f(
     >                memspace,hdferr)
      CALL h5sclose_f(
     >                fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_write_integer0:",hdferr)

      END SUBROUTINE  io_write_integer0

!----------------------------------------------------------------------
      SUBROUTINE io_write_integer1(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE

!arguments
      INTEGER(HID_T),INTENT(IN)  :: did                ! setid
      INTEGER,       INTENT(IN)  :: start(:),COUNT(:)
      INTEGER,       INTENT(IN)  :: DATA(:)            ! 1-dim 
      INTEGER(HID_T),INTENT(IN),OPTIONAL:: transprop

! locals
      INTEGER(HSIZE_T) :: dims(7)
      INTEGER(HSIZE_T) :: foffset(SIZE(start)),fcount(SIZE(count))
      INTEGER(HID_T)   :: trans,fspace,memspace
      INTEGER          :: hdferr,s,n

      dims = (/SIZE(DATA,1),0,0,0,0,0,0/)   ! write 1-dim array

      IF (.NOT.PRESENT(transprop)) THEN
         trans = gettransprop()
      ELSE
         trans = transprop
      ENDIF
      foffset = start-1
      fcount  = count
      IF (ANY(count<1)) RETURN ! write nothing
      s=1                      ! check if size of count is ok
      DO n = 1, SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'

!do I/O
      CALL h5dget_space_f(
     >                    did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(
     >                           fspace,H5S_SELECT_SET_F,
     >                           foffset,fcount,hdferr)
      CALL h5screate_simple_f(
     >                        1,dims,memspace,hdferr)
      CALL h5dwrite_f(
     >                did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,
     >                memspace,fspace,trans)
      CALL h5sclose_f(
     >                memspace,hdferr)
      CALL h5sclose_f(
     >                fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_write_integer1:",hdferr)

      END SUBROUTINE  io_write_integer1

!----------------------------------------------------------------------

      SUBROUTINE io_write_integer2(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      INTEGER,INTENT(IN)    ::DATA(:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop
      !locals
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n
      IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims=(/SIZE(DATA,1),SIZE(DATA,2),0,0,0,0,0/)   !
      !check if size of count is ok!
      IF (ANY(count<1)) RETURN !read nothing
      s=1
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(2,dims,memspace,hdferr)
      CALL h5dwrite_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_write_integer2:",hdferr)

      END SUBROUTINE

      SUBROUTINE io_write_integer3(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      INTEGER,INTENT(IN)    ::DATA(:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop
      !locals
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n
      IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),0,0,0,0/)   !
      !check if size of count is ok!
      IF (ANY(count<1)) RETURN !read nothing
      s=1
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(3,dims,memspace,hdferr)
      CALL h5dwrite_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_write_integer3:",hdferr)
      END SUBROUTINE

      SUBROUTINE io_write_integer4(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      INTEGER,INTENT(IN)    ::DATA(:,:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop
      !locals
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n
      IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4),0,0,0/)
      !check if size of count is ok!
      IF (ANY(count<1)) RETURN !read nothing
      s=1
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(4,dims,memspace,hdferr)
      CALL h5dwrite_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_write_integer4:",hdferr)
      END SUBROUTINE

      SUBROUTINE io_write_integer5(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      INTEGER,INTENT(IN)    ::DATA(:,:,:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop
      !locals
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n
      IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4) !
     $     ,SIZE(DATA,5),0,0/) 
      !check if size of count is ok!
      IF (ANY(count<1)) RETURN !read nothing
      s=1
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(5,dims,memspace,hdferr)
      CALL h5dwrite_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_write_integer5:",hdferr)
      END SUBROUTINE

      SUBROUTINE io_write_integer6(did,start,count,DATA,transprop)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      INTEGER,INTENT(IN)    ::DATA(:,:,:,:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop
      !locals
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fcount(SIZE(count))      
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n
      IF (.NOT.PRESENT(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fcount=count     !
      dims(:)=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4) !
     $     ,SIZE(DATA,5),SIZE(DATA,6),0/) 
      !check if size of count is ok!
      IF (ANY(count<1)) RETURN !read nothing
      s=1
      DO n=1,SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'
      !DO IO
      CALL h5dget_space_f(did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fcount
     +     ,hdferr)
      CALL h5screate_simple_f(6,dims,memspace,hdferr)
      CALL h5dwrite_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace
     +     ,fspace,trans)
      CALL h5sclose_f(memspace,hdferr)
      CALL h5sclose_f(fspace,hdferr)
      CALL cleartransprop(trans)
      CALL io_check("io_write_integer6:",hdferr)
      END SUBROUTINE
c*****************************************************************
c
c     The following subroutines READ or WRITE COMPLEX values 
c     from hdf-file
c
c*****************************************************************

      SUBROUTINE io_read_complex0(did,start,count,DATA,trans)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      COMPLEX(ckind),INTENT(OUT)    ::DATA
      INTEGER(HID_T),INTENT(IN),OPTIONAL::trans
      !locals
      INTEGER::foffset(SIZE(start))      
      REAL          ::a,b
      foffset=start  !
      !DO 2 calls to read real values
      WHERE(start<0) 
         foffset=1
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_read_real0(did,foffset,count,a,trans)
      ELSE
         CALL io_read_real0(did,foffset,count,a)
      ENDIF
      WHERE(start<0) 
         foffset=2
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_read_real0(did,foffset,count,b,trans)
      ELSE
         CALL io_read_real0(did,foffset,count,b)
      ENDIF
      DATA=CMPLX(a,b)
      END SUBROUTINE
      SUBROUTINE io_read_complex1(did,start,count,DATA,trans)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      COMPLEX(ckind),INTENT(OUT)    ::DATA(:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::trans
      !locals
      INTEGER::foffset(SIZE(start))      
      REAL,ALLOCATABLE::a(:),b(:)
      ALLOCATE(A(SIZE(DATA,1)),B(SIZE(DATA,1)))
      foffset=start  !
      !DO 2 calls to read real values
      WHERE(start<0) 
         foffset=1
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_read_real1(did,foffset,count,a,trans)
      ELSE
         CALL io_read_real1(did,foffset,count,a)
      ENDIF
      WHERE(start<0) 
         foffset=2
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_read_real1(did,foffset,count,b,trans)
      ELSE
         CALL io_read_real1(did,foffset,count,b)
      ENDIF
      DATA=CMPLX(a,b)
      DEALLOCATE(a,b)
      END SUBROUTINE
!----------------------------------------------------------------------

      SUBROUTINE io_read_complex2(did,start,count,DATA,trans)
c*****************************************************************
      USE hdf5
      IMPLICIT NONE
! arguments
      INTEGER(HID_T),INTENT(IN)  :: did
      INTEGER,       INTENT(IN)  :: start(:),COUNT(:)
      COMPLEX,       INTENT(OUT) :: DATA(:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL :: trans
! locals
      INTEGER           :: i,j
      INTEGER           :: foffset(SIZE(start))      
      REAL, ALLOCATABLE :: a(:,:),b(:,:)

      ALLOCATE (A(SIZE(DATA,1),SIZE(DATA,2)),
     +          B(SIZE(DATA,1),SIZE(DATA,2)))

      foffset=start  !
      !DO 2 calls to read real values
      WHERE(start<0) 
         foffset=1
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_read_real2(did,foffset,count,a,trans)
      ELSE
         CALL io_read_real2(did,foffset,count,a)
      ENDIF
      WHERE(start<0) 
         foffset=2
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_read_real2(did,foffset,count,b,trans)
      ELSE
         CALL io_read_real2(did,foffset,count,b)
      ENDIF
      DATA=CMPLX(a,b)
      DEALLOCATE(a,b)
      END SUBROUTINE
      SUBROUTINE io_read_complex3(did,start,count,DATA,trans)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      COMPLEX(ckind),INTENT(OUT)    ::DATA(:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::trans
      !locals
      INTEGER::foffset(SIZE(start))      
      REAL,ALLOCATABLE::a(:,:,:),b(:,:,:)
      ALLOCATE(A(SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3)),B(SIZE(DATA,1)
     +     ,SIZE(DATA,2),SIZE(DATA,3)))
      foffset=start  !
      !DO 2 calls to read real values
      WHERE(start<0) 
         foffset=1
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_read_real3(did,foffset,count,a,trans)
      ELSE
         CALL io_read_real3(did,foffset,count,a)
      ENDIF
      WHERE(start<0) 
         foffset=2
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_read_real3(did,foffset,count,b,trans)
      ELSE
         CALL io_read_real3(did,foffset,count,b)
      ENDIF
      DATA=CMPLX(a,b)
      DEALLOCATE(a,b)
      END SUBROUTINE
      SUBROUTINE io_read_complex4(did,start,count,DATA,trans)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      COMPLEX(ckind),INTENT(OUT)    ::DATA(:,:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::trans
      !locals
      INTEGER::foffset(SIZE(start))      
      REAL,ALLOCATABLE::a(:,:,:,:),b(:,:,:,:)
      ALLOCATE(A(SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4)),
     $     B(SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4)))
      foffset=start  !
      !DO 2 calls to read real values
      WHERE(start<0) 
         foffset=1
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_read_real4(did,foffset,count,a,trans)
      ELSE
         CALL io_read_real4(did,foffset,count,a)
      ENDIF
      WHERE(start<0) 
         foffset=2
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_read_real4(did,foffset,count,b,trans)
      ELSE
         CALL io_read_real4(did,foffset,count,b)
      ENDIF
      DATA=CMPLX(a,b)
      DEALLOCATE(a,b)
      END SUBROUTINE
      SUBROUTINE io_read_complex5(did,start,count,DATA,trans)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      COMPLEX(ckind),INTENT(OUT)    ::DATA(:,:,:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::trans
      !locals
      INTEGER::foffset(SIZE(start))      
      REAL,ALLOCATABLE::a(:,:,:,:,:),b(:,:,:,:,:)
      ALLOCATE(A(SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4)
     $     ,SIZE(DATA,5)),B(SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3)
     $     ,SIZE(DATA,4),SIZE(DATA,5)))
      foffset=start  !
      !DO 2 calls to read real values
      WHERE(start<0) 
         foffset=1
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_read_real5(did,foffset,count,a,trans)
      ELSE
         CALL io_read_real5(did,foffset,count,a)
      ENDIF
      WHERE(start<0) 
         foffset=2
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_read_real5(did,foffset,count,b,trans)
      ELSE
         CALL io_read_real5(did,foffset,count,b)
      ENDIF
      DATA=CMPLX(a,b)
      DEALLOCATE(a,b)
      END SUBROUTINE
      SUBROUTINE io_write_complex0(did,start,count,DATA,trans)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      COMPLEX(ckind),INTENT(IN)    ::DATA
      INTEGER(HID_T),INTENT(IN),OPTIONAL::trans
      !locals
      INTEGER::foffset(SIZE(start))      
      foffset=start  !
      !DO 2 calls to read real values
      WHERE(start<0) 
         foffset=1
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_write_real0(did,foffset,count,REAL(DATA),trans)
      ELSE
         CALL io_write_real0(did,foffset,count,REAL(DATA))
      ENDIF
      WHERE(start<0) 
         foffset=2
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_write_real0(did,foffset,count,AIMAG(DATA),trans)
      ELSE
         CALL io_write_real0(did,foffset,count,AIMAG(DATA))
      ENDIF
      END SUBROUTINE
      SUBROUTINE io_write_complex1(did,start,count,DATA,trans)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      COMPLEX(ckind),INTENT(IN)    ::DATA(:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::trans
      !locals
      INTEGER::foffset(SIZE(start))      
      foffset=start  !
      !DO 2 calls to read real values
      WHERE(start<0) 
         foffset=1
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_write_real1(did,foffset,count,REAL(DATA),trans)
      ELSE
         CALL io_write_real1(did,foffset,count,REAL(DATA))
      ENDIF
      WHERE(start<0) 
         foffset=2
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_write_real1(did,foffset,count,AIMAG(DATA),trans)
      ELSE
         CALL io_write_real1(did,foffset,count,AIMAG(DATA))
      ENDIF
      END SUBROUTINE
      SUBROUTINE io_write_complex2(did,start,count,DATA,trans)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      COMPLEX(ckind),INTENT(IN)    ::DATA(:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::trans
      !locals
      INTEGER::foffset(SIZE(start))      
      foffset=start  !
      !DO 2 calls to read real values
      WHERE(start<0) 
         foffset=1
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_write_real2(did,foffset,count,REAL(DATA),trans)
      ELSE
         CALL io_write_real2(did,foffset,count,REAL(DATA))
      ENDIF
      WHERE(start<0) 
         foffset=2
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_write_real2(did,foffset,count,AIMAG(DATA),trans)
      ELSE
         CALL io_write_real2(did,foffset,count,AIMAG(DATA))
      ENDIF
      END SUBROUTINE
      SUBROUTINE io_write_complex3(did,start,count,DATA,trans)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      COMPLEX(ckind),INTENT(IN)    ::DATA(:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::trans
      !locals
      INTEGER::foffset(SIZE(start))      
      foffset=start  !
      !DO 2 calls to read real values
      WHERE(start<0) 
         foffset=1
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_write_real3(did,foffset,count,REAL(DATA),trans)
      ELSE
         CALL io_write_real3(did,foffset,count,REAL(DATA))
      ENDIF
      WHERE(start<0) 
         foffset=2
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_write_real3(did,foffset,count,AIMAG(DATA),trans)
      ELSE
         CALL io_write_real3(did,foffset,count,AIMAG(DATA))
      ENDIF
      END SUBROUTINE
      SUBROUTINE io_write_complex4(did,start,count,DATA,trans)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      COMPLEX(ckind),INTENT(IN)    ::DATA(:,:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::trans
      !locals
      INTEGER::foffset(SIZE(start))      
      foffset=start  !
      !DO 2 calls to read real values
      WHERE(start<0) 
         foffset=1
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_write_real4(did,foffset,count,REAL(DATA),trans)
      ELSE
         CALL io_write_real4(did,foffset,count,REAL(DATA))
      ENDIF
      WHERE(start<0) 
         foffset=2
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_write_real4(did,foffset,count,AIMAG(DATA),trans)
      ELSE
         CALL io_write_real4(did,foffset,count,AIMAG(DATA))
      ENDIF
      END SUBROUTINE
      SUBROUTINE io_write_complex5(did,start,count,DATA,trans)
c*****************************************************************
      USE hdf5
      
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)  ::did
      INTEGER,INTENT(IN)         ::start(:),
     +     COUNT(:)
      COMPLEX(ckind),INTENT(IN)    ::DATA(:,:,:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL::trans
      !locals
      INTEGER::foffset(SIZE(start))      
      foffset=start  !
      !DO 2 calls to read real values
      WHERE(start<0) 
         foffset=1
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_write_real5(did,foffset,count,REAL(DATA),trans)
      ELSE
         CALL io_write_real5(did,foffset,count,REAL(DATA))
      ENDIF
      WHERE(start<0) 
         foffset=2
      ENDWHERE
      IF (PRESENT(trans)) THEN
         CALL io_write_real5(did,foffset,count,AIMAG(DATA),trans)
      ELSE
         CALL io_write_real5(did,foffset,count,AIMAG(DATA))
      ENDIF
      END SUBROUTINE

      FUNCTION gettransprop()RESULT(trans)
c**********************************************************************
c  local FUNCTION to get default transfer-property
c**********************************************************************
      USE hdf5
      IMPLICIT NONE
      INTEGER(HID_T)::trans
      INTEGER::hdferr
#ifdef CPP_MPI
      INCLUDE 'mpif.h'
      CALL h5pcreate_f(H5P_DATASET_XFER_F, trans, hdferr)
      CALL h5pset_dxpl_mpio_f(trans,H5FD_MPIO_INDEPENDENT_F,hdferr)
#else
      trans=H5P_DEFAULT_f
#endif
      END FUNCTION
!----------------------------------------------------------------------

      SUBROUTINE io_write_real1s(
     >                           did,start,count,DATA,stride,transprop)
c*****************************************************************
      USE hdf5
      IMPLICIT NONE

!arguments
      INTEGER(HID_T), INTENT(IN)  :: did
      INTEGER,        INTENT(IN)  :: start(:),COUNT(:),stride(:)
      REAL(rkind),    INTENT(IN)  :: DATA(:)
      INTEGER(HID_T), INTENT(IN), OPTIONAL :: transprop
!locals
      INTEGER(HSIZE_t) :: dims(7)
      INTEGER(HSIZE_t) :: foffset(SIZE(start)),fcount(SIZE(count))
      INTEGER(HSIZE_t) :: fstride(SIZE(stride))
      INTEGER(HID_t)   :: trans,fspace,memspace
      INTEGER          :: hdferr,s,n

      dims = (/SIZE(DATA,1),0,0,0,0,0,0/)   ! write 1-dim array

      IF (.NOT.PRESENT(transprop)) THEN
         trans = gettransprop()
      ELSE
         trans = transprop
      ENDIF
      foffset = start-1
      fcount  = count
      fstride = stride
      IF (ANY(count<1)) RETURN ! write nothing
      s=1                      ! check if size of count is ok
      DO n = 1, SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'

!do I/O
      CALL h5dget_space_f(
     >                    did,             ! dataset_id
     <                    fspace,hdferr)   ! dataspace_id & error
      CALL h5sselect_hyperslab_f(
     >                           fspace,H5S_SELECT_SET_F,  ! dataspace_id & operation (SET)
     >                           foffset,fcount,           ! starting point & # of blocks 
     >                           hdferr,fstride)           ! error (out) & stride (optional)
      CALL h5screate_simple_f(
     >                        1,dims,            ! rank & dimensions of dataspace
     <                        memspace,hdferr)   ! memoryspace identifier & error
      CALL h5dwrite_f(
     >                did,H5T_NATIVE_DOUBLE,     ! dataset_id, datatype_id
     >                DATA,dims,hdferr,          ! data & dimensions, error (out)
     >                memspace,fspace,trans)     ! memoryspace_id, file-dataspace_id and
                                                 ! Transfer property list id (all optional)
      CALL h5sclose_f(
     >                memspace,hdferr)
      CALL h5sclose_f(
     >                fspace,hdferr)
      CALL cleartransprop(trans)

      END SUBROUTINE io_write_real1s

!----------------------------------------------------------------------

      SUBROUTINE io_write_real2s(
     >                           did,start,count,DATA,stride,transprop)
c*****************************************************************
      USE hdf5
      IMPLICIT NONE

!arguments
      INTEGER(HID_T), INTENT(IN)  :: did
      INTEGER,        INTENT(IN)  :: start(:),COUNT(:),stride(:)
      REAL(rkind),    INTENT(IN)  :: DATA(:,:)
      INTEGER(HID_T), INTENT(IN), OPTIONAL :: transprop
!locals
      INTEGER(HSIZE_t) :: dims(7)
      INTEGER(HSIZE_t) :: foffset(SIZE(start)),fcount(SIZE(count))
      INTEGER(HSIZE_t) :: fstride(SIZE(stride))
      INTEGER(HID_t)   :: trans,fspace,memspace
      INTEGER          :: hdferr,s,n

      dims = (/SIZE(DATA,1),SIZE(DATA,2),0,0,0,0,0/)   ! write 2-dim array

      IF (.NOT.PRESENT(transprop)) THEN
         trans = gettransprop()
      ELSE
         trans = transprop
      ENDIF
      foffset = start-1
      fcount  = count
      fstride = stride
      IF (ANY(count<1)) RETURN ! write nothing
      s=1                      ! check if size of count is ok
      DO n = 1, SIZE(count)
         IF (COUNT(n)>0) s=s*COUNT(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) STOP 'Missmatch of sizes'

!     write(*,*) 'foffset',foffset
!     write(*,*) 'fcount', fcount
!     write(*,*) 'fstride',fstride
!do I/O
      CALL h5dget_space_f(
     >                    did,             ! dataset_id
     <                    fspace,hdferr)   ! dataspace_id & error
      CALL h5sselect_hyperslab_f(
     >                           fspace,H5S_SELECT_SET_F,  ! dataspace_id & operation (SET)
     >                           foffset,fcount,           ! starting point & # of blocks
     >                           hdferr,fstride)           ! error (out) & stride (optional)
      CALL h5screate_simple_f(
     >                        2,dims,            ! rank & dimensions of dataspace
     <                        memspace,hdferr)   ! memoryspace identifier & error
      CALL h5dwrite_f(
     >                did,H5T_NATIVE_DOUBLE,     ! dataset_id, datatype_id
     >                DATA,dims,hdferr,          ! data & dimensions, error (out)
     >                memspace,fspace,trans)     ! memoryspace_id, file-dataspace_id and
                                                 ! Transfer property list id (all optional)
      CALL h5sclose_f(
     >                memspace,hdferr)
      CALL h5sclose_f(
     >                fspace,hdferr)
      CALL cleartransprop(trans)

      END SUBROUTINE io_write_real2s

!----------------------------------------------------------------------
      SUBROUTINE io_write_real3s(did,start,count,data,stride,transprop)
c*****************************************************************
      USE hdf5
      IMPLICIT NONE

      INTEGER(HID_T),INTENT(IN) :: did
      INTEGER,       INTENT(IN) :: start(:),count(:),stride(:)
      REAL,          INTENT(IN) :: data(:,:,:)
      INTEGER(HID_T),INTENT(IN),OPTIONAL :: transprop
!  locals
      INTEGER(HSIZE_t) :: dims(7),foffset(size(start))
      INTEGER(HSIZE_t) :: fcount(size(count))
      INTEGER(HSIZE_t) :: fstride(SIZE(stride))
      INTEGER(HID_t)::trans,fspace,memspace
      INTEGER       ::hdferr,s,n

      IF (.not.present(transprop)) THEN
         trans=gettransprop()
      ELSE
         trans=transprop
      ENDIF
      foffset=start-1  !
      fstride = stride
      fcount=count     !
      dims=(/size(data,1),size(data,2),size(data,3),0,0,0,0/)   !
! check if size of count is ok!
      if (any(count<1)) return !read nothing
      s=1
      DO n=1,size(count)
         IF (count(n)>0) s=s*count(n)
      ENDDO
      IF (s.ne.size(data)) STOP 'Missmatch of sizes'
! DO IO
      CALL h5dget_space_f(
     >                    did,fspace,hdferr)
      CALL h5sselect_hyperslab_f(
     >                           fspace,H5S_SELECT_SET_F,
     >                           foffset,fcount,
     >                           hdferr,fstride)
      CALL h5screate_simple_f(
     >                        3,dims,
     >                        memspace,hdferr)
      CALL h5dwrite_f(
     >                did,H5T_NATIVE_DOUBLE,
     >                data,dims,hdferr,
     >                memspace,fspace,trans)
      CALL h5sclose_f(
     >                memspace,hdferr)
      CALL h5sclose_f(
     >                fspace,hdferr)
      CALL cleartransprop(trans)

      END SUBROUTINE io_write_real3s

!----------------------------------------------------------------------

      SUBROUTINE cleartransprop(trans)
c**********************************************************************
c  local FUNCTION to get default transfer-property
c**********************************************************************
      USE hdf5
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(INOUT)::trans
#ifdef CPP_MPI
      INTEGER::hdferr
      INCLUDE 'mpif.h'
      IF (trans==H5P_DEFAULT_f) RETURN
      CALL h5pclose_f(trans,hdferr)
#else
      RETURN
#endif
      END SUBROUTINE


      SUBROUTINE io_check(text,err)
c***************************************************************************
c      SUBROUTINE to check IO for error
c 
c                                  Daniel Wortmann, Juelich, 2002
c***************************************************************************
      USE hdf5
      IMPLICIT NONE    
      INTEGER,INTENT(INOUT)::err
      CHARACTER*(*)        ::text
      IF (err.GE.0) RETURN
      
      CALL h5eprint_f(err)
      WRITE(*,*) 'IO-Error detected in: ',text
      STOP 'IO-Error in hdf_tools'

      END SUBROUTINE



      SUBROUTINE checklib()
c***************************************************************************
c      SUBROUTINE to check the basic library functions
c 
c                                  Daniel Wortmann, Juelich, 2002
c***************************************************************************
      USE hdf5
      IMPLICIT NONE    
      INTEGER::hdferr,testdata
      REAL::testr
      INTEGER(HID_T)::fid,gid
      

      CALL h5fcreate_f("hdftest_tmp.hdf",H5F_ACC_TRUNC_F,fid,hdferr)
      CALL h5gcreate_f(fid,"testgroup",gid,hdferr)
      CALL io_write_att(gid,"testint",testdata)
      CALL io_write_att(gid,"testreal",testr)
      CALL h5gclose_f(gid,hdferr)
      CALL h5fclose_f(fid,hdferr)
      
      CALL h5fopen_f("hdftest_tmp.hdf",H5F_ACC_RDONLY_F,fid,hdferr)
      CALL h5gopen_f(fid,"testgroup",gid,hdferr)
      CALL io_read_att(gid,"testint",testdata)
      CALL io_read_att(gid,"testreal",testr)
      CALL h5gclose_f(gid,hdferr)
      CALL h5fclose_f(fid,hdferr)

      OPEN(99,file="hdftest_tmp.hdf")
      CLOSE(99)

      WRITE(*,*) 'HDF library was initialized'

      END SUBROUTINE

      END MODULE

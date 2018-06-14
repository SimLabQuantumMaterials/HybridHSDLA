cmake_minimum_required(VERSION 2.8)
project(Flapw)
include(CTest)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Load MPI and HDF5 packages
find_package(HDF5 REQUIRED COMPONENTS C HL)
include_directories(${HDF5_INCLUDE_DIRS})
# add HDF5 preprocessor and explicitly skip C++ bindings
# of openmpi, which are pulled in through #include <mpi.h>
# inside an included hdf5 header file and cause a problem
# with missing linkage to MPI::Op::Init(...), a static initialization
# triggered through some cmake/hdf5/cluster configuration "smartness".
# we don't use MPI anyway (yet), and if we would, it would be the C bindings.
add_definitions(${HDF5_DEFINITIONS} -DOMPI_SKIP_MPICXX)
set(LIB_hdf5
        ${HDF5_LIBRARIES}
        ${HDF5_C_LIBRARIES}
        ${HDF5_HL_LIBRARIES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -pedantic -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Wall")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O1 -Wall")

if(CMAKE_BUILD_TYPE MATCHES Release)
     # disable assert() and armadillo bound checks
    add_definitions(-DNDEBUG -DARMA_NO_DEBUG)
endif(CMAKE_BUILD_TYPE MATCHES Release)

# Target libs/system
if(TARGET_IMPLEMENTATION MATCHES CPU_STRIPPED)
    add_definitions(-DSTRIPPED)
elseif(TARGET_IMPLEMENTATION MATCHES HYBRID)
    find_package(CUDA)
    link_directories(${CUDA_TOOLKIT_ROOT_DIR}/lib64) 
    add_definitions(-DSTRIPPED -DGPU -DOOC -DMKL)
elseif(TARGET_IMPLEMENTATION MATCHES DUMP_FLAPW_DATA)
    add_definitions(-DSTRIPPED -DDUMP_FLAPW_DATA)
elseif(TARGET_IMPLEMENTATION MATCHES HS_CUDA)
    find_package(CUDA)
    link_directories(${CUDA_TOOLKIT_ROOT_DIR}/lib64) 
    add_definitions(-DSTRIPPED -DGPU -DHS_CUDA)
elseif(TARGET_IMPLEMENTATION MATCHES HS_MIC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -restrict -fPIC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -restrict -fPIC")
    add_definitions(-DSTRIPPED -DHS_HSTREAMS)
    include_directories(${PROJECT_SOURCE_DIR}/src/hetero-mic-gpu/)
    include_directories(${PROJECT_SOURCE_DIR}/src/hetero-mic-gpu/hetero-streams/include)
endif(TARGET_IMPLEMENTATION MATCHES CPU_STRIPPED)

include_directories(${PROJECT_SOURCE_DIR}/src/include)
include_directories(${PROJECT_SOURCE_DIR}/src/include/acc)
include_directories(${PROJECT_SOURCE_DIR}/external/easyloggingpp/src)
include_directories(${PROJECT_SOURCE_DIR}/external/armadillo/include)
include_directories(${PROJECT_SOURCE_DIR}/external/UnitTest++)
add_subdirectory(${PROJECT_SOURCE_DIR}/external/armadillo)
add_subdirectory(${PROJECT_SOURCE_DIR}/external/compat)
add_subdirectory(${PROJECT_SOURCE_DIR}/external/UnitTest++)
add_subdirectory(${PROJECT_SOURCE_DIR}/src)
# unit test framework
add_custom_target(build_and_test ${CMAKE_CTEST_COMMAND})
add_subdirectory(${PROJECT_SOURCE_DIR}/tests
                    EXCLUDE_FROM_ALL)

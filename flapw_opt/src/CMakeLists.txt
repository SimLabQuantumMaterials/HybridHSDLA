# put binaries in top level build folder
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# configure paths to samples directory
set(SAMPLES_DIR ${CMAKE_SOURCE_DIR}/tests/sample_input)
configure_file(include/configured_input.h.in include/configured_input.h)
# configure_file() places the output in the *corresponding* binary dir, 
# so we need to add it to our include path
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)
include_directories(${CUDA_INCLUDE_DIRS}) #/rwthfs/rz/SW/cuda/7.5.18/RHEL_7/cuda/include)
include_directories(${MKLROOT}/include)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

# move core flapw logic into lib
add_library(FLAPW STATIC types_helper.cpp tensor_AB_match.cpp special_funcs.cpp
                         FleurData.cpp binmat.cpp AtomData.cpp HDF5Store.cpp)
# compare smat
add_executable(compare_smat compare_smat.cpp binmat.cpp)

if(TARGET_IMPLEMENTATION MATCHES BLASXT)
    set(SPECIFIC_OBJS h_and_s_stripped_v2.cpp wrappers.cpp timing.c)
    set(SPECIFIC_LIBS ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})
elseif(TARGET_IMPLEMENTATION MATCHES BLAS)
    set(SPECIFIC_OBJS "")
    #set(SPECIFIC_LIBS "")
    set(SPECIFIC_LIBS "-mkl")
elseif(TARGET_IMPLEMENTATION MATCHES CPU_STRIPPED)
    set(SPECIFIC_OBJS h_and_s_stripped_v2.cpp timing.c wrappers.cpp)
    set(SPECIFIC_LIBS "-mkl")
elseif(TARGET_IMPLEMENTATION MATCHES HYBRID)
    set(SPECIFIC_OBJS h_and_s_stripped_v2.cpp wrappers.cpp timing.c 
                      hybrid/zherk_hyb.c hybrid/zher2k_hyb.c hybrid/zherkx_hyb.c hybrid/utils.c)
    set(SPECIFIC_LIBS ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})
elseif(TARGET_IMPLEMENTATION MATCHES DUMP_FLAPW_DATA)
    add_definitions(-DSTRIPPED -DDUMP_FLAPW_DATA)
elseif(TARGET_IMPLEMENTATION MATCHES HS_CUDA)
    set(SPECIFIC_OBJS h_and_s_stripped_v2.cpp wrappers.cpp timing.c
                      hetero-mic-gpu/source.cpp)
    set(SPECIFIC_LIBS ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})
elseif(TARGET_IMPLEMENTATION MATCHES HS_MIC)
    set(SPECIFIC_OBJS h_and_s_stripped_v2.cpp wrappers.cpp timing.c
                      hetero-mic-gpu/source.cpp)
    link_directories(${PROJECT_SOURCE_DIR}/src/hetero-mic-gpu/hetero-streams/bin/host) 
    set(SPECIFIC_LIBS "-lhstreams_source -lrt")

    add_library(main_mic SHARED ${PROJECT_SOURCE_DIR}/src/hetero-mic-gpu/sink.c)
    set_target_properties(main_mic PROPERTIES PREFIX "")
    set_target_properties(main_mic PROPERTIES LINK_FLAGS -mmic)
    target_compile_options(main_mic PRIVATE -mmic)
endif(TARGET_IMPLEMENTATION MATCHES BLASXT)

# h and s
add_executable(h_and_s h_and_s.cpp logging.cpp ${SPECIFIC_OBJS})
target_link_libraries(h_and_s FLAPW)
target_link_libraries(h_and_s ${LIB_armadillo} flerp ${LIB_hdf5})
target_link_libraries(h_and_s ${SPECIFIC_LIBS})

install(TARGETS h_and_s
        CONFIGURATIONS Release
        DESTINATION bin)


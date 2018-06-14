# Automatically add an executable with given name and sources,
# linking in the Flapw-Libraries and create a CTest target
function(utest name sources)
    set(exename ${name}-Unit)
    add_executable(${exename} ${sources})
    target_link_libraries(${exename}
                            UnitTest++
                            FLAPW
                            ${LIB_armadillo}
                            flerp
                            ${LIB_hdf5})
    add_test(${name} ${exename})
    add_dependencies(build_and_test ${exename})
endfunction(utest)

utest(FleurData FleurData_test.cpp)
# TODO: mostly, these tests are only smoke tests
# and do not test against predefined test data (which would be feasible)
utest(FleurWrap FleurWrap_test.cpp)
utest(Binmat binmat_test.cpp)
# TODO: this one is mostly a placeholder
utest(Config conf_test.cpp)

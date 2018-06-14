set(unittestpp_SRCS
    AssertException.cpp
    Test.cpp
    Checks.cpp
    TestRunner.cpp
    TestResults.cpp
    TestReporter.cpp
    TestReporterStdout.cpp
    ReportAssert.cpp
    TestList.cpp
    TimeConstraint.cpp
    TestDetails.cpp
    MemoryOutStream.cpp
    DeferredTestReporter.cpp
    DeferredTestResult.cpp
    XmlTestReporter.cpp
    CurrentTest.cpp)

if (${WIN32})
    set(unittestpp_SRCS ${unittestpp_SRCS}
        Win32/TimeHelpers.cpp)
else()
    set(unittestpp_SRCS ${unittestpp_SRCS}
        Posix/SignalTranslator.cpp
        Posix/TimeHelpers.cpp)
endif()

add_library(UnitTest++ ${unittestpp_SRCS})
